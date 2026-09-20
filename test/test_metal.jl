#= Metal hardware tests.

   Metal must not be a dependency of Luna, so it goes into a *separate* environment
   stacked on the package (which is what the CI job builds), never into Project.toml.
   That environment needs Luna developed into it plus the packages this file imports
   directly:

       julia --project=<env> -e '
         using Pkg
         Pkg.develop(path=".")
         Pkg.add(["Metal", "Test", "GPUArraysCore", "Adapt"])'
       julia --project=<env> -t 1 -e 'using Luna, Metal; include("test/test_metal.jl")'

   Without Metal loaded, or on a machine where it is not functional, everything here is
   skipped.

   This is the only reliable detector of a stray Float64: Metal refuses Float64 arrays
   outright and its kernel compiler rejects any `double` which survives optimisation, so
   a Float64 struct field or scalar which JLArrays and a Float32 `Array` silently promote
   fails here. Every branch which adds or changes a kernel runs this file. =#

import Test: @test, @testset, @test_throws
import Luna
import Luna: Utils, Output, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear,
             NonlinearRHS, PhysData, RK45, DeviceSpec, HostSpec
import LinearAlgebra
import GPUArraysCore
import Adapt

have_metal = try
    Metal = Base.require(Base.PkgId(
        Base.UUID("dde4c033-4e86-420c-a63e-0dd931031962"), "Metal"))
    Base.invokelatest(getproperty(Metal, :functional))
catch
    false
end

if !have_metal
    @warn "Metal is not loaded or not functional; test_metal.jl is skipped."
else

const Metal = Base.require(Base.PkgId(
    Base.UUID("dde4c033-4e86-420c-a63e-0dd931031962"), "Metal"))
const MtlArray = getproperty(Metal, :MtlArray)
const MetalSpec = DeviceSpec(MtlArray, Float32)

# Scalar indexing on a device array is a silent, catastrophic slowdown; nothing Luna runs
# per step may do it.
GPUArraysCore.allowscalar(false)

#= An output which copies what the solver hands it down to the host. The real one
   (`ScaledOutput`, which also unscales) is gpu/11's. =#
struct ToHostM{O}
    o::O
end
(h::ToHostM)(y, t, dt, yfun) = h.o(Array(y), t, dt, ti -> Array(yfun(ti)))
(h::ToHostM)(args...; kwargs...) = h.o(args...; kwargs...)
Base.getindex(h::ToHostM, k) = h.o[k]

#= The same mode-averaged Kerr propagation on whichever spec is asked for. `boundary=:none`
   because the absorbers are host code until gpu/11. =#
function metalcase(GT, spec; gas=:He, pres=1.0, energy=1e-6, flength=1e-2, λ0=800e-9,
                   precision=nothing, thg=false)
    grid = GT === Grid.RealGrid ?
        Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15) :
        Grid.EnvGrid(λ0, (300e-9, 2000e-9), 400e-15; thg)
    m = Capillary.MarcatiliMode(75e-6, gas, pres, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    #= Precomputed: `PhysData.density` goes through CoolProp and costs more per call
       than the whole right-hand side. =#
    ρ = PhysData.density(gas, pres)
    dens = z -> ρ
    resp = if GT === Grid.RealGrid
        (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    elseif thg
        (Nonlinear.Kerr_env_thg(PhysData.γ3_gas(gas), grid.ω0, grid.to),)
    else
        (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),)
    end
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0;
                                                    (GT === Grid.EnvGrid ? (; thg) : ())...)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=spec, precision)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    output = Utils.isdevice(Eω) ? ToHostM(out) : out
    Luna.run(Eω, grid, linop, transform, FT, output;
             zmax=flength, boundary=:none, init_dz=flength/20, rtol=1e-8)
    out, transform
end

#= A pressure gradient: a z-dependent operator closure and `constβ=false`, which is the
   only case exercising `NormModeAvg`'s `HostMirror` branch and `RK45.make_prop!`'s
   host-buffer branch -- the two pieces which still upload from the host on every stage
   until gpu/23. Fixed steps, so the runs differ only in arithmetic. =#
function metalgradientcase(spec; gas=:Ar, pin=1.0, pout=0.0, flength=1e-2, λ0=800e-9)
    grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15)
    coren, densityfun = Capillary.gradient(gas, flength, pin, pout)
    m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun! = LinearOps.make_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=1e-6)
    Eω, transform, FT = Luna.setup(grid, densityfun, resp, inputs, βfun!, aeff;
                                   device=spec)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    output = Utils.isdevice(Eω) ? ToHostM(out) : out
    dz = flength/20
    Luna.run(Eω, grid, linop, transform, FT, output;
             zmax=flength, boundary=:none, init_dz=dz, min_dz=dz, max_dz=dz)
    out, transform
end

@testset "Metal registration" begin
    @test haskey(Luna.DEVICES, :metal)
    @test Luna.DEVICES[:metal].spec === MetalSpec
    #= Loading Metal sets settings["device"] = :auto, and :auto then picks the GPU. An
       explicit set_device(:cpu) must still win. =#
    @test Luna.resolve_device(:auto) === MetalSpec
    @test Luna.resolve_device(:metal) === MetalSpec
    old = get(Luna.settings, "device", nothing)
    try
        Luna.set_device(:cpu)
        @test Luna.device() === HostSpec()
        Luna.set_device(:auto)
        @test Luna.device() === MetalSpec
    finally
        isnothing(old) ? delete!(Luna.settings, "device") :
                         (Luna.settings["device"] = old)
    end
    @test Luna.device_synchronize(MetalSpec) === nothing
    ms = Luna.device_memory_status(MetalSpec)
    @test ms isa Tuple{Integer, Integer}
end

@testset "allocation and transfer on Metal" begin
    x = Luna.alloc(MetalSpec, ComplexF32, (8,))
    @test x isa MtlArray{ComplexF32, 1}
    v = rand(5)
    dv = Luna.todevice(MetalSpec, v)
    @test dv isa MtlArray{Float32, 1}
    @test Array(dv) ≈ v
    b = BitVector([true, false, true])
    @test Luna.todevice(MetalSpec, b) isa MtlArray{Bool, 1}
    @test Luna.assert_resident(MetalSpec, x, dv) === nothing
    @test_throws ErrorException Luna.assert_resident(MetalSpec, v)
    # Float64 never reaches the device
    @test_throws Exception MtlArray(rand(4))
end

#= One broadcast per kernel-touching struct, on the device. This is the stray-Float64
   smoke test: if any field of the struct, or any literal in the kernel body, is a
   Float64 which survives optimisation, the Metal compiler refuses to build the kernel
   and this throws. =#
@testset "no stray Float64 in the kernels" begin
    n = 64
    ρ = PhysData.density(:He, 1.0)
    sc = Luna.UnitScaling(1024.0, PhysData.ε_0)

    # Kerr, real field, scalar and vector
    kf = Nonlinear.rescale(Nonlinear.Kerr_field(PhysData.γ3_gas(:He)), MetalSpec, sc)
    @test kf isa Nonlinear.KerrField{Float32}
    E = Luna.todevice(MetalSpec, rand(n))
    out = Luna.alloc(MetalSpec, Float32, (n,))
    kf(out, E, ρ)
    @test all(isfinite, Array(out))
    Ev = Luna.todevice(MetalSpec, rand(n, 2))
    outv = Luna.alloc(MetalSpec, Float32, (n, 2))
    kf(outv, Ev, ρ)
    @test all(isfinite, Array(outv))

    # Kerr, envelope, scalar and vector
    ke = Nonlinear.rescale(Nonlinear.Kerr_env(PhysData.γ3_gas(:He)), MetalSpec, sc)
    @test ke isa Nonlinear.KerrEnv{Float32}
    Ec = Luna.todevice(MetalSpec, rand(ComplexF64, n))
    outc = Luna.alloc(MetalSpec, ComplexF32, (n,))
    ke(outc, Ec, ρ)
    @test all(isfinite, Array(outc))
    Ecv = Luna.todevice(MetalSpec, rand(ComplexF64, n, 2))
    outcv = Luna.alloc(MetalSpec, ComplexF32, (n, 2))
    ke(outcv, Ecv, ρ)
    @test all(isfinite, Array(outcv))

    # Kerr with THG: carries an array, so it also has an Adapt rule
    t = collect(range(0, 1e-13, length=n))
    kt = Nonlinear.rescale(Nonlinear.Kerr_env_thg(PhysData.γ3_gas(:He), 2.35e15, t),
                           MetalSpec, sc)
    @test kt.C isa MtlArray{ComplexF32, 1}
    @test kt.γ3 isa Float32
    fill!(outc, 0)
    kt(outc, Ec, ρ)
    @test all(isfinite, Array(outc))
    @test Adapt.adapt(MtlArray, kt).C isa MtlArray

    # The RK45 per-step kernels
    ks = ntuple(_ -> Luna.todevice(MetalSpec, rand(ComplexF64, n)), 7)
    y = Luna.todevice(MetalSpec, rand(ComplexF64, n))
    yn = similar(y)
    RK45.combine!(yn, y, ks, 1.3e-4, RK45.B[3], 3)
    RK45.combine!(yn, y, ks, 1.3e-4, RK45.b5, 7)
    RK45.combine!(yn, y, ks, 1.3e-4, RK45.b4, 7)
    yerr = similar(y)
    RK45.errorestimate!(yerr, ks, 1.3e-4)
    @test all(isfinite, Array(yerr))
    for nrm in (RK45.weaknorm, RK45.maxnorm, RK45.maxnorm_ratio, RK45.normnorm)
        v = nrm(yerr, y, yn, 1e-6, 1e-10)
        @test isfinite(v)
    end

    # The constant-operator propagator
    linop = Luna.todevice(MetalSpec, -1e-3im .* rand(n))
    prop! = RK45.make_prop!(linop, y)
    prop!(y, 0.0, 1e-4)
    @test all(isfinite, Array(y))
end

@testset "mode-averaged Kerr on Metal" begin
    #= Metal against the CPU at the *same* precision and the same scaling: what is being
       tested is the device path, not single precision. The remaining difference is the
       order of the FFT and of the reductions. =#
    for GT in (Grid.RealGrid, Grid.EnvGrid)
        href, htr = metalcase(GT, DeviceSpec(Array, Float32))
        dref, dtr = metalcase(GT, MetalSpec)

        @test dtr.Eto isa MtlArray
        @test dtr.Eωo isa MtlArray{ComplexF32}
        @test dtr.gv.ω isa MtlArray{Float32}
        @test dtr.gv.sidx isa MtlArray{Bool}
        @test dtr.norm!.pre isa MtlArray{ComplexF32}
        @test dtr.scaling.Eref == htr.scaling.Eref

        @test size(dref["Eω"]) == size(href["Eω"])
        for idx in axes(href["Eω"], 2)
            h = href["Eω"][:, idx]
            d = dref["Eω"][:, idx]
            @test maximum(abs, d .- h)/maximum(abs, h) < 1e-4
        end
    end
end

@testset "Metal against the Float64 CPU path" begin
    #= The dynamic-range case: helium at 0.3 bar, where the unscaled Kerr coefficient is
       below the smallest Float32 subnormal and Metal flushes subnormals to zero. This is
       what the unit scaling exists for. =#
    @test PhysData.density(:He, 0.3)*PhysData.ε_0*PhysData.γ3_gas(:He) < floatmin(Float32)

    href, _ = metalcase(Grid.RealGrid, HostSpec(); pres=0.3)
    dref, dtr = metalcase(Grid.RealGrid, MetalSpec; pres=0.3)
    Eref = dtr.scaling.Eref
    for idx in axes(href["Eω"], 2)
        h = href["Eω"][:, idx]
        d = ComplexF64.(dref["Eω"][:, idx]) .* Eref
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-4
    end
end

@testset "pressure gradient on Metal" begin
    href, htr = metalgradientcase(DeviceSpec(Array, Float32))
    dref, dtr = metalgradientcase(MetalSpec)

    # The z-dependent branch: β is staged on the host and uploaded per evaluation
    @test !isnothing(dtr.norm!.β)
    @test dtr.norm!.β.host isa Vector{Float64}
    @test dtr.norm!.β.stage isa Vector{Float32}
    @test dtr.norm!.β.dev isa MtlArray{Float32, 1}

    for idx in axes(href["Eω"], 2)
        h = href["Eω"][:, idx]
        d = dref["Eω"][:, idx]
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-4
    end
end

#= Finding 1 of review round 1: loading Metal sets settings["device"] = :auto for the
   whole process, and the simple interface is not device-capable in this branch. It must
   therefore give the same answer whatever the setting says, which is what
   `Interface` passing `device=Luna.HostSpec()` guarantees until gpu/11 plumbs the
   keywords through. =#
@testset "the simple interface stays on the CPU" begin
    old = get(Luna.settings, "device", nothing)
    capargs = (125e-6, 1e-3, :He, 1.0)
    capkw = (; λ0=800e-9, energy=1e-9, τfwhm=10e-15, λlims=(300e-9, 2e-6),
             trange=400e-15, saveN=3, plasma=false, shotnoise=false,
             PPT_options=Dict(:cache => false))
    gnlsekw = (; λ0=835e-9, τfwhm=100e-15, power=1e3, pulseshape=:sech,
               λlims=(450e-9, 2e-6), trange=1e-12, saveN=3, raman=false,
               shotnoise=false)
    try
        Luna.set_device(:cpu)
        ocap = Luna.prop_capillary(capargs...; capkw...)
        ognlse = Luna.prop_gnlse(0.1, 1e-3, [0.0, 0.0, -1e-26]; gnlsekw...)

        Luna.set_device(:auto)
        @test Luna.device() === MetalSpec # the GPU really is selected globally
        dcap = Luna.prop_capillary(capargs...; capkw...)
        dgnlse = Luna.prop_gnlse(0.1, 1e-3, [0.0, 0.0, -1e-26]; gnlsekw...)

        # Same answer, on the host, in double precision
        @test eltype(dcap["Eω"]) === ComplexF64
        @test dcap["Eω"] == ocap["Eω"]
        @test eltype(dgnlse["Eω"]) === ComplexF64
        @test dgnlse["Eω"] == ognlse["Eω"]
    finally
        isnothing(old) ? delete!(Luna.settings, "device") :
                         (Luna.settings["device"] = old)
    end
end

@testset "Metal refuses what it cannot run" begin
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, :He, 1.0, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    dens = z -> PhysData.density(:He, 1.0)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, 800e-9)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=MetalSpec)
    out = ToHostM(Output.MemoryOutput(0, 1e-3, 3, Output.nostats))
    # The absorbers are host scalar code until gpu/11
    @test_throws ErrorException Luna.run(Eω, grid, linop, transform, FT, out;
                                         zmax=1e-3, boundary=:rate)
end

end # have_metal
