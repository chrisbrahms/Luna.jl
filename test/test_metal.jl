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

import Test: @test, @testset, @test_throws, @test_logs
import Luna
import Luna: Utils, Output, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear,
             NonlinearRHS, PhysData, RK45, Stats, Boundaries, DeviceSpec, HostSpec,
             UnitScaling, UNIT_SCALING, Maths, Ionisation, Raman
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

#= A user-written columnwise response: a plain closure over the contract Luna has always
   had, with no knowledge of devices, precision or units. On a device run
   `Nonlinear.HostResponse` is what has to make it work, and this is the only place that
   is checked on real hardware. =#
usercubic(s) = let s = s
    (out, E, ρ) -> (out .+= (ρ*s) .* E.^3)
end

#= The plasma response on Metal. Two rates: the analytic ADK formula, whose nine
   constants are struct fields a kernel reads, and a tabulated rate with the interface of
   a cached PPT rate, whose kernel indexes a spline. Both are the stray-Float64 case the
   plan names (GPU_PLAN.md section 4.1): a Float64 field read inside a kernel never
   compiles on Metal, so if `device_rate` missed one, this file is where it shows.

   The table is built here rather than pre-calculated: `IonRatePPTAccel(E, rate)` is the
   constructor the cache calls, the axis is uniform, and this takes milliseconds. =#
metal_adkrate() = Ionisation.IonRateADK(:Ar)

function metal_tablerate()
    E = collect(range(1e9, 3e11, length=1024))
    Ionisation.IonRatePPTAccel(E, metal_adkrate().(E))
end

#= A field which ionises, on a short grid: a 10 fs pulse at 800 nm, peak field `E0`. =#
function metal_plasmafield(nt=512, E0=6e10, twidth=60e-15)
    t = collect(range(-twidth, twidth, length=nt))
    t, @. E0*exp(-t^2/(2*(10e-15/1.66)^2))*cos(2π*PhysData.c/800e-9*t)
end

#= The same mode-averaged Kerr propagation on whichever spec is asked for. `Luna.run`
   wraps the output in `ScaledOutput` itself now (gpu/11), so this test file never needs
   its own host-copy wrapper the way it did under gpu/10 -- `out`, the plain
   `Output.MemoryOutput` this builds, already holds host, physical-unit data when the
   propagation returns. `boundary=:none` by default; the `:rate` case (`RateAbsorber`,
   `Boundaries.jl`) has its own testset below, since it is the whole point of this
   branch. =#
function metalcase(GT, spec; gas=:He, pres=1.0, energy=1e-6, flength=1e-2, λ0=800e-9,
                   precision=nothing, thg=false, boundary=:none, stats=false, fixed=false,
                   extraresp=(), plasma=false, raman=false, nothg=false)
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
        #= `thg=false` on a `RealGrid` is `prop_capillary`'s no-THG Kerr response, which
           needs the analytic signal of the block and is therefore batched. =#
        (nothg ? Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas), length(grid.to)) :
                 Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    elseif thg
        (Nonlinear.Kerr_env_thg(PhysData.γ3_gas(gas), grid.ω0, grid.to),)
    else
        (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),)
    end
    #= The default response set of a field-resolved `prop_capillary` call in a
       non-Raman gas is Kerr and plasma, which is what gpu/13 has to make run on a
       device. The rate is tabulated here rather than pre-calculated; what is under test
       is the response, not the PPT series. =#
    if plasma
        resp = (resp...,
                Nonlinear.PlasmaCumtrapz(grid.to, zeros(length(grid.to)),
                                         metal_tablerate(),
                                         PhysData.ionisation_potential(gas)))
    end
    #= The default response set of a molecular gas: Kerr and the Raman polarisation,
       which is what gpu/14 has to make run on a device. =#
    if raman
        rr = Raman.raman_response(grid.to, gas)
        resp = (resp..., GT === Grid.RealGrid ?
                         Nonlinear.RamanPolarField(grid.to, rr; thg=!nothg) :
                         Nonlinear.RamanPolarEnv(grid.to, rr))
    end
    resp = (resp..., extraresp...)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0;
                                                    (GT === Grid.EnvGrid ? (; thg) : ())...)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=spec, precision)
    #= Stats.jl is host-only: its EnvGrid plan_analytic plans an FFTW transform
       directly on a copy of the given Eω, so construction needs a host-shaped
       template, not the device state itself (found here, on real hardware). =#
    shost = Utils.isdevice(Eω) ? Luna.tohost(Eω) : Eω
    statsfun = stats ? Stats.default(grid, shost, m, linop, transform; gas) : Output.nostats
    out = Output.MemoryOutput(0, flength, 5, statsfun)
    #= `fixed=true`: min_dz == max_dz == init_dz bypasses the step-size controller
       (RK45.steplims!), so every difference between two runs is attributable to the
       arithmetic rather than to a different sequence of steps -- the same reasoning as
       the regression gate's own `:fixed` mode. Review round 1, finding 1: needed to keep
       a tight (1e-4) tolerance meaningful once the state includes boundaries/statistics
       and the run is long/strong enough for the Kerr effect to be visible. =#
    dz = flength/20
    if fixed
        Luna.run(Eω, grid, linop, transform, FT, out;
                 zmax=flength, boundary, init_dz=dz, min_dz=dz, max_dz=dz)
    else
        Luna.run(Eω, grid, linop, transform, FT, out;
                 zmax=flength, boundary, init_dz=dz, rtol=1e-8)
    end
    out, transform
end

#= A pressure gradient: a z-dependent operator closure and `constβ=false`, which is the
   only case exercising `NormModeAvg`'s `HostMirror` branch and `RK45.make_prop!`'s
   host-buffer branch -- the two pieces which still upload from the host on every stage
   until gpu/23. Fixed steps, so the runs differ only in arithmetic. =#
function metalgradientcase(spec; gas=:Ar, pin=1.0, pout=0.0, energy=1e-6, flength=1e-2,
                           λ0=800e-9, boundary=:none, stats=false, fixed=true)
    grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15)
    coren, densityfun = Capillary.gradient(gas, flength, pin, pout)
    m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun! = LinearOps.make_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, densityfun, resp, inputs, βfun!, aeff;
                                   device=spec)
    #= Stats.jl is host-only: its EnvGrid plan_analytic plans an FFTW transform
       directly on a copy of the given Eω, so construction needs a host-shaped
       template, not the device state itself (found here, on real hardware). =#
    shost = Utils.isdevice(Eω) ? Luna.tohost(Eω) : Eω
    statsfun = stats ? Stats.default(grid, shost, m, linop, transform; gas) : Output.nostats
    out = Output.MemoryOutput(0, flength, 5, statsfun)
    dz = flength/20
    if fixed
        Luna.run(Eω, grid, linop, transform, FT, out;
                 zmax=flength, boundary, init_dz=dz, min_dz=dz, max_dz=dz)
    else
        Luna.run(Eω, grid, linop, transform, FT, out;
                 zmax=flength, boundary, init_dz=dz, rtol=1e-8)
    end
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

    #= The responses go through `Et_to_Pt!`, which is how a transform reaches them: the
       fused broadcast is the kernel, and its scalar coefficients are what must not be
       Float64. The response structs themselves keep their physical `Float64` constants
       and never enter a kernel (gpu/12). =#
    kf = Nonlinear.rescale(Nonlinear.Kerr_field(PhysData.γ3_gas(:He)), MetalSpec, sc)
    @test kf isa Nonlinear.KerrField{Float64}
    E = Luna.todevice(MetalSpec, rand(n))
    out = Luna.alloc(MetalSpec, Float32, (n,))
    @test Nonlinear.pointwise_kernel(kf, out, ρ, sc)(1f0) isa Float32
    NonlinearRHS.Et_to_Pt!(out, E, (kf,), ρ; scaling=sc)
    @test all(isfinite, Array(out))
    @test !all(iszero, Array(out)) # the scaled coefficient is not flushed to zero
    Ev = Luna.todevice(MetalSpec, rand(n, 2))
    outv = Luna.alloc(MetalSpec, Float32, (n, 2))
    NonlinearRHS.Et_to_Pt!(outv, Ev, (kf,), ρ; scaling=sc)
    @test all(isfinite, Array(outv))
    @test !all(iszero, Array(outv))

    # Kerr, envelope, scalar and vector
    ke = Nonlinear.rescale(Nonlinear.Kerr_env(PhysData.γ3_gas(:He)), MetalSpec, sc)
    Ec = Luna.todevice(MetalSpec, rand(ComplexF64, n))
    outc = Luna.alloc(MetalSpec, ComplexF32, (n,))
    NonlinearRHS.Et_to_Pt!(outc, Ec, (ke,), ρ; scaling=sc)
    @test all(isfinite, Array(outc))
    Ecv = Luna.todevice(MetalSpec, rand(ComplexF64, n, 2))
    outcv = Luna.alloc(MetalSpec, ComplexF32, (n, 2))
    NonlinearRHS.Et_to_Pt!(outcv, Ecv, (ke,), ρ; scaling=sc)
    @test all(isfinite, Array(outcv))
    @test !all(iszero, Array(outcv))

    # Two pointwise responses in one fused broadcast
    fill!(outv, 0)
    NonlinearRHS.Et_to_Pt!(outv, Ev, (kf, Nonlinear.Kerr_field(2PhysData.γ3_gas(:He))),
                           ρ; scaling=sc)
    @test all(isfinite, Array(outv))

    # Kerr with THG: carries an array, so it also has an Adapt rule
    t = collect(range(0, 1e-13, length=n))
    kt = Nonlinear.rescale(Nonlinear.Kerr_env_thg(PhysData.γ3_gas(:He), 2.35e15, t),
                           MetalSpec, sc)
    @test kt.C isa MtlArray{ComplexF32, 1}
    NonlinearRHS.Et_to_Pt!(outc, Ec, (kt,), ρ; scaling=sc)
    @test all(isfinite, Array(outc))
    @test !all(iszero, Array(outc))
    @test Adapt.adapt(MtlArray, kt).C isa MtlArray

    #= The batched responses of gpu/14: the Raman convolution and the no-THG Kerr. Their
       kernels are broadcasts and planned FFTs over the block, and the scalars they carry
       (the frequency-domain factor and the density) go through `Luna.scalar`. The
       frequency-domain response function is the array a kernel broadcasts against, and
       the only part of the kernel machinery which leaves the host. =#
    tr = collect(range(-100e-15, 100e-15, length=n))
    ρn = PhysData.density(:N2, 1.0)
    Er = Luna.todevice(MetalSpec, rand(n) .- 0.5)
    outr = Luna.alloc(MetalSpec, Float32, (n,))
    #= Its own scaling, with an `E_ref` of the size `Luna.unitscaling` picks for a real
       pulse (`8.6e9` V/m). The Raman polarisation in scaled units carries `E_ref^2`, so
       with the `E_ref = 1024` of the rest of this testset -- a field of a kilovolt per
       metre, which nothing here propagates -- the frequency-domain product is genuinely
       around 1e-44 and underflows in Float32. See PR_14-raman.md's audit. =#
    scr = Luna.UnitScaling(exp2(33), PhysData.ε_0)
    for thg in (true, false)
        Rd = Nonlinear.rescale(
            Nonlinear.RamanPolarField(tr, Raman.raman_response(tr, :N2); thg),
            MetalSpec, scr, Er)
        @test Rd.hω isa MtlArray{ComplexF32, 1}
        @test Rd.E2 isa MtlArray{Float32}
        @test Luna.all_resident(MetalSpec, Nonlinear.resident_arrays(Rd)...)
        fill!(outr, 0)
        NonlinearRHS.Et_to_Pt!(outr, Er, (Rd,), ρn; scaling=scr)
        @test all(isfinite, Array(outr))
        @test !all(iszero, Array(outr))
        # the scaled frequency-domain factor is a normal Float32, not a flushed zero
        hfac = Luna.scalar(outr, Nonlinear.coefficients(Rd, ρn, scr)[1])
        @test hfac isa Float32
        @test floatmin(Float32) < abs(hfac) < floatmax(Float32)
        @test floatmin(Float32) < maximum(abs, Array(Rd.hω)) < floatmax(Float32)
    end
    Red = Nonlinear.rescale(Nonlinear.RamanPolarEnv(tr, Raman.raman_response(tr, :N2)),
                            MetalSpec, scr, Ec)
    fill!(outc, 0)
    NonlinearRHS.Et_to_Pt!(outc, Ec, (Red,), ρn; scaling=scr)
    @test all(isfinite, Array(outc))
    @test !all(iszero, Array(outc))

    knd = Nonlinear.rescale(Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(:He), n),
                            MetalSpec, sc, Er)
    @test knd.an.mask isa MtlArray{Float32, 1}
    @test knd.an.c1 isa MtlArray{ComplexF32}
    fill!(outr, 0)
    NonlinearRHS.Et_to_Pt!(outr, Er, (knd,), ρ; scaling=sc)
    @test all(isfinite, Array(outr))
    @test !all(iszero, Array(outr))

    #= The ionisation rates: nine Float64 constants in the ADK struct, and a spline
       whose knots, values, coefficients and index function are all Float64 as built.
       `device_rate` is what converts them, and this is the only place which proves it
       missed none -- a Float64 struct field read inside a kernel never compiles here. =#
    Ei = Luna.todevice(MetalSpec, collect(range(1e9, 2e11, length=n)))
    ri = Luna.alloc(MetalSpec, Float32, (n,))
    for ir in (Ionisation.IonRateADK(:Ar), metal_tablerate())
        ird = Ionisation.device_rate(ir, MetalSpec)
        Ionisation.ionrate!(ri, ird, Ei)
        @test all(isfinite, Array(ri))
        @test !all(iszero, Array(ri))
        #= The kernel's return type, on a *host* Float32 copy of the same rate: the
           device one indexes device arrays and cannot be called from here. =#
        irh = Ionisation.device_rate(ir, DeviceSpec(Array, Float32))
        @test Ionisation.ratekernel(irh, 1f0)(1f10) isa Float32
        # ... and with a unit scaling, where the kernel reconstructs the physical field
        fill!(ri, 0)
        Ionisation.ionrate!(ri, ird, Ei ./ 1f10, 1e10)
        @test all(isfinite, Array(ri))
        @test !all(iszero, Array(ri))
    end
    @test Adapt.adapt(MtlArray, Ionisation.device_rate(metal_tablerate(),
                                                       MetalSpec)).spline.x isa MtlArray

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

    # gpu/11: the boundary kernels. αt/αxy mirrored to Float32 MtlArray, no stray Float64.
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    zmax = 1e-2
    αt = Boundaries.temporal_rate(grid, zmax)
    Et = Luna.todevice(MetalSpec, zeros(Float64, length(grid.t)))
    FT = Utils.plan_ft(Et, 1)
    Utils.plan_ift(FT)
    Eωb = Luna.todevice(MetalSpec, rand(ComplexF64, length(grid.ω)))
    ra = Boundaries.RateAbsorber(αt, Et, FT, (args...; kwargs...) -> nothing, 0.0)
    ra(Eωb, zmax/20, zmax/20, nothing)
    @test all(isfinite, Array(Eωb))
    la = Boundaries.LegacyAbsorber(grid, Et, FT, (args...; kwargs...) -> nothing)
    la(Eωb, zmax/20, zmax/20, nothing)
    @test all(isfinite, Array(Eωb))
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
    # ScaledOutput unscales on the way into the output now: no manual `* Eref` here
    for idx in axes(href["Eω"], 2)
        h = href["Eω"][:, idx]
        d = ComplexF64.(dref["Eω"][:, idx])
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

#= The rms spectral width of a save, and the ratio of the last save's to the first's, as
   evidence that a comparison actually exercises the Kerr nonlinearity rather than mostly
   linear dispersion (review round 1, finding 1: the original exit-criterion case, at the
   brief's own parameters, broadens by only ×1.000 -- not visibly nonlinear). A crude,
   grid-index-based second moment is enough for a ratio between two saves of the same run. =#
function rmswidth(Eω)
    p = abs2.(Eω)
    s = sum(p)
    idx = 1:length(p)
    μ = sum(idx .* p)/s
    sqrt(sum(@. (idx - μ)^2 * p)/s)
end
broadening(out) = rmswidth(out["Eω"][:, end])/rmswidth(out["Eω"][:, 1])

#= gpu/11's exit criteria: RateAbsorber and the default statistics on Metal, low-level
   interface, both constant and z-dependent operators, against genuine CPU Float32 *and*
   Float64 references (not Metal against itself -- review round 1, finding 1) and at the
   brief's own fibre length (review round 1, finding 9: `1e-2` was a tenth of it). Fixed
   steps (`fixed=true`, `min_dz == max_dz`) throughout, so the comparison is of arithmetic
   and not of a different step sequence; a genuinely nonlinear case (He at 5 bar, 300 µJ,
   which broadens the spectrum by about ×2, `metalcase`'s `energy`/`pres`) is included
   alongside the brief's own (weakly nonlinear) parameters, so the boundaries and the
   statistics are compared with the Kerr response actually doing something. Eω unscaled
   and on the host already (`Luna.run`'s `ScaledOutput`); the energy statistic is computed
   from a host copy on both paths and should agree closely even though it is itself
   derived from Eω. =#
@testset "boundaries and default statistics on Metal" begin
    for GT in (Grid.RealGrid, Grid.EnvGrid), (energy, pres) in ((100e-9, 1.0), (300e-6, 5.0))
        h32, _ = metalcase(GT, DeviceSpec(Array, Float32);
                           energy, pres, flength=0.1, boundary=:rate, stats=true, fixed=true)
        h64, _ = metalcase(GT, HostSpec();
                           energy, pres, flength=0.1, boundary=:rate, stats=true, fixed=true)
        dm, _ = metalcase(GT, MetalSpec;
                          energy, pres, flength=0.1, boundary=:rate, stats=true, fixed=true)
        if GT === Grid.RealGrid
            @test broadening(h64) > 1.5 || (energy, pres) == (100e-9, 1.0)
        end
        for idx in axes(h64["Eω"], 2)
            @test maximum(abs, dm["Eω"][:, idx] .- h32["Eω"][:, idx]) /
                  maximum(abs, h32["Eω"][:, idx]) < 1e-4
            @test maximum(abs, dm["Eω"][:, idx] .- h64["Eω"][:, idx]) /
                  maximum(abs, h64["Eω"][:, idx]) < 1e-4
        end
        @test isapprox(dm["stats"]["energy"], h32["stats"]["energy"]; rtol=1e-3)
        @test length(dm["stats"]["z"]) == length(h32["stats"]["z"])
    end

    #= gradient, the brief's own (weak) parameters, fixed steps. A pressure gradient
       redistributes the nonlinearity along z as the pressure falls, so it needs finer
       step resolution than a constant-pressure run of the same energy to stay stable: at
       20 fixed steps over the full 10 cm, both a strongly nonlinear gradient (300 µJ at
       5 bar -- the same energy/pressure the constant-pressure "visible Kerr" case above
       resolves cleanly) and the weaker one `metalgradientcase` defaults to previously
       used produced `NaN`/`Inf` here; that is a resolution problem with the *fixed* step
       count for that case, not a boundary or device defect -- the adaptive comparison
       just below, at the same strongly nonlinear parameters and with the controller free
       to refine the step size, is well behaved. =#
    hgrad32, _ = metalgradientcase(DeviceSpec(Array, Float32);
                                   energy=100e-9, pin=1.0, flength=0.1,
                                   boundary=:rate, stats=true)
    hgrad64, _ = metalgradientcase(HostSpec();
                                   energy=100e-9, pin=1.0, flength=0.1,
                                   boundary=:rate, stats=true)
    dgrad, _ = metalgradientcase(MetalSpec;
                                 energy=100e-9, pin=1.0, flength=0.1,
                                 boundary=:rate, stats=true)
    @test all(isfinite, hgrad32["Eω"]) && all(isfinite, dgrad["Eω"])
    for idx in axes(hgrad64["Eω"], 2)
        @test maximum(abs, dgrad["Eω"][:, idx] .- hgrad32["Eω"][:, idx]) /
              maximum(abs, hgrad32["Eω"][:, idx]) < 1e-4
        @test maximum(abs, dgrad["Eω"][:, idx] .- hgrad64["Eω"][:, idx]) /
              maximum(abs, hgrad64["Eω"][:, idx]) < 1e-4
    end
    @test isapprox(dgrad["stats"]["energy"], hgrad32["stats"]["energy"]; rtol=1e-3)

    #= A visibly nonlinear gradient (300 µJ, 5 bar), adaptive (the step-size controller
       genuinely active, `fixed=false`, well-resolved there -- see the comment above) with
       a documented looser tolerance: review round 1 measured 1.86e-4 for a strongly
       nonlinear adaptive gradient case (above the 1e-4 the other comparisons here use)
       because Float32 rounding perturbs the controller's accept/reject decisions, and the
       two runs can end up at slightly different step counts as well as slightly
       different arithmetic -- not a device defect, since every fixed-step comparison in
       this testset agrees to 1e-4. =#
    hgrad_a, _ = metalgradientcase(DeviceSpec(Array, Float32);
                                   energy=300e-6, pin=5.0, flength=0.1,
                                   boundary=:rate, stats=true, fixed=false)
    dgrad_a, _ = metalgradientcase(MetalSpec;
                                   energy=300e-6, pin=5.0, flength=0.1,
                                   boundary=:rate, stats=true, fixed=false)
    @test all(isfinite, hgrad_a["Eω"]) && all(isfinite, dgrad_a["Eω"])
    for idx in axes(hgrad_a["Eω"], 2)
        @test maximum(abs, dgrad_a["Eω"][:, idx] .- hgrad_a["Eω"][:, idx]) /
              maximum(abs, hgrad_a["Eω"][:, idx]) < 5e-4
    end
end

#= gpu/11's actual exit criterion: `prop_capillary` itself, unmodified apart from the new
   keywords, runs on the GPU with the boundaries and the default statistics -- constant
   and gradient pressure, `:auto` (what a plain `using Metal` sets) as well as an explicit
   `device=MetalSpec`.

   Review round 1, finding 1: `href` must be built with an *explicit* `device`, never
   left to resolve through the sentinel -- `test_metal.jl` runs with Metal loaded, so
   `Luna.settings["device"]` is `:auto` throughout the file except where a testset
   deliberately overrides it, and `precision=Float32` alone (the original code here) does
   not select the CPU: it resolves through the same `:auto` and lands on Metal, so the
   comparison was Metal against itself. Both a genuine CPU Float32 and a genuine CPU
   Float64 reference are used below. Finding 9: the fibre length is the brief's own
   (`0.1` m, not a tenth of it). At these parameters (the brief's own: 100 nJ, Kerr only)
   the spectrum barely broadens (`broadening` above is ~1.00), so this is a weakly
   nonlinear case by construction; `test_metal.jl`'s "boundaries and default statistics on
   Metal" carries the case with the Kerr effect clearly visible (×2 broadening) and fixed
   steps, which is where the tight 1e-4 tolerance is actually exercised by the
   nonlinearity. The constant-pressure comparisons here measure in the 1e-6-1e-5 range in
   practice (adaptive stepping, but the nonlinearity is too weak for that to matter at
   this energy); the gradient one needs a separate, looser, documented tolerance against
   the Float64 reference specifically -- see the comment at that assertion. =#
@testset "vector pointwise responses on Metal" begin
    n = 512
    ρ = PhysData.density(:He, 1.0)
    γ3 = PhysData.γ3_gas(:He)
    sc = Luna.UnitScaling(1024.0, PhysData.ε_0)
    hspec = DeviceSpec(Array, Float32)
    for (resp, T) in ((Nonlinear.Kerr_field(γ3), Float32),
                      (Nonlinear.Kerr_env(γ3), ComplexF32))
        Eh = T <: Complex ? randn(ComplexF32, n, 2) : randn(Float32, n, 2)
        Ph = zeros(T, n, 2)
        NonlinearRHS.Et_to_Pt!(Ph, Eh, (Nonlinear.rescale(resp, hspec, sc),), ρ;
                               scaling=sc)
        Ed = Luna.todevice(MetalSpec, Eh)
        Pd = Luna.alloc(MetalSpec, T, (n, 2))
        NonlinearRHS.Et_to_Pt!(Pd, Ed, (Nonlinear.rescale(resp, MetalSpec, sc),), ρ;
                               scaling=sc)
        @test maximum(abs, Array(Pd) .- Ph)/maximum(abs, Ph) < 1e-5
    end
end

#= The hackability fallback on real hardware: a user-written columnwise closure, which
   knows nothing about devices or units, run through `HostResponse` in a Metal
   propagation. Correct but slow by construction. =#
@testset "a user closure response through HostResponse on Metal" begin
    cw = usercubic(PhysData.ε_0*PhysData.γ3_gas(:He)/10)
    href32, htr32 = metalcase(Grid.RealGrid, DeviceSpec(Array, Float32);
                              extraresp=(cw,), fixed=true)
    href64, _ = metalcase(Grid.RealGrid, HostSpec(); extraresp=(cw,), fixed=true)
    mref, mtr = metalcase(Grid.RealGrid, MetalSpec; extraresp=(cw,), fixed=true)

    @test htr32.resp[2] isa Nonlinear.HostResponse
    @test mtr.resp[2] isa Nonlinear.HostResponse
    @test Nonlinear.kind(mtr.resp[2]) isa Nonlinear.Batched
    @test mtr.resp[2].resp === cw
    @test eltype(mref["Eω"]) === ComplexF32

    for (ref, tol) in ((href32, 1e-4), (href64, 1e-4))
        d = 0.0
        for idx in axes(ref["Eω"], 2)
            h = ComplexF64.(ref["Eω"][:, idx])
            m = ComplexF64.(mref["Eω"][:, idx])
            d = max(d, maximum(abs, m .- h)/maximum(abs, h))
        end
        @test d < tol
    end

    # The closure really contributes
    plain, _ = metalcase(Grid.RealGrid, HostSpec(); fixed=true)
    @test maximum(abs, href64["Eω"][:, end] .- plain["Eω"][:, end])/
          maximum(abs, plain["Eω"][:, end]) > 1e-6

    # It says so, once, at setup
    @test_logs (:info,) match_mode=:any Nonlinear.rescale(
        cw, MetalSpec, UNIT_SCALING, Luna.alloc(MetalSpec, Float32, (16,)))
end

@testset "plasma on Metal" begin
    t, E = metal_plasmafield()
    Ev = hcat(E, 0.6 .* circshift(E, 7))
    ionpot = PhysData.ionisation_potential(:Ar)
    ρ = PhysData.density(:Ar, 1.0)
    Eref = exp2(round(Int, log2(maximum(abs, E))))
    sc = Luna.UnitScaling(Eref, PhysData.ε_0)
    cpu32 = DeviceSpec(Array, Float32)

    for (nm, ir) in (("ADK", metal_adkrate()), ("table", metal_tablerate())),
        Eh in (E, Ev)
        p = Nonlinear.PlasmaCumtrapz(t, Eh, ir, ionpot)

        # host Float64, physical units: the reference
        P64 = zeros(size(Eh)); p(P64, Eh, ρ)

        # host Float32, scaled
        Es = Float32.(Eh ./ Eref)
        ph = Nonlinear.rescale(p, cpu32, sc, Es)
        Ph = zeros(Float32, size(Eh))
        Nonlinear.batched!(ph, Ph, Es, ρ, sc)

        # Metal Float32, scaled
        Ed = Luna.todevice(MetalSpec, Eh ./ Eref)
        pd = Nonlinear.rescale(p, MetalSpec, sc, Ed)
        @test pd.ratedev isa Ionisation.AbstractIonRate
        @test pd.J isa MtlArray{Float32}
        @test Luna.all_resident(MetalSpec, Nonlinear.resident_arrays(pd)...)
        Pd = Luna.alloc(MetalSpec, Float32, size(Eh))
        Nonlinear.batched!(pd, Pd, Ed, ρ, sc)

        Pdh = Array(Pd)
        @test all(isfinite, Pdh)
        @test !all(iszero, Pdh) # the field really ionises: not two zeros agreeing
        # Metal against the CPU at the same precision: the device path, not Float32
        @test maximum(abs, Pdh .- Ph)/maximum(abs, Ph) < 1e-4
        # ... and the scaled Float32 answer against the physical Float64 one
        phys = Pdh .* Float32(PhysData.ε_0*Eref)
        @test maximum(abs, phys .- P64)/maximum(abs, P64) < 1e-3
    end

    #= Dynamic range, case 1: helium at 0.3 bar and a field far below the ionisation
       threshold. The rate underflows to zero in Float32 (and is ~1e-300 in Float64), and
       what matters is that the response produces zeros rather than NaNs: an underflow
       inside `exp` followed by a division by a zero field is exactly where the `ifelse`
       loss term could poison the whole block. =#
    tl, El = metal_plasmafield(512, 1e9)
    ρhe = PhysData.density(:He, 0.3)
    ipe = PhysData.ionisation_potential(:He)
    pl = Nonlinear.PlasmaCumtrapz(tl, El, Ionisation.IonRateADK(:He), ipe)
    ErefL = exp2(round(Int, log2(maximum(abs, El))))
    scl = Luna.UnitScaling(ErefL, PhysData.ε_0)
    Edl = Luna.todevice(MetalSpec, El ./ ErefL)
    pdl = Nonlinear.rescale(pl, MetalSpec, scl, Edl)
    Pdl = Luna.alloc(MetalSpec, Float32, size(El))
    Nonlinear.batched!(pdl, Pdl, Edl, ρhe, scl)
    @test all(isfinite, Array(Pdl))
    @test all(iszero, Array(Pdl))

    #= Dynamic range, case 2: argon at a field near barrier suppression, where the
       ionisation fraction saturates at 1. `1 - exp(-x)` with a large `x`, a rate of
       ~1e16 s^-1 accumulated over the whole time axis, and a loss term divided by the
       field: the largest intermediates the response produces anywhere in Luna's
       parameter range. =#
    Ebs = Ionisation.barrier_suppression(ionpot, 1.0)
    th, Eh2 = metal_plasmafield(512, Ebs)
    ph2 = Nonlinear.PlasmaCumtrapz(th, Eh2, metal_adkrate(), ionpot)
    P64h = zeros(size(Eh2)); ph2(P64h, Eh2, ρ)
    ErefH = exp2(round(Int, log2(maximum(abs, Eh2))))
    sch = Luna.UnitScaling(ErefH, PhysData.ε_0)
    Es = Float32.(Eh2 ./ ErefH)
    phh = Nonlinear.rescale(ph2, cpu32, sch, Es)
    Phh = zeros(Float32, size(Eh2)); Nonlinear.batched!(phh, Phh, Es, ρ, sch)
    Edh = Luna.todevice(MetalSpec, Eh2 ./ ErefH)
    pdh = Nonlinear.rescale(ph2, MetalSpec, sch, Edh)
    Pdh2 = Luna.alloc(MetalSpec, Float32, size(Eh2))
    Nonlinear.batched!(pdh, Pdh2, Edh, ρ, sch)
    @test all(isfinite, Array(Pdh2))
    @test maximum(abs, Array(Pdh2) .- Phh)/maximum(abs, Phh) < 1e-4
    physh = Array(Pdh2) .* Float32(PhysData.ε_0*ErefH)
    @test maximum(abs, physh .- P64h)/maximum(abs, P64h) < 1e-2
    #= The case really is the extreme one: a rate above 1e14 1/s and a tenth of the gas
       ionised by the end of the pulse. =#
    @test maximum(Array(pdh.rate)) > 1e14
    @test maximum(Array(pdh.fraction)) > 0.1f0
end

#= The Raman polarisation on Metal: nitrogen and hydrogen -- the gas with the largest
   Raman gain Luna is used with -- in the three forms the response takes. Compared with
   the CPU at the same precision, which is the device path rather than single precision,
   and with the physical Float64 answer. =#
@testset "Raman and the no-THG Kerr on Metal" begin
    t, E = metal_plasmafield(512, 1e10)
    Eenv = complex.(@. 1e10*exp(-t^2/(2*(10e-15/1.66)^2)))
    Eref = exp2(round(Int, log2(maximum(abs, E))))
    sc = Luna.UnitScaling(Eref, PhysData.ε_0)
    cpu32 = DeviceSpec(Array, Float32)
    for gas in (:N2, :H2)
        ρ = PhysData.density(gas, 1.0)
        cases = (("field, THG",
                  () -> Nonlinear.RamanPolarField(t, Raman.raman_response(t, gas)), E),
                 ("field, no THG",
                  () -> Nonlinear.RamanPolarField(t, Raman.raman_response(t, gas);
                                                  thg=false), E),
                 ("envelope",
                  () -> Nonlinear.RamanPolarEnv(t, Raman.raman_response(t, gas)), Eenv))
        for (nm, make, Eh) in cases
            # host Float64, physical units: the reference
            R64 = Nonlinear.rescale(make(), HostSpec(), UNIT_SCALING, Eh)
            P64 = zeros(eltype(Eh), size(Eh)); R64(P64, Eh, ρ)

            # host Float32, scaled
            Es = Luna.todevice(cpu32, Eh ./ Eref)
            Rh = Nonlinear.rescale(make(), cpu32, sc, Es)
            Ph = zeros(eltype(Es), size(Eh))
            Nonlinear.batched!(Rh, Ph, Es, ρ, sc)

            # Metal Float32, scaled
            Ed = Luna.todevice(MetalSpec, Eh ./ Eref)
            Rd = Nonlinear.rescale(make(), MetalSpec, sc, Ed)
            @test Rd.hω isa MtlArray{ComplexF32, 1}
            @test Rd.E2 isa MtlArray
            @test Luna.all_resident(MetalSpec, Nonlinear.resident_arrays(Rd)...)
            Pd = Luna.alloc(MetalSpec, eltype(Es), size(Eh))
            Nonlinear.batched!(Rd, Pd, Ed, ρ, sc)

            Pdh = Array(Pd)
            @test all(isfinite, Pdh)
            @test !all(iszero, Pdh) # not two zeros agreeing
            @test maximum(abs, Pdh .- Ph)/maximum(abs, Ph) < 1e-4
            phys = Pdh .* Float32(PhysData.ε_0*Eref)
            @test maximum(abs, phys .- P64)/maximum(abs, P64) < 1e-3
        end
    end

    #= The no-THG Kerr response, whose analytic signal is the same transform, on the same
       field. The device comparison in the propagation test below cannot separate it from
       the plain Kerr response -- the third-harmonic term is smaller than the Float32
       tolerance there -- but here it can: the two responses differ by 30 % of the peak. =#
    γ3 = PhysData.γ3_gas(:He)
    ρ = PhysData.density(:He, 1.0)
    kn = Nonlinear.Kerr_field_nothg(γ3, length(E))
    P64 = zeros(size(E)); kn(P64, E, ρ)
    Es = Luna.todevice(cpu32, E ./ Eref)
    knh = Nonlinear.rescale(kn, cpu32, sc, Es)
    Ph = zeros(Float32, size(E)); Nonlinear.batched!(knh, Ph, Es, ρ, sc)
    Ed = Luna.todevice(MetalSpec, E ./ Eref)
    knd = Nonlinear.rescale(kn, MetalSpec, sc, Ed)
    Pd = Luna.alloc(MetalSpec, Float32, size(E))
    Nonlinear.batched!(knd, Pd, Ed, ρ, sc)
    Pdh = Array(Pd)
    @test all(isfinite, Pdh)
    @test !all(iszero, Pdh)
    @test maximum(abs, Pdh .- Ph)/maximum(abs, Ph) < 1e-4
    @test maximum(abs, Pdh .* Float32(PhysData.ε_0*Eref) .- P64)/maximum(abs, P64) < 1e-3
    # ... and it is not the plain Kerr response
    Pk = zeros(size(E)); Nonlinear.Kerr_field(γ3)(Pk, E, ρ)
    @test maximum(abs, Pk .- P64)/maximum(abs, P64) > 0.1
end

@testset "prop_capillary on Metal" begin
    capkw = (; λ0=800e-9, energy=100e-9, τfwhm=10e-15, λlims=(300e-9, 2000e-9),
             trange=400e-15, saveN=5, plasma=false, raman=false, shotnoise=false)

    h32 = Luna.prop_capillary(125e-6, 0.1, :He, 1.0; capkw..., device=DeviceSpec(Array, Float32))
    h64 = Luna.prop_capillary(125e-6, 0.1, :He, 1.0; capkw..., device=HostSpec())
    dm = Luna.prop_capillary(125e-6, 0.1, :He, 1.0; capkw..., device=MetalSpec)
    for idx in axes(h64["Eω"], 2)
        @test maximum(abs, dm["Eω"][:, idx] .- h32["Eω"][:, idx]) /
              maximum(abs, h32["Eω"][:, idx]) < 1e-4
        @test maximum(abs, dm["Eω"][:, idx] .- h64["Eω"][:, idx]) /
              maximum(abs, h64["Eω"][:, idx]) < 1e-4
    end
    @test isapprox(dm["stats"]["energy"], h32["stats"]["energy"]; rtol=1e-3)

    #= a pressure gradient, same parameters. Metal vs CPU Float32 (same precision, the
       most direct test of the device path) stays at the 1e-4 tolerance; Metal vs CPU
       Float64 gets a documented, looser one. A pressure gradient makes the adaptive
       controller's decisions more sensitive to Float32 rounding than a constant-pressure
       run, `prop_capillary` has no fixed-step option, and Float32 vs Float64 alone (no
       device involved at all) already sits close to this size for a gradient -- review
       round 1's own measurement of the same kind of case was 2.74e-4 for CPU Float32 vs
       CPU Float64. Measured here: Metal vs CPU Float32 agrees to ~1e-5; Metal vs CPU
       Float64 reaches ~2e-4, i.e. it is the Float32/Float64 gap, not the device, that
       sets the size of this one. =#
    hgrad32 = Luna.prop_capillary(125e-6, 0.1, :He, (1.0, 0.0);
                                  capkw..., device=DeviceSpec(Array, Float32))
    hgrad64 = Luna.prop_capillary(125e-6, 0.1, :He, (1.0, 0.0); capkw..., device=HostSpec())
    dgrad = Luna.prop_capillary(125e-6, 0.1, :He, (1.0, 0.0); capkw..., device=MetalSpec)
    for idx in axes(hgrad64["Eω"], 2)
        @test maximum(abs, dgrad["Eω"][:, idx] .- hgrad32["Eω"][:, idx]) /
              maximum(abs, hgrad32["Eω"][:, idx]) < 1e-4
        @test maximum(abs, dgrad["Eω"][:, idx] .- hgrad64["Eω"][:, idx]) /
              maximum(abs, hgrad64["Eω"][:, idx]) < 3e-4
    end
    @test isapprox(dgrad["stats"]["energy"], hgrad32["stats"]["energy"]; rtol=1e-3)

    #= gpu/13's exit condition: Kerr *and* plasma, the default physics of a
       field-resolved `prop_capillary` call in a non-Raman gas, on Metal end to end.
       Argon at an intensity which ionises, and fixed steps, so the only difference
       between the runs is the arithmetic. =#
    p32, _ = metalcase(Grid.RealGrid, DeviceSpec(Array, Float32);
                       gas=:Ar, energy=150e-6, plasma=true, fixed=true)
    p64, _ = metalcase(Grid.RealGrid, HostSpec(); gas=:Ar, energy=150e-6, plasma=true,
                       fixed=true)
    pdm, ptr = metalcase(Grid.RealGrid, MetalSpec; gas=:Ar, energy=150e-6, plasma=true,
                         fixed=true)
    @test ptr.resp[2] isa Nonlinear.PlasmaCumtrapz
    @test ptr.resp[2].J isa MtlArray{Float32}
    @test ptr.resp[2].ratedev.spline.x isa MtlArray{Float32}
    for idx in axes(p64["Eω"], 2)
        @test maximum(abs, pdm["Eω"][:, idx] .- p32["Eω"][:, idx]) /
              maximum(abs, p32["Eω"][:, idx]) < 1e-4
        @test maximum(abs, pdm["Eω"][:, idx] .- p64["Eω"][:, idx]) /
              maximum(abs, p64["Eω"][:, idx]) < 1e-3
    end
    #= The plasma really contributes: a Kerr-only run of the same case differs by far
       more than the tolerances above. =#
    nop64, _ = metalcase(Grid.RealGrid, HostSpec(); gas=:Ar, energy=150e-6, fixed=true)
    @test maximum(abs, p64["Eω"][:, end] .- nop64["Eω"][:, end]) /
          maximum(abs, nop64["Eω"][:, end]) > 1e-2

    #= gpu/14's exit condition: Kerr *and* the Raman polarisation, the default physics of
       a field-resolved `prop_capillary` call in a molecular gas, on Metal end to end.
       Nitrogen and hydrogen, hydrogen being the largest Raman gain Luna is used with.
       Fixed steps, so the only difference between the runs is the arithmetic. =#
    for gas in (:N2, :H2)
        rkw = (; gas, energy=50e-6, raman=true, fixed=true)
        r32, _ = metalcase(Grid.RealGrid, DeviceSpec(Array, Float32); rkw...)
        r64, _ = metalcase(Grid.RealGrid, HostSpec(); rkw...)
        rdm, rtr = metalcase(Grid.RealGrid, MetalSpec; rkw...)
        @test rtr.resp[2] isa Nonlinear.RamanPolarField
        @test rtr.resp[2].hω isa MtlArray{ComplexF32, 1}
        @test rtr.resp[2].E2 isa MtlArray{Float32}
        for idx in axes(r64["Eω"], 2)
            @test maximum(abs, rdm["Eω"][:, idx] .- r32["Eω"][:, idx]) /
                  maximum(abs, r32["Eω"][:, idx]) < 1e-4
            @test maximum(abs, rdm["Eω"][:, idx] .- r64["Eω"][:, idx]) /
                  maximum(abs, r64["Eω"][:, idx]) < 1e-3
        end
        #= The Raman term really contributes: a Kerr-only run of the same case differs by
           far more than the tolerances above. =#
        nor64, _ = metalcase(Grid.RealGrid, HostSpec(); gas, energy=50e-6, fixed=true)
        @test maximum(abs, r64["Eω"][:, end] .- nor64["Eω"][:, end]) /
              maximum(abs, nor64["Eω"][:, end]) > 1e-2
    end

    #= The no-THG Kerr response, which `prop_capillary(...; thg=false)` selects on a
       `RealGrid`: batched, because removing the third harmonic needs the analytic signal
       of the whole column. =#
    tkw = (; energy=100e-6, nothg=true, fixed=true)
    t32, _ = metalcase(Grid.RealGrid, DeviceSpec(Array, Float32); tkw...)
    t64, _ = metalcase(Grid.RealGrid, HostSpec(); tkw...)
    tdm, ttr = metalcase(Grid.RealGrid, MetalSpec; tkw...)
    @test ttr.resp[1] isa Nonlinear.KerrFieldNoTHG
    @test ttr.resp[1].an.c1 isa MtlArray{ComplexF32}
    @test ttr.resp[1].an.mask isa MtlArray{Float32, 1}
    for idx in axes(t64["Eω"], 2)
        @test maximum(abs, tdm["Eω"][:, idx] .- t32["Eω"][:, idx]) /
              maximum(abs, t32["Eω"][:, idx]) < 1e-4
        @test maximum(abs, tdm["Eω"][:, idx] .- t64["Eω"][:, idx]) /
              maximum(abs, t64["Eω"][:, idx]) < 1e-3
    end
    #= Removing THG really changes the answer, so this is not a comparison of two
       ordinary Kerr runs. The threshold is small because the third-harmonic term itself
       is: helium at 100 µJ over a centimetre gives 5.4e-5, and no parameters in this
       file's range make it larger. That is a comparison of two *Float64* runs, where
       5e-5 is enormous, but it is below the 1e-4 the device comparison above allows --
       so what rules out the wrong response on the device is the structural check on
       `ttr.resp[1]`, and the bit-exact host test of `AnalyticSignal` against
       `Maths.plan_hilbert` in `test_device.jl`, rather than this number. =#
    k64, _ = metalcase(Grid.RealGrid, HostSpec(); energy=100e-6, fixed=true)
    @test maximum(abs, t64["Eω"][:, end] .- k64["Eω"][:, end]) /
          maximum(abs, k64["Eω"][:, end]) > 1e-5

    #= ... and through the simple interface, with a real cached PPT rate, which is what
       a user gets from `prop_capillary(...; plasma=true)`. Argon at 0.1 bar and 300 µJ,
       which ionises about 1 % of the gas: review round 1, finding 8 -- the first version
       of this used the regression matrix's helium parameters, where the rate is exactly
       zero, so the spline kernel was only ever exercised at zero. Built through
       `prop_capillary_args` so that the steps can be fixed, which `prop_capillary`
       itself has no keyword for. =#
    plkw = (; λ0=800e-9, energy=300e-6, τfwhm=10e-15, λlims=(300e-9, 2000e-9),
            trange=400e-15, saveN=3, plasma=true, raman=false, shotnoise=false)
    function ppt_prop(spec)
        Eω, grid, linop, tr, FT, o = Luna.Interface.prop_capillary_args(
            125e-6, 1e-2, :Ar, 0.1; plkw..., device=spec)
        h = 1e-2/20
        Luna.run(Eω, grid, linop, tr, FT, o;
                 zmax=1e-2, init_dz=h, min_dz=h, max_dz=h, status_period=1e6)
        o, tr
    end
    ph32, ptr32 = ppt_prop(DeviceSpec(Array, Float32))
    pdm2, ptrm = ppt_prop(MetalSpec)
    @test ptrm.resp[2].ratedev isa Ionisation.IonRatePPTAccel
    @test ptrm.resp[2].ratedev.spline.x isa MtlArray{Float32}
    @test eltype(pdm2["Eω"]) === ComplexF32
    for idx in axes(ph32["Eω"], 2)
        @test maximum(abs, pdm2["Eω"][:, idx] .- ph32["Eω"][:, idx]) /
              maximum(abs, ph32["Eω"][:, idx]) < 1e-4
    end
    #= The cached rate really fires: the electron density is a fraction of a percent of
       the gas, not zero. `Stats` computes it on the host from the saved field, so this
       is the device run's own output. =#
    @test maximum(pdm2["stats"]["electrondensity"])/PhysData.density(:Ar, 0.1) > 1e-3

    # `:auto` (what loading Metal sets) resolves to the same device as the explicit spec
    old = get(Luna.settings, "device", nothing)
    try
        Luna.set_device(:auto)
        @test Luna.device() === MetalSpec
        aref = Luna.prop_capillary(125e-6, 0.1, :He, 1.0; capkw...)
        # Metal defaults to Float32; the saved field is unscaled but stays that precision
        @test eltype(aref["Eω"]) === ComplexF32
        @test aref["Eω"] == dm["Eω"] # the same device, so the same answer as above
    finally
        isnothing(old) ? delete!(Luna.settings, "device") :
                         (Luna.settings["device"] = old)
    end
end

#= `Luna.set_device(:cpu)` opts out, whatever `settings["device"]` is otherwise: this is
   exit criterion 3. `prop_gnlse` and multimode/radial `prop_capillary` are not
   device-capable (`_cpu_only!`, `Interface.jl`) and keep giving the CPU, Float64 answer
   under `:auto` too -- refusing only when the caller explicitly asks for something else. =#
@testset "Luna.set_device(:cpu) opts out" begin
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
        @test eltype(ocap["Eω"]) === ComplexF64
        @test eltype(ognlse["Eω"]) === ComplexF64

        Luna.set_device(:auto)
        @test Luna.device() === MetalSpec # the GPU really is selected globally
        # prop_gnlse is not device-capable and always refuses anything but the CPU
        @test_throws ErrorException Luna.prop_gnlse(0.1, 1e-3, [0.0, 0.0, -1e-26];
                                                     gnlsekw..., device=MetalSpec)
        # ... but is untouched when device is left at its (CPU-resolving) default
        dgnlse = Luna.prop_gnlse(0.1, 1e-3, [0.0, 0.0, -1e-26]; gnlsekw...)
        @test dgnlse["Eω"] == ognlse["Eω"]

        # explicitly asking for the CPU under :auto still gives the CPU
        dcap = Luna.prop_capillary(capargs...; capkw..., device=:cpu)
        @test dcap["Eω"] == ocap["Eω"]

        #= Multimode propagation is not device-capable and must stay on the CPU by
           default under :auto too -- it must not turn a working run into an error just
           because a GPU package happens to be loaded (the bug an earlier version of
           this branch had: `device`'s default resolved through `:auto` even for paths
           that can never honour it). =#
        om = Luna.prop_capillary(capargs...; capkw..., modes=4)
        @test size(om["Eω"], 2) == 4
        # ... but an explicit device request for multimode still errors
        @test_throws ErrorException Luna.prop_capillary(capargs...; capkw..., modes=4,
                                                         device=MetalSpec)
    finally
        isnothing(old) ? delete!(Luna.settings, "device") :
                         (Luna.settings["device"] = old)
    end
end

@testset "Metal refuses what it cannot run" begin
    #= A columnwise response is wrapped in a `HostResponse` rather than refused
       (gpu/12); what is refused is one which claims a device kernel for arrays nothing
       has converted. =#
    wrapped = Nonlinear.rescale((out, E, ρ) -> nothing, MetalSpec, UNIT_SCALING,
                                Luna.alloc(MetalSpec, Float32, (16,)))
    @test wrapped isa Nonlinear.HostResponse
    # Its buffers exist at construction, the device one on the device
    @test wrapped.Eh isa Vector{Float64}
    @test wrapped.stage isa Vector{Float32}
    @test wrapped.Pd isa MtlArray{Float32, 1}
    # An unwrapped columnwise response applied to a device array is refused, not run
    Pd = Luna.alloc(MetalSpec, Float32, (16,))
    Ed = Luna.todevice(MetalSpec, randn(16))
    @test_throws ErrorException NonlinearRHS.Et_to_Pt!(
        Pd, Ed, ((out, E, ρ) -> nothing,), 1.0)

    # multimode propagation is not device-capable through the simple interface either
    @test_throws ErrorException Luna.prop_capillary(
        125e-6, 1e-3, :He, 1.0; λ0=800e-9, energy=1e-9, τfwhm=10e-15,
        λlims=(300e-9, 2e-6), trange=400e-15, saveN=3, plasma=false, shotnoise=false,
        modes=4, device=MetalSpec)
end

end # have_metal
