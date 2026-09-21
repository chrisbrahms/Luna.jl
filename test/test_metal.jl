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

#= A user-written *two-component* columnwise response: a toy χ⁽²⁾ crystal written the way a
   user would write one, with scalar indexing, in physical SI units and Float64. Nothing
   about it can run on a Metal array. =#
userchi2(d) = let d = d
    (out, E, ρ) -> begin
        for i in axes(E, 1)
            @inbounds out[i, 1] += d*2*E[i, 1]*E[i, 2]
            @inbounds out[i, 2] += d*(E[i, 1]^2 - E[i, 2]^2)
        end
        out
    end
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
                   extraresp=(), plasma=false, raman=false, nothg=false,
                   stats_device=:auto)
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
    #= The state itself: `Stats.default` builds its buffers and plans its transform for
       the array type and precision `Eω` has, so the default statistics run on the
       device with the state left where it is (gpu/24). =#
    statsfun = stats ?
        Stats.default(grid, Eω, m, linop, transform; gas, stats_device) : Output.nostats
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

#= How many times a piece of host code was called, so that "no host work per stage" can
   be checked by counting rather than by inspecting types (the same wrapper
   `test_device.jl` uses). `densityfun` is called exactly once per right-hand side, so it
   counts the stages. =#
mutable struct CountCalls{F}
    f::F
    n::Int
end
CountCalls(f) = CountCalls(f, 0)
(c::CountCalls)(args...; kwargs...) = (c.n += 1; c.f(args...; kwargs...))

#= A pressure gradient: a z-dependent operator closure and `constβ=false`, which is the
   only case exercising `NormModeAvg`'s `HostMirror` branch and `RK45.make_prop!`'s
   host-buffer branch -- the two pieces which upload from the host on every stage unless
   `tabulate_linop=true` replaces them with tables (gpu/23). Fixed steps, so the runs
   differ only in arithmetic. The third element of the return value counts the host calls
   the run made. =#
function metalgradientcase(spec; gas=:Ar, pin=1.0, pout=0.0, energy=1e-6, flength=1e-2,
                           λ0=800e-9, boundary=:none, stats=false, fixed=true,
                           tabulate_linop=false, stats_device=:auto,
                           linop_tol=LinearOps.DEFAULT_LINOP_TOL, nsteps=20, rtol=1e-8)
    grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15)
    coren, densityfun0 = Capillary.gradient(gas, flength, pin, pout)
    m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
    aeff = CountCalls(z -> Modes.Aeff(m, z=z))
    densityfun = CountCalls(densityfun0)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop0, βfun0! = LinearOps.make_linop(grid, m, λ0)
    linop = CountCalls(linop0)
    βfun! = CountCalls(βfun0!)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, densityfun, resp, inputs, βfun!, aeff;
                                   device=spec)
    #= The state itself: `Stats.default` builds its buffers and plans its transform for
       the array type and precision `Eω` has, so the default statistics run on the
       device with the state left where it is (gpu/24). =#
    statsfun = stats ?
        Stats.default(grid, Eω, m, linop, transform; gas, stats_device) : Output.nostats
    out = Output.MemoryOutput(0, flength, 5, statsfun)
    #= `max_dz` does not depend on the step count, so that two runs with different step
       counts tabulate over the same interval and build the same tables. =#
    maxdz = flength/20
    dz = flength/nsteps
    if fixed
        Luna.run(Eω, grid, linop, transform, FT, out;
                 zmax=flength, boundary, init_dz=dz, min_dz=dz, max_dz=maxdz, rtol,
                 tabulate_linop, linop_tol)
    else
        Luna.run(Eω, grid, linop, transform, FT, out;
                 zmax=flength, boundary, init_dz=dz, rtol=1e-8,
                 tabulate_linop, linop_tol)
    end
    out, transform, (; linop=linop.n, β=βfun!.n, aeff=aeff.n, rhs=densityfun.n)
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

#= The χ⁽²⁾ responses on Metal. The free-space transforms which use them are host-only
   until Group E, so this is the block path: one fused pair of component broadcasts over an
   `(nt, 2, ncols)` device array. This is the only place a stray Float64 in the crystal
   matrices or in the carrier phase would show up -- Metal's kernel compiler rejects any
   `double` which survives optimisation, and `MtlArray{Float64}` does not exist. =#
@testset "χ⁽²⁾ responses on Metal" begin
    nt, ncols = 256, 3
    θ, ϕ = deg2rad(29.2), deg2rad(30)
    sc = Luna.UnitScaling(1024.0, PhysData.ε_0)
    hspec = DeviceSpec(Array, Float32)
    to = collect(range(0, 1e-13, length=nt))
    for (c, T) in ((Nonlinear.Chi2Field(θ, ϕ, PhysData.χ2(:BBO)), Float32),
                   (Nonlinear.Chi2Env(θ, ϕ, PhysData.χ2(:BBO), PhysData.wlfreq(800e-9),
                                      to), ComplexF32))
        Eh = randn(T, nt, 2, ncols)
        Ph = zeros(T, nt, 2, ncols)
        NonlinearRHS.Et_to_Pt!(Ph, Eh, (Nonlinear.rescale(c, hspec, sc),), 1.0; scaling=sc)
        Ed = Luna.todevice(MetalSpec, Eh)
        rd = Nonlinear.rescale(c, MetalSpec, sc)
        Pd = Luna.alloc(MetalSpec, T, (nt, 2, ncols))
        NonlinearRHS.Et_to_Pt!(Pd, Ed, (rd,), 1.0; scaling=sc)
        # the crystal matrices and the carrier phase reach the kernel in Float32
        @test eltype(rd.χ2_toLab) === Float32
        @test eltype(rd.toCrystal) === Float32
        if c isa Nonlinear.Chi2Env
            @test eltype(rd.C) === ComplexF32
        end
        @test maximum(abs, Ph) > 0
        @test maximum(abs, Array(Pd) .- Ph)/maximum(abs, Ph) < 1e-5
    end
end

#= The χ⁽²⁾ case of the hackability fallback on real hardware: a two-component response a
   user wrote as a Float64 closure with scalar indexing, applied to a Metal block through
   `HostResponse` alongside a χ⁽²⁾ response which has a kernel. =#
@testset "a user χ⁽²⁾ closure through HostResponse on Metal" begin
    nt, ncols = 256, 2
    sc = Luna.UnitScaling(1024.0, PhysData.ε_0)
    hspec = DeviceSpec(Array, Float32)
    #= `d` is of the order of ε₀χ⁽²⁾ for BBO, so that the closure and the response with a
       kernel contribute comparably. =#
    cw = userchi2(1e-23)
    c = Nonlinear.Chi2Field(deg2rad(29.2), deg2rad(30), PhysData.χ2(:BBO))
    Eh = randn(Float32, nt, 2, ncols)
    Ph = zeros(Float32, nt, 2, ncols)
    NonlinearRHS.Et_to_Pt!(Ph, Eh,
                           (Nonlinear.rescale(c, hspec, sc),
                            Nonlinear.rescale(cw, hspec, sc, Eh)), 1.0; scaling=sc)

    Ed = Luna.todevice(MetalSpec, Eh)
    hr = @test_logs (:info,) match_mode=:any Nonlinear.rescale(cw, MetalSpec, sc, Ed)
    @test hr isa Nonlinear.HostResponse
    @test Nonlinear.kind(hr) isa Nonlinear.Batched
    Pd = Luna.alloc(MetalSpec, Float32, (nt, 2, ncols))
    NonlinearRHS.Et_to_Pt!(Pd, Ed, (Nonlinear.rescale(c, MetalSpec, sc), hr), 1.0;
                           scaling=sc)
    @test maximum(abs, Array(Pd) .- Ph)/maximum(abs, Ph) < 1e-5

    # the closure really contributes
    Pc = zeros(Float32, nt, 2, ncols)
    NonlinearRHS.Et_to_Pt!(Pc, Eh, (Nonlinear.rescale(c, hspec, sc),), 1.0; scaling=sc)
    @test maximum(abs, Pc .- Ph)/maximum(abs, Ph) > 1e-2
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

#= gpu/23's exit criterion: a pressure gradient on Metal with `tabulate_linop=true` and
   no host work inside the propagation. Three things are checked -- that the tables reach
   the kernels as Float32/ComplexF32 MtlArrays (the only place a stray Float64 in them
   would show, since Metal refuses one), that the tabulated Metal run agrees with a
   tabulated CPU Float32 run, and that nothing on the host is evaluated per stage, which
   is counted rather than inferred. =#
@testset "tabulated operator on Metal" begin
    flength = 1e-2
    h32, htr, _ = metalgradientcase(DeviceSpec(Array, Float32); flength,
                                    tabulate_linop=true)
    dm, dtr, dcount = metalgradientcase(MetalSpec; flength, tabulate_linop=true)

    #= `Luna.run` tabulates into a transform of its own and does not modify the caller's,
       so the tables are inspected by building the same thing here. =#
    dtab = NonlinearRHS.tabulate(dtr, 0.0, 1.05flength, 1e-6, dtr.Eωo)
    @test dtab.norm!.β isa LinearOps.TabulatedVector
    @test dtab.norm!.β.f isa MtlArray{Float32, 2}
    @test dtab.norm!.β.buf isa MtlArray{Float32, 1}
    @test dtab.norm!.β(0.3flength) isa MtlArray{Float32, 1}
    @test dtab.aeff isa LinearOps.TabulatedScalar

    # and the operator's own tables, which the propagator broadcasts over
    linop0, _ = LinearOps.make_linop(Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15),
                                     Capillary.MarcatiliMode(
                                         75e-6, first(Capillary.gradient(:Ar, flength, 1.0, 0.0)),
                                         loss=false),
                                     800e-9)
    proto = Luna.alloc(MetalSpec, ComplexF32, (length(dtr.grid.ω),))
    tab = LinearOps.TabulatedLinop(linop0, proto, 0.0, 1.05flength; tol=1e-6, quiet=true)
    @test tab.Φ isa MtlArray{ComplexF32, 2}
    @test tab.dΦ isa MtlArray{ComplexF32, 2}
    @test tab.secant isa MtlArray{ComplexF32, 1}
    # the readback and the propagator are broadcasts over those, with no Float64 anywhere
    out = similar(proto)
    LinearOps.phase!(out, tab, 0.4flength)
    @test all(isfinite, Array(out))
    y = Luna.todevice(MetalSpec, ones(ComplexF64, length(dtr.grid.ω)))
    RK45.make_prop!(tab, y)(y, 0.4flength, 0.41flength)
    @test all(isfinite, Array(y))
    @test maximum(abs, Array(y)) > 0

    for idx in axes(h32["Eω"], 2)
        @test maximum(abs, dm["Eω"][:, idx] .- h32["Eω"][:, idx]) /
              maximum(abs, h32["Eω"][:, idx]) < 1e-4
    end

    #= No per-stage host work: two runs with different step counts but the same `max_dz`,
       so the tables span the same interval and cost the same number of evaluations.
       Anything the stepper evaluated on the host would scale with the stage count, which
       the right-hand side count shows really did change. =#
    _, _, fine = metalgradientcase(MetalSpec; flength, tabulate_linop=true,
                                   nsteps=80, rtol=1e-13)
    @test fine.rhs > 1.5*dcount.rhs
    @test fine.linop == dcount.linop
    @test fine.β == dcount.β
    @test fine.aeff == dcount.aeff
    # ... where the untabulated path evaluates and uploads all three at every stage
    _, _, ecoarse = metalgradientcase(MetalSpec; flength, nsteps=20)
    _, _, efine = metalgradientcase(MetalSpec; flength, nsteps=80, rtol=1e-13)
    @test efine.linop > 1.5*ecoarse.linop
    @test efine.β > 1.5*ecoarse.β
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

    #= gpu/23: the same gradient with `tabulate_linop=true`, which is the keyword a user
       passes to run a pressure-graded capillary on the GPU with nothing left on the host
       inside the step. Metal against CPU Float32 at the same tolerance as the untabulated
       comparison above: the two runs use the same discretisation as each other, so this
       is the device path and nothing else. The tabulated answer is *not* compared with
       the untabulated one -- it is a different discretisation of the linear step, and
       much better resolved; that difference is measured in `test_device.jl`. =#
    tgrad32 = Luna.prop_capillary(125e-6, 0.1, :He, (1.0, 0.0);
                                  capkw..., device=DeviceSpec(Array, Float32),
                                  tabulate_linop=true)
    tgrad = Luna.prop_capillary(125e-6, 0.1, :He, (1.0, 0.0);
                                capkw..., device=MetalSpec, tabulate_linop=true)
    for idx in axes(tgrad32["Eω"], 2)
        @test maximum(abs, tgrad["Eω"][:, idx] .- tgrad32["Eω"][:, idx]) /
              maximum(abs, tgrad32["Eω"][:, idx]) < 1e-4
    end
    @test isapprox(tgrad["stats"]["energy"], tgrad32["stats"]["energy"]; rtol=1e-3)

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

#= gpu/24's exit criterion on real hardware: every default statistic computed on the
   device, and no per-step copy of the state to the host to do it.

   The comparison is `prop_capillary`'s own setup on Metal against the same on the CPU in
   Float32 -- the same precision, so what is being compared is the device arithmetic and
   the device branches of `Stats.jl` against the host branches, not Float32 against
   Float64. The steps are fixed (built through `prop_capillary_args`, which
   `prop_capillary` itself has no keyword for), because in Float32 the adaptive
   controller's accept/reject decisions differ between the two and the runs then record
   their statistics at different `z`: the same *number* of steps, but not the same ones,
   which makes every statistic incomparable for a reason that has nothing to do with the
   statistics. The regression gate excludes `stats/z` and `stats/dz` from its adaptive
   comparison for the same reason (GPU_PLAN.md section 11).

   The copy count is the same instrument gpu/11's `stats_period` test uses:
   `ScaledOutput.ybuf` is the only place a device-to-host copy of the state lands and
   `_tohost_unscale!` is the only thing which writes it, so a sentinel which survives
   several calls is the count. =#
@testset "default statistics on Metal" begin
    capkw = (; λ0=800e-9, energy=300e-6, τfwhm=10e-15, λlims=(300e-9, 2000e-9),
             trange=400e-15, saveN=5, plasma=false, raman=false, shotnoise=false)
    function statsprop(spec; flength=0.1, stats_device=:auto)
        Eω, grid, linop, tr, FT, o = Luna.Interface.prop_capillary_args(
            125e-6, flength, :He, 5.0; capkw..., device=spec,
            stats_kwargs=Dict{Symbol, Any}(:stats_device => stats_device))
        h = flength/20
        Luna.run(Eω, grid, linop, tr, FT, o;
                 zmax=flength, init_dz=h, min_dz=h, max_dz=h, status_period=1e6)
        o
    end
    h32 = statsprop(DeviceSpec(Array, Float32))
    dm = statsprop(MetalSpec; stats_device=:device)
    @test length(dm["stats"]["z"]) == length(h32["stats"]["z"])
    @test sort(collect(keys(dm["stats"]))) == sort(collect(keys(h32["stats"])))
    for key in sort(collect(keys(h32["stats"])))
        h = h32["stats"][key]
        d = dm["stats"][key]
        scale = maximum(abs, h)
        err = scale > 0 ? maximum(abs, d .- h)/scale : maximum(abs, d .- h)
        @test (key, err <= 1e-3) == (key, true)
    end
    # the Kerr effect is visible over this propagation, so the comparison has content
    @test maximum(h32["stats"]["fwhm_t_min"])/minimum(h32["stats"]["fwhm_t_min"]) > 1.05

    #= The same set built for a Metal state, wrapped the way `Luna.run` wraps it, and
       called: nothing is copied down. =#
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(125e-6, :He, 5.0, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    ρ = PhysData.density(:He, 5.0)
    dens = z -> ρ
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, 800e-9)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=10e-15, energy=300e-6)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=MetalSpec)
    sf = Stats.default(grid, Eω, m, linop, transform; gas=:He, stats_device=:device)
    @test Stats.device_capable(sf)
    @test isempty(Stats.host_statistics(sf))
    out = Output.MemoryOutput(0, 1.0, 2, sf)
    so = Luna.ScaledOutput(out, Eω, Luna.runscaling(transform).Eref)
    @test so.devstats
    fill!(so.ybuf, 7)
    snap = copy(so.ybuf)
    for t in (0.0, 0.1, 0.2)
        so(Eω, t, 0.05, _ -> Eω)
    end
    @test so.ybuf == snap   # the state never reached the host
    @test !so.warned[]
    @test length(out["stats"]["z"]) == 3
    # the statistics really ran on the device, and in physical units
    @test out["stats"]["energy"][1] ≈ 300e-6 rtol=1e-2
    @test all(isfinite, out["stats"]["peakpower"])

    #= A user statistic has no device form, so the copy comes back and the warning names
       it. This is the fallback path on real hardware. =#
    uf = (d, Eω, Et, z, dz) -> d["mine"] = maximum(abs2, Et)
    sfu = Stats.default(grid, Eω, m, linop, transform;
                        gas=:He, userfuns=Any[uf], stats_device=:device)
    @test !Stats.device_capable(sfu)
    @test Stats.host_statistics(sfu) == ["userfuns[1]"]
    outu = Output.MemoryOutput(0, 1.0, 2, sfu)
    sou = Luna.ScaledOutput(outu, Eω, Luna.runscaling(transform).Eref)
    @test !sou.devstats
    #= ... and it is built for the host copy it will be handed, not for the device state
       it was constructed from: a device buffer with a host field in the same broadcast
       is what Metal refuses and JLArrays does not. =#
    @test sfu.Et isa Array
    @test_logs (:warn, r"have no device form") match_mode=:any begin
        sou(Eω, 0.0, 0.05, _ -> Eω)
    end
    @test sou.warned[]
    @test haskey(outu["stats"], "mine")

    #= `:auto` chooses the host path for this state: it is a single column well below
       `Stats.STATS_DEVICE_MINLEN`, where six device-to-host round trips cost more than
       the copy (review round 1, finding 3). =#
    sfa = Stats.default(grid, Eω, m, linop, transform; gas=:He)
    @test !Stats.device_capable(sfa)
    @test isempty(Stats.host_statistics(sfa))   # capable, but not the chosen path
    soa = Luna.ScaledOutput(Output.MemoryOutput(0, 1.0, 2, sfa), Eω,
                            Luna.runscaling(transform).Eref)
    @test !soa.devstats
end

#= Review round 1, finding 1: `prop_capillary(…; filepath=…)` builds an `HDF5Output` with
   a resume cache, which writes the per-step `y` to the file *and* passes it to its
   statistics function. On this branch's first version that made `ScaledOutput` hand a
   host array to a set built for the device, and Metal threw. Both paths are covered:
   `:auto` (the host path, which is what a user gets) and `:device` (the fixed one). =#
@testset "HDF5 file output with statistics on Metal" begin
    capkw = (; λ0=800e-9, energy=100e-9, τfwhm=10e-15, λlims=(300e-9, 2000e-9),
             trange=400e-15, saveN=5, plasma=false, raman=false, shotnoise=false)
    mem = Luna.prop_capillary(125e-6, 1e-2, :He, 1.0; capkw..., device=MetalSpec)
    mktempdir() do dir
        for (nm, sd) in (("auto", :auto), ("device", :device))
            fp = joinpath(dir, "metal_$nm.h5")
            out = Luna.prop_capillary(
                125e-6, 1e-2, :He, 1.0; capkw..., device=MetalSpec, filepath=fp,
                stats_kwargs=Dict{Symbol, Any}(:stats_device => sd))
            @test isfile(fp)
            @test eltype(out["Eω"]) === ComplexF32
            @test length(out["stats"]["z"]) == length(mem["stats"]["z"])
            for key in sort(collect(keys(mem["stats"])))
                h = mem["stats"][key]
                d = out["stats"][key]
                scale = maximum(abs, h)
                err = scale > 0 ? maximum(abs, d .- h)/scale : maximum(abs, d .- h)
                @test (nm, key, err <= 1e-3) == (nm, key, true)
            end
        end
    end
end

#= A multimode propagation on the fixed transverse quadrature rule
   (`NonlinearRHS.TransModalFixed`), which is the transverse integral with a device path:
   the adaptive cubature driver is host scalar code returning `Vector{Float64}`.

   Fixed steps, so that the only difference between two runs is the arithmetic. `nr=32`
   rather than the default 64 keeps the block small; the transverse integral of a smooth
   HE1m set converges long before that (test_device.jl measures it against the adaptive
   rule). =#
function metalmodalcase(spec; nmodes=4, components=:y, gas=:Ar, pres=0.1, energy=50e-6,
                        flength=2e-3, λ0=800e-9, plasma=false, nr=32, nθ=16, full=false,
                        precision=nothing, boundary=:none)
    grid = Grid.RealGrid(λ0, (200e-9, 3000e-9), 400e-15)
    modes = Tuple(Capillary.MarcatiliMode(75e-6, gas, pres; n=1, m=mi, loss=false)
                  for mi in 1:nmodes)
    ρ = PhysData.density(gas, pres)
    dens = z -> ρ
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    if plasma
        resp = (resp...,
                Nonlinear.PlasmaCumtrapz(grid.to, zeros(length(grid.to)),
                                         metal_tablerate(),
                                         PhysData.ionisation_potential(gas)))
    end
    linop = LinearOps.make_const_linop(grid, modes, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, modes, components;
                                   modal_integral=:fixed, nr, nθ, full,
                                   device=spec, precision)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    h = flength/10
    Luna.run(Eω, grid, linop, transform, FT, out;
             zmax=flength, boundary, min_dz=h, max_dz=h, init_dz=h)
    out, transform
end

#= Normalised to the strongest mode rather than per mode: in a four-mode capillary run
   the fourth mode carries ~1e-10 of the energy of the first, so its own relative
   difference measures the first mode's rounding against the fourth mode's amplitude. =#
function modaldiff(a, b)
    nrm = maximum(abs, a[:, :, end])
    maximum(abs, ComplexF64.(b[:, :, end]) .- ComplexF64.(a[:, :, end]))/nrm
end

@testset "multimode on Metal" begin
    @testset "plasma=$plasma" for plasma in (false, true)
        m32, mtr = metalmodalcase(MetalSpec; plasma)
        h32, _ = metalmodalcase(DeviceSpec(Array, Float32); plasma)
        h64, htr = metalmodalcase(HostSpec(); plasma)

        # everything the kernels touch is a Float32 device array
        @test mtr isa NonlinearRHS.TransModalFixed
        @test mtr.S isa MtlArray{Float32}
        @test mtr.Wp isa MtlArray{Float32}
        @test mtr.Wd isa MtlArray{Float32}
        @test mtr.Emt isa MtlArray{Float32}
        @test mtr.Emωo isa MtlArray{ComplexF32}
        @test mtr.Pmt isa MtlArray{Float32}
        @test mtr.block.Et isa MtlArray{Float32}
        @test mtr.block.Et2 isa MtlArray{Float32} # the reshape a GEMM needs
        @test mtr.block.Pt isa MtlArray{Float32}
        @test mtr.gv.ω isa MtlArray{Float32}
        @test mtr.norm!.pre isa MtlArray{ComplexF32}
        plasma && @test mtr.block.resp[2].J isa MtlArray{Float32}
        plasma && @test mtr.block.resp[2].ratedev.spline.x isa MtlArray{Float32}
        # the host Float64 run is unscaled; the Float32 ones are not
        @test Luna.isunity(htr.scaling)
        @test !Luna.isunity(mtr.scaling)

        @test size(m32["Eω"]) == size(h64["Eω"])
        @test eltype(m32["Eω"]) === ComplexF32
        # Metal against the same arithmetic on the CPU, and against Float64
        @test modaldiff(h32["Eω"], m32["Eω"]) < 1e-4
        @test modaldiff(h64["Eω"], m32["Eω"]) < 1e-3
    end

    #= Two polarisation components and the full 2-D rule: the θ nodes and the vector form
       of the Kerr response on a device block. =#
    mxy, mxytr = metalmodalcase(MetalSpec; nmodes=2, components=:xy, full=true, nr=16,
                                nθ=8)
    hxy32, _ = metalmodalcase(DeviceSpec(Array, Float32); nmodes=2, components=:xy,
                              full=true, nr=16, nθ=8)
    hxy64, _ = metalmodalcase(HostSpec(); nmodes=2, components=:xy, full=true, nr=16,
                              nθ=8)
    @test size(mxytr.block.Et) == (length(mxytr.grid.to), 2, 16*8)
    @test modaldiff(hxy32["Eω"], mxy["Eω"]) < 1e-4
    @test modaldiff(hxy64["Eω"], mxy["Eω"]) < 1e-3

    #= The absorbing boundaries, which `prop_capillary` turns on by default, on a
       multimode device state. =#
    mrate, _ = metalmodalcase(MetalSpec; boundary=:rate)
    hrate32, _ = metalmodalcase(DeviceSpec(Array, Float32); boundary=:rate)
    @test modaldiff(hrate32["Eω"], mrate["Eω"]) < 1e-4

    #= `modal_integral=:adaptive` -- the default -- on a device says what to do about it
       rather than failing somewhere inside Cubature. =#
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    ms = (Capillary.MarcatiliMode(75e-6, :He, 1.0, loss=false),)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6)
    err = try
        Luna.setup(grid, z -> 1.0, resp, inputs, ms, :y; device=MetalSpec)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("modal_integral=:fixed", err.msg)
end

#= The same thing through the simple interface, which is how a user gets it: the default
   response set of a field-resolved `prop_capillary` call in a non-Raman gas (Kerr, and
   Kerr plus plasma), four modes, the rate-based absorbing boundaries and the default
   statistics. `prop_capillary_args` plus `Luna.run` rather than `prop_capillary`, so
   that the step sequence can be fixed and the comparison is of arithmetic only. =#
@testset "prop_capillary multimode on Metal" begin
    args = (125e-6, 2e-3, :Ar, 0.1)
    kwargs = (λ0=800e-9, energy=50e-6, τfwhm=20e-15, trange=400e-15,
              λlims=(200e-9, 3000e-9), shotnoise=false, saveN=3, modes=4,
              modal_integral=:fixed, modal_nr=32)
    function runfixed(; kw...)
        Eω, grid, linop, transform, FT, output =
            Luna.Interface.prop_capillary_args(args...; kwargs..., kw...)
        h = args[2]/10
        Luna.run(Eω, grid, linop, transform, FT, output;
                 zmax=args[2], min_dz=h, max_dz=h, init_dz=h)
        output, transform
    end
    @testset "plasma=$plasma" for plasma in (false, true)
        m32, mtr = runfixed(; plasma, device=MetalSpec)
        h32, _ = runfixed(; plasma, device=DeviceSpec(Array, Float32))
        h64, _ = runfixed(; plasma, device=HostSpec())
        @test mtr isa NonlinearRHS.TransModalFixed
        @test mtr.block.Et isa MtlArray{Float32}
        @test eltype(m32["Eω"]) === ComplexF32
        @test modaldiff(h32["Eω"], m32["Eω"]) < 1e-4
        @test modaldiff(h64["Eω"], m32["Eω"]) < 1e-3
        # the statistics survive the device-to-host boundary
        @test isapprox(m32["stats"]["energy"][1, end], h32["stats"]["energy"][1, end];
                       rtol=1e-3)
        #= `Stats.mode_reconstruction_error` is the adaptive transform's; the fixed rule
           collects the other default statistics and not that one. =#
        @test !haskey(m32["stats"], "mode_reconstruction_error")
        @test haskey(m32["stats"], "fwhm_r")
    end

    # the adaptive default is refused on a device, naming the fix
    err = try
        Luna.Interface.prop_capillary_args(args...; kwargs..., modal_integral=:adaptive,
                                           device=MetalSpec)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("modal_integral=:fixed", err.msg)
end

#= `Luna.set_device(:cpu)` opts out, whatever `settings["device"]` is otherwise: this is
   exit criterion 3. `prop_gnlse`, radial `prop_capillary` and multimode `prop_capillary`
   with the *adaptive* transverse integral (the default) are not device-capable, and keep
   giving the CPU, Float64 answer under `:auto` too -- refusing only when the caller
   explicitly asks for something else. Multimode with `modal_integral=:fixed` is
   device-capable and does follow `:auto`, which is checked below. =#
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

        #= Multimode propagation with the adaptive transverse integral (the default) is
           not device-capable and must stay on the CPU by default under :auto too -- it
           must not turn a working run into an error just because a GPU package happens
           to be loaded (the bug an earlier version of this branch had: `device`'s
           default resolved through `:auto` even for paths that can never honour it). =#
        om = Luna.prop_capillary(capargs...; capkw..., modes=4)
        @test size(om["Eω"], 2) == 4
        @test eltype(om["Eω"]) === ComplexF64

        #= ... while multimode with `modal_integral=:fixed` *is* device-capable, so the
           sentinel resolves it to the GPU under :auto with no `device` keyword at all.
           This is the only hardware test of that branch of `Interface.prop_capillary_args`;
           everything else passes `device=MetalSpec` explicitly. =#
        omf = Luna.prop_capillary(capargs...; capkw..., modes=4,
                                  modal_integral=:fixed, modal_nr=16)
        @test size(omf["Eω"], 2) == 4
        @test eltype(omf["Eω"]) === ComplexF32
        # ... and `device=:cpu` still opts that out
        omfc = Luna.prop_capillary(capargs...; capkw..., modes=4, modal_integral=:fixed,
                                   modal_nr=16, device=:cpu)
        @test eltype(omfc["Eω"]) === ComplexF64

        # ... but an explicit device request for the adaptive rule still errors
        @test_throws ErrorException Luna.prop_capillary(capargs...; capkw..., modes=4,
                                                         device=MetalSpec)
    finally
        isnothing(old) ? delete!(Luna.settings, "device") :
                         (Luna.settings["device"] = old)
    end
end

#= ---------------------------------------------------------------- the radial transform =#

#= A small radially symmetric free-space propagation on Metal. `boundary=:rate` so that
   `Boundaries.RadialCollar`, the k-space absorber and the evanescent source taper run
   too; fixed steps, so two runs differ only in their arithmetic; `boundary_N` small
   enough that the absorber's reference length does not cap `max_dz` and undo that. =#
function metalradialcase(GT, spec; gas=:Ar, pres=1.0, energy=1e-6, flength=2e-3,
                         λ0=800e-9, R=1e-3, N=24, w0=200e-6, plasma=false, raman=false,
                         precision=nothing, boundary=:rate)
    grid = GT === Grid.RealGrid ?
        Grid.RealGrid(λ0, (400e-9, 2000e-9), 100e-15) :
        Grid.EnvGrid(λ0, (400e-9, 2000e-9), 100e-15)
    rg = Grid.RadialGrid(R, N)
    nfunλ = PhysData.ref_index_fun(gas, pres)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    linop = LinearOps.make_const_linop(grid, rg, nfun)
    ρ = PhysData.density(gas, pres)
    dens = z -> ρ
    resp = GT === Grid.RealGrid ?
        Any[Nonlinear.Kerr_field(PhysData.γ3_gas(gas))] :
        Any[Nonlinear.Kerr_env(PhysData.γ3_gas(gas))]
    plasma && push!(resp, Nonlinear.PlasmaCumtrapz(
        grid.to, zeros(length(grid.to)), metal_tablerate(),
        PhysData.ionisation_potential(gas)))
    if raman
        rr = Raman.raman_response(grid.to, gas)
        push!(resp, GT === Grid.RealGrid ? Nonlinear.RamanPolarField(grid.to, rr) :
                                           Nonlinear.RamanPolarEnv(grid.to, rr))
    end
    normfun = NonlinearRHS.const_norm_radial(grid, rg, nfun)
    inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy, w0, propz=-flength)
    Eω, transform, FT = Luna.setup(grid, rg, dens, normfun, Tuple(resp), inputs;
                                   device=spec, precision)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    h = flength/8
    Luna.run(Eω, grid, linop, transform, FT, out;
             zmax=flength, boundary, boundary_N=4, min_dz=h, max_dz=h, init_dz=h)
    out, transform
end

#= Per save, normalised by the largest `|Eω|` in that save: an elementwise relative
   difference is meaningless in the k-channels the evanescent taper has emptied. Named
   apart from `test_device.jl`'s `radialdiff`, which is the same idea without the
   `ComplexF64` conversion a `ComplexF32` saved field needs. =#
function metalradialdiff(a, b)
    A, B = a["Eω"], b["Eω"]
    size(A) == size(B) || return Inf
    worst = 0.0
    for isave in axes(A, ndims(A))
        h = ComplexF64.(selectdim(A, ndims(A), isave))
        d = ComplexF64.(selectdim(B, ndims(B), isave))
        m = maximum(abs, h)
        m == 0 && continue
        worst = max(worst, maximum(abs, d .- h)/m)
    end
    worst
end

#= The Hankel step is one GEMM on the block reshaped to `(nto*npol, nr)`. On Metal that
   matters twice over: MPSGraph's matmul needs plain zero-offset operands of equal element
   type, which a `view` is not (it becomes an `MtlMatrixOperand` and falls back to a
   scalar kernel for complex), and `ComplexF32 x ComplexF32` is the least-exercised of the
   paths Luna uses (GPU_PLAN.md section 8). Checked against the host product here rather
   than only through a propagation. =#
@testset "the Hankel GEMM on Metal" begin
    rg = Grid.RadialGrid(1e-3, 32)
    for np in (1, 2), TT in (Float32, ComplexF32)
        A = TT <: Complex ? complex.(randn(Float32, 64, np, rg.N),
                                     randn(Float32, 64, np, rg.N)) :
                            randn(Float32, 64, np, rg.N)
        Tm = convert(Matrix{TT}, rg.Tfwd)
        href = similar(A)
        Grid.radial_matmul!(href, A, Tm)
        dA = MtlArray(A)
        dT = MtlArray(Tm)
        # The operands the GEMM actually sees: plain matrices, no view, no offset
        @test reshape(dA, :, rg.N) isa MtlArray{TT, 2}
        @test eltype(dT) === eltype(dA)
        dout = similar(dA)
        Grid.radial_matmul!(dout, dA, dT)
        @test maximum(abs, Array(dout) .- href)/maximum(abs, href) < 1e-5
        # out === A: `radial_matmul!` copies, which works on a device too
        Grid.radial_matmul!(dA, dA, dT)
        @test maximum(abs, Array(dA) .- href)/maximum(abs, href) < 1e-5
    end
end

#= The stray-Float64 smoke test for the radial pieces: the transform's own buffers and
   mirrors, the free-space normalisation (whose kernel carries `c`, `μ₀`, `κmax` and `ℓ`
   as captured scalars) and the transverse collar. Metal refuses a `Float64` array and its
   kernel compiler rejects any `double` which survives optimisation, so an element type
   here which is not `Float32`/`ComplexF32`/`Bool` is the failure this file exists for. =#
@testset "no stray Float64 in the radial kernels" begin
    #= A fine transverse grid over a small aperture, so that the largest k⊥ on the grid
       exceeds k(ω) at the long-wavelength end and the evanescent branch of the kernel --
       the one which takes the taper -- is actually reached. =#
    grid = Grid.RealGrid(800e-9, (400e-9, 4000e-9), 100e-15)
    rg = Grid.RadialGrid(100e-6, 64)
    nfunλ = PhysData.ref_index_fun(:Ar, 1.0)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    nrm = NonlinearRHS.const_norm_radial(grid, rg, nfun; spec=MetalSpec)
    @test nrm.out isa MtlArray{ComplexF32, 3}
    @test nrm.kperp2m isa MtlArray{Float32}
    @test nrm.kwinm isa MtlArray{Float32}
    @test nrm.sidxm isa MtlArray{Bool}
    @test nrm.ωm isa MtlArray{Float32}
    @test nrm.nm.stage isa Vector{Float32}
    @test nrm.nm.dev isa MtlArray{Float32, 1}
    out = nrm(0.0) # compiles and runs the fill kernel on the GPU
    @test out isa MtlArray{ComplexF32, 3}
    @test all(isfinite, Array(out))
    # the physics, on both sides of cutoff, against the host at the same parameters
    hnrm = NonlinearRHS.const_norm_radial(grid, rg, nfun)
    hout = hnrm(0.0)
    @test maximum(abs, ComplexF64.(Array(out)) .- hout)/maximum(abs, hout) < 1e-5
    @test any(x -> imag(x) != 0, hout) # there really are evanescent channels here

    #= The taper branch of the kernel: `reflength!` sets `ℓ`, `κmax` and the k-window and
       invalidates the mirror, which the next fill rebuilds. =#
    kwin = Boundaries.kprofile(rg, 0.1)
    before = copy(Array(out))
    NonlinearRHS.reflength!(nrm, 1e-3; κmax=1e4, kwin)
    NonlinearRHS.reflength!(hnrm, 1e-3; κmax=1e4, kwin)
    @test !nrm.mirrored
    out2 = nrm(0.0)
    @test nrm.mirrored
    hout2 = hnrm(0.0)
    @test maximum(abs, ComplexF64.(Array(out2)) .- hout2)/maximum(abs, hout2) < 1e-5
    @test maximum(abs, Array(out2) .- before) > 0 # the taper did something

    # the transverse collar, built the way `Boundaries.setup` builds it
    Et = Luna.alloc(MetalSpec, Float32, (length(grid.t), 1, rg.N))
    αr = Boundaries.rate(Boundaries.rprofile(rg, 0.1), 1e-3)
    collar = Boundaries.spatialcollar(rg, αr, grid, Et)
    @test collar isa Boundaries.RadialCollar
    @test collar.Tfwd isa MtlArray{ComplexF32, 2}
    @test collar.Tbwd isa MtlArray{ComplexF32, 2}
    @test collar.αr isa MtlArray{Float32, 1}
    @test collar.weight isa MtlArray{Float32, 1}
    @test collar.fac isa MtlArray{Float32, 1}
    @test collar.buf isa MtlArray{ComplexF32, 3}
    Eωd = Luna.todevice(MetalSpec, rand(ComplexF64, length(grid.ω), 1, rg.N))
    Boundaries.apply_kspace!(collar, Eωd, 1e-4)
    @test all(isfinite, Array(Eωd))
    @test collar.reference[] > 0
    @test isfinite(collar.removed[])
end

@testset "radial propagation on Metal" begin
    #= Metal against the CPU at the same precision and the same scaling: what is under
       test is the device path, not single precision. =#
    for GT in (Grid.RealGrid, Grid.EnvGrid)
        href, htr = metalradialcase(GT, DeviceSpec(Array, Float32))
        dref, dtr = metalradialcase(GT, MetalSpec)

        @test dtr.Eto_r isa MtlArray
        @test dtr.Eto_k isa MtlArray
        @test dtr.Eωo isa MtlArray{ComplexF32, 3}
        @test dtr.Tfwd isa MtlArray{eltype(dtr.Eto_r), 2}
        @test dtr.Tbwd isa MtlArray{eltype(dtr.Eto_r), 2}
        @test dtr.prefac isa MtlArray{ComplexF32, 1}
        @test dtr.gv.towin isa MtlArray{Float32}
        @test dtr.normfun.out isa MtlArray{ComplexF32, 3}
        @test dtr.scaling.Eref == htr.scaling.Eref
        @test eltype(dref["Eω"]) === ComplexF32

        @test size(dref["Eω"]) == size(href["Eω"])
        @test metalradialdiff(href, dref) < 1e-4
    end
end

@testset "radial Kerr on Metal against the Float64 CPU path" begin
    href, _ = metalradialcase(Grid.RealGrid, HostSpec())
    dref, _ = metalradialcase(Grid.RealGrid, MetalSpec)
    @test metalradialdiff(href, dref) < 1e-4
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

    #= multimode propagation with the adaptive transverse integral is not device-capable
       through the simple interface either (with `modal_integral=:fixed` it is, which the
       "prop_capillary multimode on Metal" testset above covers) =#
    @test_throws ErrorException Luna.prop_capillary(
        125e-6, 1e-3, :He, 1.0; λ0=800e-9, energy=1e-9, τfwhm=10e-15,
        λlims=(300e-9, 2e-6), trange=400e-15, saveN=3, plasma=false, shotnoise=false,
        modes=4, device=MetalSpec)
end

#= ------------------------------------------------- the Cartesian free-space transforms =#

#= The multi-axis FFT plans `TransFree`/`TransFree2D` need, against FFTW on the same random
   data. `(1, 3)` and `(1, 3, 4)` -- the time axis and the one or two transverse axes, with
   the polarisation axis skipped -- are what Luna plans; Metal.jl's own tests cover `(1, 3)`
   and `(1, 4)` but not `(1, 3, 4)` (GPU_PLAN.md section 2), so the three-axis real and
   complex regions are checked here before anything is built on them. Forward and inverse,
   one and two polarisation components, powers of two and not.

   This runs first in the free-space part of the file: if a region were unsupported or
   wrong, the transforms below would have to be built out of two plans instead of one. =#
@testset "multi-axis FFT plans on Metal" begin
    shapes = (((64, 1, 16, 8), (1, 3, 4)),
              ((64, 2, 16, 8), (1, 3, 4)),
              ((96, 2, 15, 7), (1, 3, 4)), # transverse axes not powers of two
              ((64, 1, 32), (1, 3)),
              ((64, 2, 32), (1, 3)),
              ((385, 1, 24), (1, 3)))      # odd time axis: real output length matters
    for (sz, region) in shapes, TT in (Float32, ComplexF32)
        A = TT <: Complex ? complex.(randn(Float32, sz), randn(Float32, sz)) :
                            randn(Float32, sz)
        Ah = convert(Array{TT <: Complex ? ComplexF64 : Float64}, A)
        # the host plan Luna would make for the same buffer, in double precision
        FTh = Utils.plan_ft(Ah, region)
        href = FTh * Ah

        dA = MtlArray(A)
        FT = Utils.plan_ft(dA, region)
        dout = similar(dA, ComplexF32, size(href))
        LinearAlgebra.mul!(dout, FT, dA) # plain arrays, exactly the planned shape
        @test size(dout) == size(href)
        @test maximum(abs, ComplexF64.(Array(dout)) .- href)/maximum(abs, href) < 1e-5

        #= The inverse plan Luna holds explicitly: an `AbstractFFTs.ScaledPlan` whose
           `scale` is 1/N over the *output* lengths, which for a real transform is not the
           input's. `to_time!` folds that scale into the oversampling copy and applies
           `Utils.iplan(IFT)`, so both halves are checked separately here. =#
        IFT = Utils.plan_ift(FT)
        @test Utils.iscale(IFT) ≈ 1/prod(sz[i] for i in region)
        back = similar(dA)
        LinearAlgebra.mul!(back, Utils.iplan(IFT), copy(dout))
        B = Array(back) .* Float32(Utils.iscale(IFT))
        @test maximum(abs, convert(Array{eltype(Ah)}, B) .- Ah)/maximum(abs, Ah) < 1e-5
    end
end

#= The χ⁽²⁾ response matching the grid, built from the grid's own carrier frequency and
   oversampled time axis (which `Chi2Env` needs for the carrier phase). =#
metalchi2(::Grid.RealGrid, θ, ϕ) = Nonlinear.Chi2Field(θ, ϕ, PhysData.χ2(:BBO))
metalchi2(grid::Grid.EnvGrid, θ, ϕ) = Nonlinear.Chi2Env(θ, ϕ, PhysData.χ2(:BBO),
                                                        grid.ω0, grid.to)

#= Type I SHG in BBO on a `Grid.Free2DGrid`, field-resolved or envelope: the 2-D Cartesian
   transform with the two-component χ⁽²⁾ response and the crystal-optics normalisation
   (whose fill is host root-finding, staged through `ohost` and uploaded). `boundary=:rate`
   so the k-space absorber, the evanescent source taper and `Boundaries.CartesianCollar`
   run too; fixed steps, so two runs differ only in their arithmetic. =#
const METAL_BBO_θ = deg2rad(29.2)
const METAL_BBO_ϕ = deg2rad(30)

function metalfree2dcase(GT, spec; λ0=800e-9, τfwhm=30e-15, w0=20e-6, energy=10e-9,
                         thickness=30e-6, Nx=2^5, precision=nothing, boundary=:rate)
    grid = GT === Grid.RealGrid ?
        Grid.RealGrid(λ0, (250e-9, 2e-6), 120e-15) :
        Grid.EnvGrid(λ0, (250e-9, 2e-6), 120e-15; thg=true)
    xgrid = Grid.Free2DGrid(4w0, Nx)
    nfuns = PhysData.ref_index_fun_xy(:BBO, METAL_BBO_θ)
    linop = LinearOps.make_const_linop(grid, xgrid, nfuns)
    normfun = NonlinearRHS.const_norm_free2D(grid, xgrid, nfuns)
    densityfun = z -> 1 # unity density: this is a solid
    resp = (metalchi2(grid, METAL_BBO_θ, METAL_BBO_ϕ),)
    inputs = Fields.GaussGaussField(;λ0, τfwhm, energy=energy/(sqrt(π/2)*w0), w0)
    Eω, transform, FT = Luna.setup(grid, xgrid, densityfun, normfun, resp, inputs;
                                   device=spec, precision)
    out = Output.MemoryOutput(0, thickness, 3, Output.nostats)
    h = thickness/8
    Luna.run(Eω, grid, linop, transform, FT, out;
             zmax=thickness, boundary, boundary_N=4, min_dz=h, max_dz=h, init_dz=h)
    out, transform
end

#= A small 3-D free-space envelope Kerr propagation on a `Grid.FreeGrid`: the block is
   `(nto, npol, Nx, Ny)` and the transform is one region-(1,3,4) FFT each way. =#
function metalfree3dcase(spec; gas=:Ar, pres=1.0, λ0=800e-9, energy=1e-9, flength=2e-3,
                         R=1e-3, Nx=8, Ny=6, w0=200e-6, precision=nothing, boundary=:rate)
    grid = Grid.EnvGrid(λ0, (400e-9, 2000e-9), 100e-15)
    xygrid = Grid.FreeGrid(R, Nx, R, Ny)
    nfunλ = PhysData.ref_index_fun(gas, pres)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    linop = LinearOps.make_const_linop(grid, xygrid, nfun)
    ρ = PhysData.density(gas, pres)
    dens = z -> ρ
    resp = (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),)
    normfun = NonlinearRHS.const_norm_free(grid, xygrid, nfun)
    inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy, w0, propz=-flength)
    Eω, transform, FT = Luna.setup(grid, xygrid, dens, normfun, resp, inputs;
                                   device=spec, precision)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    h = flength/8
    Luna.run(Eω, grid, linop, transform, FT, out;
             zmax=flength, boundary, boundary_N=4, min_dz=h, max_dz=h, init_dz=h)
    out, transform
end

#= The stray-Float64 smoke test for the Cartesian free-space pieces: the transform's own
   buffers and mirrors, both normalisation fills (the isotropic broadcast and the
   crystal-optics host fill with its `ohost` staging buffer) and the transverse collar.
   Metal refuses a Float64 array and its kernel compiler rejects any `double` which
   survives optimisation, so an element type here which is not Float32/ComplexF32/Bool is
   the failure this file exists for. =#
@testset "no stray Float64 in the Cartesian free-space kernels" begin
    #= A fine transverse grid over a small box, so that the largest k⊥ on the grid exceeds
       k(ω) at the long-wavelength end and the evanescent branch of the kernel -- the one
       which takes the taper -- is actually reached. =#
    grid = Grid.RealGrid(800e-9, (400e-9, 4000e-9), 100e-15)
    xygrid = Grid.FreeGrid(10e-6, 32, 10e-6, 16)
    nfunλ = PhysData.ref_index_fun(:Ar, 1.0)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    nrm = NonlinearRHS.const_norm_free(grid, xygrid, nfun; spec=MetalSpec)
    @test nrm.out isa MtlArray{ComplexF32, 4}
    @test nrm.kperp2m isa MtlArray{Float32, 4}
    @test nrm.kwinm isa MtlArray{Float32, 4}
    @test nrm.sidxm isa MtlArray{Bool}
    @test nrm.ωm isa MtlArray{Float32}
    @test nrm.nm.dev isa MtlArray{Float32, 1}
    out = nrm(0.0) # compiles and runs the fill kernel on the GPU
    @test all(isfinite, Array(out))
    hnrm = NonlinearRHS.const_norm_free(grid, xygrid, nfun)
    hout = hnrm(0.0)
    @test maximum(abs, ComplexF64.(Array(out)) .- hout)/maximum(abs, hout) < 1e-5
    @test any(x -> imag(x) != 0, hout) # there really are evanescent channels here

    #= The taper branch: `reflength!` sets ℓ, κmax and the k-window and invalidates the
       mirror. The window is clamped at its floor exactly as `Boundaries.setup` clamps it
       -- the raw profile is exactly zero at the Nyquist wavevector of an FFT grid, and
       the normalisation divides by it. =#
    kwin = max.(Boundaries.kprofile(xygrid, 0.1), exp(-Boundaries.MAX_αℓ/2))
    before = copy(Array(out))
    NonlinearRHS.reflength!(nrm, 1e-3; κmax=1e4, kwin)
    NonlinearRHS.reflength!(hnrm, 1e-3; κmax=1e4, kwin)
    @test !nrm.mirrored
    out2 = nrm(0.0)
    @test nrm.mirrored
    hout2 = hnrm(0.0)
    @test maximum(abs, ComplexF64.(Array(out2)) .- hout2)/maximum(abs, hout2) < 1e-5
    @test maximum(abs, Array(out2) .- before) > 0 # the taper did something

    #= The crystal-optics fill: host root-finding per (ω, kx), staged through `ohost` and
       uploaded. `ohost` is a host ComplexF64 buffer by design -- it is never broadcast
       against the state -- and what reaches the device is `out`. =#
    bgrid = Grid.RealGrid(800e-9, (250e-9, 2e-6), 120e-15)
    xgrid = Grid.Free2DGrid(80e-6, 2^5)
    nfuns = PhysData.ref_index_fun_xy(:BBO, METAL_BBO_θ)
    cnrm = NonlinearRHS.const_norm_free2D(bgrid, xgrid, nfuns; spec=MetalSpec)
    @test cnrm.ohost isa Array{ComplexF64, 3}
    cout = cnrm(0.0)
    @test cout isa MtlArray{ComplexF32, 3}
    @test size(cout, 2) == 2 # crystal optics is a two-polarisation normalisation
    hcnrm = NonlinearRHS.const_norm_free2D(bgrid, xgrid, nfuns)
    hcout = hcnrm(0.0)
    @test maximum(abs, ComplexF64.(Array(cout)) .- hcout)/maximum(abs, hcout) < 1e-5

    # the transverse collar, built the way `Boundaries.setup` builds it, in both shapes
    for sg in (xgrid, xygrid)
        xyshape = sg isa Grid.Free2DGrid ? (length(sg.x),) : (length(sg.x), length(sg.y))
        Et = Luna.alloc(MetalSpec, Float32, (length(grid.t), 2, xyshape...))
        αr = Boundaries.rate(Boundaries.rprofile(sg, 0.1), 1e-3)
        collar = Boundaries.spatialcollar(sg, αr, grid, Et)
        @test collar isa Boundaries.CartesianCollar
        @test collar.αxy isa MtlArray{Float32}
        @test collar.fac isa MtlArray{Float32}
        @test size(collar.αxy) == xyshape
        copyto!(Et, randn(Float32, size(Et)))
        Boundaries.apply_realspace!(collar, Et, 1e-4)
        @test all(isfinite, Array(Et))
        @test collar.reference[] > 0
        @test collar.removed[] > 0
    end
end

#= The exit test of the branch on hardware: 2-D and 3-D Cartesian free space, Kerr and
   χ⁽²⁾, `boundary=:rate`, through the low-level interface. Metal against the CPU at the
   *same* precision and the same scaling -- an explicit `DeviceSpec(Array, Float32)`, not
   the sentinel, which with Metal loaded would resolve to the GPU and compare Metal with
   itself. =#
@testset "2-D free-space χ⁽²⁾ on Metal" begin
    for GT in (Grid.RealGrid, Grid.EnvGrid)
        href, htr = metalfree2dcase(GT, DeviceSpec(Array, Float32))
        dref, dtr = metalfree2dcase(GT, MetalSpec)

        @test dtr isa NonlinearRHS.TransFree2D
        @test dtr.Eto isa MtlArray
        @test dtr.Eωo isa MtlArray{ComplexF32, 3}
        @test dtr.Pωo === dtr.Eωo # one field-sized buffer fewer
        @test dtr.prefac isa MtlArray{ComplexF32, 1}
        @test dtr.gv.towin isa MtlArray{Float32}
        @test dtr.normfun.out isa MtlArray{ComplexF32, 3}
        @test dtr.scaling.Eref == htr.scaling.Eref
        @test eltype(dref["Eω"]) === ComplexF32

        # both polarisation components carry field: the second harmonic is on x
        for ip in 1:2
            @test maximum(abs, href["Eω"][:, ip, :, end]) > 0
        end
        @test size(dref["Eω"]) == size(href["Eω"])
        @test metalradialdiff(href, dref) < 1e-4
    end
end

@testset "3-D free-space Kerr on Metal" begin
    href, htr = metalfree3dcase(DeviceSpec(Array, Float32))
    dref, dtr = metalfree3dcase(MetalSpec)

    @test dtr isa NonlinearRHS.TransFree
    @test dtr.Eto isa MtlArray{ComplexF32, 4}
    @test dtr.Eωo isa MtlArray{ComplexF32, 4}
    @test dtr.Pωo === dtr.Eωo
    @test dtr.normfun.out isa MtlArray{ComplexF32, 4}
    @test ndims(dref["Eω"]) == 5 # (ω, pol, kx, ky, z)
    @test size(dref["Eω"]) == size(href["Eω"])
    @test metalradialdiff(href, dref) < 1e-4
end

@testset "free space on Metal against the Float64 CPU path" begin
    href, _ = metalfree2dcase(Grid.RealGrid, HostSpec())
    dref, _ = metalfree2dcase(Grid.RealGrid, MetalSpec)
    @test metalradialdiff(href, dref) < 1e-4

    h3, _ = metalfree3dcase(HostSpec())
    d3, _ = metalfree3dcase(MetalSpec)
    @test metalradialdiff(h3, d3) < 1e-4
end

end # have_metal
