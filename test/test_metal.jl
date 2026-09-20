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
             UnitScaling, UNIT_SCALING
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
                   extraresp=())
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
