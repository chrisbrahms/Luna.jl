#= Tests for Luna's device (precision and array type) model.

   Luna's per-step code is made of broadcasts, planned FFTs applied with `mul!` and
   reductions, so it can be exercised without GPU hardware:

   * the `Utils.backend` trait is purely type-level, so a dummy `AbstractGPUArray`
     subtype is enough to test its dispatch;
   * `JLArrays.JLArray` is a host-backed `AbstractGPUArray` which takes the same generic
     GPUArrays code path a real device does and, with `allowscalar(false)`, enforces the
     same no-scalar-indexing contract. It is a test-only dependency, so these tests are
     skipped if it cannot be loaded.

   `Pkg.test()` runs this file with everything it needs (`JLArrays` is in
   `[extras]`/`[targets]`; the rest are Luna's own dependencies). To run it on its own,
   the environment needs Luna developed into it plus the packages this file imports
   directly:

       julia --project=<env> -e '
         using Pkg
         Pkg.develop(path=".")
         Pkg.add(["JLArrays", "Test", "GPUArraysCore", "Adapt", "AbstractFFTs", "FFTW"])'
       julia --project=<env> -t 1 -e 'using Luna; include("test/test_device.jl")'


   What this file cannot catch, and what `test/test_metal.jl` must: JLArrays interprets
   its kernels on the host, so a stray `Float64` in a struct field or a mixed host/device
   broadcast passes here and would still fail on real hardware. =#

import Test: @test, @testset, @test_throws, @test_logs, @inferred
import Luna
import Luna: Utils, Output, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear,
             NonlinearRHS, PhysData, RK45, Stats
import Luna: DeviceSpec, HostSpec, UnitScaling, UNIT_SCALING
import GPUArraysCore
import AbstractFFTs
import FFTW
import Adapt
import LinearAlgebra
import LinearAlgebra: mul!

# A minimal device-array type: enough to test the `backend` trait, which never touches
# the data.
struct DummyGPUArray{T, N} <: GPUArraysCore.AbstractGPUArray{T, N}
    a::Array{T, N}
end
Base.size(x::DummyGPUArray) = size(x.a)

@testset "backend trait" begin
    host = zeros(ComplexF64, 4, 3, 2)
    dev = DummyGPUArray(host)

    @test Utils.backend(host) isa Utils.CPUBackend
    @test Utils.backend(dev) isa Utils.DeviceBackend
    @test !Utils.isdevice(host)
    @test Utils.isdevice(dev)

    # Types dispatch identically to instances
    @test Utils.backend(typeof(host)) isa Utils.CPUBackend
    @test Utils.backend(typeof(dev)) isa Utils.DeviceBackend

    # Wrappers report their parent's backend
    @test Utils.backend(view(host, :, 1, 1)) isa Utils.CPUBackend
    @test Utils.backend(view(dev, :, 1, 1)) isa Utils.DeviceBackend
    @test Utils.backend(reshape(host, 24)) isa Utils.CPUBackend
    @test Utils.backend(reshape(dev, 24)) isa Utils.DeviceBackend
    @test Utils.backend(PermutedDimsArray(host, (3, 2, 1))) isa Utils.CPUBackend
    @test Utils.backend(PermutedDimsArray(dev, (3, 2, 1))) isa Utils.DeviceBackend

    # Unknown array types count as host: reading a host wrapper as CPU is a no-op,
    # whereas the reverse would break it
    @test Utils.backend(1:10) isa Utils.CPUBackend
    @test Utils.backend(reshape(1:24, 4, 6)) isa Utils.CPUBackend

    # Resolved at compile time, so it costs nothing where it is used
    @test @inferred(Utils.backend(host)) isa Utils.CPUBackend
    @test @inferred(Utils.backend(dev)) isa Utils.DeviceBackend
end

@testset "device spec and settings" begin
    s = HostSpec()
    @test Luna.arraytype(s) === Array
    @test Luna.realtype(s) === Float64
    @test !Luna.isdevicespec(s)
    @test Luna.realtype(Luna.withprecision(s, Float32)) === Float32
    @test Luna.withprecision(s, nothing) === s

    # No GPU package loaded here, so every route resolves to the host
    @test Luna.resolve_device(:cpu) === HostSpec()
    @test Luna.resolve_device(:auto) === HostSpec()
    @test Luna.resolve_device(s) === s
    @test_throws ErrorException Luna.resolve_device(:metal)
    @test_throws ErrorException Luna.resolve_device("gpu")

    #= The initial state is no key at all, which `device()` reads as :cpu. That matters:
       a GPU extension sets the key to :auto only when it is absent, so an explicit
       set_device(:cpu) is never overridden. =#
    had = haskey(Luna.settings, "device")
    old = get(Luna.settings, "device", nothing)
    try
        delete!(Luna.settings, "device")
        @test Luna.device() === HostSpec()
        Luna.set_device(:cpu)
        @test Luna.settings["device"] === :cpu
        @test Luna.device() === HostSpec()
        @test_throws ErrorException Luna.set_device(:nonsense)
        @test Luna.settings["device"] === :cpu # a failed set changes nothing
    finally
        had ? (Luna.settings["device"] = old) : delete!(Luna.settings, "device")
    end

    #= The *unresolved* request is what `Luna.setup` defaults to, so that `log_device`
       can tell `:auto` which found no GPU from an explicit `:cpu`. =#
    had = haskey(Luna.settings, "device")
    old = get(Luna.settings, "device", nothing)
    try
        delete!(Luna.settings, "device")
        @test Luna.device_request() === :cpu
        Luna.settings["device"] = :auto
        @test Luna.device_request() === :auto
        @test Luna.device() === HostSpec() # no GPU package loaded here
        #= GPU_PLAN.md section 3 requires this message: it is what a `Scans` worker which
           only did `using Luna` sees, and section 8 relies on it as the safety net
           against `:auto` picking a GPU unexpectedly. =#
        @test_logs (:info, r"no GPU package is loaded") match_mode=:any begin
            Luna.log_device(Luna.resolve_device(:auto), :auto)
        end
        @test_logs (:info, r"precision") match_mode=:any begin
            Luna.log_device(HostSpec(), :cpu)
        end
    finally
        had ? (Luna.settings["device"] = old) : delete!(Luna.settings, "device")
    end

    # The vendor hooks are no-ops without a registered device
    @test Luna.device_synchronize(HostSpec()) === nothing
    @test Luna.device_reclaim(HostSpec()) === nothing
    @test Luna.device_memory_status(HostSpec()) === nothing
end

@testset "allocation and transfer" begin
    s = HostSpec()
    x = Luna.alloc(s, ComplexF64, (4, 2))
    @test x isa Matrix{ComplexF64}
    @test all(iszero, x)
    @test Luna.alloc(s, Float64, 3) isa Vector{Float64}

    # On the host in Float64 nothing is copied or converted
    v = rand(5)
    @test Luna.todevice(s, v) === v
    b = BitVector([true, false, true, true, false])
    @test Luna.todevice(s, b) === b
    @test Luna.todevice(s, nothing) === nothing
    @test Luna.tohost(v) === v

    s32 = DeviceSpec(Array, Float32)
    @test Luna.todevice(s32, v) isa Vector{Float32}
    @test Luna.todevice(s32, complex(v)) isa Vector{ComplexF32}
    @test Luna.todevice(s32, b) === b # masks keep their element type

    @test Luna.scalar(zeros(Float32, 2), 1.5) === 1.5f0
    @test Luna.scalar(zeros(ComplexF32, 2), 1.5) === 1.5f0
    @test Luna.scalar(zeros(ComplexF64, 2), 1.5) === 1.5

    # upload_like: the linear operator in the state's array type and precision
    y = zeros(ComplexF64, 6)
    l = rand(ComplexF64, 6)
    @test Luna.upload_like(y, l) === l # host Float64: untouched
    y32 = zeros(ComplexF32, 6)
    l32 = Luna.upload_like(y32, l)
    @test l32 isa Vector{ComplexF32}
    @test l32 ≈ l
    @test Luna.upload_like(y, z -> z) isa Function # a closure operator passes through

    @test Luna.mask_like(y, b) === b
end

@testset "residency assertions" begin
    s = HostSpec()
    @test Luna.assert_resident(s, zeros(4), zeros(ComplexF64, 4), BitVector([true])) ===
          nothing
    @test Luna.assert_resident(s, nothing) === nothing
    # Wrong precision
    @test_throws ErrorException Luna.assert_resident(s, zeros(Float32, 4))
    # Not an array at all
    @test_throws ErrorException Luna.assert_resident(s, 1.0)
    # A device array where a host one is required
    @test_throws ErrorException Luna.assert_resident(s, DummyGPUArray(zeros(4)))
    # A host array where a device one is required
    s32 = DeviceSpec(Array, Float32)
    @test_throws ErrorException Luna.assert_resident(s32, zeros(Float64, 4))
end

@testset "unit scaling" begin
    @test Luna.isunity(UNIT_SCALING)
    # Float64 never scales, whatever the field
    @test Luna.unitscaling(Float64, [1e7, 2e7], PhysData.ε_0) === UNIT_SCALING
    @test Luna.unitscaling(Float64, () -> error("must not be evaluated"),
                           PhysData.ε_0) === UNIT_SCALING
    # Float32 takes the peak of the time-domain field, rounded to a power of two
    sc = Luna.unitscaling(Float32, [300.0, -100.0], PhysData.ε_0)
    @test sc.Eref == 256.0
    @test log2(sc.Eref) == round(log2(sc.Eref)) # exact scaling in every linear operation
    @test sc.Pref == PhysData.ε_0
    @test !Luna.isunity(sc)
    # Degenerate input falls back to 1
    @test Luna.unitscaling(Float32, [0.0, 0.0], PhysData.ε_0).Eref == 1.0
end

@testset "FFT planner dispatch" begin
    #= The host planner takes Luna's FFTW flags, the device planner must not (device FFT
       libraries reject them), and `plan_ift` splits the inverse into an unnormalised
       plan and its 1/N so that the factor can be folded into the oversampling copy. =#
    xr = zeros(Float64, 16)
    pr = Utils.plan_ft(xr, 1)
    @test pr isa FFTW.rFFTWPlan
    ipr = Utils.plan_ift(pr)
    @test Utils.iscale(ipr) == 1/16
    @test Utils.iplan(ipr) !== ipr

    xc = zeros(ComplexF64, 16)
    pc = Utils.plan_ft(xc, 1)
    @test pc isa FFTW.cFFTWPlan
    @test Utils.iscale(Utils.plan_ift(pc)) == 1/16

    #= A forward plan where an inverse one belongs is refused rather than transforming
       the wrong way (a method error on a real grid, silently wrong output on an envelope
       grid). Dispatch decides, so the check costs nothing per step. =#
    @test_throws ErrorException Utils.iscale(pr)
    @test_throws ErrorException Utils.iplan(pr)
    @test_throws ErrorException NonlinearRHS.to_time!(zeros(16), rand(ComplexF64, 9),
                                                      zeros(ComplexF64, 9), pc)

    #= The whole point of folding: with a power-of-two length the two routes agree
       bitwise, which is why the default CPU path did not move. =#
    v = rand(ComplexF64, 9)
    a = similar(xr)
    b = similar(xr)
    buf = zeros(ComplexF64, 9)
    NonlinearRHS.to_time!(a, v, buf, ipr)
    fill!(buf, 0); NonlinearRHS.copy_scale!(buf, v, 9, 1.0)
    LinearAlgebra.ldiv!(b, pr, buf)
    @test a == b
end

#= One small mode-averaged Kerr propagation, built once and run on whichever spec is
   asked for. `Luna.run` wraps `out` in `ScaledOutput` itself, whenever the state is on a
   device or the run is scaled, so the test only ever sees the plain `Output.MemoryOutput`
   it built -- already on the host, already unscaled. `boundary=:none` by default: the
   `:rate` case, which exercises `Boundaries.RateAbsorber`, has its own testset below. =#
function kerrcase(GT, spec; gas=:He, pres=1.0, energy=1e-6, flength=1e-2, λ0=800e-9,
                  precision=nothing, boundary=:none, stats=false)
    grid = GT === Grid.RealGrid ?
        Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15) :
        Grid.EnvGrid(λ0, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, gas, pres, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    #= Precomputed: `PhysData.density` goes through CoolProp and costs more per call
       than the whole right-hand side. =#
    ρ = PhysData.density(gas, pres)
    dens = z -> ρ
    resp = GT === Grid.RealGrid ?
        (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),) :
        (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=spec, precision)
    #= Stats.jl is host-only: its EnvGrid plan_analytic builds an FFTW plan directly on
       a copy of the given Eω, so construction needs a host-shaped template, not the
       (possibly device) state itself. See Interface.jl's prop_capillary_args. =#
    shost = Utils.isdevice(Eω) ? Luna.tohost(Eω) : Eω
    statsfun = stats ? Stats.default(grid, shost, m, linop, transform; gas) : Output.nostats
    out = Output.MemoryOutput(0, flength, 3, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, out;
             zmax=flength, boundary, init_dz=flength/20, rtol=1e-8)
    out, transform
end


#= A pressure gradient. `LinearOps.make_linop` gives a z-dependent operator closure and
   `constβ` is left at its default of false, so this is the only case which exercises
   `NormModeAvg`'s `HostMirror` branch and `RK45.make_prop!`'s host-buffer branch -- the
   two pieces of GPU_PLAN.md section 4.5 layer 1 which upload from the host on every stage
   until gpu/23 tabulates them. Fixed steps, so that the only difference between the runs
   is the arithmetic. =#
function gradientcase(spec; gas=:Ar, pin=1.0, pout=0.0, flength=1e-2, λ0=800e-9,
                      precision=nothing)
    grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15)
    coren, densityfun = Capillary.gradient(gas, flength, pin, pout)
    m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun! = LinearOps.make_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=1e-6)
    Eω, transform, FT = Luna.setup(grid, densityfun, resp, inputs, βfun!, aeff;
                                   device=spec, precision)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    dz = flength/20
    Luna.run(Eω, grid, linop, transform, FT, out;
             zmax=flength, boundary=:none, init_dz=dz, min_dz=dz, max_dz=dz)
    out, transform
end

@testset "constβ is checked, not trusted" begin
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    coren, densityfun = Capillary.gradient(:Ar, 1e-2, 1.0, 0.0)
    m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    _, βfun! = LinearOps.make_linop(grid, m, 800e-9)
    #= `constβ=true` with a z-dependent βfun! would freeze β at z = 0 and give a silently
       wrong propagation, so setup evaluates it twice and refuses. =#
    @test_throws ErrorException NonlinearRHS.norm_mode_average(grid, βfun!, aeff;
                                                               constβ=true)
    # The constant operator's βfun! passes, which is what Interface relies on
    mc = Capillary.MarcatiliMode(75e-6, :Ar, 1.0, loss=false)
    _, βc!, _, _ = LinearOps.make_const_linop(grid, mc, 800e-9)
    @test NonlinearRHS.norm_mode_average(grid, βc!, z -> Modes.Aeff(mc, z=z);
                                         constβ=true) isa NonlinearRHS.NormModeAvg
end

# --- The device path proper, skipped without JLArrays -------------------------
have_jlarrays = try
    @eval import JLArrays
    true
catch
    false
end

if !have_jlarrays
    @warn "JLArrays is not available; the JLArray device tests are skipped. Run through "*
          "`Pkg.test()` or add JLArrays to the environment."
else

#= AbstractFFTs plan shims for JLArray. JLArrays provides no FFTs, so supply the minimum
   the device path needs, backed by FFTW on a host copy. Deliberately no `ldiv!` method:
   Luna must reach the transform through the explicit inverse plan, which is the dispatch
   a real device plan takes. `pinv` is part of the `Plan` contract (`inv` memoises into
   it), and `plan_inv` returns a `ScaledPlan` so that `Utils.iscale` finds the 1/N. =#
mutable struct JLPlan{T, N, P} <: AbstractFFTs.Plan{T}
    hp::P
    sz::NTuple{N, Int}
    dims::Any
    pinv::AbstractFFTs.ScaledPlan
    JLPlan{T, N, P}(hp, sz, dims) where {T, N, P} = new{T, N, P}(hp, sz, dims)
end
JLPlan(hp, sz::NTuple{N, Int}, dims, T=ComplexF64) where {N} =
    JLPlan{T, N, typeof(hp)}(hp, sz, dims)
Base.size(p::JLPlan) = p.sz
Base.eltype(::JLPlan{T}) where {T} = T
AbstractFFTs.plan_fft(x::JLArrays.JLArray{ComplexF64}, dims) =
    JLPlan(FFTW.plan_fft(Array(x), dims), size(x), dims)
AbstractFFTs.plan_rfft(x::JLArrays.JLArray{Float64}, dims) =
    JLPlan(FFTW.plan_rfft(Array(x), dims), size(x), dims, Float64)
AbstractFFTs.plan_inv(p::JLPlan) =
    AbstractFFTs.ScaledPlan(JLPlan(inv(p.hp).p, p.sz, p.dims, ComplexF64),
                            AbstractFFTs.normalization(Float64, p.sz, p.dims))
Base.:*(p::JLPlan, x::JLArrays.JLArray) = JLArrays.JLArray(p.hp * Array(x))
LinearAlgebra.mul!(y::JLArrays.JLArray, p::JLPlan, x::JLArrays.JLArray) =
    (copyto!(y, p.hp * Array(x)); y)

const JLArray = JLArrays.JLArray
const JLSpec = DeviceSpec(JLArray, Float64)

#= `allowscalar(false)` writes a task-local key and a process-global default. The
   task-local key is restored at the end of the file, so that one test file does not
   change what the rest of `Pkg.test()` sees; the process-global default is
   `ScalarDisallowed` in a non-interactive run anyway, which is what this sets it to. =#
const SCALAR_WAS = get(task_local_storage(), :ScalarIndexing, nothing)
GPUArraysCore.allowscalar(false)

@testset "JLArray basics" begin
    @test Utils.backend(JLArray(zeros(4))) isa Utils.DeviceBackend
    @test Utils.backend(view(JLArray(zeros(4)), 1:2)) isa Utils.DeviceBackend

    x = Luna.alloc(JLSpec, ComplexF64, (8,))
    @test x isa JLArray{ComplexF64, 1}
    @test all(iszero, Array(x))

    v = rand(5)
    dv = Luna.todevice(JLSpec, v)
    @test dv isa JLArray{Float64, 1}
    @test Array(dv) == v
    @test Luna.tohost(dv) == v

    # A BitArray becomes a Bool device array, since a mask has to be broadcastable there
    b = BitVector([true, false, true])
    db = Luna.todevice(JLSpec, b)
    @test db isa JLArray{Bool, 1}
    @test Array(db) == Array(b)

    @test Luna.assert_resident(JLSpec, x, dv, db) === nothing
    @test_throws ErrorException Luna.assert_resident(JLSpec, v)

    # upload_like moves and converts in one go
    y = Luna.alloc(JLSpec, ComplexF64, (6,))
    l = rand(ComplexF64, 6)
    dl = Luna.upload_like(y, l)
    @test dl isa JLArray{ComplexF64, 1}
    @test Array(dl) == l
    @test Luna.mask_like(y, b) isa JLArray{Bool, 1}
end

@testset "RK45 kernels on JLArray" begin
    n = 64
    ks = ntuple(_ -> rand(ComplexF64, n), 7)
    y = rand(ComplexF64, n)
    yn = similar(y)
    dks = map(JLArray, ks)
    dy = JLArray(y)
    dyn = similar(dy)
    dt = 1.3e-4

    for (ii, b) in ((1, RK45.B[1]), (3, RK45.B[3]), (6, RK45.B[6]),
                    (7, RK45.b5), (7, RK45.b4))
        RK45.combine!(yn, y, ks, dt, b, ii)
        RK45.combine!(dyn, dy, dks, dt, b, ii)
        @test Array(dyn) == yn # no scalar indexing, and the same arithmetic
    end

    yerr = similar(y)
    dyerr = similar(dy)
    RK45.errorestimate!(yerr, ks, dt)
    RK45.errorestimate!(dyerr, dks, dt)
    @test Array(dyerr) == yerr

    #= The norms are tree reductions on a device, so they agree only to rounding; the
       error metric is a scalar the controller compares to 1, so that is plenty. =#
    for nrm in (RK45.weaknorm, RK45.maxnorm, RK45.maxnorm_ratio, RK45.normnorm)
        h = nrm(yerr, y, yn, 1e-6, 1e-10)
        d = nrm(dyerr, dy, dyn, 1e-6, 1e-10)
        @test isapprox(h, d, rtol=1e-12)
    end
    # The max-based norms have no summation order to differ in, so they are exact
    @test RK45.maxnorm(dyerr, dy, dyn, 1e-6, 1e-10) ==
          RK45.maxnorm(yerr, y, yn, 1e-6, 1e-10)
end

@testset "mode-averaged Kerr on JLArray" begin
    for GT in (Grid.RealGrid, Grid.EnvGrid)
        href, htr = kerrcase(GT, HostSpec())
        dref, dtr = kerrcase(GT, JLSpec)

        # The transform really is resident on the device
        @test dtr.Eto isa JLArray
        @test dtr.Eωo isa JLArray
        @test dtr.gv.ω isa JLArray
        @test dtr.gv.sidx isa JLArray{Bool}
        @test dtr.norm!.pre isa JLArray
        # ... and the host one still aliases the grid's own vectors
        @test htr.gv.ω === htr.grid.ω
        @test htr.gv.sidx === htr.grid.sidx

        @test size(dref["Eω"]) == size(href["Eω"])
        @test dref["z"] ≈ href["z"]
        for idx in axes(href["Eω"], 2)
            h = href["Eω"][:, idx]
            d = dref["Eω"][:, idx]
            @test maximum(abs, d .- h)/maximum(abs, h) < 1e-10
        end
    end
end

@testset "pressure gradient on JLArray" begin
    href, htr = gradientcase(HostSpec())
    dref, dtr = gradientcase(JLSpec)

    # The z-dependent branch: β is mirrored rather than folded into `pre`
    @test !isnothing(dtr.norm!.β)
    @test dtr.norm!.β.host isa Vector{Float64}
    @test dtr.norm!.β.dev isa JLArray{Float64, 1}
    @test dtr.norm!.β.stage isa Vector{Float64}
    # ... and on the host the mirror is the host buffer itself, so `upload!` does nothing
    @test htr.norm!.β.dev === htr.norm!.β.host
    @test isnothing(htr.norm!.β.stage)

    @test dref["z"] ≈ href["z"]
    for idx in axes(href["Eω"], 2)
        h = href["Eω"][:, idx]
        d = dref["Eω"][:, idx]
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-10
    end
end

#= gpu/11's exit condition for this test file: RateAbsorber and the default statistics,
   both broadcasts and reductions now (Boundaries.jl, Stats.jl is untouched but the field
   it is called with is host, unscaled data whatever device the state lives on), run on a
   device end to end through `Luna.run`, with no explicit host-copy wrapper from the
   caller -- `ScaledOutput` is automatic. =#
@testset "boundaries and default statistics on JLArray" begin
    for GT in (Grid.RealGrid, Grid.EnvGrid)
        href, htr = kerrcase(GT, HostSpec(); boundary=:rate, stats=true)
        dref, dtr = kerrcase(GT, JLSpec; boundary=:rate, stats=true)

        @test dref["z"] ≈ href["z"]
        for idx in axes(href["Eω"], 2)
            h = href["Eω"][:, idx]
            d = dref["Eω"][:, idx]
            @test maximum(abs, d .- h)/maximum(abs, h) < 1e-10
        end
        # The statistics agree too: they were computed from a host copy on both paths
        @test dref["stats"]["energy"] ≈ href["stats"]["energy"] rtol=1e-8
        @test length(dref["stats"]["z"]) == length(href["stats"]["z"])
    end
end

#= Review round 1, finding 3: `stats_period` has to actually skip the device-to-host copy
   on a step whose statistics `PeriodicStats` is going to discard, not only skip the
   statistics arithmetic. `ScaledOutput.ybuf` is the buffer that copy lands in, and it is
   never touched except by `_tohost_unscale!`, so a step which does not fire leaves it
   bit-for-bit as the previous fire left it -- constructing `ScaledOutput` directly (as
   `Luna.run` does internally) makes that directly observable, without needing to
   instrument `copyto!` or time anything. =#
@testset "stats_period skips the device-to-host copy" begin
    calls = Ref(0)
    statsfun(y, t, dt) = (calls[] += 1; Dict("s" => sum(abs2, y)))
    periodic = Output.maybe_periodic(statsfun, 3) # fires on the 1st, 4th, 7th call
    out = Output.MemoryOutput(0, 1.0, 2, periodic)
    y1 = JLArray(ComplexF32[1, 2, 3])
    so = Luna.ScaledOutput(out, y1, 1.0)

    so(y1, 0.0, 0.1, _ -> y1) # call 1: fires
    @test calls[] == 1
    snap = copy(Array(so.ybuf))
    @test snap == ComplexF32[1, 2, 3]

    # calls 2 and 3 do not fire: ybuf must be untouched, whatever y they are given
    so(JLArray(ComplexF32[9, 9, 9]), 0.1, 0.1, _ -> y1)
    @test calls[] == 1
    @test Array(so.ybuf) == snap
    so(JLArray(ComplexF32[8, 8, 8]), 0.2, 0.1, _ -> y1)
    @test calls[] == 1
    @test Array(so.ybuf) == snap

    # call 4 fires again: ybuf now reflects it
    y4 = JLArray(ComplexF32[7, 7, 7])
    so(y4, 0.3, 0.1, _ -> y4)
    @test calls[] == 2
    @test Array(so.ybuf) == ComplexF32[7, 7, 7]

    # the one-time host-statistics warning still fires exactly once, on the first copy
    calls2 = Ref(0)
    statsfun2(y, t, dt) = (calls2[] += 1; Dict("s" => 0.0))
    out2 = Output.MemoryOutput(0, 1.0, 2, Output.maybe_periodic(statsfun2, 5))
    so2 = Luna.ScaledOutput(out2, y1, 1.0)
    @test_logs (:warn, r"Per-step statistics run on the host") match_mode=:any begin
        so2(y1, 0.0, 0.1, _ -> y1) # call 1: fires, first device copy -> warns
    end
    @test_logs begin # calls 2-4 do not fire and do not copy, so no repeat warning either
        for t in (0.1, 0.2, 0.3)
            so2(JLArray(ComplexF32[1, 1, 1]), t, 0.1, _ -> y1)
        end
    end
end

@testset "a device run refuses host-only machinery" begin
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, :He, 1.0, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    dens = z -> PhysData.density(:He, 1.0)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, 800e-9)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6)

    # A normalisation built for the host cannot be used for a device run
    hostnorm = NonlinearRHS.norm_mode_average(grid, βfun!, aeff)
    @test_throws ErrorException Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                           norm! = hostnorm, device=JLSpec)
    #= Nor can one built for the right array type but the wrong units: that would be
       wrong by the factor Pref = ε₀ with nothing to say so. =#
    s32 = DeviceSpec(Array, Float32)
    wrongunits = NonlinearRHS.norm_mode_average(grid, βfun!, aeff; spec=s32)
    @test_throws ErrorException NonlinearRHS.check_norm(
        wrongunits, s32, Luna.UnitScaling(1024.0, PhysData.ε_0))
    rightunits = NonlinearRHS.norm_mode_average(
        grid, βfun!, aeff; spec=s32, scaling=Luna.UnitScaling(1024.0, PhysData.ε_0))
    @test NonlinearRHS.check_norm(
        rightunits, s32, Luna.UnitScaling(1024.0, PhysData.ε_0)) === nothing
    # A response with no `rescale` method is refused in a scaled run
    @test_throws ErrorException Nonlinear.rescale(
        (out, E, ρ) -> nothing, DeviceSpec(Array, Float32), UNIT_SCALING)
end

isnothing(SCALAR_WAS) ? delete!(task_local_storage(), :ScalarIndexing) :
                        task_local_storage(:ScalarIndexing, SCALAR_WAS)

end # have_jlarrays

#= The reduced-precision path on the CPU. This is the scaling layer under a magnifying
   glass: Float32 on the host keeps subnormals, so an unscaled run would *work* here and
   fail only on Metal. The check is therefore on the numbers, not on whether it runs. =#
@testset "Float32 on the CPU" begin
    #= Helium at 0.3 bar is the dynamic-range case: ρ ε₀ γ₃ is 3.3e-39 there, below the
       smallest Float32 subnormal (1.2e-38 normal, 1.4e-45 subnormal), and Metal flushes
       subnormals to zero. With the scaling it is ~1e-34, comfortably normal. =#
    ρ = PhysData.density(:He, 0.3)
    raw = ρ*PhysData.ε_0*PhysData.γ3_gas(:He)
    @test raw < floatmin(Float32) # unscaled, this coefficient does not exist in Float32

    href, _ = kerrcase(Grid.RealGrid, HostSpec(); pres=0.3)
    f32, tr32 = kerrcase(Grid.RealGrid, DeviceSpec(Array, Float32); pres=0.3)

    #= The state and every buffer are single precision. =#
    @test tr32.Eto isa Vector{Float32}
    @test tr32.Eωo isa Vector{ComplexF32}
    @test tr32.gv.ω isa Vector{Float32}
    @test tr32.scaling.Eref > 1
    @test log2(tr32.scaling.Eref) == round(log2(tr32.scaling.Eref))
    @test tr32.scaling.Pref == PhysData.ε_0
    # The scaled coefficient is a normal Float32 with room to spare
    @test floatmin(Float32) < abs(tr32.resp[1].γ3*ρ*PhysData.ε_0) < floatmax(Float32)

    #= `ScaledOutput` unscales on the way into the output and `Output.MemoryOutput` now
       allocates with `eltype(y)`, so the saved field is `ComplexF32`, in physical units,
       directly comparable with the `Float64` run -- no manual `* Eref` here any more.
       The tolerance is measured, not aspirational: see PR_10-device-model.md. =#
    @test eltype(f32["Eω"]) === ComplexF32
    for idx in axes(href["Eω"], 2)
        h = href["Eω"][:, idx]
        d = ComplexF64.(f32["Eω"][:, idx])
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-5
    end
end

