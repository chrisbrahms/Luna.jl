#= Tests for Luna's device (precision and array type) model.

   Luna's per-step code is made of broadcasts, planned FFTs applied with `mul!` and
   reductions, so it can be exercised without GPU hardware:

   * the `Utils.backend` trait is purely type-level, so a dummy `AbstractGPUArray`
     subtype is enough to test its dispatch;
   * `JLArrays.JLArray` is a host-backed `AbstractGPUArray` which takes the same generic
     GPUArrays code path a real device does and, with `allowscalar(false)`, enforces the
     same no-scalar-indexing contract. It is a test-only dependency, so these tests are
     skipped if it cannot be loaded.

   What this file cannot catch, and what `test/test_metal.jl` must: JLArrays interprets
   its kernels on the host, so a stray `Float64` in a struct field or a mixed host/device
   broadcast passes here and would still fail on real hardware. =#

import Test: @test, @testset, @test_throws, @inferred
import Luna
import Luna: Utils, Output, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear,
             NonlinearRHS, PhysData, RK45
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

    # A plan which is already normalised reports a factor of 1
    @test Utils.iscale(pr) == 1
    @test Utils.iplan(pr) === pr

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

#= An output which copies whatever the solver hands it down to the host. The real one
   (`ScaledOutput`, which also unscales) is gpu/11's; this is the minimum this branch
   needs in order to save a device propagation. =#
struct ToHost{O}
    o::O
end
(h::ToHost)(y, t, dt, yfun) = h.o(Array(y), t, dt, ti -> Array(yfun(ti)))
(h::ToHost)(args...; kwargs...) = h.o(args...; kwargs...)
Base.getindex(h::ToHost, k) = h.o[k]

#= One small mode-averaged Kerr propagation, built once and run on whichever spec is
   asked for. `boundary=:none` because the absorbers are host code until gpu/11. =#
function kerrcase(GT, spec; gas=:He, pres=1.0, energy=1e-6, flength=1e-2, λ0=800e-9,
                  precision=nothing)
    grid = GT === Grid.RealGrid ?
        Grid.RealGrid(λ0, (300e-9, 2000e-9), 400e-15) :
        Grid.EnvGrid(λ0, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, gas, pres, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    dens = z -> PhysData.density(gas, pres)
    resp = GT === Grid.RealGrid ?
        (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),) :
        (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=energy)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=spec, precision)
    out = Output.MemoryOutput(0, flength, 3, Output.nostats)
    output = Utils.isdevice(Eω) ? ToHost(out) : out
    Luna.run(Eω, grid, linop, transform, FT, output;
             zmax=flength, boundary=:none, init_dz=flength/20, rtol=1e-8)
    out, transform
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

@testset "a device run refuses host-only machinery" begin
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, :He, 1.0, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    dens = z -> PhysData.density(:He, 1.0)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, 800e-9)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6)
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, device=JLSpec)
    out = ToHost(Output.MemoryOutput(0, 1e-3, 3, Output.nostats))
    # The absorbers are host scalar code, so they are refused rather than run slowly
    @test_throws ErrorException Luna.run(Eω, grid, linop, transform, FT, out;
                                         zmax=1e-3, boundary=:rate)

    # A normalisation built for the host cannot be used for a device run
    hostnorm = NonlinearRHS.norm_mode_average(grid, βfun!, aeff)
    @test_throws ErrorException Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                           norm! = hostnorm, device=JLSpec)
    # A response with no `rescale` method is refused in a scaled run
    @test_throws ErrorException Nonlinear.rescale(
        (out, E, ρ) -> nothing, DeviceSpec(Array, Float32), UNIT_SCALING)
end

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

    #= The state and every buffer are single precision. `Output.MemoryOutput` still
       allocates ComplexF64 and widens on the way in (Output.jl:44); making the output
       carry `eltype(y)` is gpu/11's, together with the unscaling wrapper. =#
    @test tr32.Eto isa Vector{Float32}
    @test tr32.Eωo isa Vector{ComplexF32}
    @test tr32.gv.ω isa Vector{Float32}
    @test tr32.scaling.Eref > 1
    @test log2(tr32.scaling.Eref) == round(log2(tr32.scaling.Eref))
    @test tr32.scaling.Pref == PhysData.ε_0
    # The scaled coefficient is a normal Float32 with room to spare
    @test floatmin(Float32) < abs(tr32.resp[1].γ3*ρ*PhysData.ε_0) < floatmax(Float32)

    #= The state is scaled, and unscaling is gpu/11's job (the output wrapper), so
       multiply by Eref here. The tolerance is measured, not aspirational: see
       PR_10-device-model.md. =#
    Eref = tr32.scaling.Eref
    for idx in axes(href["Eω"], 2)
        h = href["Eω"][:, idx]
        d = ComplexF64.(f32["Eω"][:, idx]) .* Eref
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-5
    end
end

