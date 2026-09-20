#= Device (GPU) and precision support.

   Luna has no GPU dependency. Everything here is generic over an array type and a real
   element type: the per-step code elsewhere in Luna consists of broadcasts, planned FFTs
   applied with `mul!`, and reductions, so it compiles for any array following the
   GPUArrays interface (`MtlArray`, `CuArray`, `JLArray`, ...). The operations which
   genuinely cannot be written generically -- is there a device, block until it is idle,
   return pooled memory, how much is free -- are function-valued hooks installed by a
   package extension (`ext/LunaMetalExt.jl`, `ext/LunaCUDAExt.jl`).

   These definitions live directly in the `Luna` module (the file is `include`d without a
   `module` wrapper) so that extensions can extend them as `Luna.register_device!` etc.
   It is included immediately after `Utils`, which provides the backend trait. =#

import Adapt
import .Utils: Backend, CPUBackend, DeviceBackend, backend, isdevice

#=================================================#
#===============  DEVICE SPEC  ===================#
#=================================================#

"""
    DeviceSpec{A, T}()

Where a propagation runs and in what precision. `A` is the array type (`Array`,
`MtlArray`, `CuArray`, `JLArray`) and `T` the *real* element type (`Float64` or
`Float32`); the state itself is `Complex{T}`.

A `DeviceSpec` is a singleton type, so it carries no data and everything derived from it
(`alloc`, `todevice`, the transforms' buffer types) is resolved at compile time.

Construct with `DeviceSpec(Array, Float64)` or obtain the one a run will use from
[`Luna.device`](@ref).
"""
struct DeviceSpec{A, T} end

DeviceSpec(::Type{A}, ::Type{T}) where {A, T} = DeviceSpec{A, T}()

"The array type of a [`DeviceSpec`](@ref)."
arraytype(::DeviceSpec{A, T}) where {A, T} = A

"The real element type of a [`DeviceSpec`](@ref); the state is `Complex{realtype(spec)}`."
realtype(::DeviceSpec{A, T}) where {A, T} = T

"`true` if the spec's array type is not `Array`, i.e. the state lives on a GPU."
isdevicespec(::DeviceSpec{A, T}) where {A, T} = !(A === Array)

"Same device, different precision."
withprecision(::DeviceSpec{A, T}, ::Type{S}) where {A, T, S} = DeviceSpec{A, S}()
withprecision(s::DeviceSpec, ::Nothing) = s

Base.show(io::IO, ::DeviceSpec{A, T}) where {A, T} = print(io, "DeviceSpec($A, $T)")

"The default spec: host arrays in double precision, which is what Luna has always done."
const HostSpec = DeviceSpec{Array, Float64}

#=================================================#
#==========  DEVICE REGISTRY AND HOOKS  ==========#
#=================================================#

"""
    DeviceHooks

The vendor operations for one registered device, held as function-valued fields.

They are fields rather than methods because an extension cannot *overwrite* an
identically signed stub defined in the main package -- Julia rejects that during
precompilation, which leaves the extension uncompilable. They are called through
`Base.invokelatest` (the GPU package may have been loaded after the calling frame was
compiled) and only at setup or teardown, never per step, so the dynamic dispatch costs
nothing.
"""
struct DeviceHooks{S<:DeviceSpec}
    spec::S
    functional::Any
    synchronize::Any
    reclaim::Any
    memory_status::Any
end

"Registered devices, keyed by name (`:metal`, `:cuda`). Filled by the extensions."
const DEVICES = Dict{Symbol, DeviceHooks}()

"""
    register_device!(name, spec; functional, synchronize, reclaim, memory_status)

Register a GPU backend under `name` (`:metal`, `:cuda`). Called from the `__init__` of a
package extension; user code never calls it.

`spec` is the [`DeviceSpec`](@ref) this backend runs, and the four keyword arguments are
the vendor operations (see [`DeviceHooks`](@ref)), each a zero-argument callable:
`functional()` says whether this process can actually use the device, `synchronize()`
blocks until queued work has finished, `reclaim()` returns pooled memory to the driver,
and `memory_status()` returns `(free, total)` in bytes or `nothing`.

Registering also sets `Luna.settings["device"] = :auto` **if and only if the key is
absent**, so that `using Metal` before a run switches Luna to the GPU while an explicit
`Luna.set_device(:cpu)` is never overridden.
"""
function register_device!(name::Symbol, spec::DeviceSpec;
                          functional = () -> true,
                          synchronize = () -> nothing,
                          reclaim = () -> nothing,
                          memory_status = () -> nothing)
    DEVICES[name] = DeviceHooks(spec, functional, synchronize, reclaim, memory_status)
    haskey(settings, "device") || (settings["device"] = :auto)
    name
end

_functional(h::DeviceHooks) = Base.invokelatest(h.functional)::Bool

"The registered device names, in a fixed order, so that `:auto` is deterministic."
devicenames() = sort!(collect(keys(DEVICES)))

"""
    device_functional(name) -> Bool

Whether the GPU registered under `name` (see [`register_device!`](@ref)) is registered at
all and usable in this process.
"""
device_functional(name::Symbol) = haskey(DEVICES, name) && _functional(DEVICES[name])

"""
    Luna.device() -> DeviceSpec

The [`DeviceSpec`](@ref) the next run will use, resolved from `Luna.settings["device"]`.

The setting may be

- absent (the initial state): the CPU in double precision;
- `:cpu`: the same, explicitly, and not overridden when a GPU package is loaded later;
- `:auto`: the registered GPU if one is registered and functional, else the CPU;
- `:metal` or `:cuda`: that GPU, erroring if it is not registered or not functional;
- a `DeviceSpec`: used as it stands.

See [`Luna.set_device`](@ref).
"""
device() = resolve_device(device_request())

"""
    Luna.device_request()

`Luna.settings["device"]` as the user set it, or `:cpu` if the key is absent -- the
unresolved form of [`Luna.device`](@ref). This is what `Luna.setup` takes as the default
of its `device` keyword, so that [`log_device`](@ref) can distinguish `:auto` which found
no GPU from an explicit `:cpu`.
"""
device_request() = get(settings, "device", :cpu)

"Resolve a device specification (see [`Luna.device`](@ref)) to a [`DeviceSpec`](@ref)."
resolve_device(s::DeviceSpec) = s

function resolve_device(s::Symbol)
    s === :cpu && return HostSpec()
    if s === :auto
        for name in devicenames()
            _functional(DEVICES[name]) && return DEVICES[name].spec
        end
        return HostSpec()
    end
    haskey(DEVICES, s) || error(
        "no device is registered under `:$s`. Load the GPU package first "*
        "(`using Metal` or `using CUDA`), which registers it through Luna's package "*
        "extension. Registered: $(isempty(DEVICES) ? "none" : join(devicenames(), ", ")).")
    _functional(DEVICES[s]) || error(
        "`:$s` is registered but not functional in this process on $(gethostname()).")
    DEVICES[s].spec
end

resolve_device(x) = error(
    "`$x` is not a device specification: expected `:cpu`, `:auto`, `:metal`, `:cuda` "*
    "or a `Luna.DeviceSpec`.")

"""
    Luna.set_device(x)

Set which device Luna runs on. `x` is `:cpu`, `:auto`, `:metal`, `:cuda` or a
[`DeviceSpec`](@ref); see [`Luna.device`](@ref) for what each means. The value is
validated immediately, so an unavailable device is reported here rather than at setup.

Like the FFTW settings this changes `Luna.settings` in the calling process only, so a
`Scans` worker needs `@everywhere using Metal` (or `@everywhere Luna.set_device(...)`).
"""
function set_device(x)
    resolve_device(x) # validate now, not at setup
    settings["device"] = x
end

"""
    device_synchronize([spec])

Block until all queued work on the device has finished. A no-op on the CPU. Only needed
for timing: Luna's own code paths are synchronised by the data dependencies between
kernels.
"""
device_synchronize(spec = device()) = _hookcall(spec, :synchronize)

"""
    device_reclaim([spec])

Return cached device memory to the driver. A no-op on the CPU. GPU array libraries keep a
memory pool which garbage collection alone does not release, so call this between
independent propagations in a long-lived process.
"""
device_reclaim(spec = device()) = _hookcall(spec, :reclaim)

"""
    device_memory_status([spec]) -> (free, total) in bytes, or `nothing`

Free and total memory on the device, or `nothing` on the CPU.
"""
device_memory_status(spec = device()) = _hookcall(spec, :memory_status)

function _hookcall(spec, field::Symbol)
    for name in devicenames()
        h = DEVICES[name]
        h.spec === spec && return Base.invokelatest(getfield(h, field))
    end
    nothing
end

#=================================================#
#==========  ALLOCATION AND TRANSFER  ============#
#=================================================#

"""
    log_device(spec, request)

Log, once per `Luna.setup`, which device and precision the run will use, and warn when
`:auto` was asked for but no GPU package is loaded in this process (which is what a
`Scans` worker that only did `using Luna` sees).
"""
function log_device(spec::DeviceSpec, request)
    if request === :auto && !isdevicespec(spec)
        Logging.@info(
            "`:auto` requested but no GPU package is loaded in this process; running on "*
            "the CPU. `Scans` workers load Luna on their own, so a scan needs "*
            "`@everywhere using Metal` (or CUDA).")
    end
    Logging.@info("Propagating on $(arraytype(spec)) in $(realtype(spec)) precision.")
    nothing
end

"""
    alloc(spec, T, dims)

A zero-filled array of element type `T`, shape `dims` and the array type of `spec`.
Equivalent to `zeros(T, dims)` on the host.

Field-sized buffers are allocated this way rather than built on the host and copied, so
that a device run never materialises a host array it does not need.
"""
alloc(spec::DeviceSpec, ::Type{T}, dims::Tuple) where {T} =
    fill!(arraytype(spec){T}(undef, dims...), zero(T))
alloc(spec::DeviceSpec, ::Type{T}, dims::Integer...) where {T} = alloc(spec, T, dims)

"""
    todevice(spec, x)

Move the host array `x` to `spec`'s array type, converting `Float64` to `realtype(spec)`
and `ComplexF64` to `Complex{realtype(spec)}`. Boolean masks (including `BitArray`s) keep
their element type.

On `DeviceSpec(Array, Float64)` this returns `x` itself: the host double-precision path
allocates nothing and the grid mirrors alias the grid's own vectors.
"""
todevice(spec::DeviceSpec, x::AbstractArray) =
    _adapt(arraytype(spec), _convertprec(realtype(spec), x))

todevice(spec::DeviceSpec, x::AbstractArray{Bool}) =
    isdevicespec(spec) ? _adapt(arraytype(spec), Array{Bool}(x)) : x

todevice(::DeviceSpec, ::Nothing) = nothing

#= The array type's own constructor rather than `Adapt.adapt`: not every GPU array
   package defines an `Adapt.adapt_storage` rule for its bare array type (JLArrays does
   not), and `Adapt.adapt` silently returns the host array when none exists, which would
   put a host array into a device kernel. The constructor is the documented way to move
   an array and every backend has it. The `Array` method still has to copy a device array
   down, so that `todevice(HostSpec(), x)` means what its name says whatever `x` is. =#
_adapt(::Type{Array}, x::Array) = x
_adapt(::Type{Array}, x::AbstractArray) = isdevice(x) ? Array(x) : x
_adapt(::Type{A}, x::AbstractArray) where {A} = A(x)

_convertprec(::Type{T}, x::AbstractArray{T}) where {T<:AbstractFloat} = x
_convertprec(::Type{T}, x::AbstractArray{<:AbstractFloat}) where {T} = convert(Array{T}, x)
_convertprec(::Type{T}, x::AbstractArray{Complex{T}}) where {T<:AbstractFloat} = x
_convertprec(::Type{T}, x::AbstractArray{<:Complex}) where {T} = convert(Array{Complex{T}}, x)

"""
    tohost(x)

A host `Array` holding the same data as `x`. Returns `x` unchanged if it is already one.
"""
tohost(x::Array) = x
tohost(x::AbstractArray) = Array(x)

"""
    scalar(x, s)

Convert the `Float64` scalar `s` to `real(eltype(x))`, so that a broadcast over `x` never
sees a `Float64`. Metal's kernel compiler rejects any `double` which survives
optimisation, and its kernel adaptor does not convert scalars.
"""
scalar(x::AbstractArray, s) = convert(real(eltype(x)), s)

"""
    upload_like(y, x)

`x` in the array type and precision of `y`: a `ComplexF64` host array becomes a
`Complex{real(eltype(y))}` array of `y`'s array type. Returns `x` itself when `y` is a
host array of the same precision, so the default CPU path is untouched.

Used for the linear operator, which is always built on the host in `Float64` (through
`Boundaries.setup`) and uploaded once by [`Luna.run`](@ref).
"""
function upload_like(y::AbstractArray, x::AbstractArray)
    ET = _matcheltype(eltype(y), eltype(x))
    (!isdevice(y) && eltype(x) === ET) && return x
    out = similar(y, ET, size(x))
    copyto!(out, eltype(x) === ET ? x : convert(Array{ET}, x))
    out
end

upload_like(y::AbstractArray, x) = x # a closure operator: nothing to upload

"""
    mask_like(y, m)

The boolean mask `m` on the array type of `y`. Returns `m` itself when `y` is a host
array, so a `BitArray` stays one and the host path allocates nothing.
"""
mask_like(y::AbstractArray, m::AbstractArray{Bool}) =
    isdevice(y) ? copyto!(similar(y, Bool, size(m)), Array{Bool}(m)) : m

_matcheltype(::Type{Complex{T}}, ::Type{<:Complex}) where {T} = Complex{T}
_matcheltype(::Type{Complex{T}}, ::Type{<:Real}) where {T} = T
_matcheltype(::Type{T}, ::Type{<:Complex}) where {T<:Real} = Complex{T}
_matcheltype(::Type{T}, ::Type{<:Real}) where {T<:Real} = T

#=================================================#
#===========  RESIDENCY ASSERTIONS  ==============#
#=================================================#

"""
    assert_resident(spec, arrays...)

Check that every array lives on `spec`'s array type with element type `realtype(spec)`,
`Complex{realtype(spec)}` or `Bool`, and error otherwise. `nothing` entries are skipped.

Every transform constructor calls this on its buffers, mirrors and coefficient arrays.
A mixed host/device broadcast is not caught by `JLArrays` and is a silent, very slow
fallback on some backends, so residency is checked structurally at construction rather
than left to the first kernel.
"""
function assert_resident(spec::DeviceSpec, xs...)
    for (i, x) in enumerate(xs)
        isnothing(x) && continue
        _resident(spec, x) || error(
            "argument $i of assert_resident is a $(typeof(x)), which does not live on "*
            "$(spec). Every buffer, grid mirror and coefficient array of a transform "*
            "must be allocated with `Luna.alloc` or moved with `Luna.todevice`.")
    end
    nothing
end

"""
    all_resident(spec, arrays...) -> Bool

The predicate form of [`assert_resident`](@ref): `true` if every array (skipping
`nothing`) lives on `spec`'s array type and precision. Used where a failure is reported
with a different message.
"""
all_resident(spec::DeviceSpec, xs...) =
    all(x -> isnothing(x) || _resident(spec, x), xs)

function _resident(spec::DeviceSpec, x::AbstractArray)
    T = realtype(spec)
    eltype(x) === T || eltype(x) === Complex{T} || eltype(x) === Bool || return false
    isdevicespec(spec) || return !isdevice(x)
    Base.typename(typeof(x)).wrapper === arraytype(spec)
end

_resident(::DeviceSpec, ::Any) = false

#=================================================#
#===============  UNIT SCALING  ==================#
#=================================================#

"""
    UnitScaling(Eref, Pref)

The units the propagating state and the nonlinear polarisation are expressed in:
the state is `e = E/Eref` and the polarisation buffer holds `p = P/(Pref*Eref)`.

SI-unit coefficients do not fit in `Float32`: the combined Kerr coefficient
`ρ ε₀ γ₃` is 3e-39 for helium at 0.3 bar, below the smallest `Float32` subnormal, and
Metal flushes subnormals to zero. Scaling moves every coefficient a response sees into
the middle of the exponent range; the powers of `Eref` and `Pref` are combined with the
physical constants on the host, in `Float64`, and only the result is converted to the
device precision.

`Eref == Pref == 1` (see [`unitscaling`](@ref)) for every `Float64` run, which makes the
scaled and unscaled arithmetic identical, so the default CPU path is unchanged.
"""
struct UnitScaling
    Eref::Float64
    Pref::Float64
end

"The identity scaling, used by every `Float64` run."
const UNIT_SCALING = UnitScaling(1.0, 1.0)

"`true` if this scaling is the identity, i.e. the state is in physical units."
isunity(s::UnitScaling) = (s.Eref == 1.0) && (s.Pref == 1.0)

"""
    polscale(scaling, degree)

The factor a nonlinear response of polynomial degree `degree` in the electric field
carries when the state is expressed in the units of `scaling`.

A response contributes ``P = c E^n`` in physical units. With ``E = E_{ref} e`` and
``p = P/(P_{ref} E_{ref})`` the same response contributes
``p = c E_{ref}^{n-1}/P_{ref}\\,e^n``, so the coefficient its kernel needs is
`c * polscale(scaling, n)`. The Kerr responses are cubic (`n = 3`); a ``χ^{(2)}``
response is quadratic (`n = 2`).

Exactly `1` for every `Float64` run, where `Eref == Pref == 1`, so multiplying by it
changes nothing on the default CPU path.

A response with no single polynomial degree (an ionisation rate, the plasma current)
scales its intermediates individually instead; see the developer guide.
"""
polscale(s::UnitScaling, degree::Integer) = s.Eref^(degree-1)/s.Pref

Base.show(io::IO, s::UnitScaling) = print(io, "UnitScaling(Eref=$(s.Eref), Pref=$(s.Pref))")

"""
    unitscaling(T, Et, Pref)

The [`UnitScaling`](@ref) for a run in real precision `T` with time-domain initial state
`Et`.

For `T === Float64` this is the identity: with a scaled state the per-step statistics,
the HDF5 cache, user `stepfun`s and the interpolant would all see non-physical data, and
`Float64` has the range to do without it.

For `Float32`, `Eref` is the peak of `abs.(Et)` rounded to a power of two -- so that
dividing by it is exact in every linear operation -- and `Pref` is the constant the
polarisation is measured in (`ε₀`, supplied by the caller). A zero or non-finite peak
falls back to 1.
"""
unitscaling(::Type{Float64}, Et, Pref) = UNIT_SCALING

#= `Et` may be given as a thunk, so that a Float64 run -- which never needs it -- does not
   pay for the inverse transform of the input field. =#
unitscaling(::Type{Float64}, Et::Function, Pref) = UNIT_SCALING
unitscaling(::Type{T}, Et::Function, Pref) where {T} = unitscaling(T, Et(), Pref)

function unitscaling(::Type{T}, Et, Pref) where {T}
    m = Float64(maximum(abs, Et))
    Eref = (isfinite(m) && m > 0) ? exp2(round(Int, log2(m))) : 1.0
    UnitScaling(Eref, Pref)
end

#=================================================#
#=============  GRID VECTOR MIRRORS  =============#
#=================================================#

"""
    GridVectors

The grid vectors which appear in per-step kernels alongside the propagating field: the
angular frequency axis, the spectral and temporal apodisation windows, and `sidx` as a
`Bool` mask.

`Grid.RealGrid`/`Grid.EnvGrid` hold host `Vector{Float64}`s, and a host vector cannot be
broadcast against a device array. Rather than making the grid types themselves
device-aware -- which would put device arrays into the metadata written to output files
-- each transform keeps a `GridVectors` mirror, adapted once to its own array type and
precision at construction. On `DeviceSpec(Array, Float64)` the mirror aliases the grid's
own vectors: no copy, no extra memory, and the CPU path is unchanged.
"""
struct GridVectors{V, M}
    ω::V
    ωwin::V
    twin::V
    towin::V
    sidx::M
end

"""
    gridvectors(grid, spec)

[`GridVectors`](@ref) for `grid`, on `spec`'s array type and precision.
"""
gridvectors(grid, spec::DeviceSpec) = GridVectors(
    todevice(spec, grid.ω), todevice(spec, grid.ωwin),
    todevice(spec, grid.twin), todevice(spec, grid.towin),
    todevice(spec, grid.sidx))

Adapt.adapt_structure(to, g::GridVectors) = GridVectors(
    Adapt.adapt(to, g.ω), Adapt.adapt(to, g.ωwin), Adapt.adapt(to, g.twin),
    Adapt.adapt(to, g.towin), Adapt.adapt(to, g.sidx))

#=================================================#
#===============  HOST MIRRORS  ==================#
#=================================================#

"""
    HostMirror(spec, n)

A length-`n` vector which is filled on the host in `Float64` and then made available on
`spec`'s array type and precision by [`upload!`](@ref).

This is the interim mechanism for quantities which are still evaluated by host scalar
code on every right-hand side -- the propagation constant `β` of a tapered or
pressure-graded waveguide, and a user-supplied `linop!` -- until `gpu/23` tabulates them.
On `DeviceSpec(Array, Float64)` the device array *is* the host buffer and `upload!` does
nothing, so the CPU path pays neither a copy nor a conversion.

# Fields
- `host`: the `Float64` buffer host code writes into
- `stage`: a host buffer in the device precision, or `nothing` when none is needed
- `dev`: the array the kernels broadcast against
"""
struct HostMirror{D, S}
    host::Vector{Float64}
    stage::S
    dev::D
end

function HostMirror(spec::DeviceSpec, n::Integer)
    host = zeros(Float64, n)
    T = realtype(spec)
    if arraytype(spec) === Array
        return HostMirror(host, nothing, T === Float64 ? host : zeros(T, n))
    end
    HostMirror(host, zeros(T, n), alloc(spec, T, (n,)))
end

"""
    upload!(m::HostMirror)

Make the contents of `m.host` visible to the kernels and return `m.dev`.
"""
function upload!(m::HostMirror{D, Nothing}) where {D}
    m.dev === m.host || copyto!(m.dev, m.host)
    m.dev
end

function upload!(m::HostMirror)
    copyto!(m.stage, m.host)
    copyto!(m.dev, m.stage)
    m.dev
end

#=================================================#
#===========  OUTPUT (DEVICE BOUNDARY)  ==========#
#=================================================#

#= `Output.jl` stays device-unaware -- it knows how to save an array and a dictionary,
   nothing about where the array lives or what units it is in. This section is what
   `Luna.run` wraps an output handler in when it needs to. It comes after Output.jl in
   Luna.jl's include list so that it can dispatch on `Output.MemoryOutput`/
   `Output.HDF5Output` and extend `Output.willsave`/`Output.check_cache`/`Output.hasdata`
   for the new wrapper type. =#

"""
    needs_host_y(o, t) -> Bool

Whether the output handler `o` is, *on the step about to be reported at coordinate `t`*,
going to inspect the per-step solution `y` (as opposed to only the interpolated values it
requests through `yfun`), and therefore needs it on the host, in physical units.

True whenever `o`'s statistics function is not `Output.nostats` and, if it is an
[`Output.PeriodicStats`](@ref), will actually fire this step ([`Output.willfire`](@ref)):
a statistic which cannot be evaluated where the state lives has to have it copied down on
every step whose statistics are really computed -- but not on a step `PeriodicStats` is
going to skip, which is the whole point of `stats_period` on a device. Conservatively
`true` for an output whose statistics function cannot be inspected this way.

Whether the copy is needed at all is the second question, answered by
[`stats_device_capable`](@ref): a statistics set which runs on the state as the stepper
holds it needs no copy on any step.
"""
needs_host_y(o, t) = true
needs_host_y(o::Output.MemoryOutput, t) = _stats_will_run(o.statsfun, t)
needs_host_y(o::Output.HDF5Output, t) = _stats_will_run(o.statsfun, t)

_stats_will_run(f, t) = f !== Output.nostats
_stats_will_run(p::Output.PeriodicStats, t) = Output.willfire(p, t)

"""
    stats_device_capable(x) -> Bool

Whether the statistics `x` stands for can be evaluated on the propagating state as the
stepper holds it: a device array, in the scaled units of [`UnitScaling`](@ref). `x` is
either an output handler or a statistics function.

The default sets [`Stats.default`](@ref Luna.Stats.default) builds answer this through
[`Stats.device_capable`](@ref Luna.Stats.device_capable); anything else -- a user-written
closure, a statistics function from somewhere else -- is `false`, which is what makes
[`ScaledOutput`](@ref) copy the state to the host every step. `Output.nostats` is
trivially capable: it reads nothing.
"""
stats_device_capable(x) = false
stats_device_capable(o::Output.MemoryOutput) = stats_device_capable(o.statsfun)
stats_device_capable(o::Output.HDF5Output) = stats_device_capable(o.statsfun)
stats_device_capable(p::Output.PeriodicStats) = stats_device_capable(p.f)
stats_device_capable(::typeof(Output.nostats)) = true

"""
    stats_host_list(x) -> Vector{String}

The names of the statistics in `x` which are not [`stats_device_capable`](@ref), for the
one-time warning [`ScaledOutput`](@ref) emits. `x` is either an output handler or a
statistics function. Extended by
[`Stats.host_statistics`](@ref Luna.Stats.host_statistics) for the sets `Stats.default`
builds; the fallback is a single entry naming the type of the statistics function.
"""
stats_host_list(x) = stats_device_capable(x) ? String[] : [string(nameof(typeof(x)))]
stats_host_list(o::Output.MemoryOutput) = stats_host_list(o.statsfun)
stats_host_list(o::Output.HDF5Output) = stats_host_list(o.statsfun)
stats_host_list(p::Output.PeriodicStats) = stats_host_list(p.f)

"""
    needs_host_cache(o) -> Bool

Whether `o` is an `Output.HDF5Output` with `cache=true`: it writes the raw per-step `y`
into the file's resume cache (not only the interpolated saves), and that write needs a
host array in physical units. Gated by [`Output.willsave`](@ref) in `ScaledOutput`, since
the cache is only written on a save step.
"""
needs_host_cache(o) = false
needs_host_cache(o::Output.HDF5Output) = o.cache

"""
    ScaledOutput(o, y)

Wrap the output handler `o` so that it receives host arrays in physical units, whatever
units and array type the propagating state `y` is in. Constructed by [`Luna.run`](@ref)
whenever the state is on a device or the run is scaled (`E_ref != 1`, i.e. every `Float32`
run, host or device); `Output.jl` itself never sees a device array or a scaled one.

Two reusable host buffers, in `y`'s element type (so a `Float32` run saves `Float32`):

- `ybuf` holds the unscaled, host copy of the per-step `y`, used for statistics and for an
  `HDF5Output`'s resume cache. Only filled when [`needs_host_y`](@ref) (re-evaluated every
  step, so `Output.PeriodicStats` skips the copy on a step it is not going to fire on) or
  ([`needs_host_cache`](@ref) and [`Output.willsave`](@ref)) says it is needed this step
  -- so a device run with `Output.nostats` and no HDF5 cache never pays for it, and an
  `HDF5Output` with caching pays only on a save step.
- `ibuf` holds the unscaled, host copy of a **saved** field. It is filled lazily, inside
  the closure `o` calls as `yfun`, so it costs nothing on a step which does not save
  (`o`'s own save condition decides whether to call it at all) and it is a different
  buffer from `ybuf` because the two can be needed in the same call with different
  contents (an `HDF5Output` reads its cache value `y` -- the step endpoint -- after
  writing possibly several interpolated saves at earlier `yfun(ts)`).

Multiplying by `E_ref` happens on the host, after the copy, in `y`'s own precision, and is
skipped entirely when `E_ref == 1` (every `Float64` run): the stepper's own arrays are
never modified, only these two buffers.
"""
mutable struct ScaledOutput{O, A<:AbstractArray}
    o::O
    Eref::Float64
    needcache::Bool     # `o` is an HDF5Output with a resumable cache
    ybuf::A
    ibuf::A
    warned::Base.RefValue{Bool}
end

function ScaledOutput(o, y::AbstractArray, Eref::Real)
    A = Array{eltype(y), ndims(y)}
    ScaledOutput{typeof(o), A}(o, Float64(Eref), needs_host_cache(o),
                               A(undef, size(y)), A(undef, size(y)), Ref(false))
end

#= Device-to-host copy plus, when the run is scaled, the unscaling multiply, into the
   reusable buffer `buf`. A plain host-to-host copy on the default CPU path (`Eref == 1`,
   `y` already an `Array`) still has to happen: `buf` must not alias the stepper's own
   arrays or the interpolant's, which are reused every call. =#
function _tohost_unscale!(buf::AbstractArray, y::AbstractArray, Eref::Float64)
    copyto!(buf, y)
    Eref == 1.0 && return buf
    buf .*= Eref
    buf
end

function _warn_host_stats!(so::ScaledOutput, y)
    so.warned[] && return nothing
    isdevice(y) || return nothing
    so.warned[] = true
    Logging.@warn(
        "Per-step statistics run on the host: the propagating field is copied to the "*
        "device every accepted step to compute them. Pass `stats_period` (or "*
        "`Output.PeriodicStats`) to run them less often, or `Output.nostats` to disable "*
        "them. Device statistics are `gpu/24`'s. (Reported once.)")
    nothing
end

function (so::ScaledOutput)(y, t, dt, yfun)
    needy = needs_host_y(so.o, t)
    needcache = so.needcache && Output.willsave(so.o, y, t, dt)
    if needy || needcache
        yh = _tohost_unscale!(so.ybuf, y, so.Eref)
        needy && _warn_host_stats!(so, y)
    else
        yh = y # nothing inspects it this step: pass the state through untouched, no copy
    end
    so.o(yh, t, dt, ts -> _tohost_unscale!(so.ibuf, yfun(ts), so.Eref))
end

# Metadata and any other call (e.g. `output(dict; group=...)`) pass straight through.
(so::ScaledOutput)(args...; kwargs...) = so.o(args...; kwargs...)

Base.getindex(so::ScaledOutput, args...) = getindex(so.o, args...)
Base.haskey(so::ScaledOutput, key) = haskey(so.o, key)
Output.hasdata(so::ScaledOutput, key) = Output.hasdata(so.o, key)
Output.willsave(so::ScaledOutput, y, t, dt) = Output.willsave(so.o, y, t, dt)
Output.check_cache(so::ScaledOutput, y, t, dt) = Output.check_cache(so.o, y, t, dt)
