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
device() = resolve_device(get(settings, "device", :cpu))

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

_adapt(::Type{Array}, x) = x
_adapt(::Type{A}, x) where {A} = Adapt.adapt(A, x)

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
