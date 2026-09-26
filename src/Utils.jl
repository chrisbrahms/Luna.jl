module Utils
import Dates
import FFTW
import AbstractFFTs
import GPUArraysCore
import Logging
import LibGit2
import FileWatching.Pidfile: mkpidlock
import HDF5
import Luna: settings
import Printf: @sprintf
import Scratch: get_scratch!, clear_scratchspaces!
import Luna

subzero = '\u2080'
subscript(digit::Char) = string(Char(codepoint(subzero)+parse(Int, digit)))
subscript(num::AbstractString) = prod([subscript(chi) for chi in num])
subscript(num::Int) = num >= 0 ? subscript(string(num)) : "₋"*subscript(string(abs(num)))

unsubscript(digit::Char) = string(codepoint(digit)-codepoint(subzero))
unsubscript(num::AbstractString) = prod([unsubscript(chi) for chi in num])

function git_commit()
    try
        repo = LibGit2.GitRepo(lunadir())
        commit = string(LibGit2.GitHash(LibGit2.head(repo)))
        LibGit2.isdirty(repo) && (commit *= " (dirty)")
        return commit
    catch
        "unavailable (Luna is not checkout out for development)"
    end
end

function git_branch()
    try
        repo = LibGit2.GitRepo(lunadir())
        n = string(LibGit2.name(LibGit2.head(repo)))
        branch = split(n, "/")[end]
        return branch
    catch
        "unavailable (Luna is not checkout out for development)"
    end
end

srcdir() = dirname(@__FILE__)

lunadir() = dirname(srcdir())

datadir() = joinpath(srcdir(), "data")

cachedir() = get_scratch!(Luna, "lunacache")

clear_cache() = clear_scratchspaces!(Luna)

function sourcecode()
    src = dirname(@__FILE__)
    luna = dirname(src)
    out = "#= Date: $(Dates.now())\n"
    out *= "git branch: $(git_branch())\n"
    out *= "git commit: $(git_commit())\n"
    out *= "hostname: $(gethostname())\n"
    out *= "=#"
    for folder in (src, luna)
        for obj in readdir(folder)
            if isfile(joinpath(folder, obj))
                if split(obj, ".")[end] in ("md", "jl", "txt", "toml") # avoid binary files
                    out *= "\n" * "#" * "="^8 * obj * "="^8 * "#" * "\n"^2
                    open(joinpath(folder, obj), "r") do file
                        out *= read(file, String)
                    end
                end
            end
        end
    end
    return out
end

function FFTWthreads()
    if Threads.nthreads() == 1
        1
    else
        settings["fftw_threads"] == 0 ? 4*Threads.nthreads() : settings["fftw_threads"]
    end
end

"""
    loadFFTwisdom()

Import accumulated FFTW wisdom from the cache file in `cachedir()`, unless the wisdom
cache is disabled (see [`Luna.set_fftw_wisdom`](@ref)), in which case only the FFTW thread
count is re-asserted and nothing is read.

`Luna.setup` calls this immediately before planning, so the `FFTW.set_num_threads` here is
what makes the thread count Luna plans with independent of anything else in the process
that may have called `FFTW.set_num_threads`. It therefore happens whether or not the
wisdom cache is enabled.
"""
function loadFFTwisdom()
    FFTW.set_num_threads(FFTWthreads())
    settings["fftw_wisdom"] || return
    fpath = joinpath(cachedir(), "FFTWcache_$(FFTWthreads())threads")
    lockpath = joinpath(cachedir(), "FFTWlock")
    isdir(cachedir()) || mkpath(cachedir())
    if isfile(fpath)
        Logging.@info("Found FFTW wisdom at $fpath")
        mkpidlock(lockpath; stale_age=600) do
            try
                FFTW.import_wisdom(fpath)
            catch
                Logging.@info("FFTW wisdom at $fpath is incompatible and has been deleted.")
                rm(fpath)
            end
        end
    else
        Logging.@info("No FFTW wisdom found")
    end
end

"""
    saveFFTwisdom()

Write the FFTW wisdom accumulated in this process to the cache file in `cachedir()`,
unless the wisdom cache is disabled (see [`Luna.set_fftw_wisdom`](@ref)), in which case
this does nothing.
"""
function saveFFTwisdom()
    settings["fftw_wisdom"] || return
    fpath = joinpath(cachedir(), "FFTWcache_$(FFTWthreads())threads")
    lockpath = joinpath(cachedir(), "FFTWlock")
    mkpidlock(lockpath; stale_age=600) do
        isfile(fpath) && rm(fpath)
        isdir(cachedir()) || mkpath(cachedir())
        FFTW.export_wisdom(fpath)
    end
    Logging.@info("FFTW wisdom saved to $fpath")
end

function save_dict_h5(fpath, d; force=false, rmold=false)
    if isfile(fpath) && rmold
        rm(fpath)
    end

    function dict2h5(k::AbstractString, v, parent)
        if HDF5.haskey(parent, k) && !force
            error("Dataset $k exists in $fpath. Set force=true to overwrite.")
        end
        parent[k] = v
    end

    function dict2h5(k::AbstractString, v::BitArray, parent)
        if HDF5.haskey(parent, k) && !force
            error("Dataset $k exists in $fpath. Set force=true to overwrite.")
        end
        parent[k] = Array{Bool, 1}(v)
    end

    function dict2h5(k::AbstractString, v::Nothing, parent)
        if HDF5.haskey(parent, k) && !force
            error("Dataset $k exists in $fpath. Set force=true to overwrite.")
        end
        parent[k] = Float64[]
    end

    function dict2h5(k::AbstractString, v::AbstractDict, parent)
        if !HDF5.haskey(parent, k)
            subparent = HDF5.create_group(parent, k)
        else
            subparent = parent[k]
        end
        for (kk, vv) in pairs(v)
            dict2h5(kk, vv, subparent)
        end
    end

    HDF5.h5open(fpath, "cw") do file
        for (k, v) in pairs(d)
            dict2h5(k, v, file)
        end
    end
end

function save_dict_h5(fpath, t::NamedTuple; kwargs...)
    d = Dict{String, Any}()
    for (k, v) in pairs(t)
        d[string(k)] = v
    end
    save_dict_h5(fpath, d; kwargs...)
end

function load_dict_h5(fpath)
    isfile(fpath) || error("Error loading file $fpath: file does not exist")

    function h52dict(x::HDF5.Dataset)
        return read(x)
    end

    function h52dict(x::Union{HDF5.Group, HDF5.File})
        dd = Dict{String, Any}()
        for n in keys(x)
            dd[n] = h52dict(x[n])
        end
        return dd
    end

    d = HDF5.h5open(fpath) do file
        h52dict(file)
    end
end

#=================================================#
#===============  ARRAY BACKENDS   ===============#
#=================================================#

"""
    Backend

Trait distinguishing host arrays ([`CPUBackend`](@ref)) from device (GPU) arrays
([`DeviceBackend`](@ref)). A device array is anything subtyping
`GPUArraysCore.AbstractGPUArray` (`MtlArray`, `CuArray`, `JLArray`, ...), so the trait
needs no GPU package.

It is used only where the two genuinely differ: FFT planning (FFTW flags and wisdom
against the generic `AbstractFFTs` planners), the host-copy fallbacks, residency checks,
and whether `threaded` may share a broadcast out over threads. It is never used to
select a kernel -- there is one implementation of every per-step operation and it runs on
both.

Query it with [`backend`](@ref).
"""
abstract type Backend end

"Trait value for host (CPU) arrays. See [`Backend`](@ref)."
struct CPUBackend <: Backend end

"Trait value for device (GPU) arrays. See [`Backend`](@ref)."
struct DeviceBackend <: Backend end

"""
    backend(x)

[`CPUBackend()`](@ref CPUBackend) or [`DeviceBackend()`](@ref DeviceBackend) for the
array (or array type) `x`. Wrappers (`SubArray`, `ReshapedArray`, `PermutedDimsArray`)
report the backend of their parent, so a view or a reshape of a device array is
identified correctly.

Anything unrecognised is treated as a host array: mis-classifying an exotic host wrapper
as CPU is a no-op, whereas the reverse would break it.
"""
backend(x) = backend(typeof(x))
backend(::Type) = CPUBackend()
backend(::Type{<:GPUArraysCore.AbstractGPUArray}) = DeviceBackend()
backend(::Type{<:SubArray{T, N, P}}) where {T, N, P} = backend(P)
backend(::Type{<:Base.ReshapedArray{T, N, P}}) where {T, N, P} = backend(P)
backend(::Type{<:PermutedDimsArray{T, N, A, B, P}}) where {T, N, A, B, P} = backend(P)

"Whether `x` lives on a device (GPU). See [`backend`](@ref)."
isdevice(x) = backend(x) isa DeviceBackend

#=================================================#
#============  THREADED BROADCASTS   =============#
#=================================================#

#= Per-step elementwise work (the propagator, the ionisation rate, the plasma terms, the
   fused pointwise responses, the stage combines) is written as broadcasts, which Base
   runs on one thread. `threaded(dest)` marks a destination so that the same broadcast
   expression is shared out over Julia's threads on the host; on a device, or when it
   would not pay, it is exactly the plain broadcast. Each element is computed by the same
   scalar code whichever task runs it, so the result is bit-identical to the serial
   broadcast for any number of threads. A kernel-launch framework (KernelAbstractions)
   was measured for this and was no faster than `Threads.@spawn` over contiguous chunks,
   so Luna takes no extra dependency for it. =#

"""
    THREAD_MINLEN

Default length below which `threaded` runs a broadcast serially. Spawning and
joining the tasks costs ≈15–25 µs on an M1 Pro, and more with 8 threads, where equal
chunks also wait for the slowest (efficiency) core; a cheap broadcast (a few arithmetic
operations per element) only recovers that above ≈2¹⁷ elements.
"""
const THREAD_MINLEN = 1 << 17

"""
    THREAD_MINLEN_HEAVY

`threaded` length threshold for broadcasts dominated by a transcendental function
(`exp` of a complex number, a spline evaluation), which pay from ≈2¹⁵ elements. A
mode-averaged state (a few thousand frequency samples) stays below both thresholds, where
threading was measured to slow a run down.
"""
const THREAD_MINLEN_HEAVY = 1 << 15

"""
    Threaded(dest, minlen)

A broadcast destination wrapper; see `threaded`.
"""
struct Threaded{A}
    dest::A
    minlen::Int
end

"""
    threaded(dest; minlen=THREAD_MINLEN)

Mark `dest` as the destination of a broadcast to be shared out over threads:

```julia
@. \$(Utils.threaded(y; minlen=Utils.THREAD_MINLEN_HEAVY)) = y * exp(linop*dt)
```

(the `\$` keeps `@.` from dotting the call). The broadcast runs on up to
`Threads.nthreads()` tasks, each a slab of `dest` along its last non-singleton dimension (at
least `minlen ÷ 4` elements each) when all of these hold, and as the ordinary
broadcast otherwise:

- `dest` is a host array ([`backend`](@ref)),
- `length(dest) >= minlen`,
- Julia was started with more than one thread,
- `Luna.settings["threaded_broadcasts"]` is `true` (the default; see
  [`Luna.set_threaded_broadcasts`](@ref)),
- the caller is not already inside a threaded region (see `serial_region`).

The results do not depend on which path is taken. The expression must be elementwise,
which a broadcast is; an argument which aliases `dest` without being `dest` itself is
copied first, as Base does.
"""
threaded(dest; minlen=THREAD_MINLEN) = Threaded(dest, minlen)

#= Nesting guard: a caller which already shares its work out over threads (the plasma
   response's column loop) runs its inner broadcasts serially. A process-wide counter
   rather than a task-local flag, because tasks spawned inside the region do not inherit
   task-local storage; two unrelated concurrent runs in one process at worst make each
   other's broadcasts serial, which changes the speed and not the result. =#
const _THREADED_DEPTH = Threads.Atomic{Int}(0)

"""
    serial_region(f)

Call `f()` with `threaded` broadcasts turned into plain ones, for code which is
itself running on several threads.
"""
function serial_region(f)
    Threads.atomic_add!(_THREADED_DEPTH, 1)
    try
        return f()
    finally
        Threads.atomic_sub!(_THREADED_DEPTH, 1)
    end
end

_threadable(dest, minlen) =
    (length(dest) >= minlen) && (Threads.nthreads() > 1) && !isdevice(dest) &&
    (_THREADED_DEPTH[] == 0) && (settings["threaded_broadcasts"] === true)

@inline Base.Broadcast.materialize!(t::Threaded, x) = Base.Broadcast.materialize!(
    t, Base.Broadcast.instantiate(
        Base.Broadcast.Broadcasted(identity, (x,), axes(t.dest))))

function Base.Broadcast.materialize!(t::Threaded, bc::Base.Broadcast.Broadcasted)
    dest = t.dest
    _threadable(dest, t.minlen) || return Base.Broadcast.materialize!(dest, bc)
    #= What Base's `materialize!`/`copyto!` do before their loop: fix the axes to
       `dest`'s (checking that the arguments broadcast to them), drop the style, and
       unalias/extrude the arguments. =#
    bc′ = Base.Broadcast.preprocess(
        dest, Base.Broadcast.instantiate(_nostyle(bc, axes(dest))))
    #= Chunks are slabs along the last non-singleton dimension, each iterated as a
       `CartesianIndices` block the way Base iterates the whole array: a linear index
       per element would cost an integer division per dimension, which is more than a
       cheap broadcast's arithmetic. =#
    ax = axes(bc′)
    #= Only one dimension is split, so the task count is capped by its length: a
       (2^17, 2) array gets two tasks. Luna's multi-column blocks have their columns (or
       polarisation, mode or transverse axes) last and long enough for this. =#
    d = something(findlast(a -> length(a) > 1, ax), 1)
    r = ax[d]
    # at least minlen/4 elements per task, so a short array does not pay for idle tasks
    nt = clamp(length(dest) ÷ max(t.minlen ÷ 4, 1), 1, min(Threads.nthreads(), length(r)))
    chunk = cld(length(r), nt)
    Threads.atomic_add!(_THREADED_DEPTH, 1)
    try
        @sync for k in 1:nt
            lo = first(r) + (k - 1)*chunk
            hi = min(lo + chunk - 1, last(r))
            lo <= hi || continue
            R = CartesianIndices(ntuple(i -> i == d ? (lo:hi) : ax[i], length(ax)))
            Threads.@spawn _bcchunk!(dest, bc′, R)
        end
    finally
        Threads.atomic_sub!(_THREADED_DEPTH, 1)
    end
    dest
end

# The style-free `Broadcasted` Base's `copyto!` loops over, with the axes of `dest`.
@static if VERSION >= v"1.10"
    _nostyle(bc, ax) = Base.Broadcast.Broadcasted(nothing, bc.f, bc.args, ax)
else
    _nostyle(bc, ax) = Base.Broadcast.Broadcasted{Nothing}(bc.f, bc.args, ax)
end

@noinline function _bcchunk!(dest, bc, R)
    @inbounds @simd for I in R
        dest[I] = bc[I]
    end
    nothing
end

#=================================================#
#===============  FFT PLANNING   =================#
#=================================================#

"""
    plan_ft(x, dims)

Plan the forward time-to-frequency transform of an array like `x` along `dims`: a
real-to-complex transform if `x` is real (field-resolved grids) and a complex-to-complex
one if it is complex (envelope grids).

On the host this is FFTW with Luna's configured planning flags, so the wisdom logic of
`Utils.loadFFTwisdom`/`Utils.saveFFTwisdom` applies. On a device it is the generic
`AbstractFFTs` planner, which device FFT libraries implement and which takes no flags.

Device plans work on plain arrays of exactly the planned shape -- Metal's in particular
reject views -- so `x` must be the buffer the transform will actually be applied to (or
one just like it).
"""
plan_ft(x, dims) = _plan_ft(backend(x), x, dims)
_plan_ft(::CPUBackend, x::AbstractArray{<:Real}, dims) =
    FFTW.plan_rfft(x, dims, flags=settings["fftw_flag"])
_plan_ft(::CPUBackend, x::AbstractArray{<:Complex}, dims) =
    FFTW.plan_fft(x, dims, flags=settings["fftw_flag"])
_plan_ft(::DeviceBackend, x::AbstractArray{<:Real}, dims) = AbstractFFTs.plan_rfft(x, dims)
_plan_ft(::DeviceBackend, x::AbstractArray{<:Complex}, dims) = AbstractFFTs.plan_fft(x, dims)

"""
    plan_ift(FT)

The explicit inverse of the forward plan `FT`, as an `AbstractFFTs.ScaledPlan` carrying
the unnormalised backward plan and the `1/N` factor.

Luna holds the inverse plan rather than calling `ldiv!(y, FT, x)`, on every backend, for
two reasons: `ldiv!` is a multiply followed by a separate scaling pass, and the `1/N` can
instead be folded into the scale factor the oversampling copy already applies
([`iscale`](@ref), [`iplan`](@ref)), which removes that pass. The result differs from
`ldiv!` at rounding level.

For a real-to-complex `FT` the inverse is a `brfft`, which **overwrites its input**.
"""
plan_ift(FT) = inv(FT)

"""
    iplan(IFT)

The unnormalised backward plan held by an inverse plan. Split from its normalisation
factor ([`iscale`](@ref)) so that the factor can be folded into the scale of the
oversampling copy; see [`plan_ift`](@ref).
"""
iplan(p::AbstractFFTs.ScaledPlan) = p.p
iplan(p) = _notinverse(p)

"""
    iscale(IFT)

The normalisation factor of an inverse plan, `1/N`, separated from the plan itself
([`iplan`](@ref)).
"""
iscale(p::AbstractFFTs.ScaledPlan) = p.scale
iscale(p) = _notinverse(p)

#= Every inverse plan Luna makes is a `ScaledPlan`: that is what `inv` returns for an
   FFTW plan, for a Metal or CUDA plan, and for the JLArray shim in the tests. Anything
   else reaching here is a forward plan passed where an inverse one belongs, which would
   otherwise transform the wrong way -- a method error on a real grid, silently wrong
   output on an envelope grid. Dispatch decides, so the check costs nothing. =#
_notinverse(p) = error(
    "an inverse plan is required here, not $(typeof(p)). Build it with "*
    "`Utils.plan_ift(FT)`, or use the transform's own `IFT` field.")

function format_elapsed(ms::Dates.Millisecond)
    stot = Dates.value(ms)/1000 # total seconds
    seconds = stot % 60
    stot -= seconds
    mtot = stot ÷ 60
    minutes = mtot % 60
    mtot -= minutes
    hours = mtot ÷ 60
    out = @sprintf("%.3f seconds", seconds)
    minstr = abs(minutes) == 1 ? "minute" : "minutes"
    hrstr = abs(hours) == 1 ? "hour" : "hours"
    if abs(hours) > 0
        out = @sprintf("%d %s, ", minutes, minstr) * out
        out = @sprintf("%d %s, ", hours, hrstr) * out
    elseif abs(minutes) > 0
        out = @sprintf("%d %s, ", minutes, minstr) * out
    end
    out
end

end
