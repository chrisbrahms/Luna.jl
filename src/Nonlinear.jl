module Nonlinear
import Luna
import Luna.PhysData: ε_0, e_ratio
import Luna: Maths, Utils
import Adapt
import FFTW
import Logging
import LinearAlgebra: mul!, ldiv!
import Rotations: RotZY, RotYZ, RotMatrix, RotMatrix3
import StaticArrays: SMatrix, SVector, MArray

#=================================================#
#=============  RESPONSE PROTOCOL  ===============#
#=================================================#

"""
    ResponseKind

How a nonlinear response is evaluated by
[`NonlinearRHS.Et_to_Pt!`](@ref Luna.NonlinearRHS.Et_to_Pt!). One of

| kind | contract |
| --- | --- |
| [`Pointwise`](@ref) | the contribution at one element of the field block depends only on that element |
| [`VectorPointwise`](@ref) | the contribution at one time sample couples the two polarisation components of that sample |
| [`Batched`](@ref) | the response is called once with the whole `(nt, npol, ncols)` block |
| [`Columnwise`](@ref) | the response is called once per column, on the host |

[`kind`](@ref) returns it, `Columnwise()` by default.
"""
abstract type ResponseKind end

"""
    Pointwise()

The contribution of the response at one element of the time-domain field block depends
only on the field at that element (and on per-sample coefficient arrays aligned with the
time axis, e.g. a carrier phase). Every pointwise response of a transform is evaluated in
**one** fused broadcast per right-hand side, with no intermediate buffer.

A response of this kind provides [`pointwise_kernel`](@ref), or
[`pointwise_expr`](@ref) if it needs per-sample arrays of its own.
"""
struct Pointwise <: ResponseKind end

"""
    VectorPointwise()

The contribution of the response at one time sample depends only on the two polarisation
components of the field at that sample. Applies only to a two-component field block; the
same response reports [`Pointwise`](@ref) (or `Columnwise`) for a scalar one.

A response of this kind provides [`vector_kernel`](@ref), or [`vector_expr`](@ref) if it
needs per-sample arrays of its own; either returns the two lab-frame components as an
`SVector{2}`.

The two components are written by two broadcasts over
`view(E, :, 1, ..)`/`view(E, :, 2, ..)`, one per output component, fused across all the
pointwise responses of the transform. Luna's buffers are `(nt, npol, ncols)`, so the
polarisation index is the *slow* axis and a single broadcast writing an `SVector{2}`
through a `reinterpret`ed view is not available; see the developer guide.
"""
struct VectorPointwise <: ResponseKind end

"""
    Batched()

The response is called once per right-hand side with the whole `(nt, npol, ncols)` block
as `resp!(out, E, ρ)`, accumulating into `out`. It owns whatever full-size buffers it
needs, in the array type and precision of the run. [`HostResponse`](@ref) is one.
"""
struct Batched <: ResponseKind end

"""
    Columnwise()

Luna's historical contract, and the default: the response is called as `resp!(out, E, ρ)`
once per column of the block, accumulating into `out`, with host arrays in physical SI
units and `Float64`. Any callable `(out, E, ρ)` works, which is what makes an ad hoc
response possible.

A columnwise response cannot run on a device or in reduced precision as it stands; on
such a run [`rescale`](@ref) wraps it in a [`HostResponse`](@ref), which copies the block
to the host and back at every right-hand side.
"""
struct Columnwise <: ResponseKind end

"""
    kind(response) -> ResponseKind
    kind(response, npol) -> ResponseKind

How this response is evaluated; see [`ResponseKind`](@ref). `Columnwise()` by default, so
an ad hoc response keeps working without declaring anything.

The two-argument form is the kind for a field block with `npol` polarisation components,
given as `Val(1)`/`Val(2)` (an `Integer` is accepted and converted). It defaults to the
one-argument form, and exists for responses whose physics couples the components, such as
the vector Kerr forms:

    Nonlinear.kind(::KerrField) = Nonlinear.Pointwise()
    Nonlinear.kind(::KerrField, ::Val{2}) = Nonlinear.VectorPointwise()

The dispatcher resolves `npol` from the field block before it groups the responses, so a
response is free to report different kinds for the scalar and vector cases, including a
`Columnwise()` for one and a fused kind for the other.
"""
kind(r) = Columnwise()
kind(r, ::Val) = kind(r)
kind(r, npol::Integer) = kind(r, Val(Int(npol)))

"""
    coefficients(response, ρ, scaling)

The combined scalar coefficient(s) of `response` at number density `ρ`, in the units of
`scaling` (see [`Luna.UnitScaling`](@ref)).

This is the one place where a response's physical constants, the density and the powers
of `E_ref`/`P_ref` are put together. It runs on the host, in `Float64`; the kernel
converts the result to the run's real element type with
[`Luna.scalar`](@ref Luna.scalar), so no `Float64` reaches a device. Use
[`Luna.polscale`](@ref Luna.polscale) for the power of `E_ref` a response of a given
polynomial degree carries.

Returns a single scalar for a response with one coefficient, or a tuple for one with
several.
"""
function coefficients end

"""
    pointwise_kernel(response, E, ρ, scaling) -> f

A callable `f(e)` returning the contribution of a [`Pointwise`](@ref) response to the
nonlinear polarisation at one element of the field block `E`, at density `ρ`.

`E` is passed so that the kernel can convert its coefficients to `real(eltype(E))` with
[`Luna.scalar`](@ref Luna.scalar); it must not be indexed. `f` must capture nothing but
`isbits` scalars, since it is the body of a broadcast which may run on a GPU.

This is the simple form of [`pointwise_expr`](@ref), which is what the dispatcher
actually calls.
"""
function pointwise_kernel end

"""
    pointwise_expr(response, E, ρ, scaling) -> Broadcasted

The lazy elementwise expression for the contribution of a [`Pointwise`](@ref) response to
the nonlinear polarisation of field block `E` at density `ρ`, with the same axes as `E`.

The default builds it from [`pointwise_kernel`](@ref). Override it for a response which
needs per-sample arrays of its own (a carrier phase, a tabulated coefficient): return
`Base.broadcasted(f, E, arr, ...)`, where every array is aligned with the time axis and
already in the run's array type and precision (see [`rescale`](@ref) and
[`resident_arrays`](@ref)).

The dispatcher sums the expressions of all the pointwise responses of a transform and
materialises them as one broadcast, so no response allocates or writes a buffer.
"""
pointwise_expr(r, E, ρ, scaling) =
    Base.broadcasted(pointwise_kernel(r, E, ρ, scaling), E)

"""
    vector_kernel(response, E, ρ, scaling) -> f

A callable `f(ex, ey)` returning the two lab-frame polarisation components of the
contribution of a [`VectorPointwise`](@ref) response at one time sample, as an
`SVector{2}`. `E` is passed for its element type only, as in
[`pointwise_kernel`](@ref).
"""
function vector_kernel end

"""
    vector_expr(response, Ex, Ey, ρ, scaling) -> Broadcasted

The lazy expression for a [`VectorPointwise`](@ref) response over the two polarisation
component views `Ex`, `Ey`, with `SVector{2}` elements. The default builds it from
[`vector_kernel`](@ref); override it for a response which needs per-sample arrays of its
own, as for [`pointwise_expr`](@ref).
"""
vector_expr(r, Ex, Ey, ρ, scaling) =
    Base.broadcasted(vector_kernel(r, Ex, ρ, scaling), Ex, Ey)

"""
    batched!(response, out, E, ρ, scaling)

Evaluate a [`Batched`](@ref) response on the whole field block `E`, accumulating into
`out`. `scaling` is the [`Luna.UnitScaling`](@ref) the block and `out` are expressed in.

The default ignores the scaling and calls `response(out, E, ρ)`, which is right for a
response whose coefficients already carry it — every response built by
[`rescale`](@ref)`(r, spec, scaling, Et)`, including [`HostResponse`](@ref). Override it
only for a response which would rather combine the scaling per call than at construction.

This is the batched counterpart of [`coefficients`](@ref): a batched response owns its
buffers and its own loop, so the dispatcher can give it nothing but the block, the
density and the units.
"""
batched!(r, out, E, ρ, scaling) = r(out, E, ρ)

"""
    rescale(response, spec, scaling)
    rescale(response, spec, scaling, Et)

The same nonlinear response, with every array it carries converted to `spec`'s array type
and precision (see [`Luna.DeviceSpec`](@ref)) and with whatever state depends on the shape
of the field block allocated, ready for a run whose state is expressed in the units of
`scaling` (see [`Luna.UnitScaling`](@ref)).

Scalar coefficients are *not* converted here: they are combined with the density and the
unit scaling once per right-hand side by [`coefficients`](@ref), in `Float64`, and only
the result is converted. A response struct therefore never enters a kernel — only the
scalars its kernel captures and the arrays this function converted.

**Which form to implement.** A transform calls the **four**-argument form on each of its
responses at construction, passing `Et`, a prototype of the time-domain field block the
response will be called with: `size(Et)` is `(nt, npol, ncols...)`, `eltype(Et)` says
whether the field is real or complex and in what precision, and `similar(Et)` allocates a
buffer of the run's array type. A [`Batched`](@ref) response, which owns full-size
buffers, implements that form, so that its buffers exist at construction and the
transform's [`Luna.assert_resident`](@ref) check can see them (GPU_PLAN.md §4.2 rule 5).
Everything else implements the **three**-argument form, which the four-argument fallback
delegates to.

**What the fallbacks do.** Both forms pass the response through unchanged for an unscaled
`Float64` run on host arrays, which is every run on the default CPU path, so an ad hoc
response written as a closure keeps working. Otherwise:

- a [`Columnwise`](@ref) response is wrapped in a [`HostResponse`](@ref), which runs it on
  the host. This needs `Et`, so it is only available from the four-argument form;
- a [`Pointwise`](@ref) or [`VectorPointwise`](@ref) response with no arrays of its own
  ([`resident_arrays`](@ref) empty) is passed through, which is right for every response
  whose coefficients are scalars: `coefficients` combines them per call in `Float64`;
- a [`Batched`](@ref) response is **not** passed through, because nothing else would give
  it the scaling or allocate its buffers in the run's array type. It needs its own
  `rescale` method unless the run is unscaled and on the host;
- a response which carries an array and has no `rescale` method is an error, since nothing
  would have converted it.

A response which declares a device kind must also list every non-`isbits` array it carries
in [`resident_arrays`](@ref); the fallbacks check this structurally and name the field
which is missing.
"""
function rescale(r, spec, scaling)
    _isdefaultrun(spec, scaling) && return r
    _rescale_fallback(kind(r), r, spec, scaling)
end

rescale(r, spec, scaling, Et) = _rescale_dims(kind(r), r, spec, scaling, Et)

#= Only a columnwise response needs the prototype in the generic path: its wrapper is
   built here rather than by the response, which by definition has no `rescale` method of
   its own. Everything else implements the three-argument form, or overrides the
   four-argument one. =#
_rescale_dims(::ResponseKind, r, spec, scaling, Et) = rescale(r, spec, scaling)

function _rescale_dims(::Columnwise, r, spec, scaling, Et)
    _isdefaultrun(spec, scaling) && return r
    HostResponse(r, spec, scaling, Et)
end

"The default CPU path: host arrays, double precision, physical units."
_isdefaultrun(spec, scaling) =
    !Luna.isdevicespec(spec) && Luna.realtype(spec) === Float64 && Luna.isunity(scaling)

_rescale_fallback(::Columnwise, r, spec, scaling) = error(
    "the nonlinear response $(nameof(typeof(r))) is columnwise, so it has to be wrapped "*
    "in a `Nonlinear.HostResponse` to run on $(Luna.arraytype(spec))/"*
    "$(Luna.realtype(spec)) in $(scaling), and that needs the shape of the field block. "*
    "Call `Nonlinear.rescale(response, spec, scaling, Et)` with a prototype of the block.")

#= A response with a device kernel and no arrays of its own -- which is every response
   whose coefficients are scalars -- needs no conversion: `coefficients` combines them in
   Float64 per call and the kernel converts the result. One which carries an array has to
   say how it moves. =#
function _rescale_fallback(k::Union{Pointwise, VectorPointwise}, r, spec, scaling)
    _check_listed_arrays(r, k)
    isempty(resident_arrays(r)) && return r
    _no_rescale_error(k, r, spec)
end

#= A batched response is different: it is handed the block and the density and nothing
   else, so the only place its coefficients can meet `E_ref`/`P_ref` is its own `rescale`
   method, and the only place its buffers can be allocated in the run's array type is the
   same method. Passing it through would run physical-unit coefficients against a scaled
   state -- silently wrong rather than an error. =#
function _rescale_fallback(k::Batched, r, spec, scaling)
    _check_listed_arrays(r, k)
    (Luna.isunity(scaling) && !Luna.isdevicespec(spec)) && return r
    error("the nonlinear response $(nameof(typeof(r))) declares `Nonlinear.kind` = "*
          "Batched() but has no `Nonlinear.rescale` method. A batched response is called "*
          "with the block, the density and the units and nothing else, so it needs a "*
          "`rescale(r, spec, scaling, Et)` method which combines its coefficients with "*
          "`scaling` (see `Luna.polscale`) and allocates its buffers with `similar(Et)`. "*
          "Without one it would run physical-unit coefficients against a state in "*
          "$(scaling) on $(Luna.arraytype(spec))/$(Luna.realtype(spec)).")
end

_no_rescale_error(k, r, spec) = error(
    "the nonlinear response $(nameof(typeof(r))) declares `Nonlinear.kind` = $(k) and "*
    "carries $(length(resident_arrays(r))) array(s), but has no `Nonlinear.rescale` "*
    "method, so they cannot be converted to $(Luna.arraytype(spec))/"*
    "$(Luna.realtype(spec)). Give it a `rescale` method, or declare `Nonlinear.kind` = "*
    "Columnwise() to run it on the host.")

#= Structural check, once per setup: a response which declares a device kind may not keep
   an array on the host, so every non-`isbits` `AbstractArray` field of it has to be one
   `resident_arrays` names (a `StaticArrays` matrix or a `Rotations` matrix is `isbits`
   and travels inside the struct, so it is exempt). Catches the array an author forgot to
   list, which the pass-through would otherwise hand to a kernel as a host array --
   a kernel-compilation error far from the cause on Metal, and nothing at all on
   JLArrays. =#
function _check_listed_arrays(r, k)
    listed = resident_arrays(r)
    for f in fieldnames(typeof(r))
        x = getfield(r, f)
        (x isa AbstractArray && !isbits(x)) || continue
        any(a -> a === x, listed) && continue
        error("the nonlinear response $(nameof(typeof(r))) declares `Nonlinear.kind` = "*
              "$(k) and carries the array `$(f)::$(typeof(x))`, which "*
              "`Nonlinear.resident_arrays` does not list. Every array a device kernel "*
              "broadcasts against has to be listed, so that `rescale` converts it and "*
              "the transform's residency assertion checks it.")
    end
    nothing
end

"""
    rescale_responses(responses, spec, scaling, Et)

[`rescale`](@ref) applied to a transform's whole response collection, including the
tuple-of-tuples a gas mixture is (each inner tuple is rescaled element by element and
stays a tuple, so that it still lines up with the density it is paired with).

`Et` is the prototype of the time-domain field block, as for the four-argument
[`rescale`](@ref). Every transform calls this at construction, so that a response which
owns block-sized buffers ([`Batched`](@ref)) has them in the right shape, array type and
precision before the transform asserts residency.

The collection type is preserved, because whether the responses are a tuple is what
[`NonlinearRHS.Et_to_Pt!`](@ref Luna.NonlinearRHS.Et_to_Pt!) uses to decide between the
grouped dispatch and the historical per-response loop.

On the default CPU path — host arrays, `Float64`, physical units — every fallback
returns the response unchanged, so this is a no-op except for a response which
implements the four-argument form because it has something to size to the block.
"""
rescale_responses(responses, spec, scaling, Et) =
    map(r -> _rescale_each(r, spec, scaling, Et), responses)

_rescale_each(r::Tuple, spec, scaling, Et) =
    map(x -> _rescale_each(x, spec, scaling, Et), r)
_rescale_each(r, spec, scaling, Et) = rescale(r, spec, scaling, Et)

"""
    resident_arrays_all(responses) -> Tuple

Every array [`resident_arrays`](@ref) names, over a whole response collection and
through the tuple-of-tuples a gas mixture is. This is what a transform hands to
[`Luna.assert_resident`](@ref) along with its own buffers.
"""
resident_arrays_all(responses) =
    reduce((a, r) -> (a..., _resident_each(r)...), Tuple(responses); init=())

_resident_each(r::Tuple) = resident_arrays_all(r)
_resident_each(r) = resident_arrays(r)

"""
    resident_arrays(response) -> Tuple

The arrays a nonlinear response carries which a kernel broadcasts against, and which
therefore have to live on the run's array type and precision. A transform passes them to
[`Luna.assert_resident`](@ref) at construction, together with its own buffers and grid
mirrors.

Empty by default, which is right for a response whose coefficients are all scalars and
for anything which does not run on a device at all.
"""
resident_arrays(r) = ()

"""
    device_capable(response) -> Bool

Whether `response` has a kernel of its own which runs in the array type and precision of
the run, i.e. whether its [`kind`](@ref) is anything but [`Columnwise`](@ref).

A columnwise response still *works* on a device, through the [`HostResponse`](@ref)
fallback, but at the cost of a device-to-host copy of the whole block at every right-hand
side. `device_capable` is therefore about performance, not possibility.

`Interface.jl` uses it to decide, for an *unspecified* `device` request, whether a
mode-averaged `prop_capillary` call can follow `Luna.settings["device"]` (every response
it built is device-capable) or must stay on the CPU regardless of the global setting (at
least one is not, e.g. plasma or Raman) -- so that loading a GPU package does not turn a
default, field-resolved `prop_capillary` call into a slow one.
"""
device_capable(r) = !(kind(r) isa Columnwise)

#=================================================#
#=============  THE HOST FALLBACK  ===============#
#=================================================#

"""
    HostResponse(response, spec, scaling, Et)

A [`Columnwise`](@ref) response made to work in a run whose state is not host `Float64`:
a [`Batched`](@ref) wrapper which, at every right-hand side, copies the whole field block
to a host `Float64`/`ComplexF64` buffer in physical SI units, calls `response` on it
column by column exactly as the host path does, converts the result back into the units
and precision of the run, and adds it to the output.

`Et` is a prototype of the block (see [`rescale`](@ref)): every buffer is allocated here,
at construction, so the per-call code is fully typed. A device run needs four of them —
`Eh` and `Ph` on the host in physical `Float64`, `stage` on the host in the run's element
type, and `Pd` in the run's array type — because `copyto!` between a host array and a
device array does not convert the precision. A scaled host run needs only `Eh` and `Ph`
and leaves `stage`/`Pd` as `nothing`.

This is the fallback which keeps Luna hackable (GPU_PLAN.md §3): a user-written
`resp!(out, E, ρ)` closure, and every response which has not yet been given a device
kernel, works on a GPU or in `Float32` without being rewritten. It is correct and slow —
two copies and a host evaluation of the response per right-hand side, which on a GPU also
serialises the step — and its constructor logs one line saying so.

Built by [`rescale`](@ref); there is no reason to build one directly.
"""
struct HostResponse{R, H, S, D}
    resp::R
    Eref::Float64 # the field the state is measured in
    invfac::Float64 # 1/(Pref*Eref): the polarisation the buffer is measured in
    Eh::H # host field buffer, physical units, Float64/ComplexF64
    Ph::H # host polarisation buffer, physical units
    stage::S # host buffer in the run's element type, or nothing on a host run
    Pd::D # buffer in the run's array type and element type, or nothing on a host run
end

function HostResponse(resp, spec, scaling, Et)
    Logging.@info(
        "The nonlinear response $(nameof(typeof(resp))) runs on the host: the field "*
        "block is copied to the host in physical units at every right-hand side, the "*
        "response evaluated column by column in Float64, and the result copied back. "*
        "Correct but slow; give it a `Nonlinear.kind` and a kernel to run it in place.")
    HT = _hosteltype(eltype(Et))
    Eh = zeros(HT, size(Et))
    Ph = zeros(HT, size(Et))
    if Luna.isdevicespec(spec)
        stage = zeros(eltype(Et), size(Et))
        Pd = fill!(similar(Et), zero(eltype(Et)))
    else
        stage = nothing
        Pd = nothing
    end
    HostResponse(resp, scaling.Eref, 1/(scaling.Pref*scaling.Eref), Eh, Ph, stage, Pd)
end

kind(::HostResponse) = Batched()

#= `Pd` is the one buffer which has to be on the device; the others are host by design.
   Listed so that the transform's residency assertion covers it. =#
resident_arrays(h::HostResponse) = (h.Pd,)

"The element type a host copy of a block of element type `T` is held in."
_hosteltype(::Type{T}) where {T<:Real} = Float64
_hosteltype(::Type{Complex{T}}) where {T<:Real} = ComplexF64

function (h::HostResponse)(out, E, ρ)
    _tohost!(h.Eh, E, h.stage, h.Eref)
    fill!(h.Ph, 0)
    _hostcolumns!(h.Ph, h.Eh, h.resp, ρ)
    _toout!(out, h.Ph, h.Pd, h.stage, h.invfac)
    out
end

# A host run needs no staging copy: the broadcast converts the precision directly.
_tohost!(Eh, E, ::Nothing, Eref) = (@. Eh = Eref * E)
_tohost!(Eh, E, stage, Eref) = (copyto!(stage, E); @. Eh = Eref * stage)

_toout!(out, Ph, Pd, ::Nothing, invfac) = (@. out += invfac * Ph)
function _toout!(out, Ph, Pd, stage, invfac)
    @. stage = invfac * Ph
    copyto!(Pd, stage)
    @. out += Pd
    out
end

#= The columnwise contract: a 1-D or `(nt, npol)` block is one column, anything bigger is
   sliced along every axis past the polarisation one, exactly as `Et_to_Pt!`'s `idcs`
   loop does. =#
function _hostcolumns!(Ph::AbstractArray{<:Any, N}, Eh, resp, ρ) where {N}
    N <= 2 && return (resp(Ph, Eh, ρ); Ph)
    for i in CartesianIndices(size(Eh)[3:end])
        resp(view(Ph, :, :, i), view(Eh, :, :, i), ρ)
    end
    Ph
end

#= The Kerr responses are structs rather than closures so that they can carry an `Adapt`
   rule for their arrays and take the protocol's methods. The constructors
   `Kerr_field(γ3)` etc. keep their signatures and return the structs. =#

"""
    KerrField(γ3)

Kerr response for a real (field-resolved) field; built by `Kerr_field(γ3)`. `γ3` is the
third-order hyperpolarisability in physical units, e.g.
[`PhysData.γ3_gas`](@ref Luna.PhysData.γ3_gas).

[`Pointwise`](@ref) for a scalar field, [`VectorPointwise`](@ref) for a two-component one.
"""
struct KerrField{T}
    γ3::T
end

"Kerr response for real field"
Kerr_field(γ3) = KerrField(γ3)

kind(::KerrField) = Pointwise()
kind(::KerrField, ::Val{2}) = VectorPointwise()

coefficients(k::KerrField, ρ, scaling) = ρ*ε_0*k.γ3*Luna.polscale(scaling, 3)

#= The arithmetic, once. Everything which evaluates this response -- the fused broadcast,
   the columnwise call operator below -- goes through these two, so there is one kernel
   body per response and the two paths cannot drift apart. =#
_kerrfield(fac) = e -> fac*e^3
_kerrfield(fac, ::Val{2}) = (ex, ey) -> (s = fac*(ex^2 + ey^2); SVector(s*ex, s*ey))

pointwise_kernel(k::KerrField, E, ρ, scaling) =
    _kerrfield(Luna.scalar(E, coefficients(k, ρ, scaling)))

vector_kernel(k::KerrField, E, ρ, scaling) =
    _kerrfield(Luna.scalar(E, coefficients(k, ρ, scaling)), Val(2))

#= The columnwise contract, in physical units: what a low-level caller, or a response
   collection which is not a tuple, still reaches. =#
function (k::KerrField)(out, E, ρ)
    fac = Luna.scalar(E, coefficients(k, ρ, Luna.UNIT_SCALING))
    if size(E, 2) == 1
        KerrScalar!(out, E, fac)
    else
        KerrVector!(out, E, fac)
    end
end

"Accumulate the scalar Kerr polarisation for a combined coefficient `fac`."
function KerrScalar!(out, E, fac)
    f = _kerrfield(fac)
    @. out += f(E)
end

#= One broadcast per polarisation component, over views of the two columns: Luna's
   buffers put the polarisation index on the slow axis, so a single broadcast cannot
   write both components at once. =#
@doc (@doc KerrScalar!)
function KerrVector!(out, E, fac)
    f = _kerrfield(fac, Val(2))
    Ex = view(E, :, 1)
    Ey = view(E, :, 2)
    ox = view(out, :, 1)
    oy = view(out, :, 2)
    @. ox += first(f(Ex, Ey))
    @. oy += last(f(Ex, Ey))
    out
end

"Kerr response for real field but without THG"
function Kerr_field_nothg(γ3, n)
    E = Array{Float64}(undef, n)
    hilbert = Maths.plan_hilbert(E)
    Kerr = let γ3 = γ3, hilbert = hilbert
        function Kerr(out, E, ρ)
            out .+= ρ*3/4*ε_0*γ3.*abs2.(hilbert(E)).*E
        end
    end
end

"""
    KerrEnv(γ3)

Kerr response for an envelope field without THG; built by `Kerr_env(γ3)`. See
[`KerrField`](@ref) for `γ3` and for the kinds.
"""
struct KerrEnv{T}
    γ3::T
end

"Kerr response for envelope"
Kerr_env(γ3) = KerrEnv(γ3)

kind(::KerrEnv) = Pointwise()
kind(::KerrEnv, ::Val{2}) = VectorPointwise()

#= The 3/4 is part of the coefficient rather than of the kernel body, where its Float64
   literal would promote a Float32 kernel. The value is unchanged: the broadcast evaluated
   the same product of scalars for every element. =#
coefficients(k::KerrEnv, ρ, scaling) = 3/4*(ρ*ε_0*k.γ3)*Luna.polscale(scaling, 3)

_kerrenv(fac) = e -> fac*abs2(e)*e

function _kerrenv(fac, ::Val{2}, ::Type{R}) where {R}
    c23 = convert(R, 2/3)
    c13 = convert(R, 1/3)
    (ex, ey) -> SVector(fac*((abs2(ex) + c23*abs2(ey))*ex + c13*conj(ex)*ey^2),
                        fac*((abs2(ey) + c23*abs2(ex))*ey + c13*conj(ey)*ex^2))
end

pointwise_kernel(k::KerrEnv, E, ρ, scaling) =
    _kerrenv(Luna.scalar(E, coefficients(k, ρ, scaling)))

vector_kernel(k::KerrEnv, E, ρ, scaling) =
    _kerrenv(Luna.scalar(E, coefficients(k, ρ, scaling)), Val(2), real(eltype(E)))

function (k::KerrEnv)(out, E, ρ)
    fac = Luna.scalar(E, coefficients(k, ρ, Luna.UNIT_SCALING))
    if size(E, 2) == 1
        KerrScalarEnv!(out, E, fac)
    else
        KerrVectorEnv!(out, E, fac)
    end
end

"`fac` includes the factor 3/4; see [`KerrEnv`](@ref)."
function KerrScalarEnv!(out, E, fac)
    f = _kerrenv(fac)
    @. out += f(E)
end

@doc (@doc KerrScalarEnv!)
function KerrVectorEnv!(out, E, fac)
    f = _kerrenv(fac, Val(2), real(eltype(E)))
    Ex = view(E, :, 1)
    Ey = view(E, :, 2)
    ox = view(out, :, 1)
    oy = view(out, :, 2)
    @. ox += first(f(Ex, Ey))
    @. oy += last(f(Ex, Ey))
    out
end

"""
    KerrEnvTHG(γ3, C)

Kerr response for an envelope field including THG; built by `Kerr_env_thg(γ3, ω0, t)`.
`C` is the carrier factor `exp(2iω₀t)` on the (oversampled) time grid. See Eq. 4, Genty
et al., Opt. Express 15 5382 (2007).

[`Pointwise`](@ref) whatever the number of polarisation components: the expression is
elementwise in the whole block, with `C` broadcast along the time axis. Because it
carries a per-sample array it overrides [`pointwise_expr`](@ref) rather than supplying a
[`pointwise_kernel`](@ref), and `C` is a [`resident_arrays`](@ref) entry.
"""
struct KerrEnvTHG{T, V}
    γ3::T
    C::V
end

"Kerr response for envelope but with THG"
Kerr_env_thg(γ3, ω0, t) = KerrEnvTHG(γ3, exp.(2im*ω0.*t))

kind(::KerrEnvTHG) = Pointwise()

coefficients(k::KerrEnvTHG, ρ, scaling) = ρ*ε_0*k.γ3/4*Luna.polscale(scaling, 3)

_kerrthg(fac) = (e, c) -> fac*(3*abs2(e) + c*e^2)*e

#= The carrier array is a second broadcast argument rather than something the kernel
   captures, so that it is aligned with the time axis whatever the block's shape. =#
pointwise_expr(k::KerrEnvTHG, E, ρ, scaling) =
    Base.broadcasted(_kerrthg(Luna.scalar(E, coefficients(k, ρ, scaling))), E, k.C)

function (k::KerrEnvTHG)(out, E, ρ)
    f = _kerrthg(Luna.scalar(E, coefficients(k, ρ, Luna.UNIT_SCALING)))
    C = k.C
    @. out += f(E, C)
end

rescale(k::KerrEnvTHG, spec, scaling) = KerrEnvTHG(k.γ3, Luna.todevice(spec, k.C))

Adapt.adapt_structure(to, k::KerrEnvTHG) = KerrEnvTHG(k.γ3, Adapt.adapt(to, k.C))

resident_arrays(k::KerrEnvTHG) = (k.C,)

struct Chi2Field{χT}
    χ2::χT
    toCrystal::RotMatrix3{Float64}
    toLab::RotMatrix3{Float64}
    χ2_toLab::χT # combined matrix to multiply by toLab * χ2
    El::Vector{Float64} # field in the lab frame
    Ec::Vector{Float64} # field in the crystal frame
    Enl::Vector{Float64} # field products in the crystal frame
    Pl::Vector{Float64} # polarisation in the lab frame
end

"""
    Chi2Field(θ, ϕ, χ2)

Construct a second-order nonlinear polarisation response for real, two-component
electric fields in the lab frame.

`θ` and `ϕ` (radians) define the crystal orientation relative to the lab frame.
`χ2` must be a 3×6 second-order susceptibility tensor in contracted notation,
with column order `[xx, yy, zz, yz, xz, xy]` (where mixed terms are multiplied by 2
in [`field_products!`](@ref)).

The returned callable adds \$ε_0 P_{NL}\$ to `out` when invoked as
`response(out, E, ρ)`. Note that the density `ρ` is ignored.
"""
function Chi2Field(θ, ϕ, χ2)
    toCrystal = RotMatrix(RotZY(-ϕ, -θ)) # RotMatrix converts to static matrix
    toLab = RotMatrix(RotYZ(θ, ϕ))
    sv3 = zeros(3)
    sv6 = zeros(6)
    χ2 = SMatrix{3, 6}(χ2) # just χ2
    χ2_toLab = SMatrix{3, 6}(toLab * χ2) # χ2 and coordinate transform in one step
    Chi2Field(χ2, toCrystal, toLab, χ2_toLab, sv3, copy(sv3), sv6, copy(sv3))
end

function (c::Chi2Field)(out, E, ρ)
    for i in axes(E, 1)
        @inbounds c.El[1] = E[i, 1]
        @inbounds c.El[2] = E[i, 2]
        # note c.El[3] (Ez in the lab frame) is always zero here
        mul!(c.Ec, c.toCrystal, c.El) # transform to crystal frame
        @inbounds field_products!(c.Enl, c.Ec) # calculate nonlinear products
        mul!(c.Pl, c.χ2_toLab, c.Enl) # multiply by χ2 tensor and transform to lab frame
        @inbounds out[i, 1] += ε_0*c.Pl[1]
        @inbounds out[i, 2] += ε_0*c.Pl[2]
    end
end


"""
    field_products!(Enl, Ec)

Fill the contracted second-order field-product vector `Enl` from crystal-frame
field components `Ec`.

The output ordering is
`[Ex^2, Ey^2, Ez^2, 2EyEz, 2ExEz, 2ExEy]`, matching the 3x6 `χ2` tensor
column order `[xx, yy, zz, yz, xz, xy]` used by [`Chi2Field`](@ref).

Both `Enl` and `Ec` are mutated/read in place and are expected to have length 6
and 3, respectively.
"""
function field_products!(Enl, Ec)
    Enl[1] = Ec[1]^2
    Enl[2] = Ec[2]^2
    Enl[3] = Ec[3]^2
    Enl[4] = 2*Ec[2]*Ec[3]
    Enl[5] = 2*Ec[1]*Ec[3]
    Enl[6] = 2*Ec[1]*Ec[2]
end

struct Chi2Env{χT}
    χ2::χT
    toCrystal::RotMatrix3{Float64}
    toLab::RotMatrix3{Float64}
    χ2_toLab::χT # combined matrix to multiply by toLab * χ2
    C::Vector{ComplexF64} # carrier phase exp(iω0t) on the (oversampled) time grid
    Al::Vector{ComplexF64} # envelope in the lab frame
    Ac::Vector{ComplexF64} # envelope in the crystal frame
    Anl::Vector{ComplexF64} # combined SFG+DFG envelope products in the crystal frame
    Pl::Vector{ComplexF64} # polarisation envelope in the lab frame
end

"""
    Chi2Env(θ, ϕ, χ2, ω0, t)

Construct a second-order nonlinear polarisation response for complex envelope,
two-component electric fields in the lab frame. Envelope counterpart of
[`Chi2Field`](@ref).

`θ` and `ϕ` (radians) define the crystal orientation relative to the lab frame.
`χ2` must be a 3×6 second-order susceptibility tensor in contracted notation,
with column order `[xx, yy, zz, yz, xz, xy]`. `ω0` is the carrier frequency and `t`
the time axis on which the response is evaluated—for propagation simulations these
must be `grid.ω0` and `grid.to` (the oversampled time axis) of a `Grid.EnvGrid`.

The response includes both the sum-frequency term (``ω + ω → 2ω``, carrying the phase
factor ``e^{+iω_0t}``) and the difference-frequency term (``2ω - ω → ω``, back-conversion,
carrying ``e^{-iω_0t}``), see [`env_products!`](@ref). Optical-rectification content near
zero absolute frequency lies outside the frequency window and is removed by the grid
apodisation. Note that the grid must contain the second harmonic—use
`Grid.EnvGrid(...; thg=true)` or wavelength limits reaching below `λ0/2`.

!!! warning
    The linear operator must use a reference frame which is transparent to carrier-mixing
    nonlinearities, i.e. one whose subtracted phase is strictly linear in the absolute
    frequency (a pure time shift), like the crystal operators
    `LinearOps.make_const_linop(grid, xgrid, nfuns::Tuple)`. Envelope operators which
    subtract the carrier phase `β0` at `grid.ω0` (`thg=false`) introduce a spurious phase
    mismatch `β0 - β1ω_0` into χ⁽²⁾ processes.

The returned callable adds \$ε_0 P_{NL}\$ to `out` when invoked as
`response(out, E, ρ)`. Note that the density `ρ` is ignored.
"""
function Chi2Env(θ, ϕ, χ2, ω0, t)
    toCrystal = RotMatrix(RotZY(-ϕ, -θ)) # RotMatrix converts to static matrix
    toLab = RotMatrix(RotYZ(θ, ϕ))
    χ2 = SMatrix{3, 6}(χ2) # just χ2
    χ2_toLab = SMatrix{3, 6}(toLab * χ2) # χ2 and coordinate transform in one step
    C = exp.(1im*ω0.*t)
    sv3 = zeros(ComplexF64, 3)
    sv6 = zeros(ComplexF64, 6)
    Chi2Env(χ2, toCrystal, toLab, χ2_toLab, C, sv3, copy(sv3), sv6, copy(sv3))
end

function (c::Chi2Env)(out, E, ρ)
    size(E, 2) == 2 || error("Chi2Env requires a two-component (Nt×2) envelope field")
    length(c.C) == size(E, 1) || error(
        "Chi2Env carrier phase array does not match the field length. "
        * "The response must be constructed with the oversampled time axis `grid.to`.")
    for i in axes(E, 1)
        cp = 0.5*c.C[i] # ½exp(iω0t): sum-frequency (ω + ω → 2ω)
        cm = conj(c.C[i]) # exp(-iω0t): difference-frequency (2ω - ω → ω)
        @inbounds c.Al[1] = E[i, 1]
        @inbounds c.Al[2] = E[i, 2]
        # note c.Al[3] (Az in the lab frame) is always zero here
        mul!(c.Ac, c.toCrystal, c.Al) # transform to crystal frame
        @inbounds env_products!(c.Anl, c.Ac, cp, cm) # calculate nonlinear products
        mul!(c.Pl, c.χ2_toLab, c.Anl) # multiply by χ2 tensor and transform to lab frame
        @inbounds out[i, 1] += ε_0*c.Pl[1]
        @inbounds out[i, 2] += ε_0*c.Pl[2]
    end
end

"""
    env_products!(Anl, Ac, cp, cm)

Fill the contracted second-order envelope-product vector `Anl` from crystal-frame
envelope components `Ac`. Envelope counterpart of [`field_products!`](@ref).

With the real field given by ``E_j = \\mathrm{Re}[A_j e^{iω_0t}]``, the envelope of each
second-order product is

``(E_jE_k)_{env} = \\frac{1}{2}A_jA_k e^{+iω_0t} + \\frac{1}{2}(A_jA_k^* + A_j^*A_k)e^{-iω_0t}``

where the first (sum-frequency) and second (difference-frequency) terms enter here via
`cp` = ``\\frac{1}{2}e^{+iω_0t}`` and `cm` = ``e^{-iω_0t}`` at the current time sample.

The output ordering is `[xx, yy, zz, yz, xz, xy]` with mixed terms multiplied by 2,
matching the 3x6 `χ2` tensor column order used by [`Chi2Env`](@ref).

Both `Anl` and `Ac` are mutated/read in place and are expected to have length 6
and 3, respectively.
"""
function env_products!(Anl, Ac, cp, cm)
    Anl[1] = cp*Ac[1]^2 + cm*abs2(Ac[1])
    Anl[2] = cp*Ac[2]^2 + cm*abs2(Ac[2])
    Anl[3] = cp*Ac[3]^2 + cm*abs2(Ac[3])
    Anl[4] = 2*(cp*Ac[2]*Ac[3] + cm*real(Ac[2]*conj(Ac[3])))
    Anl[5] = 2*(cp*Ac[1]*Ac[3] + cm*real(Ac[1]*conj(Ac[3])))
    Anl[6] = 2*(cp*Ac[1]*Ac[2] + cm*real(Ac[1]*conj(Ac[2])))
end

"""
    PlasmaCumtrapz(t, E, ratefunc, ionpot; preionfrac=0.0)

Cumulative-trapezoid plasma polarisation response, adapted from M. Geissler, G. Tempea,
A. Scrinzi, M. Schnürer, F. Krausz, and T. Brabec, Physical Review Letters 83, 2930
(1999).

`t` is the (oversampled) time grid, `E` a prototype of the time-domain field block the
response will be applied to, `ratefunc` an ionisation rate (see
`Ionisation`) and `ionpot` the ionisation potential in Joules.
`preionfrac` is a pre-ionised fraction, which is not a well founded physical model.

[`Batched`](@ref): the response is handed the whole `(nt, npol, ncols...)` block and
evaluates it with whole-array operations — one broadcast for the rate, one prefix scan
plus one broadcast for each of the three cumulative integrals
(`Maths.cumtrapz_scan!`), and one `ifelse` broadcast
for the ionisation-loss term. That is the same code on the host and on a device
(GPU_PLAN.md §4.2). On the host the columns are shared out over threads when there are
enough of them; each column's arithmetic is the same either way, and the same as if the
block were passed one column at a time.

The buffers are sized for the block, so a transform reallocates them for the block it
will pass by calling [`rescale`](@ref), which is also where the rate function is
converted to the run's precision and array type.

Its result differs from a serial `Maths.cumtrapz!` at rounding level: the scan
accumulates the integrand and corrects, where the loop accumulates the trapezoid
increments.

# Fields
- `ratefunc`: the rate as given, in physical units on the host. `Stats` evaluates it.
- `ratedev`: the same rate in the run's precision and array type; `=== ratefunc` on the
  default host path.
- `rate`, `fraction`, `Em`: buffers with one polarisation component (`Em`, the field
  magnitude, only for a two-component field; `nothing` otherwise).
- `J`, `P`: block-sized buffers. `P` holds the phase-modulation term before it holds the
  polarisation; the two are never needed at the same time.
"""
struct PlasmaCumtrapz{R, RD, tType, EType, mType}
    ratefunc::R # the ionisation rate function, as given (host, physical units)
    ratedev::RD # the same rate in the run's precision and array type
    ionpot::Float64 # the ionisation potential (for the ionisation loss term)
    rate::tType # buffer to hold the rate
    fraction::tType # buffer to hold the ionisation fraction
    Em::mType # buffer for the field magnitude (two-component field), or nothing
    J::EType # buffer to hold the plasma current
    P::EType # buffer to hold the phase modulation, then the plasma polarisation
    δt::Float64 # the time step
    preionfrac::Float64 # the pre-ionisation fraction
    scaling::Luna.UnitScaling # units the field block and the output are expressed in
end

function PlasmaCumtrapz(t, E, ratefunc, ionpot; preionfrac=0.0)
    !(0.0 <= preionfrac <= 1.0) && throw(DomainError(preionfrac, "preionfrac must be between 0 and 1"))
    if preionfrac > 0.0
        @warn("Using preionfrac > 0.0 is not a well founded physical model. Use only after careful consideration.")
    end
    rate, fraction, Em, J, P = _plasmabuffers(E)
    PlasmaCumtrapz(ratefunc, ratefunc, ionpot, rate, fraction, Em, J, P,
                   t[2]-t[1], preionfrac, Luna.UNIT_SCALING)
end

#= The rate, the ionisation fraction and the field magnitude are the same for both
   polarisation components, so they are held with a singleton second axis and broadcast
   against the block. A one-dimensional block (mode-averaged) has no such axis. =#
_singletonpol(dims::Tuple{}) = ()
_singletonpol(dims::Tuple{Int}) = dims
_singletonpol(dims::Tuple) = (dims[1], 1, dims[3:end]...)

function _plasmabuffers(E)
    RT = real(eltype(E))
    pdims = _singletonpol(size(E))
    rate = fill!(similar(E, RT, pdims), zero(RT))
    fraction = similar(rate)
    Em = _npol(E) == 2 ? similar(rate) : nothing
    J = fill!(similar(E), zero(eltype(E)))
    P = similar(J)
    (rate, fraction, Em, J, P)
end

"The number of polarisation components of a field block."
_npol(E::AbstractArray) = ndims(E) < 2 ? 1 : size(E, 2)

"The number of columns (transverse points) of a field block."
_ncols(E::AbstractArray) = ndims(E) < 3 ? 1 : prod(size(E)[3:end])

kind(::PlasmaCumtrapz) = Batched()
kind(::PlasmaCumtrapz, ::Val) = Batched()

resident_arrays(p::PlasmaCumtrapz) =
    (p.rate, p.fraction, p.Em, p.J, p.P,
     Luna.Ionisation.resident_arrays(p.ratedev)...)

"""
    coefficients(p::PlasmaCumtrapz, ρ, scaling)

The four scalars the plasma kernels need at number density `ρ`, in the units of
`scaling`: `(Eref, cphase, closs, cout)`.

The response is not polynomial in the field, so [`Luna.polscale`](@ref) does not apply
and the scaling is carried by each stage separately. With `e = E/E_ref` and the output
in `P/(P_ref E_ref)`:

| | |
| --- | --- |
| `Eref` | the rate is evaluated at the physical field `Eref*e` |
| `cphase` | `e_ratio*Eref`, so that `fraction*cphase*e` is the physical phase term |
| `closs` | `ionpot/Eref`; the loss term divides by the physical field (twice, for a two-component field, and multiplies by it once) |
| `cout` | `ρ/(P_ref E_ref)`, applied to the physical polarisation |

All four are `Float64` here and converted once, by the kernel, with
[`Luna.scalar`](@ref). Every one of them reduces to the physical constant when
`E_ref == P_ref == 1`, which is every `Float64` run.
"""
function coefficients(p::PlasmaCumtrapz, ρ, scaling)
    Eref = scaling.Eref
    (Eref, e_ratio*Eref, p.ionpot/Eref, ρ/(scaling.Pref*Eref))
end

"""
    rescale(p::PlasmaCumtrapz, spec, scaling, Et)

A copy of the plasma response with buffers sized for the block `Et`, in its array type
and precision, and with the ionisation rate converted to match
([`Ionisation.device_rate`](@ref Luna.Ionisation.device_rate)).

Every transform calls this on its responses at construction, so the buffers always match
the block the response is handed — which for a radial or free-space transform is the
whole transverse grid, not the single column the constructor was given a prototype of.

The rate as given is kept in `ratefunc` (the host, physical-unit object `Stats`
evaluates) whatever the run.
"""
function rescale(p::PlasmaCumtrapz, spec, scaling, Et)
    ratedev = Luna.Ionisation.device_rate(p.ratefunc, spec)
    rate, fraction, Em, J, P = _plasmabuffers(Et)
    Luna.assert_resident(spec, rate, fraction, Em, J, P)
    PlasmaCumtrapz(p.ratefunc, ratedev, p.ionpot, rate, fraction, Em, J, P,
                   p.δt, p.preionfrac, scaling)
end

#= Below this many elements the block is not worth the task overhead; above it, one task
   per column. A column costs an ionisation-rate evaluation and three prefix scans over
   the time axis, so the grain is large even for one column. =#
const PLASMA_THREAD_MINLEN = 1 << 14

"Whether to share this block's columns out over threads."
_plasma_threaded(E) =
    (Threads.nthreads() > 1) && !Utils.isdevice(E) &&
    (_ncols(E) > 1) && (length(E) >= PLASMA_THREAD_MINLEN)

#= The dispatcher hands a batched response the units the block is in, which is what a
   transform is running in; `p.scaling` is only what `rescale` was told, and is what a
   direct call falls back on. =#
batched!(p::PlasmaCumtrapz, out, E, ρ, scaling) =
    _plasma_run!(p, out, E, coefficients(p, ρ, scaling))

(p::PlasmaCumtrapz)(out, Et, ρ) =
    _plasma_run!(p, out, Et, coefficients(p, ρ, p.scaling))

function _plasma_run!(p::PlasmaCumtrapz, out, Et, c)
    size(out) == size(Et) || throw(DimensionMismatch(
        "PlasmaCumtrapz: output block is $(size(out)), field block is $(size(Et))"))
    size(Et) == size(p.J) || error(
        "PlasmaCumtrapz was given a $(join(size(Et), "x")) field block but its buffers "*
        "are $(join(size(p.J), "x")). A batched response is handed the whole block, so "*
        "its buffers have to match it: call "*
        "`Nonlinear.rescale(response, spec, scaling, Et)` with a prototype of the block "*
        "(every transform does this at construction).")
    _npol(Et) in (1, 2) || error(
        "PlasmaCumtrapz: a field block has one or two polarisation components along "*
        "dimension 2, not $(_npol(Et)).")
    if _plasma_threaded(Et)
        nc = _ncols(Et)
        #= Host arrays only, so the reshapes are free. Each task gets one column of
           every buffer; the columns are independent, so the arithmetic in each is the
           same as it would be in one call over the whole block. =#
        o3, E3, rate3, frac3, Em3, J3, P3 =
            map(x -> _as3d(x, nc), (out, Et, p.rate, p.fraction, p.Em, p.J, p.P))
        Threads.@threads :dynamic for i in 1:nc
            _plasma_block!(_col(o3, i), _col(E3, i), _col(rate3, i), _col(frac3, i),
                           _col(Em3, i), _col(J3, i), _col(P3, i),
                           p.ratedev, p.δt, p.preionfrac, c)
        end
    else
        _plasma_block!(out, Et, p.rate, p.fraction, p.Em, p.J, p.P,
                       p.ratedev, p.δt, p.preionfrac, c)
    end
    out
end

_as3d(x::AbstractArray, nc) = reshape(x, size(x, 1), _npol(x), nc)
_as3d(::Nothing, nc) = nothing

_col(x::AbstractArray, i) = view(x, :, :, i:i)
_col(::Nothing, i) = nothing

#= The whole of the plasma response, for a block of any shape and on any array type: one
   broadcast for the rate, three (scan + broadcast) cumulative integrals, one `ifelse`
   broadcast for the ionisation loss and one for the output. `P` is the phase-modulation
   buffer up to the last scan and the polarisation afterwards.

   The ionisation rate, the fraction and the field magnitude have a singleton
   polarisation axis and broadcast against both components. =#
function _plasma_block!(out, E, rate, fraction, Em, J, P, ratefunc, δt, preionfrac, c)
    Eref, cphase, closs, cout = c
    _ratefield!(rate, Em, E, ratefunc, Eref)
    Maths.cumtrapz_scan!(fraction, rate, δt)
    pf = Luna.scalar(E, preionfrac)
    @. fraction = pf + 1 - exp(-fraction)
    cp = Luna.scalar(E, cphase)
    @. P = fraction * cp * E
    Maths.cumtrapz_scan!(J, P, δt)
    _plasma_loss!(J, E, Em, rate, fraction, Luna.scalar(E, closs))
    Maths.cumtrapz_scan!(P, J, δt)
    co = Luna.scalar(E, cout)
    @. out += co * P
    out
end

# Scalar field: the rate is a function of the field itself.
_ratefield!(rate, ::Nothing, E, ratefunc, Eref) =
    Luna.Ionisation.ionrate!(rate, ratefunc, E, Eref)

#= Two-component field: the magnitude of the field drives the ionisation, and the plasma
   polarisation is then solved component by component. See C Tailliez et al 2020 New J.
   Phys. 22 103038. =#
function _ratefield!(rate, Em, E, ratefunc, Eref)
    Ex = view(E, :, 1:1, ntuple(_ -> Colon(), ndims(E)-2)...)
    Ey = view(E, :, 2:2, ntuple(_ -> Colon(), ndims(E)-2)...)
    @. Em = hypot(Ex, Ey)
    Luna.Ionisation.ionrate!(rate, ratefunc, Em, Eref)
end

#= The ionisation-loss term, as one `ifelse` broadcast rather than the branch of a loop.
   Both arms are evaluated, so the division happens where the field is zero too; `ifelse`
   is a select, not a branch, so the infinity it produces there is discarded rather than
   propagated.

   The guard is on the denominator, which for a two-component field is `Em^2` and not
   `Em`. In the wings of a pulse `Em` is many orders below the peak -- `exp(-50)` is
   1e-22 -- and its square underflows in Float32 where `Em` itself does not. A device
   flushes that subnormal to zero, and since the rate there is zero too the term becomes
   `0/0`, one NaN of which poisons the whole column through the scans which follow. In
   Float64 the two conditions differ only below 1e-162 V/m, which is not a field. =#
_plasma_loss!(J, E, ::Nothing, rate, fraction, closs) =
    @. J += ifelse(abs(E) > 0, closs*rate*(1-fraction)/E, zero(E))

function _plasma_loss!(J, E, Em, rate, fraction, closs)
    @. J += ifelse(Em^2 > 0, closs*rate*(1-fraction)/Em^2*E, zero(E))
end

"Raman polarisation response type"
abstract type RamanPolar end

"Raman polarisation response type for a carrier resolved field"
struct RamanPolarField{TR, Tt, Thv, Tω, Tv, FTt, HTt} <: RamanPolar
    r::TR # Raman response
    h::Tt # doubled buffer to hold response + padding
    ht::Thv # buffer to hold time domain response
    hω::Tω # the frequency domain Raman response function
    Eω2::Tω # buffer to hold the Fourier transform of E^2
    Pω::Tω # buffer to hold the frequency domain polarisation
    E2::Tt # buffer to hold E^2
    E2v::Tv # view into first half of E2
    P::Tt # buffer to hold the time domain polarisation
    Pout::Tt # buffer to hold the output portion of the time domain polarisation
    FT::FTt # Fourier transform plan
    HT::HTt # Hilbert transform
    thg::Bool # do we include third harmonic generation
    dt::Float64 # time step for scaling
end

"Raman polarisation response type for an envelope"
struct RamanPolarEnv{TR, Tt, Thv, Tω, Tv, FTt} <: RamanPolar
    r::TR # Raman response
    h::Tt # doubled buffer to hold response + padding
    ht::Thv # buffer to hold time domain response
    hω::Tω # the frequency domain Raman response function
    Eω2::Tω # buffer to hold the Fourier transform of E^2
    Pω::Tω # buffer to hold the frequency domain polarisation
    E2::Tω # buffer to hold E^2
    E2v::Tv # view into first half of E2
    P::Tω # buffer to hold the time domain polarisation
    Pout::Tω # buffer to hold the output portion of the time domain polarisation
    FT::FTt # Fourier transform plan
    dt::Float64 # time step for scaling
end

"""
    RamanPolarField(t, ht; thg=true)

Construct Raman polarisation response for a field on time grid `t`
using response function `r`. If `thg=false` then exclude the third
harmonic generation component of the response.
"""
function RamanPolarField(t, r; thg=true)
    h = zeros(length(t)*2) # note double grid size, see explanation below
    ht = view(h, 1:length(t))
    Utils.loadFFTwisdom()
    FT = FFTW.plan_rfft(h, 1, flags=Luna.settings["fftw_flag"])
    inv(FT)
    Utils.saveFFTwisdom()
    hω = FT * h
    Eω2 = similar(hω)
    Pω = similar(hω)
    E2 = similar(h)
    E2v = view(E2, 1:length(t))
    P = similar(h)
    Pout = similar(t)
    HT = Maths.plan_hilbert(Pout)
    fill!(E2, 0.0)
    RamanPolarField(r, h, ht, hω, Eω2, Pω, E2, E2v, P, Pout, FT, HT, thg, t[2] - t[1])
end

"""
    RamanPolarEnv(t, ht)

Construct Raman polarisation response for an envelope on time grid `t`
using response function `r`.
"""
function RamanPolarEnv(t, r)
    h = zeros(length(t)*2) # note double grid size, see explanation below
    ht = view(h, 1:length(t))
    Utils.loadFFTwisdom()
    FT = FFTW.plan_fft(h, 1, flags=Luna.settings["fftw_flag"])
    inv(FT)
    Utils.saveFFTwisdom()
    hω = FT * h
    Eω2 = similar(hω)
    Pω = similar(hω)
    E2 = similar(hω)
    P = similar(hω)
    Pout = Array{ComplexF64,}(undef,size(t))
    E2v = view(E2, 1:length(t))
    fill!(E2, 0.0)
    RamanPolarEnv(r, h, ht, hω, Eω2, Pω, E2, E2v, P, Pout, FT, t[2] - t[1])
end

"Square the field or envelope"
function sqr!(R::RamanPolarField, E)
    if !R.thg
        # see documentation for factor of 1/2 here
        R.E2v .= 1/2 .* abs2.(R.HT(E))
    else
        R.E2v .= E.^2
    end
end

function sqr!(R::RamanPolarEnv, E)
    # see documentation for factor of 1/2 here
    R.E2v .= 1/2 .* abs2.(E)
end

"Calculate Raman polarisation for field/envelope Et"
function (R::RamanPolar)(out, Et, ρ)
    # get the field as a 1D Array
    n = size(Et, 1)
    if ndims(Et) > 1
        if size(Et, 2) == 1 # handle scalar case but within modal simulation
            E = reshape(Et, n)
        else
            # handle vector case
            error("vector Raman not yet implemented")
        end
    else
        E = Et # handle straight scalar case
    end

    # square the field or envelope in first half
    # corresponding to the field/envelope grid size
    sqr!(R, E)

    # update frequency domain response function `hω`.
    # we fill only up to the first half of h (using the view ht)
    # i.e. only the part corresponding to the original time grid
    # note that the response function time 0 is put into the first element of the response array
    # this ensures that causality is maintained, and no artificial delay between the field and
    # the start of the response function occurs, at each convolution point.
    R.r(R.ht, ρ)
    R.hω .= R.FT * R.h

    # convolution by multiplication in frequency domain
    # The double grid gives us accurate full convolution between the full field grid
    # and full response function. It is unnecessary for highly damped responses, like
    # in glass. But for gases with very long decay times it prevents artefacts due to
    # truncation of the response function. There is likely a more efficient way. But
    # this is safe, until we come up with one.
    # we scale to correct for missing dt*dt*df from IFFT(FFT*FFT)
    # the ifft already scales by 1/n = dt*df, so we need an additional dt
    R.Eω2 .= R.FT * R.E2
    @. R.Pω = R.hω * R.Eω2 * R.dt
    R.P .= R.FT \ R.Pω

    # calculate full polarisation, extracting only the valid
    # grid region, which is the first length(E) part.
    for i = 1:length(E)
        R.Pout[i] = ρ*E[i]*R.P[i]
    end

    # copy to output in dimensions requested
    if ndims(Et) > 1
        out .+= reshape(R.Pout, size(Et))
    else
        out .+= R.Pout
    end
end

end
