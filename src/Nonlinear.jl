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

#=================================================#
#===========  SECOND-ORDER RESPONSES  ============#
#=================================================#

#= The χ⁽²⁾ responses are [`VectorPointwise`](@ref): the contraction at one time sample
   couples the two lab-frame polarisation components and nothing else. The matrices the
   kernel contracts with are `StaticArrays`/`Rotations` matrices, which are `isbits` and
   travel inside the closure a broadcast compiles, so no buffer and no host array is
   involved -- that is what makes them run wherever the block does. =#

"""
    Chi2Field(θ, ϕ, χ2)

Second-order nonlinear polarisation response for real, two-component (field-resolved)
electric fields in the lab frame.

`θ` and `ϕ` (radians) define the crystal orientation relative to the lab frame.
`χ2` must be a 3×6 second-order susceptibility tensor in contracted notation,
with column order `[xx, yy, zz, yz, xz, xy]` (where mixed terms are multiplied by 2
in [`field_products!`](@ref)).

The returned callable adds \$ε_0 P_{NL}\$ to `out` when invoked as
`response(out, E, ρ)`. Note that the density `ρ` is ignored.

[`VectorPointwise`](@ref): it needs a two-component field block and errors on a
one-component one.
"""
struct Chi2Field{T}
    χ2::SMatrix{3, 6, T, 18}
    toCrystal::RotMatrix3{T}
    toLab::RotMatrix3{T}
    χ2_toLab::SMatrix{3, 6, T, 18} # combined matrix toLab * χ2
end

Chi2Field(θ, ϕ, χ2) = Chi2Field(_chi2matrices(θ, ϕ, χ2)...)

#= The crystal matrices, promoted to one real type so that the struct has a single
   parameter and the kernel a single conversion. `toLab*χ2` is formed before the
   conversion, as it always was, so that the default Float64 path is unchanged. =#
function _chi2matrices(θ, ϕ, χ2)
    toCrystal = RotMatrix(RotZY(-ϕ, -θ)) # RotMatrix converts to static matrix
    toLab = RotMatrix(RotYZ(θ, ϕ))
    χ2 = SMatrix{3, 6}(χ2) # just χ2
    χ2_toLab = SMatrix{3, 6}(toLab * χ2) # χ2 and coordinate transform in one step
    T = promote_type(eltype(χ2), eltype(toCrystal), eltype(toLab), eltype(χ2_toLab))
    (SMatrix{3, 6, T}(χ2), RotMatrix3{T}(toCrystal), RotMatrix3{T}(toLab),
     SMatrix{3, 6, T}(χ2_toLab))
end

kind(::Chi2Field) = VectorPointwise()

coefficients(::Chi2Field, ρ, scaling) = ε_0*Luna.polscale(scaling, 2)

#= The rotation and the contracted tensor in the real element type of the block. A
   response is built from crystal data in Float64 and `rescale` converts it for a
   reduced-precision run, but the conversion is repeated here -- it is 27 host scalars
   once per right-hand side -- so that a kernel built from a response which was never
   rescaled still carries no Float64 (GPU_PLAN.md section 4.2 rule 3). It is the identity
   when the response is already in the block's precision. =#
_chi2mats(::Type{T}, c) where {T} =
    (SMatrix{3, 3, T}(c.toCrystal), SMatrix{3, 6, T}(c.χ2_toLab))

#= The arithmetic, once. The fused component broadcasts and the columnwise call operator
   below both go through this, so there is one kernel body per response. =#
_chi2field(fac, toCrystal, χ2_toLab) = (ex, ey) -> begin
    # the third lab-frame component (Ez) is always zero
    Ec = toCrystal*SVector(ex, ey, zero(ex)) # transform to crystal frame
    Enl = _field_products(Ec) # calculate nonlinear products
    Pl = χ2_toLab*Enl # multiply by χ2 tensor and transform to lab frame
    SVector(fac*Pl[1], fac*Pl[2])
end

vector_kernel(c::Chi2Field, E, ρ, scaling) =
    _chi2field(Luna.scalar(E, coefficients(c, ρ, scaling)),
               _chi2mats(real(eltype(E)), c)...)

#= The columnwise contract, in physical units: what a direct caller still reaches. Like
   the dispatcher's vector path it evaluates the kernel once per output component; see
   `docs/src/developer/device_model.md`. =#
function (c::Chi2Field)(out, E, ρ)
    size(E, 2) == 2 || error("Chi2Field requires a two-component (Nt×2) field")
    f = _chi2field(Luna.scalar(E, coefficients(c, ρ, Luna.UNIT_SCALING)),
                   _chi2mats(real(eltype(E)), c)...)
    Ex = selectdim(E, 2, 1)
    Ey = selectdim(E, 2, 2)
    ox = selectdim(out, 2, 1)
    oy = selectdim(out, 2, 2)
    @. ox += first(f(Ex, Ey))
    @. oy += last(f(Ex, Ey))
    out
end

function rescale(c::Chi2Field, spec, scaling)
    _isdefaultrun(spec, scaling) && return c
    _chi2rescale(Luna.realtype(spec), c)
end

_chi2rescale(::Type{T}, c::Chi2Field) where {T} =
    Chi2Field(SMatrix{3, 6, T}(c.χ2), RotMatrix3{T}(c.toCrystal),
              RotMatrix3{T}(c.toLab), SMatrix{3, 6, T}(c.χ2_toLab))

#= Every field is an `isbits` static matrix, so `Adapt` has nothing to move; the rule
   exists so that adapting a container which holds a response is not a no-op by accident.
   Precision is `rescale`'s job: `Adapt` does not know the target element type. =#
Adapt.adapt_structure(to, c::Chi2Field) =
    Chi2Field(Adapt.adapt(to, c.χ2), Adapt.adapt(to, c.toCrystal),
              Adapt.adapt(to, c.toLab), Adapt.adapt(to, c.χ2_toLab))

#= The contracted second-order field products, as an `SVector{6}`: the one body, which
   both the kernel and `field_products!` use. =#
_field_products(Ec) = SVector(Ec[1]^2, Ec[2]^2, Ec[3]^2,
                              2*Ec[2]*Ec[3], 2*Ec[1]*Ec[3], 2*Ec[1]*Ec[2])

"""
    field_products!(Enl, Ec)

Fill the contracted second-order field-product vector `Enl` from crystal-frame
field components `Ec`.

The output ordering is
`[Ex^2, Ey^2, Ez^2, 2EyEz, 2ExEz, 2ExEy]`, matching the 3x6 `χ2` tensor
column order `[xx, yy, zz, yz, xz, xy]` used by [`Chi2Field`](@ref).

Both `Enl` and `Ec` are mutated/read in place and are expected to have length 6
and 3, respectively. This is the out-of-place expression the kernel uses, written into a
vector; it is not what the response itself calls.
"""
field_products!(Enl, Ec) = (Enl .= _field_products(SVector{3}(Ec)); Enl)

"""
    Chi2Env(θ, ϕ, χ2, ω0, t)

Second-order nonlinear polarisation response for complex envelope, two-component
electric fields in the lab frame. Envelope counterpart of [`Chi2Field`](@ref).

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

[`VectorPointwise`](@ref): it needs a two-component field block and errors on a
one-component one. The carrier phase is a second broadcast argument rather than something
the kernel captures, so it is aligned with the time axis whatever the block's shape.
"""
struct Chi2Env{T, cT}
    χ2::SMatrix{3, 6, T, 18}
    toCrystal::RotMatrix3{T}
    toLab::RotMatrix3{T}
    χ2_toLab::SMatrix{3, 6, T, 18} # combined matrix to multiply by toLab * χ2
    C::cT # carrier phase exp(iω0t) on the (oversampled) time grid
end

Chi2Env(θ, ϕ, χ2, ω0, t) = Chi2Env(_chi2matrices(θ, ϕ, χ2)..., exp.(1im*ω0.*t))

kind(::Chi2Env) = VectorPointwise()

coefficients(::Chi2Env, ρ, scaling) = ε_0*Luna.polscale(scaling, 2)

_chi2env(fac, toCrystal, χ2_toLab) = (ax, ay, C) -> begin
    cp = C/2 # ½exp(iω0t): sum-frequency (ω + ω → 2ω)
    cm = conj(C) # exp(-iω0t): difference-frequency (2ω - ω → ω)
    # the third lab-frame component (Az) is always zero
    Ac = toCrystal*SVector(ax, ay, zero(ax)) # transform to crystal frame
    Anl = _env_products(Ac, cp, cm) # calculate nonlinear products
    Pl = χ2_toLab*Anl # multiply by χ2 tensor and transform to lab frame
    SVector(fac*Pl[1], fac*Pl[2])
end

vector_expr(c::Chi2Env, Ex, Ey, ρ, scaling) =
    Base.broadcasted(_chi2envkernel(c, Ex, ρ, scaling), Ex, Ey, c.C)

function _chi2envkernel(c::Chi2Env, E, ρ, scaling)
    size(E, 1) == length(c.C) || error(
        "Chi2Env carrier phase array does not match the field length. "
        * "The response must be constructed with the oversampled time axis `grid.to`.")
    _chi2env(Luna.scalar(E, coefficients(c, ρ, scaling)),
             _chi2mats(real(eltype(E)), c)...)
end

function (c::Chi2Env)(out, E, ρ)
    size(E, 2) == 2 || error("Chi2Env requires a two-component (Nt×2) envelope field")
    Ex = selectdim(E, 2, 1)
    Ey = selectdim(E, 2, 2)
    f = _chi2envkernel(c, Ex, ρ, Luna.UNIT_SCALING)
    C = c.C
    ox = selectdim(out, 2, 1)
    oy = selectdim(out, 2, 2)
    @. ox += first(f(Ex, Ey, C))
    @. oy += last(f(Ex, Ey, C))
    out
end

function rescale(c::Chi2Env, spec, scaling)
    _isdefaultrun(spec, scaling) && return c
    _chi2rescale(Luna.realtype(spec), c, Luna.todevice(spec, c.C))
end

_chi2rescale(::Type{T}, c::Chi2Env, C) where {T} =
    Chi2Env(SMatrix{3, 6, T}(c.χ2), RotMatrix3{T}(c.toCrystal),
            RotMatrix3{T}(c.toLab), SMatrix{3, 6, T}(c.χ2_toLab), C)

Adapt.adapt_structure(to, c::Chi2Env) =
    Chi2Env(Adapt.adapt(to, c.χ2), Adapt.adapt(to, c.toCrystal),
            Adapt.adapt(to, c.toLab), Adapt.adapt(to, c.χ2_toLab),
            Adapt.adapt(to, c.C))

resident_arrays(c::Chi2Env) = (c.C,)

#= The contracted second-order envelope products, as an `SVector{6}`: the one body, which
   both the kernel and `env_products!` use. =#
_env_products(Ac, cp, cm) = SVector(
    cp*Ac[1]^2 + cm*abs2(Ac[1]),
    cp*Ac[2]^2 + cm*abs2(Ac[2]),
    cp*Ac[3]^2 + cm*abs2(Ac[3]),
    2*(cp*Ac[2]*Ac[3] + cm*real(Ac[2]*conj(Ac[3]))),
    2*(cp*Ac[1]*Ac[3] + cm*real(Ac[1]*conj(Ac[3]))),
    2*(cp*Ac[1]*Ac[2] + cm*real(Ac[1]*conj(Ac[2]))))

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
and 3, respectively. This is the out-of-place expression the kernel uses, written into a
vector; it is not what the response itself calls.
"""
env_products!(Anl, Ac, cp, cm) = (Anl .= _env_products(SVector{3}(Ac), cp, cm); Anl)

"Response type for cumtrapz-based plasma polarisation, adapted from:
M. Geissler, G. Tempea, A. Scrinzi, M. Schnürer, F. Krausz, and T. Brabec, Physical Review Letters 83, 2930 (1999)."
struct PlasmaCumtrapz{R, EType, tType}
    ratefunc::R # the ionization rate function
    ionpot::Float64 # the ionization potential (for calculation of ionization loss)
    rate::tType # buffer to hold the rate
    fraction::tType # buffer to hold the ionization fraction
    phase::EType # buffer to hold the plasma induced (mostly) phase modulation
    J::EType # buffer to hold the plasma current
    P::EType # buffer to hold the plasma polarisation
    δt::Float64 # the time step
    preionfrac::Float64 # the pre-ionisation fraction
end

"""
    PlasmaCumtrapz(t, E, ratefunc, ionpot)

Construct the Plasma polarisation response for a field on time grid `t`
with example electric field like `E`, an ionization rate callable
`ratefunc` and ionization potential `ionpot`.
"""
function PlasmaCumtrapz(t, E, ratefunc, ionpot; preionfrac=0.0)
    rate = similar(t)
    fraction = similar(t)
    phase = similar(E)
    J = similar(E)
    P = similar(E)
    !(0.0 <= preionfrac <= 1.0) && throw(DomainError(preionfrac, "preionfrac must be between 0 and 1"))
    if preionfrac > 0.0
        @warn("Using preionfrac > 0.0 is not a well founded physical model. Use only after careful consideration.")
    end
    return PlasmaCumtrapz(ratefunc, ionpot, rate, fraction, phase, J, P, t[2]-t[1], preionfrac)
end

"The plasma response for a scalar electric field"
function PlasmaScalar!(Plas::PlasmaCumtrapz, E)
    Plas.ratefunc(Plas.rate, E)
    Maths.cumtrapz!(Plas.fraction, Plas.rate, Plas.δt)
    @. Plas.fraction = Plas.preionfrac + 1 - exp(-Plas.fraction)
    @. Plas.phase = Plas.fraction * e_ratio * E
    Maths.cumtrapz!(Plas.J, Plas.phase, Plas.δt)
    for ii in eachindex(E)
        if abs(E[ii]) > 0
            Plas.J[ii] += Plas.ionpot * Plas.rate[ii] * (1-Plas.fraction[ii])/E[ii]
        end
    end
    Maths.cumtrapz!(Plas.P, Plas.J, Plas.δt)
end

"""
The plasma response for a vector electric field.

We take the magnitude of the electric field to calculate the ionization
rate and fraction, and then solve the plasma polarisation component-wise
for the vector field.

A similar approach was used in: C Tailliez et al 2020 New J. Phys. 22 103038.
"""
function PlasmaVector!(Plas::PlasmaCumtrapz, E)
    Ex = E[:,1]
    Ey = E[:,2]
    Em = @. hypot.(Ex, Ey)
    Plas.ratefunc(Plas.rate, Em)
    Maths.cumtrapz!(Plas.fraction, Plas.rate, Plas.δt)
    @. Plas.fraction = Plas.preionfrac + 1 - exp(-Plas.fraction)
    @. Plas.phase = Plas.fraction * e_ratio * E
    Maths.cumtrapz!(Plas.J, Plas.phase, Plas.δt)
    for ii in eachindex(Em)
        if abs(Em[ii]) > 0
            pre = Plas.ionpot * Plas.rate[ii] * (1-Plas.fraction[ii])/Em[ii]^2
            Plas.J[ii,1] += pre*Ex[ii]
            Plas.J[ii,2] += pre*Ey[ii]
        end
    end
    Maths.cumtrapz!(Plas.P, Plas.J, Plas.δt)
end

"Handle plasma polarisation routing to `PlasmaVector` or `PlasmaScalar`."
function (Plas::PlasmaCumtrapz)(out, Et, ρ)
    if ndims(Et) > 1
        if size(Et, 2) == 1 # handle scalar case but within modal simulation
            PlasmaScalar!(Plas, reshape(Et, size(Et,1)))
            out .+= ρ .* reshape(Plas.P, size(Et))
        else
            PlasmaVector!(Plas, Et) # vector case
            out .+= ρ .* Plas.P
        end
    else
        PlasmaScalar!(Plas, Et) # straight scalar case
        out .+= ρ .* Plas.P
    end
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
