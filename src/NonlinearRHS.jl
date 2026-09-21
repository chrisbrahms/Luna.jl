"""
    NonlinearRHS

The right-hand side of the propagation equation: transforms which take the frequency-domain
field `E(ω)` to the nonlinear polarisation `Pₙₗ(ω)` for each modal decomposition
`Luna` supports.

Each transform also band-limits `Pₙₗ` to the simulation band. That is part of the
definition of the equation being solved and is separate from the absorbing boundaries
applied to the field itself, which live in [`Boundaries`](@ref Luna.Boundaries).
"""
module NonlinearRHS
import FFTW
import Cubature
import Base: show
import LinearAlgebra: mul!, ldiv!
import NumericalIntegration: integrate, SimpsonEven
import Luna: PhysData, Modes, Maths, Grid, Utils, Nonlinear
import Luna: DeviceSpec, HostSpec, UnitScaling, UNIT_SCALING, GridVectors, HostMirror
import Luna: alloc, todevice, gridvectors, assert_resident, all_resident, scalar,
             upload!, realtype, arraytype, isdevicespec, isunity
import Luna.PhysData: wlfreq, c, crystal_internal_angle
import Luna: LinearOps
import Luna.LinearOps: βz, transverse_k2
import Logging
using EllipsisNotation

"""
    to_time!(Ato, Aω, Aωo, IFT)

Transform ``A(ω)`` on the normal grid to ``A(t)`` on the oversampled time grid, with the
explicit inverse plan `IFT` (see [`Utils.plan_ift`](@ref Luna.Utils.plan_ift)).

The plan's `1/N` normalisation is folded into the scale factor of the oversampling copy
rather than applied as a separate pass, which is one pass fewer over the oversampled
array. Where `1/N` is a power of two the folding is exact -- scaling by a power of two is
exact, and an exactly scaled FFT input gives an exactly scaled output -- so a transform
over the time axis alone changes nothing, Luna's time grids being powers of two. A
multi-axis transform normalises by `1/(Nt·Nx·Ny)`, and the free-space grids accept any
`Nx`, `Ny`; where those are not powers of two the folded and unfolded routes differ at
rounding level instead (measured: 5.6e-17 for a length-24 transform, 0 for 16 or 32).

Dispatches on the element type of `Ato`, not on its array type, so the same code runs on
the host and on a device. For a field-resolved (real) transform the inverse plan
overwrites `Aωo`, which `fill!` refills on the next call.
"""
function to_time!(Ato::AbstractArray{<:Real}, Aω, Aωo, IFT)
    N = size(Aω, 1)
    No = size(Aωo, 1)
    # Scale factor makes up for difference in FFT array length, and normalises the
    # inverse transform
    scale = (No-1)/(N-1) * Utils.iscale(IFT)
    fill!(Aωo, 0)
    copy_scale!(Aωo, Aω, N, scale)
    mul!(Ato, Utils.iplan(IFT), Aωo)
end

function to_time!(Ato::AbstractArray{<:Complex}, Aω, Aωo, IFT)
    N = size(Aω, 1)
    No = size(Aωo, 1)
    scale = No/N * Utils.iscale(IFT)
    fill!(Aωo, 0)
    copy_scale_both!(Aωo, Aω, N÷2, scale)
    mul!(Ato, Utils.iplan(IFT), Aωo)
end

"""
    to_freq!(Aω, Aωo, Ato, FTplan)

Transform oversampled A(t) to A(ω) on normal grid. Dispatches on the element type of
`Ato`, not on its array type.
"""
function to_freq!(Aω, Aωo, Ato::AbstractArray{<:Real}, FTplan)
    N = size(Aω, 1)
    No = size(Aωo, 1)
    scale = (N-1)/(No-1) # Scale factor makes up for difference in FFT array length
    mul!(Aωo, FTplan, Ato)
    copy_scale!(Aω, Aωo, N, scale)
end

function to_freq!(Aω, Aωo, Ato::AbstractArray{<:Complex}, FTplan)
    N = size(Aω, 1)
    No = size(Aωo, 1)
    scale = N/No # Scale factor makes up for difference in FFT array length
    mul!(Aωo, FTplan, Ato)
    copy_scale_both!(Aω, Aωo, N÷2, scale)
end

#= Views of the first and last N samples along the first axis, with every other axis
   taken whole. A broadcast over these is one kernel on any array type. =#
_front(A, N) = view(A, 1:N, ntuple(_ -> Colon(), ndims(A)-1)...)
_back(A, N) = view(A, size(A, 1)-N+1:size(A, 1), ntuple(_ -> Colon(), ndims(A)-1)...)

function _checkshape(dest, source)
    size(dest)[2:end] == size(source)[2:end] || error(
        "dest and source must be same size except along first dimension")
end

"""
    copy_scale!(dest, source, N, scale)

Copy the first `N` elements from `source` to `dest` and simultaneously multiply by
`scale`. For multi-dimensional `dest` and `source`, work along the first axis.

`scale` is converted to the real element type of `dest` so that no `Float64` enters a
device kernel; on the default `Float64` path the conversion is the identity.
"""
function copy_scale!(dest, source, N, scale)
    _checkshape(dest, source)
    s = convert(real(eltype(dest)), scale)
    _front(dest, N) .= s .* _front(source, N)
    dest
end

"""
    copy_scale_both!(dest, source, N, scale)

Copy the first and last `N` elements from `source` to the first and last `N` elements of
`dest` and simultaneously multiply by `scale`. For multi-dimensional `dest` and `source`,
work along the first axis.
"""
function copy_scale_both!(dest, source, N, scale)
    _checkshape(dest, source)
    s = convert(real(eltype(dest)), scale)
    _front(dest, N) .= s .* _front(source, N)
    _back(dest, N) .= s .* _back(source, N)
    dest
end

# Note on noise and ionization/plasma: when the modified shot-noise model is active,
# Et_to_Pt! receives the combined field (field + noise). The noise amplitude is of order
# √(ħω·Δν) ≈ 5×10⁻⁴ √W per mode — roughly 10⁻¹⁴ of typical pulse peak power. This is
# completely negligible for the highly nonlinear ionization rate and plasma
# response, so including noise in the field passed to all response functions is physically
# reasonable. The noise meaningfully affects only Kerr and Raman processes, as intended.
"""
    Et_to_Pt!(Pt, Et, responses, density; scaling=UNIT_SCALING)
    Et_to_Pt!(Pt, Et, responses, density, idcs; scaling=UNIT_SCALING)

Accumulate the nonlinear polarisation induced by the time-domain field block `Et` into
`Pt`, one term per response, in the order the responses are given.

`Et` is `(nt,)`, `(nt, npol)` or `(nt, npol, ncols...)`; `idcs`, where a transform has
more than one column, indexes the axes past the polarisation one and **must cover every
column of the block**. The pointwise and batched paths act on the whole block and ignore
it — it is there for the columnwise path, which is a loop — so a subset would leave the
other columns holding the previous step's polarisation rather than zero. Every transform
which passes one builds it as `CartesianIndices(size(Pt)[3:end])`. `density` is a number,
or a vector with one entry per gas of a mixture, in which case `responses` is a tuple of
tuples, one per gas. `scaling` is the [`Luna.UnitScaling`](@ref) the state and the
polarisation buffer are expressed in.

Each response is applied according to its
[`Nonlinear.kind`](@ref Luna.Nonlinear.kind):

- consecutive **pointwise** responses (scalar or vector, and across the gases of a
  mixture) are evaluated as **one** fused broadcast, which sums their per-sample
  contributions in the broadcast body. No intermediate buffer, one pass over the block
  for the whole group, and it runs wherever the block does.
- a **batched** response is called once with the whole block.
- a **columnwise** response is called once per column, on the host. On a device this is
  refused; [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale) wraps such a response in a
  [`Nonlinear.HostResponse`](@ref Luna.Nonlinear.HostResponse), which is batched, before
  a device transform ever sees it.

The first group written *assigns* into `Pt` instead of zero-filling it and accumulating,
which is one pass over the block fewer. Everything after it accumulates, in tuple order,
so the sequence of additions each element sees is the one the per-response loop produced
— up to the sign of an exact zero, which `0 + (-0.0)` turned into `+0.0` and the
assignment does not.

A `responses` collection which is not a tuple (or a tuple of tuples for a mixture) falls
back to the historical per-response loop, unfused.
"""
function Et_to_Pt!(Pt, Et, responses, density, idcs...; scaling=UNIT_SCALING)
    _et_to_pt!(Pt, Et, _resppairs(responses, density), responses, density, scaling,
               idcs...)
end

#= Pairing each response with the density it sees turns the flat and the gas-mixture
   cases into the same list, so a mixture's Kerr responses fuse into one broadcast with
   their per-gas coefficients summed in its body, exactly as several responses of one gas
   do. `nothing` means "not a tuple": keep the historical loop. =#
_resppairs(responses::Tuple, density::Number) = map(r -> (r, density), responses)

function _resppairs(responses::Tuple{Vararg{Tuple}}, density::AbstractVector)
    length(responses) == length(density) || throw(DimensionMismatch(
        "$(length(responses)) response tuples for $(length(density)) densities"))
    _flatpairs(responses, density, 1)
end

_resppairs(responses, density) = nothing

_flatpairs(::Tuple{}, density, i) = ()
_flatpairs(rs::Tuple, density, i) =
    (map(r -> (r, density[i]), first(rs))..., _flatpairs(Base.tail(rs), density, i+1)...)

function _et_to_pt!(Pt, Et, ::Nothing, responses, density, scaling, idcs...)
    #= The legacy loop calls each response on the columnwise contract, which is physical
       SI units: it cannot carry a unit scaling. Nothing in Luna reaches it in a scaled
       run (`TransModeAvg`, the only scaled transform, always holds a tuple), but a
       low-level caller could. =#
    isunity(scaling) || error(
        "a response collection which is not a tuple is applied one response at a time on "*
        "the columnwise contract, which is in physical units, so it cannot be used in a "*
        "run with $(scaling). Pass the responses as a tuple.")
    _refuse_batched_legacy(responses)
    _et_to_pt_legacy!(Pt, Et, responses, density, idcs...)
end

function _et_to_pt!(Pt, Et, pairs::Tuple, responses, density, scaling, idcs...)
    #= The number of polarisation components is resolved to a compile-time constant
       before the responses are grouped, so that `kind(r, npol)` -- and with it which
       responses fuse -- is known to the compiler. =#
    npol = _npol(Et)
    if npol == 1
        _respgroups!(Pt, Et, pairs, scaling, Val(1), true, idcs...)
    elseif npol == 2
        _respgroups!(Pt, Et, pairs, scaling, Val(2), true, idcs...)
    else
        error("a nonlinear response block has 1 or 2 polarisation components, got $npol")
    end
    Pt
end

_npol(Et::AbstractArray{<:Any, 1}) = 1
_npol(Et::AbstractArray) = size(Et, 2)

#= A `Batched` response is handed the whole block and owns buffers sized for it, which
   the legacy loop cannot give it: for the transforms which pass `idcs` the loop calls
   each response one column at a time. Caught here, with the fix, rather than later in
   the response as a shape mismatch. =#
function _refuse_batched_legacy(responses)
    for r in responses
        if r isa Tuple
            _refuse_batched_legacy(r)
            continue
        end
        Nonlinear.kind(r) isa Nonlinear.Batched && error(
            "the nonlinear response $(nameof(typeof(r))) is batched: it is called once "*
            "with the whole field block and its buffers are sized for it. A response "*
            "collection which is not a `Tuple` is applied one response at a time, and "*
            "one column at a time where the transform has several, so it cannot carry "*
            "a batched response. Pass the responses as a `Tuple`.")
    end
    nothing
end

# The historical per-response loop, for a `responses` collection which is not a tuple.
function _et_to_pt_legacy!(Pt, Et, responses, density::Number)
    fill!(Pt, 0)
    for resp! in responses
        resp!(Pt, Et, density)
    end
end

function _et_to_pt_legacy!(Pt, Et, responses, density::AbstractVector)
    fill!(Pt, 0)
    for ii in eachindex(density)
        for resp! in responses[ii]
            resp!(Pt, Et, density[ii])
        end
    end
end

function _et_to_pt_legacy!(Pt, Et, responses, density, idcs)
    for i in idcs
        _et_to_pt_legacy!(view(Pt, .., i), view(Et, .., i), responses, density)
    end
end

#= Walk the response list in order, taking the longest run of pointwise responses at a
   time (one broadcast) and everything else one at a time. `firstgroup` is true while `Pt`
   has not been written yet. =#
_respgroups!(Pt, Et, ::Tuple{}, scaling, v::Val, firstgroup, idcs...) =
    (firstgroup && fill!(Pt, 0); Pt)

function _respgroups!(Pt, Et, pairs::Tuple, scaling, v::Val, firstgroup, idcs...)
    fused, rest = _splitfused(pairs, v)
    _respgroups_step!(Pt, Et, fused, rest, scaling, v, firstgroup, idcs...)
end

function _respgroups_step!(Pt, Et, ::Tuple{}, rest::Tuple, scaling, v::Val, firstgroup,
                           idcs...)
    firstgroup && fill!(Pt, 0)
    r, ρ = rest[1]
    _apply_unfused!(Pt, Et, Nonlinear.kind(r, v), r, ρ, scaling, idcs...)
    _respgroups!(Pt, Et, Base.tail(rest), scaling, v, false, idcs...)
end

function _respgroups_step!(Pt, Et, fused::Tuple, rest::Tuple, scaling, v::Val, firstgroup,
                           idcs...)
    _fusedbroadcast!(Pt, Et, fused, scaling, v, firstgroup)
    _respgroups!(Pt, Et, rest, scaling, v, false, idcs...)
end

#= The longest prefix of `pairs` whose responses fuse, and the rest. The `Union` below is
   the single place which decides which kinds fuse; a fifth kind is added there. =#
_splitfused(pairs::Tuple, v::Val) = _splitfused(pairs, v, ())
_splitfused(::Tuple{}, ::Val, acc) = (acc, ())
_splitfused(pairs::Tuple, v::Val, acc) =
    _splitfused(Nonlinear.kind(pairs[1][1], v), pairs, v, acc)
_splitfused(::Union{Nonlinear.Pointwise, Nonlinear.VectorPointwise}, pairs, v, acc) =
    _splitfused(Base.tail(pairs), v, (acc..., pairs[1]))
_splitfused(::Nonlinear.ResponseKind, pairs, v, acc) = (acc, pairs)

function _fusedbroadcast!(Pt, Et, fused::Tuple, scaling, v::Val{1}, firstgroup)
    exprs = map(q -> _scalarexpr(q[1], Nonlinear.kind(q[1], v), Et, q[2], scaling), fused)
    _materialise!(Pt, exprs, firstgroup)
    Pt
end

_scalarexpr(r, ::Nonlinear.Pointwise, Et, ρ, scaling) =
    Nonlinear.pointwise_expr(r, Et, ρ, scaling)

#= A response which declares `VectorPointwise()` unconditionally -- the natural way to
   write a two-component-only response -- would otherwise be handed the scalar path and
   silently evaluated as if it were elementwise. =#
_scalarexpr(r, ::Nonlinear.VectorPointwise, Et, ρ, scaling) = error(
    "the nonlinear response $(nameof(typeof(r))) reports `Nonlinear.kind` = "*
    "VectorPointwise() for a field block with one polarisation component. A "*
    "vector-pointwise response couples the two components, so it needs a two-component "*
    "block: declare `Nonlinear.kind(::$(nameof(typeof(r))), ::Val{2}) = "*
    "Nonlinear.VectorPointwise()` and give the one-component case its own kind "*
    "(`Pointwise()` if the same formula applies per component, `Columnwise()` otherwise).")

#= Two components, two broadcasts. Luna's buffers are (nt, npol, ncols...), so the
   polarisation index is the slow axis: a single broadcast writing an `SVector{2}` would
   need a `reinterpret` of a contiguous leading axis of length 2, which this layout does
   not have (GPU_PLAN.md section 4.3 assumed it did). Each component broadcast is still
   fused across the whole group, which is where the saving is. =#
function _fusedbroadcast!(Pt, Et, fused::Tuple, scaling, ::Val{2}, firstgroup)
    Ex = selectdim(Et, 2, 1)
    Ey = selectdim(Et, 2, 2)
    _fusedcomponent!(selectdim(Pt, 2, 1), Et, Ex, Ey, fused, scaling, firstgroup, Val(1))
    _fusedcomponent!(selectdim(Pt, 2, 2), Et, Ex, Ey, fused, scaling, firstgroup, Val(2))
    Pt
end

function _fusedcomponent!(o, Et, Ex, Ey, fused::Tuple, scaling, firstgroup, p::Val)
    exprs = map(q -> _componentexpr(q[1], Nonlinear.kind(q[1], Val(2)),
                                    Et, Ex, Ey, q[2], scaling, p), fused)
    _materialise!(o, exprs, firstgroup)
    o
end

_componentexpr(r, ::Nonlinear.Pointwise, Et, Ex, Ey, ρ, scaling, ::Val{p}) where {p} =
    Nonlinear.pointwise_expr(r, selectdim(Et, 2, p), ρ, scaling)

_componentexpr(r, ::Nonlinear.VectorPointwise, Et, Ex, Ey, ρ, scaling, p::Val) =
    Base.broadcasted(_component(p), Nonlinear.vector_expr(r, Ex, Ey, ρ, scaling))

_component(::Val{1}) = first
_component(::Val{2}) = last

# Left-associated, in the order the responses were given.
_sumexprs(e::Tuple{Any}) = e[1]
_sumexprs(e::Tuple) = Base.broadcasted(+, _sumexprs(Base.front(e)), e[end])

#= `dest` is folded in as the *leading* term of a group which is not the first, so that
   the element-by-element sequence of additions is `((dest + t1) + t2) + ...`, exactly the
   one the per-response loop produced. Summing the group first and adding `dest` to the
   total would associate differently and move the result at rounding level. =#
_materialise!(dest, exprs::Tuple, firstgroup) =
    Base.materialize!(dest, _sumexprs(firstgroup ? exprs : (dest, exprs...)))

_apply_unfused!(Pt, Et, ::Nonlinear.Batched, r, ρ, scaling) =
    Nonlinear.batched!(r, Pt, Et, ρ, scaling)

_apply_unfused!(Pt, Et, ::Nonlinear.Batched, r, ρ, scaling, idcs) =
    Nonlinear.batched!(r, Pt, Et, ρ, scaling)

function _apply_unfused!(Pt, Et, ::Nonlinear.Columnwise, r, ρ, scaling)
    _refuse_columnwise(Pt, r, scaling)
    r(Pt, Et, ρ)
end

function _apply_unfused!(Pt, Et, ::Nonlinear.Columnwise, r, ρ, scaling, idcs)
    _refuse_columnwise(Pt, r, scaling)
    for i in idcs
        r(view(Pt, .., i), view(Et, .., i), ρ)
    end
    Pt
end

#= The columnwise contract is host arrays in physical SI units, so both departures from
   it are refused in the same place: a device block, and a scaled state (where the
   response's own coefficients would be applied to `E/Eref` and its result read as
   `P/(Pref*Eref)`). Neither is reachable through a transform, which rescales every
   response at construction, but a low-level caller can produce both. =#
function _refuse_columnwise(Pt, r, scaling)
    Utils.isdevice(Pt) && error(
        "the nonlinear response $(nameof(typeof(r))) is evaluated column by column on "*
        "the host, so it cannot be applied to a $(typeof(Pt)). `Nonlinear.rescale` wraps "*
        "a columnwise response in a `Nonlinear.HostResponse` for a device run; this one "*
        "reached the transform unwrapped.")
    isunity(scaling) || error(
        "the nonlinear response $(nameof(typeof(r))) is evaluated column by column in "*
        "physical SI units, so it cannot be applied to a state in $(scaling). "*
        "`Nonlinear.rescale` wraps a columnwise response in a `Nonlinear.HostResponse` "*
        "for a scaled run; this one reached the transform unwrapped.")
    nothing
end

"""
    TransModal

Transform E(ω) -> Pₙₗ(ω) for multimode propagation via spatial integration.

# Fields
- `Emω_noise`: modal noise field `(nω, nmodes)` for the modified shot-noise model, or
  `nothing`. When present, the noise is projected to real space at each integration point
  and combined with the field in a separate buffer (`Er_nl`) for nonlinear evaluation.
  The propagating field (`Er`) is never modified.
- `Er_noise`: preallocated buffer for the real-space time-domain noise, same shape as `Er`.
- `Er_nl`: preallocated buffer for the combined field + noise, passed to `Et_to_Pt!`.
"""
mutable struct TransModal{tsT, lT, TT, FTT, IFTT, rT, gT, dT, ddT, nT, eT, enT, enlT}
    ts::tsT
    full::Bool
    dimlimits::lT
    Emω::Array{ComplexF64,2}
    Erω::Array{ComplexF64,2}
    Erωo::Array{ComplexF64,2}
    Er::Array{TT,2}
    Pr::Array{TT,2}
    Prω::Array{ComplexF64,2}
    Prωo::Array{ComplexF64,2}
    Prmω::Array{ComplexF64,2}
    FT::FTT
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    resp::rT
    grid::gT
    densityfun::dT
    density::ddT
    norm!::nT
    ncalls::Int
    z::Float64
    rtol::Float64
    atol::Float64
    mfcn::Int
    err::Array{ComplexF64,2}
    Emω_noise::eT # modal noise field for modified shot-noise model, or nothing
    Er_noise::enT # buffer for real-space time-domain noise, or nothing
    Er_nl::enlT # buffer for field+noise passed to Et_to_Pt!, or nothing
end

function show(io::IO, t::TransModal)
    grid = "grid type: $(typeof(t.grid))"
    modes = "modes: $(t.ts.nmodes)\n"*" "^4*join([string(mi) for mi in t.ts.ms], "\n    ")
    p = t.ts.indices == 1:2 ? "x,y" : t.ts.indices == 1 ? "x" : "y"
    pol = "polarisation: $p"
    samples = "time grid size: $(length(t.grid.t)) / $(length(t.grid.to))"
    resp = "responses: "*join([string(typeof(ri)) for ri in t.resp], "\n    ")
    full = "full: $(t.full)"
    out = join(["TransModal", modes, pol, grid, samples, full, resp], "\n  ")
    print(io, out)
end

"""
    TransModal(grid, ts, FT, resp, densityfun, norm!; rtol=1e-3, atol=0.0, mfcn=300, full=false, noise_field=nothing)

Construct a `TransModal`, transform E(ω) -> Pₙₗ(ω) for modal fields.

# Arguments
- `grid::AbstractGrid` : the grid used in the simulation
- `ts::Modes.ToSpace` : pre-created `ToSpace` for conversion from modal fields to space
- `FT::FFTW.Plan` : the time-frequency Fourier transform for the oversampled time grid
- `resp` : `Tuple` of response functions
- `densityfun` : callable which returns the gas density as a function of `z`
- `norm!` : normalisation function as fctn of `z`, can be created via [`norm_modal`](@ref)
- `rtol::Float=1e-3` : relative tolerance on the `HCubature` integration
- `atol::Float=0.0` : absolute tolerance on the `HCubature` integration
- `mfcn::Int=512` : maximum number of function evaluations for one modal integration
- `full::Bool=false` : if `true`, use full 2-D mode integral, if `false`, only do radial integral
- `noise_field=nothing` : optional `(nω, nmodes)` noise field for the modified shot-noise
  model. Each mode column should contain independent noise with the one-photon-per-mode
  spectral density. Generate with [`Fields.generate_noise_field`](@ref Luna.Fields.generate_noise_field).
"""
function TransModal(tT, grid, ts::Modes.ToSpace, FT, resp, densityfun, norm!;
                    rtol=1e-3, atol=0.0, mfcn=512, full=false, noise_field=nothing,
                    spec=HostSpec(), scaling=UNIT_SCALING)
    Emω = Array{ComplexF64,2}(undef, length(grid.ω), ts.nmodes)
    Erω = Array{ComplexF64,2}(undef, length(grid.ω), ts.npol)
    Erωo = Array{ComplexF64,2}(undef, length(grid.ωo), ts.npol)
    Er = Array{tT,2}(undef, length(grid.to), ts.npol)
    Pr = Array{tT,2}(undef, length(grid.to), ts.npol)
    Prω = Array{ComplexF64,2}(undef, length(grid.ω), ts.npol)
    Prωo = Array{ComplexF64,2}(undef, length(grid.ωo), ts.npol)
    Prmω = Array{ComplexF64,2}(undef, length(grid.ω), ts.nmodes)
    IFT = Utils.plan_ift(FT)
    # For the modified shot-noise model, store the modal noise field and allocate a buffer
    # for the real-space time-domain noise. The noise is projected to space at each
    # integration point in Erω_to_Prω!, so we store it in the modal domain.
    if !isnothing(noise_field)
        Emω_noise = copy(noise_field)
        Er_noise = Array{tT,2}(undef, length(grid.to), ts.npol)
        Er_nl = Array{tT,2}(undef, length(grid.to), ts.npol)
    else
        Emω_noise = nothing
        Er_noise = nothing
        Er_nl = nothing
    end
    #= Responses are given the prototype of the block they will be called with, so that
       a batched one has its buffers in the right shape, and as a `Tuple`, because a
       batched response cannot be applied by the legacy per-response loop (see
       `_refuse_batched_legacy`). Every fallback returns the response unchanged on the
       host Float64 path, which is the only one this transform runs today; `spec` and
       `scaling` are threaded through so that Group E, which gives it a device path,
       changes the caller and not this line. =#
    resp = Nonlinear.rescale_responses(Tuple(resp), spec, scaling, Er)
    TransModal(ts, full, Modes.dimlimits(ts.ms[1]), Emω, Erω, Erωo, Er, Pr, Prω, Prωo, Prmω,
               FT, IFT, resp, grid, densityfun, densityfun(0.0), norm!, 0, 0.0, rtol, atol, mfcn,
               similar(Prmω), Emω_noise, Er_noise, Er_nl)
end

function TransModal(grid::Grid.RealGrid, args...; kwargs...)
    TransModal(Float64, grid, args...; kwargs...)
end

function TransModal(grid::Grid.EnvGrid, args...; kwargs...)
    TransModal(ComplexF64, grid, args...; kwargs...)
end

function reset!(t::TransModal, Emω::Array{ComplexF64,2}, z::Float64)
    t.Emω .= Emω
    t.ncalls = 0
    t.z = z
    t.dimlimits = Modes.dimlimits(t.ts.ms[1], z=z)
    t.density = t.densityfun(z)
end

function pointcalc!(fval, xs, t::TransModal)
    # TODO: parallelize this in Julia 1.3
    for i in 1:size(xs, 2)
        x1 = xs[1, i]
        # on or outside boundaries are zero
        if x1 <= t.dimlimits[2][1] || x1 >= t.dimlimits[3][1]
            fval[:, i] .= 0.0
            continue
        end
        if size(xs, 1) > 1 # full 2-D mode integral
            x2 = xs[2, i]
            if t.dimlimits[1] == :polar
                pre = x1
            else
                if x2 <= t.dimlimits[2][2] || x1 >= t.dimlimits[3][2]
                    fval[:, i] .= 0.0
                    continue
                end
                pre = 1.0
            end
        else
            if t.dimlimits[1] == :polar
                x2 = 0.0
                pre = 2π*x1
            else
                x2 = 0.0
                pre = 1.0
            end
        end
        x = (x1,x2)
        Erω_to_Prω!(t, x)
        t.ncalls += 1
        # now project back to each mode
        # matrix product (nω x npol) * (npol x nmodes) -> (nω x nmodes)
        mul!(t.Prmω, t.Prω, transpose(t.ts.Ems))
        fval[:, i] .= pre.*reshape(reinterpret(Float64, t.Prmω), length(t.Emω)*2)
    end
end

function Erω_to_Prω!(t, x)
    Modes.to_space!(t.Erω, t.Emω, x, t.ts, z=t.z)
    to_time!(t.Er, t.Erω, t.Erωo, t.IFT)
    # Modified shot-noise model: project noise modes to real space at this spatial point,
    # convert to oversampled time domain, and combine with field in a separate buffer (Er_nl)
    # so the propagating field (Er) is never contaminated.
    if !isnothing(t.Emω_noise)
        Modes.to_space!(t.Erω, t.Emω_noise, x, t.ts, z=t.z)
        to_time!(t.Er_noise, t.Erω, t.Erωo, t.IFT)
        @. t.Er_nl = t.Er + t.Er_noise
        Et_to_Pt!(t.Pr, t.Er_nl, t.resp, t.density)
    else
        Et_to_Pt!(t.Pr, t.Er, t.resp, t.density)
    end
    @. t.Pr *= t.grid.towin
    to_freq!(t.Prω, t.Prωo, t.Pr, t.FT)
    @. t.Prω *= t.grid.ωwin
    t.norm!(t.Prω)
end

function (t::TransModal)(nl, Eω, z)
    reset!(t, Eω, z)
    _, ll, ul = t.dimlimits
    if t.full
        val, err = Cubature.hcubature_v(
            length(Eω)*2,
            (x, fval) -> pointcalc!(fval, x, t),
            ll, ul,
            reltol=t.rtol, abstol=t.atol, maxevals=t.mfcn, error_norm=Cubature.L2)
    else
        val, err = Cubature.pcubature_v(
            length(Eω)*2,
            (x, fval) -> pointcalc!(fval, x, t),
            (ll[1],), (ul[1],),
            reltol=t.rtol, abstol=t.atol, maxevals=t.mfcn, error_norm=Cubature.L2)
    end
    t.err .= reshape(reinterpret(ComplexF64, err), size(nl))
    nl .= reshape(reinterpret(ComplexF64, val), size(nl))
end

"""
    norm_modal(grid; shock=true)

Normalisation function for modal propagation. If `shock` is `false`, the intrinsic frequency
dependence of the nonlinear response is ignored, which turns off optical shock formation/
self-steepening.
"""
function norm_modal(grid; shock=true)
    ω0 = PhysData.wlfreq(grid.referenceλ)
    withshock!(nl) = @. nl *= (-im * grid.ω/4)
    withoutshock!(nl) = @. nl *= (-im * ω0/4)
    shock ? withshock! : withoutshock!
end

"""
    TransModeAvg

Transform E(ω) -> Pₙₗ(ω) for mode-averaged single-mode propagation.

# Fields
- `Et_noise`: precomputed time-domain noise on the oversampled grid for the modified
  shot-noise model (Chen & Wise, arXiv:2410.20567), or `nothing` for the traditional model.
- `Et_nl`: preallocated buffer for the combined field + noise. When `Et_noise` is present,
  `Et_nl = Eto + Et_noise` is computed at each step and passed to `Et_to_Pt!`. The
  propagating field (`Eto`) is never modified; dispersion acts only on the physical field.
"""
struct TransModeAvg{TT, ωT, FTT, IFTT, rT, gT, gvT, dT, nT, aT, eT, nlT}
    Pto::TT
    Eto::TT
    Eωo::ωT
    Pωo::ωT
    FT::FTT
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    resp::rT
    grid::gT # host grid, for metadata and for anything not in a kernel
    gv::gvT # mirror of the grid vectors the kernels broadcast against
    densityfun::dT
    norm!::nT
    aeff::aT # function which returns effective area
    Et_noise::eT # time-domain noise for modified shot-noise model, or nothing
    Et_nl::nlT # buffer for field+noise passed to Et_to_Pt!, or nothing
    scaling::UnitScaling # units the state and the polarisation are expressed in
end

function show(io::IO, t::TransModeAvg)
    grid = "grid type: $(typeof(t.grid))"
    samples = "time grid size: $(length(t.grid.t)) / $(length(t.grid.to))"
    resp = "responses: "*join([string(typeof(ri)) for ri in t.resp], "\n    ")
    out = join(["TransModeAvg", grid, samples, resp], "\n  ")
    print(io, out)
end

"""
    TransModeAvg(TT, grid, FT, IFT, resp, densityfun, norm!, aeff; kwargs...)

Construct a `TransModeAvg` transform for mode-averaged propagation. `TT` is the
time-domain element type (`Float64`/`Float32` on a `RealGrid`, complex on an `EnvGrid`);
the two-argument forms taking the grid pick it. `FT` and `IFT` are the forward and
inverse plans for the oversampled time grid
(see [`Utils.plan_ft`](@ref Luna.Utils.plan_ft)).

# Keyword arguments
- `noise_field=nothing`: optional frequency-domain noise field (on the normal grid) for the
  modified shot-noise model. When provided, it is converted to the oversampled time grid and
  stored as `Et_noise` for injection into the nonlinear operator at every propagation step.
  Generate with [`Fields.generate_noise_field`](@ref Luna.Fields.generate_noise_field).
- `spec=HostSpec()`: the [`Luna.DeviceSpec`](@ref) the buffers, mirrors and responses live
  on.
- `scaling=UNIT_SCALING`: the [`Luna.UnitScaling`](@ref) the state and the nonlinear
  polarisation are expressed in. The responses are converted to it with
  [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale) and the noise field is divided by
  `Eref`.
"""
function TransModeAvg(TT, grid, FT, IFT, resp, densityfun, norm!, aeff;
                      noise_field=nothing, spec=HostSpec(), scaling=UNIT_SCALING)
    nωo = length(grid.ωo)
    nto = length(grid.to)
    CT = Complex{realtype(spec)}
    Eωo = alloc(spec, CT, (nωo,))
    Eto = alloc(spec, TT, (nto,))
    Pto = similar(Eto)
    Pωo = similar(Eωo)
    gv = gridvectors(grid, spec)
    # Precompute time-domain noise on the oversampled grid if noise_field is provided.
    # Uses the same ω→t conversion path as to_time!: copy_scale! into oversampled spectral
    # array, then inverse FFT. The result is constant throughout propagation.
    if !isnothing(noise_field)
        Eωo_noise = alloc(spec, CT, (nωo,))
        Et_noise = alloc(spec, TT, (nto,))
        #= The noise is a state-unit quantity, so it carries the same 1/Eref the state
           does. Scaled here, once, rather than in the per-step kernel. =#
        to_time!(Et_noise, todevice(spec, noise_field ./ scaling.Eref), Eωo_noise, IFT)
        Et_nl = alloc(spec, TT, (nto,))
    else
        Et_noise = nothing
        Et_nl = nothing
    end
    #= The four-argument form: `Eto` is the prototype of the block the responses are
       called with, so one which owns buffers can allocate them here, in the run's array
       type, in time for the residency assertion below. =#
    resp = Nonlinear.rescale_responses(Tuple(resp), spec, scaling, Eto)
    #= Every mirror the transform holds and every array its responses carry, not only the
       ones this transform's own kernels touch: the assertion is what catches a future
       mistake, so it has to cover everything. =#
    resparrays = Nonlinear.resident_arrays_all(resp)
    assert_resident(spec, Eωo, Eto, Pto, Pωo, gv.ω, gv.ωwin, gv.twin, gv.towin, gv.sidx,
                    Et_noise, Et_nl, resparrays...)
    TransModeAvg(Pto, Eto, Eωo, Pωo, FT, IFT, resp, grid, gv, densityfun, norm!, aeff,
                 Et_noise, Et_nl, scaling)
end

function TransModeAvg(grid::Grid.RealGrid, FT, IFT, resp, densityfun, norm!, aeff;
                      spec=HostSpec(), kwargs...)
    TransModeAvg(realtype(spec), grid, FT, IFT, resp, densityfun, norm!, aeff;
                 spec, kwargs...)
end

function TransModeAvg(grid::Grid.EnvGrid, FT, IFT, resp, densityfun, norm!, aeff;
                      spec=HostSpec(), kwargs...)
    TransModeAvg(Complex{realtype(spec)}, grid, FT, IFT, resp, densityfun, norm!, aeff;
                 spec, kwargs...)
end

const nlscale = sqrt(PhysData.ε_0*PhysData.c/2)

function (t::TransModeAvg)(nl, Eω, z)
    to_time!(t.Eto, Eω, t.Eωo, t.IFT)
    sc = scalar(t.Eto, nlscale*sqrt(t.aeff(z)))
    @. t.Eto /= sc
    # Modified shot-noise model: compute field+noise in a separate buffer (Et_nl) so that
    # the propagating field (Eto) is never contaminated. The noise is scaled by the same
    # normalisation factor (nlscale × √Aeff) so it enters in physical units.
    if !isnothing(t.Et_noise)
        @. t.Et_nl = t.Eto + t.Et_noise / sc
        Et_to_Pt!(t.Pto, t.Et_nl, t.resp, t.densityfun(z); scaling=t.scaling)
    else
        Et_to_Pt!(t.Pto, t.Eto, t.resp, t.densityfun(z); scaling=t.scaling)
    end
    @. t.Pto *= t.gv.towin
    to_freq!(nl, t.Pωo, t.Pto, t.FT)
    t.norm!(nl, z)
    @. nl *= t.gv.ωwin # zero outside the simulation band, where ωwin is exactly 0
end

"""
    NormModeAvg

Normalisation of the mode-averaged nonlinear polarisation, built by
[`norm_mode_average`](@ref). A callable `norm!(nl, z)` which multiplies `nl` in place by
`pre/β(z)·√aeff(z)` inside the simulation band and zeroes it outside.

A struct rather than a closure so that its arrays can live on a device and be checked
for residency, and so that the per-`z` propagation constant can be mirrored.

# Fields
- `pre`: the z-independent part of the factor, with the unit scaling already folded in,
  and divided by `β` as well when that is z-independent
- `mask`: `grid.sidx` as a `Bool` mask
- `β`: mirror of the propagation constant, or `nothing` when it is z-independent
- `βfun!`, `aeff`: the host callables for `β(z)` and `Aeff(z)`
- `scaling`: the [`Luna.UnitScaling`](@ref) `pre` was built for, so that a normalisation
  handed to a run with a different one is refused rather than silently wrong by a factor
  of `Pref`
"""
struct NormModeAvg{vT, mT, bT, fT, aT}
    pre::vT
    mask::mT
    β::bT
    βfun!::fT
    aeff::aT
    scaling::UnitScaling
end

"""
    norm_mode_average(grid, βfun!, aeff; shock=true, spec=HostSpec(),
                      scaling=UNIT_SCALING, constβ=false)

Normalisation function for mode-averaged propagation; see [`NormModeAvg`](@ref).

If `shock` is `false`, the intrinsic frequency dependence of the nonlinear response is
ignored, which turns off optical shock formation/self-steepening.

`constβ=true` declares that `βfun!` does not depend on `z` -- which is the case whenever
the linear operator is constant, i.e. for a waveguide of fixed radius at fixed pressure.
`β` is then evaluated once at construction and divided into `pre`, so that no host code
runs inside the right-hand side. With `constβ=false` (the default, and what a taper or a
pressure gradient needs) `βfun!` is called on every evaluation and its result uploaded,
unless [`tabulate`](@ref) has replaced the mirror with a table over `z`
(`tabulate_linop=true` on [`Luna.run`](@ref)).
"""
function norm_mode_average(grid, βfun!, aeff; shock=true, spec=HostSpec(),
                           scaling=UNIT_SCALING, constβ=false)
    shockterm = shock ? grid.ω.^2 : grid.ω .* PhysData.wlfreq(grid.referenceλ)
    #= `Pref` converts the polarisation buffer's units back to physical ones. It is 1 on
       every Float64 run, so `pre` is then exactly what it always was. =#
    pre = @. -im*shockterm/4 / nlscale / PhysData.c * scaling.Pref
    if constβ
        βh = zeros(Float64, length(grid.ω))
        βfun!(βh, 0.0)
        check_constβ(βfun!, βh)
        pre = pre ./ βh
        β = nothing
    else
        β = HostMirror(spec, length(grid.ω))
    end
    pre = todevice(spec, pre)
    mask = todevice(spec, grid.sidx)
    assert_resident(spec, pre, mask, isnothing(β) ? nothing : β.dev)
    NormModeAvg(pre, mask, β, βfun!, aeff, scaling)
end

"""
The distance, in metres, at which [`check_constβ`](@ref) evaluates `βfun!` a second time.
Small enough to be inside any waveguide Luna is used for, and far enough from zero that a
taper or a pressure gradient has moved.
"""
const CONSTβ_PROBE_Z = 1e-3

"""
    check_constβ(βfun!, β0)

Check the claim `constβ=true` makes, rather than trusting it: `β` is baked into the
normalisation at `z = 0`, so a `βfun!` which does depend on `z` would give a silently
wrong propagation with `β` frozen there. Evaluating it once more at
[`CONSTβ_PROBE_Z`](@ref) at setup turns that into an error message.
"""
function check_constβ(βfun!, β0)
    βz = similar(β0)
    try
        βfun!(βz, CONSTβ_PROBE_Z)
    catch e
        error("constβ=true was passed, but βfun! could not be evaluated at "*
              "z = $(CONSTβ_PROBE_Z) m, which is how that claim is checked: $e. Pass "*
              "constβ=false for a z-dependent waveguide.")
    end
    βz == β0 || error(
        "constβ=true was passed, but βfun! gives a different propagation constant at "*
        "z = $(CONSTβ_PROBE_Z) m than at z = 0 (largest difference "*
        "$(maximum(abs, βz .- β0))). It would be folded into the normalisation at z = 0 "*
        "and the propagation would be silently wrong. Pass constβ=false, which is the "*
        "default and what a taper or a pressure gradient needs.")
    nothing
end

#= β is 1 rather than 0 outside the simulation band (`LinearOps` fills it that way), so
   the division is finite everywhere and the mask decides what survives. Zeroing out of
   band rather than skipping matters: skipping would leave the raw, unnormalised
   transform of the polarisation in place, and since the linear operator is also zero out
   of band nothing downstream would remove it. =#
function (n::NormModeAvg{vT, mT, Nothing})(nl, z) where {vT, mT}
    sqrtaeff = scalar(nl, sqrt(n.aeff(z)))
    pre = n.pre
    mask = n.mask
    z0 = zero(eltype(nl))
    @. nl = ifelse(mask, nl*(pre*sqrtaeff), z0)
end

#= `β` at `z`, on the array type the kernel broadcasts against. Two sources: a
   `Luna.HostMirror`, filled by host scalar code and uploaded on every evaluation, or a
   `LinearOps.TabulatedVector`, which reads it out of a z table and never touches the host
   (`tabulate_linop=true`; the table ignores `βfun!`, which it was built from). =#
_βdev(m::HostMirror, βfun!, z) = (βfun!(m.host, z); upload!(m))
_βdev(tab, βfun!, z) = tab(z)

function (n::NormModeAvg)(nl, z)
    β = _βdev(n.β, n.βfun!, z)
    sqrtaeff = scalar(nl, sqrt(n.aeff(z)))
    pre = n.pre
    mask = n.mask
    z0 = zero(eltype(nl))
    @. nl = ifelse(mask, nl*(pre/β*sqrtaeff), z0)
end

"""
    NormModeAvgGNLSE

Normalisation of the mode-averaged nonlinear polarisation in the GNLSE form, built by
[`norm_mode_average_gnlse`](@ref). See [`NormModeAvg`](@ref).
"""
struct NormModeAvgGNLSE{vT, mT, aT}
    pre::vT
    mask::mT
    aeff::aT
    scaling::UnitScaling
end

"""
    norm_mode_average_gnlse(grid, aeff; shock=true, spec=HostSpec(),
                            scaling=UNIT_SCALING)

Normalisation function for the GNLSE form of mode-averaged propagation; see
[`norm_mode_average`](@ref) and [`NormModeAvgGNLSE`](@ref).
"""
function norm_mode_average_gnlse(grid, aeff; shock=true, spec=HostSpec(),
                                 scaling=UNIT_SCALING)
    shockterm = shock ? grid.ω.^2 : grid.ω .* PhysData.wlfreq(grid.referenceλ)
    pre = @. -im*shockterm/(2*PhysData.c^(3/2)*sqrt(2*PhysData.ε_0))/(grid.ω/PhysData.c)
    pre = pre .* scaling.Pref
    pre = todevice(spec, pre)
    mask = todevice(spec, grid.sidx)
    assert_resident(spec, pre, mask)
    NormModeAvgGNLSE(pre, mask, aeff, scaling)
end

function (n::NormModeAvgGNLSE)(nl, z)
    sqrtaeff = scalar(nl, sqrt(n.aeff(z)))
    pre = n.pre
    mask = n.mask
    z0 = zero(eltype(nl))
    @. nl = ifelse(mask, nl*(pre*sqrtaeff), z0)
end

"""
    check_norm(norm!, spec, scaling)

Check that a caller-supplied normalisation can be used with `spec` and `scaling`.

Luna's own normalisations take both at construction, so they are checked against both:
the arrays for residency and the scaling for equality. A normalisation built with the
right array type but the wrong scaling would otherwise produce a polarisation wrong by
the factor `Pref` with nothing to say so. Anything else is accepted only for an unscaled
`Float64` host run, which is the default CPU path.
"""
_normerror(n, spec, scaling) = error(
    "the normalisation $(typeof(n)) was not built for $(spec) with $(scaling). Build it "*
    "with the `spec` and `scaling` keywords of `norm_mode_average` (or let `Luna.setup` "*
    "do it), or run on the default CPU path.")

check_norm(n, spec, scaling) = (arraytype(spec) === Array && realtype(spec) === Float64 &&
                                isunity(scaling)) ? nothing : _normerror(n, spec, scaling)

function check_norm(n::NormModeAvg, spec, scaling)
    (n.scaling == scaling &&
     all_resident(spec, n.pre, n.mask, isnothing(n.β) ? nothing : n.β.dev)) ||
        _normerror(n, spec, scaling)
    nothing
end

function check_norm(n::NormModeAvgGNLSE, spec, scaling)
    (n.scaling == scaling && all_resident(spec, n.pre, n.mask)) ||
        _normerror(n, spec, scaling)
    nothing
end

"""
    TransRadial

Transform E(ω) -> Pₙₗ(ω) for radially symmetric free-space propagation.

Parametric in its buffer array type: every buffer, grid mirror and transform matrix is
allocated with [`Luna.alloc`](@ref)/[`Luna.todevice`](@ref) from the run's
[`Luna.DeviceSpec`](@ref), so the same code runs on the host, in reduced precision and on
a device.

# Fields
- `rgrid`: the transverse grid ([`Grid.RadialGrid`](@ref Luna.Grid.RadialGrid)), kept for
  metadata; the per-step code uses `Tfwd`/`Tbwd` instead.
- `gv`: mirror of the grid vectors the kernels broadcast against.
- `Tfwd`, `Tbwd`: the grid's transform matrices in the *time-domain* element type, so that
  both operands of the Hankel GEMM have the same element type (Metal's accelerated matrix
  multiply requires that).
- `prefac`: the z-independent part of the frequency-domain normalisation,
  `ωwin·(-iω)·Pref`, precombined on the host.
- `Et_noise`: precomputed time-domain noise on the oversampled real-space grid `(nto, nr)`
  for the modified shot-noise model, or `nothing`.
- `Et_nl`: preallocated buffer for the combined field + noise, passed to `Et_to_Pt!`. The
  propagating field (`Eto`) is never modified.
- `scaling`: the [`Luna.UnitScaling`](@ref) the state and the polarisation are in.
"""
struct TransRadial{TT, ωT, RGT, FTT, IFTT, nT, rT, gT, gvT, dT, iT, mT, pT, eT, nlT}
    rgrid::RGT # transverse grid (Grid.RadialGrid: space to k-space)
    FT::FTT # Fourier transform (time to frequency)
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    normfun::nT # Function which returns normalisation factor
    resp::rT # nonlinear responses (tuple of callables)
    grid::gT # host grid, for metadata and for anything not in a kernel
    gv::gvT # mirror of the grid vectors the kernels broadcast against
    densityfun::dT # callable which returns density
    Pto_r::TT # Buffer array for NL polarisation on oversampled time grid
    Pto_k::TT # Buffer array for NL polarisation on oversampled time grid
    Eto_r::TT # Buffer array for field on oversampled time grid
    Eto_k::TT # Buffer array for field on oversampled time grid
    Eωo::ωT # Buffer array for field on oversampled frequency grid
    Pωo::ωT # Buffer array for NL polarisation on oversampled frequency grid
    idcs::iT # CartesianIndices for Et_to_Pt! to iterate over
    Tfwd::mT # forward Hankel transform matrix, in the time-domain element type
    Tbwd::mT # backward Hankel transform matrix, in the time-domain element type
    prefac::pT # ωwin*(-im*ω)*Pref: the z-independent normalisation factor
    Et_noise::eT # time-domain noise for modified shot-noise model, or nothing
    Et_nl::nlT # buffer for field+noise passed to Et_to_Pt!, or nothing
    scaling::UnitScaling # units the state and the polarisation are expressed in
end

function show(io::IO, t::TransRadial)
    grid = "grid type: $(typeof(t.grid))"
    samples = "time grid size: $(length(t.grid.t)) / $(length(t.grid.to))"
    resp = "responses: "*join([string(typeof(ri)) for ri in t.resp], "\n    ")
    nr = "radial points: $(t.rgrid.N)"
    R = "aperture: $(t.rgrid.R)"
    out = join(["TransRadial", grid, samples, nr, R, resp], "\n  ")
    print(io, out)
end

"""
    TransRadial(TT, grid, rgrid, FT, responses, densityfun, normfun, pol=false; kwargs...)

Construct a `TransRadial` to calculate the reciprocal-domain nonlinear polarisation.
`rgrid` is a [`Grid.RadialGrid`](@ref Luna.Grid.RadialGrid). `TT` is the time-domain
element type (`Float64`/`Float32` on a `RealGrid`, complex on an `EnvGrid`); the
two-argument forms taking the grid pick it.

# Keyword arguments
- `noise_field=nothing`: optional `(nω, npol, nk)` frequency/k-space noise field for the
  modified shot-noise model. When provided, it is converted to the real-space time domain
  `(nto, npol, nr)` via inverse FFT and inverse Hankel transform, and stored as `Et_noise`.
  Generate with [`Fields.generate_noise_field`](@ref Luna.Fields.generate_noise_field).
- `spec=HostSpec()`: the [`Luna.DeviceSpec`](@ref) the buffers, mirrors and responses live
  on.
- `scaling=UNIT_SCALING`: the [`Luna.UnitScaling`](@ref) the state and the nonlinear
  polarisation are expressed in. The responses are converted to it with
  [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale), the noise field is divided by `Eref`
  and `Pref` is folded into `prefac`.
"""
function TransRadial(TT, grid, rgrid::Grid.RadialGrid, FT, responses, densityfun, normfun,
                     pol=false; noise_field=nothing, spec=HostSpec(),
                     scaling=UNIT_SCALING)
    np = pol ? 2 : 1
    N = rgrid.N
    CT = Complex{realtype(spec)}
    IFT = Utils.plan_ift(FT)
    Eωo = alloc(spec, CT, (length(grid.ωo), np, N))
    Eto_r = alloc(spec, TT, (length(grid.to), np, N))
    Pto_r = similar(Eto_r)
    Eto_k = similar(Eto_r)
    Pto_k = similar(Eto_r)
    Pωo = similar(Eωo)
    idcs = CartesianIndices(size(Pto_r)[3:end])
    gv = gridvectors(grid, spec)
    #= Our own copies of the grid's transform matrices in the type we multiply: a GEMM
       needs both operands in the same element type (and on devices it is required).
       `convert` first, on the host, because a real-to-complex conversion is not one
       `todevice` performs. =#
    Tfwd = todevice(spec, convert(Matrix{TT}, rgrid.Tfwd))
    Tbwd = todevice(spec, convert(Matrix{TT}, rgrid.Tbwd))
    prefac = fsprefac(grid, spec, scaling)
    #= Precompute time-domain noise in real space: ω→t via to_time!, then k→r. This is
       Grid.to_rspace! done with our own Tbwd, so that the noise passes through exactly
       the same matrix as the field does on every step. The noise is a state-unit
       quantity, so it carries the same 1/Eref the state does. =#
    if !isnothing(noise_field)
        Eωo_noise = alloc(spec, CT, (length(grid.ωo), np, N))
        Et_noise = alloc(spec, TT, (length(grid.to), np, N))
        nf = isunity(scaling) ? noise_field : noise_field ./ scaling.Eref
        to_time!(Et_noise, todevice(spec, nf), Eωo_noise, IFT)
        Grid.radial_matmul!(Et_noise, Et_noise, Tbwd)
        Et_nl = alloc(spec, TT, (length(grid.to), np, N))
    else
        Et_noise = nothing
        Et_nl = nothing
    end
    #= Responses are given the prototype of the block they will be called with, so that
       a batched one has its buffers in the right shape, and as a `Tuple`, because a
       batched response cannot be applied by the legacy per-response loop (see
       `_refuse_batched_legacy`). =#
    responses = Nonlinear.rescale_responses(Tuple(responses), spec, scaling, Eto_r)
    check_norm(normfun, spec, scaling)
    #= Every mirror the transform holds and every array its responses carry, not only the
       ones this transform's own kernels touch. =#
    resparrays = Nonlinear.resident_arrays_all(responses)
    assert_resident(spec, Eωo, Pωo, Eto_r, Eto_k, Pto_r, Pto_k, Tfwd, Tbwd, prefac,
                    gv.ω, gv.ωwin, gv.twin, gv.towin, gv.sidx, Et_noise, Et_nl,
                    resparrays...)
    TransRadial(rgrid, FT, IFT, normfun, responses, grid, gv, densityfun,
                Pto_r, Pto_k, Eto_r, Eto_k, Eωo, Pωo, idcs,
                Tfwd, Tbwd, prefac, Et_noise, Et_nl, scaling)
end

# accept a Hankel.QDHT as before, converting it (with a deprecation warning)
function TransRadial(TT::Type, grid, q::Grid.HankelTransform, args...; kwargs...)
    TransRadial(TT, grid, Grid.RadialGrid(q), args...; kwargs...)
end

function TransRadial(grid::Grid.RealGrid, args...; spec=HostSpec(), kwargs...)
    TransRadial(realtype(spec), grid, args...; spec, kwargs...)
end

function TransRadial(grid::Grid.EnvGrid, args...; spec=HostSpec(), kwargs...)
    TransRadial(Complex{realtype(spec)}, grid, args...; spec, kwargs...)
end

"""
    (t::TransRadial)(nl, Eω, z)

Calculate the reciprocal-domain (ω-k-space) nonlinear response due to the field `Eω` and
place the result in `nl`

The two Hankel steps are **one** matrix multiply each, on the block reshaped to
`(nto·npol, nr)` ([`Grid.radial_matmul!`](@ref Luna.Grid.radial_matmul!)), rather than one
per polarisation component on a view: a device's accelerated matrix multiply needs plain
zero-offset operands of equal element type, which a view is not.
"""
function (t::TransRadial)(nl, Eω, z)
    to_time!(t.Eto_k, Eω, t.Eωo, t.IFT) # transform ω -> t
    Grid.radial_matmul!(t.Eto_r, t.Eto_k, t.Tbwd) # transform Eto k -> r
    # Modified shot-noise: compute field+noise in separate buffer (Et_nl) so the
    # propagating field is never contaminated.
    # Note that if noise_field is nothing, we pass t.Eto_r straight without copying
    # to the buffer t.Et_nl first
    if !isnothing(t.Et_noise)
        @. t.Et_nl = t.Eto_r + t.Et_noise
        Et_to_Pt!(t.Pto_r, t.Et_nl, t.resp, t.densityfun(z), t.idcs; scaling=t.scaling)
    else
        Et_to_Pt!(t.Pto_r, t.Eto_r, t.resp, t.densityfun(z), t.idcs; scaling=t.scaling)
    end
    @. t.Pto_r *= t.gv.towin # apodisation
    Grid.radial_matmul!(t.Pto_k, t.Pto_r, t.Tfwd) # transform Pto r -> k
    to_freq!(nl, t.Pωo, t.Pto_k, t.FT) # transform t -> ω
    fsnorm!(nl, t.prefac, t.normfun(z))
end

"""
    fsnorm!(nl, pre, nrm)

Apply the frequency-domain normalisation of a free-space transform to `nl` in place, as
one fused broadcast over the precombined prefactor `pre` ([`fsprefac`](@ref), which
carries `ωwin`, `-iω` and `Pref`) and the normalisation array `nrm`
([`FreeSpaceNorm`](@ref)).

Written with the same association as the expression it replaces, `pre ./ (2 .* norm)`, so
the `Float64` path is unchanged; the `2` goes through [`Luna.scalar`](@ref) so that
nothing reachable from a kernel is a `Float64` or an `Int`. Shared by
[`TransRadial`](@ref), [`TransFree`](@ref) and [`TransFree2D`](@ref).
"""
function fsnorm!(nl, pre, nrm)
    two = scalar(nl, 2.0)
    @. nl *= pre/(two*nrm)
end

"""
    fsprefac(grid, spec, scaling)

The z-independent part of the frequency-domain normalisation of a free-space transform,
`ωwin·(-iω)·Pref`, precombined on the host in `Float64` as one vector (GPU_PLAN.md §4.1)
and mirrored to `spec`. Shared by [`TransRadial`](@ref), [`TransFree`](@ref) and
[`TransFree2D`](@ref), which all multiply by it once per right-hand side in
[`fsnorm!`](@ref).

`Pref` converts the polarisation buffer's units back to physical ones and is `1` on every
`Float64` run, so this is exactly the vector the per-step expression it replaces used to
rebuild on every call.
"""
fsprefac(grid, spec, scaling) =
    todevice(spec, @. grid.ωwin * (-im*grid.ω) * scaling.Pref)

#=================================================#
#==========  FREE-SPACE NORMALISATION  ===========#
#=================================================#

"""
    FreeSpaceNorm

Normalisation factor for the free-space transforms ([`TransRadial`](@ref),
[`TransFree`](@ref), [`TransFree2D`](@ref)): a callable `normfun(z)` returning an array
`(Nω, Npol, Nk...)` by which the transforms *divide* the nonlinear polarisation, so that
the source term of the UPPE is ``-i\\omega/(2\\,\\mathrm{norm})\\,P_\\mathrm{nl}``.

# Physics

The factor is ``\\beta_z/(\\mu_0\\omega)`` with ``\\beta_z = \\sqrt{k^2 - k_\\perp^2}`` the
longitudinal wavevector ([`LinearOps.βz`](@ref)), giving the standard UPPE source
``-i\\mu_0\\omega^2 P_\\mathrm{nl}/(2\\beta_z)``. Below cutoff (``k_\\perp > k``) ``\\beta_z``
is imaginary, ``-i\\kappa``, and the same expression gives the evanescent source
``\\mu_0\\omega^2 P_\\mathrm{nl}/(2\\kappa)`` — real, finite and consistent with the decay
``e^{-\\kappa z}`` the linear operator applies. Historically this branch was set to `1.0`,
which is dimensionally inconsistent and arbitrary. The factor is exactly `1.0` only at
``\\omega = 0`` and outside the simulation band, where the nonlinear polarisation is zero
anyway.

# Source taper

An evanescent component driven by a nonlinear source settles to its adiabatic amplitude
``S/\\kappa`` within a distance ``1/\\kappa``, and it does not radiate. Luna's
interaction-picture stepper, however, cannot integrate such a channel when
``\\kappa\\,\\Delta z \\gg 1``: it back-propagates the source by ``e^{+\\kappa\\Delta z}``. The
remedy (the same one [`Luna.Boundaries`](@ref) uses for the spectral window) is to taper the
*source* in those channels: the normalisation is divided by
```math
W = \\exp\\big[-\\min(\\kappa, \\kappa_\\mathrm{max})\\,\\ell\\big] \\, W_k(k_\\perp)\\,,
```
where ``\\ell`` is the reference length over which the stepper is allowed to take one
step (`max_dz ≤ ℓ`), ``\\kappa_\\mathrm{max}`` is the cap [`Luna.Boundaries`](@ref) also
applies to the decay rate of the linear operator, and ``W_k`` is the k-space absorbing
window. The amplification ``e^{\\kappa\\Delta z}`` is then never larger than ``1/W``, which
is bounded. Channels with ``\\kappa\\ell \\gg 1`` have their source removed altogether,
which is the correct limit: their adiabatic amplitude vanishes as ``1/\\kappa``. Channels
with ``\\kappa\\ell \\lesssim 1`` are integrated exactly.

`ℓ`, `κmax` and `W_k` are set by [`reflength!`](@ref), which `Luna.run` calls through
`Boundaries.setup` once the reference length is known. Until then `ℓ = 0` and there is no
taper: the factor is the pure physics, which is only usable with steps `Δz ≲ 1/κ`.

# Fields
- `grid`, `spacegrid`: the temporal and transverse grids
- `nfun`: refractive index, either `nfun(ω; z)` returning one index or a tuple (one per
  polarisation), or a tuple `(nfunx, nfuny)` of crystal-optics functions
  `nfunx(λ, δθ; z)`, `nfuny(λ; z)` (see [`Luna.PhysData.crystal_internal_angle`](@ref))
- `kperp2`, `kidcs`: squared transverse wavevector and the indices of the k axes, on the
  host (the crystal-optics fill is host scalar code)
- `out`: the normalisation array (complex, since ``\\beta_z`` is complex below cutoff), on
  the run's array type and precision
- `ℓ`, `κmax`, `kwin`: taper parameters (see above), on the host
- `constant`: if `true`, `out` is computed once and reused (the index does not depend on `z`)
- `spec`: the [`Luna.DeviceSpec`](@ref) `out` and the mirrors live on
- `kperp2m`, `kwinm`, `sidxm`, `ωm`: mirrors of the vectors the isotropic fill broadcasts
  against, in the run's array type and precision. `kwinm` is refilled from `kwin` on the
  first `fillnorm!` after [`reflength!`](@ref) (which is called *after* construction, by
  `Boundaries.setup`), which is what `mirrored` tracks.
- `nm`: host mirror of the `(Nω, Npol)` refractive-index table, evaluated by host scalar
  code and uploaded once per `fillnorm!`
- `ohost`: host staging buffer for the crystal-optics fill, or `nothing` when `out` is
  already a host array of the right precision
"""
mutable struct FreeSpaceNorm{gT, sT, nT, kT, iT, oT, wT, ST, dT, mT, xT, nmT, hT}
    grid::gT
    spacegrid::sT
    nfun::nT
    kperp2::kT
    kidcs::iT
    out::oT
    ℓ::Float64
    κmax::Float64
    kwin::wT
    constant::Bool
    filled::Bool
    spec::ST
    kperp2m::dT
    kwinm::dT
    sidxm::mT
    ωm::xT
    nm::nmT
    ohost::hT
    mirrored::Bool
end

npol(nfun::Tuple, grid) = 2 # crystal optics: (nfunx, nfuny)
npol(nfun, grid) = length(nfun(grid.ω[findfirst(grid.sidx)]; z=0)) # 1 if single index, 2 if nx, ny

function FreeSpaceNorm(grid, spacegrid, nfun; constant, spec=HostSpec())
    kperp2, kidcs = transverse_k2(spacegrid)
    np = npol(nfun, grid)
    nω = length(grid.ω)
    T = realtype(spec)
    out = alloc(spec, Complex{T}, (nω, np, size(kidcs)...))
    kwin = ones(Float64, size(kidcs))
    #= The k-space quantities are broadcast against `out`, whose first two axes are ω and
       polarisation, so they are mirrored already reshaped to `(1, 1, Nk...)`. =#
    kshape = (1, 1, size(kidcs)...)
    kperp2m = todevice(spec, reshape(convert(Array{Float64}, kperp2), kshape))
    kwinm = todevice(spec, reshape(copy(kwin), kshape))
    ohost = (arraytype(spec) === Array && T === Float64) ? nothing :
            zeros(ComplexF64, size(out))
    FreeSpaceNorm(grid, spacegrid, nfun, kperp2, kidcs, out, 0.0, Inf, kwin, constant,
                  false, spec, kperp2m, kwinm, todevice(spec, grid.sidx),
                  todevice(spec, grid.ω), HostMirror(spec, nω*np), ohost, false)
end

"""
    retarget(normfun, spec)

The same normalisation with its arrays on `spec`'s array type and precision: `normfun`
itself when they already are, and a fresh [`FreeSpaceNorm`](@ref) otherwise, carrying over
whatever [`reflength!`](@ref) has already set.

`Luna.setup` calls this because the free-space normalisations are built by the caller,
before the run's device is known, and every low-level radial script (and every example)
builds one with no `spec`. Anything which is not a `FreeSpaceNorm` cannot be retargeted
and is accepted only for the default host `Float64` path.
"""
function retarget(nf::FreeSpaceNorm, spec)
    all_resident(spec, nf.out) && return nf
    new = FreeSpaceNorm(nf.grid, nf.spacegrid, nf.nfun; constant=nf.constant, spec)
    new.ℓ = nf.ℓ
    new.κmax = nf.κmax
    new.kwin .= nf.kwin
    new
end

retarget(normfun, spec) =
    (arraytype(spec) === Array && realtype(spec) === Float64) ? normfun : error(
        "the normalisation $(typeof(normfun)) is host Float64 code, so it cannot be used "*
        "for a run on $(spec). Build it with `NonlinearRHS.norm_radial`/`norm_free`/"*
        "`norm_free2D` (or the `const_` variants), which `Luna.setup` retargets itself.")

function check_norm(nf::FreeSpaceNorm, spec, scaling)
    #= The unit scaling is not the normalisation's: `βz/(μ0 ω)` is physics, and `Pref` is
       folded into the transform's own `prefac`. Only residency is checked here. =#
    all_resident(spec, nf.out, nf.kperp2m, nf.kwinm, nf.sidxm, nf.ωm, nf.nm.dev) ||
        _normerror(nf, spec, scaling)
    nothing
end

function (nf::FreeSpaceNorm)(z)
    if !(nf.constant && nf.filled)
        fillnorm!(nf, z)
        nf.filled = true
    end
    nf.out
end

"""
    reflength!(normfun, ℓ; κmax=Inf, kwin=nothing)
    reflength!(transform, ℓ; κmax=Inf, kwin=nothing)

Set the reference length `ℓ` of the source taper of a [`FreeSpaceNorm`](@ref) (or of the
normalisation held by a free-space transform), the cap `κmax` on the evanescent decay rate
it assumes, and the k-space window profile `kwin` (an array over the k axes, or `nothing`
for none). Called by `Boundaries.setup`; see [`FreeSpaceNorm`](@ref) for the meaning.

For transforms which are not free-space, or whose normalisation is not a
[`FreeSpaceNorm`](@ref), this does nothing (with a warning in the latter case, since the
evanescent channels are then left untapered).
"""
function reflength!(nf::FreeSpaceNorm, ℓ; κmax=Inf, kwin=nothing)
    nf.ℓ = ℓ
    nf.κmax = κmax
    isnothing(kwin) ? fill!(nf.kwin, 1) : (nf.kwin .= kwin)
    nf.filled = false
    #= `Boundaries.setup` calls this after the transform exists, so the k-window mirror is
       stale from here until the next `fillnorm!` rebuilds it. =#
    nf.mirrored = false
    nf
end

function reflength!(normfun::Function, ℓ; kwargs...)
    Logging.@warn("The normalisation function is not a FreeSpaceNorm, so the evanescent " *
                  "source cannot be tapered. Use norm_radial/norm_free/norm_free2D (or the " *
                  "const_ variants) to build it.")
    nothing
end

#= The factor for one (ω, polarisation, k⊥) element: physics divided by the taper. Exactly
   at cutoff (βsq == 0) the physical factor vanishes and the source would be infinite; that
   point is measure-zero and was always returned as 1.0, so keep doing that. =#
function normfactor(nf::FreeSpaceNorm, βsq, ω, wk)
    βsq == 0 && return complex(1.0)
    W = wk
    if βsq < 0
        W *= exp(-min(sqrt(-βsq), nf.κmax)*nf.ℓ)
    end
    βz(βsq)/(PhysData.μ_0*ω)/W
end

#= `LinearOps.βz` written so that the result is `Complex{T}` for a `T` argument without an
   `im` literal (which is a `Complex{Int}` and would put a 64-bit integer into a kernel).
   `-im*sqrt(x)` is `Complex(0.0*x, -1.0*x)`, i.e. exactly this, so the Float64 path is
   unchanged. =#
_βzc(βsq::T) where {T} = βsq < 0 ? complex(zero(T), -sqrt(-βsq)) : complex(sqrt(βsq), zero(T))

#= [`normfactor`](@ref) as a broadcast kernel: the same expression, with every constant
   converted to the element type so that nothing reachable from a device kernel is a
   `Float64`. `inband` is `grid.sidx`; out of band, and at ω = 0, the factor is 1, which is
   what the loop's `out[iω, :, ii] .= 1` wrote. =#
function _normelem(n::T, ω::T, kperp2::T, wk::T, inband::Bool, c::T, μ0::T,
                   κmax::T, ℓ::T) where {T}
    unity = complex(one(T))
    (inband && ω != zero(T)) || return unity
    βsq = (n*ω/c)^2 - kperp2
    βsq == zero(T) && return unity
    W = wk
    if βsq < zero(T)
        W *= exp(-min(sqrt(-βsq), κmax)*ℓ)
    end
    _βzc(βsq)/(μ0*ω)/W
end

#= The k-space window is set by `reflength!` after the norm is built, so its mirror is
   refilled on the first `fillnorm!` which follows. Everything else the isotropic kernel
   broadcasts against is fixed at construction. =#
function _mirrors!(nf::FreeSpaceNorm)
    nf.mirrored && return nothing
    T = realtype(nf.spec)
    copyto!(nf.kwinm, T === Float64 ? nf.kwin : convert(Array{T}, nf.kwin))
    nf.mirrored = true
    nothing
end

#= The refractive index is host scalar code (a Sellmeier equation, or a user's function),
   so it is evaluated on the host into the `(Nω, Npol)` mirror and uploaded once per call.
   Out of band it is left at 1: `nfun` is not evaluated there -- the loop this replaces did
   not either, and an index function need not be defined outside the simulation band. =#
function _fillindex!(nf::FreeSpaceNorm, z)
    ω = nf.grid.ω
    np = size(nf.out, 2)
    M = reshape(nf.nm.host, length(ω), np)
    fill!(M, 1.0)
    for iω in eachindex(ω)
        (ω[iω] == 0 || !nf.grid.sidx[iω]) && continue
        ns = nf.nfun(ω[iω]; z)
        for ip in 1:np
            M[iω, ip] = real(ns[ip])
        end
    end
    reshape(upload!(nf.nm), length(ω), np)
end

# isotropic: nfun(ω; z) -> n or (nx, ny), the same k⊥ for every polarisation
function fillnorm!(nf::FreeSpaceNorm, z)
    _mirrors!(nf)
    nd = _fillindex!(nf, z)
    T = realtype(nf.spec)
    c = convert(T, PhysData.c)
    μ0 = convert(T, PhysData.μ_0)
    κmax = convert(T, nf.κmax)
    ℓ = convert(T, nf.ℓ)
    #= One broadcast over the whole `(Nω, Npol, Nk...)` array: `nd` is `(Nω, Npol)`, the
       grid mirrors are `(Nω,)` and the k-space mirrors are `(1, 1, Nk...)`. =#
    @. nf.out = _normelem(nd, nf.ωm, nf.kperp2m, nf.kwinm, nf.sidxm, c, μ0, κmax, ℓ)
    nf.out
end

#= crystal optics: nfunx(λ, δθ; z) depends on the internal angle, which depends on kx only,
   so the angle is found once per (ω, kx) and reused along ky. For Free2DGrid the k axes are
   (Nkx,) and the trailing ky index below is the (allowed) singleton 1.

   The root-finding for the internal angle is host scalar code with no kernel, so this fill
   stays on the host and its result is uploaded (GPU_PLAN.md 4.4). `ohost` is `nothing`, and
   the copy is skipped, whenever `out` is already a host `ComplexF64` array. =#
function fillnorm!(nf::FreeSpaceNorm{<:Any, <:Any, <:Tuple}, z)
    nfunx, nfuny = nf.nfun
    ω = nf.grid.ω
    out = isnothing(nf.ohost) ? nf.out : nf.ohost
    kx = nf.spacegrid.kx
    for iω in eachindex(ω)
        if ω[iω] == 0 || !nf.grid.sidx[iω]
            out[iω, :, nf.kidcs] .= 1
            continue
        end
        ny = real(nfuny(wlfreq(ω[iω]); z))
        ksq_ypol = (ny*ω[iω]/PhysData.c)^2
        for ix in eachindex(kx)
            δθ = crystal_internal_angle((λ, δθ) -> nfunx(λ, δθ; z), ω[iω], kx[ix])
            nx = real(nfunx(wlfreq(ω[iω]), δθ; z))
            ksq_xpol = (nx*ω[iω]/PhysData.c)^2
            for iy in axes(nf.kperp2, 2)
                kperp2 = nf.kperp2[ix, iy]
                wk = nf.kwin[ix, iy]
                out[iω, 1, ix, iy] = normfactor(nf, ksq_xpol - kperp2, ω[iω], wk)
                out[iω, 2, ix, iy] = normfactor(nf, ksq_ypol - kperp2, ω[iω], wk)
            end
        end
    end
    out === nf.out ||
        copyto!(nf.out, convert(Array{eltype(nf.out)}, out))
    nf.out
end

"""
    norm_radial(grid, q, nfun)
    norm_free(grid, xygrid, nfun)
    norm_free2D(grid, xgrid, nfun)

Make the normalisation factor ([`FreeSpaceNorm`](@ref)) for radial, full-3D and 2D (x-z)
free-space propagation with a `z`-dependent refractive index, recomputed on every call.

`nfun(ω; z)` takes frequency `ω` and a keyword argument `z` and returns either one index
or a tuple of indices (one per polarisation). For crystal optics (`norm_free`,
`norm_free2D` only) pass a tuple `(nfunx, nfuny)` with `nfunx(λ, δθ; z)` and `nfuny(λ; z)`,
as for [`LinearOps.make_const_linop`](@ref).
"""
norm_radial(grid, rg::Grid.RadialGrid, nfun; spec=HostSpec()) =
    FreeSpaceNorm(grid, rg, nfun; constant=false, spec)
norm_radial(grid, q::Grid.HankelTransform, nfun; kwargs...) =
    norm_radial(grid, Grid.RadialGrid(q), nfun; kwargs...)
norm_free(grid, xygrid::Grid.FreeGrid, nfun; spec=HostSpec()) =
    FreeSpaceNorm(grid, xygrid, nfun; constant=false, spec)
norm_free2D(grid, xgrid::Grid.Free2DGrid, nfun; spec=HostSpec()) =
    FreeSpaceNorm(grid, xgrid, nfun; constant=false, spec)

"""
    const_norm_radial(grid, q, nfun)
    const_norm_free(grid, xygrid, nfun)
    const_norm_free2D(grid, xgrid, nfun)

Make the normalisation factor ([`FreeSpaceNorm`](@ref)) for a `z`-independent refractive
index, computed once and reused. `nfun(λ)` takes wavelength; for crystal optics pass
`(nfunx, nfuny)` with `nfunx(λ, δθ)` and `nfuny(λ)`.
"""
const_norm_radial(grid, rg::Grid.RadialGrid, nfun; spec=HostSpec()) =
    FreeSpaceNorm(grid, rg, _zfun(nfun); constant=true, spec)
const_norm_radial(grid, q::Grid.HankelTransform, nfun; kwargs...) =
    const_norm_radial(grid, Grid.RadialGrid(q), nfun; kwargs...)
const_norm_free(grid, xygrid::Grid.FreeGrid, nfun; spec=HostSpec()) =
    FreeSpaceNorm(grid, xygrid, _zfun(nfun); constant=true, spec)
const_norm_free2D(grid, xgrid::Grid.Free2DGrid, nfun; spec=HostSpec()) =
    FreeSpaceNorm(grid, xgrid, _zfun(nfun); constant=true, spec)

# wrap a z-independent index function in the (ω; z) / (λ, δθ; z), (λ; z) forms
_zfun(nfun) = (ω; z) -> nfun(wlfreq(ω))
function _zfun(nfuns::Tuple)
    nfunx, nfuny = nfuns
    ((λ, δθ; z) -> nfunx(λ, δθ), (λ; z) -> nfuny(λ))
end

#=================================================#
#======  CARTESIAN FREE-SPACE TRANSFORMS  ========#
#=================================================#

#= `TransFree` (3-D, `(t, x, y)`) and `TransFree2D` (2-D, `(t, x)`) are the same transform
   in two dimensionalities: one multi-axis FFT takes the state from `(ω, k⊥)` to
   `(t, r⊥)` and back, with the polarisation axis skipped. They are separate types because
   `Boundaries.spacegrid` and `Luna.setup` dispatch on them and because the transform
   region differs; everything else below is written once for both. =#

"""
    freebuffers(spec, TT, ωshape, tshape, Eto)

The buffers a Cartesian free-space transform holds: `(Eωo, Eto, Pto)`, on `spec`'s array
type and precision.

**Three** field-sized arrays, not four: the oversampled frequency-domain buffer does
double duty as the field's (`Eωo`) and the nonlinear polarisation's (`Pωo`). `to_time!`
is the only reader of the field's copy and it has finished with it before `to_freq!`
writes the polarisation's, so aliasing them removes one field-sized allocation per
transform (GPU_PLAN.md §4.4). The two never appear as the input and the output of the
same FFT call, which a device plan would reject.

`Eto`, if given, is **taken over** as the oversampled time-domain buffer rather than
allocated here: `Luna.setup` has to allocate one block of exactly that shape and type to
plan the forward transform against, and handing it to the transform instead of leaving it
to the garbage collector is one oversampled block less peak memory at setup (64 MB on the
3-D example's grid). It is zero-filled here because FFTW's planning modes other than
`:estimate` write into the array they plan against; nothing reads it before `to_time!`
overwrites it, so this is only for reproducibility with the allocated path.
"""
function freebuffers(spec, ::Type{TT}, ωshape, tshape, Eto) where {TT}
    CT = Complex{realtype(spec)}
    Eωo = alloc(spec, CT, ωshape)
    if isnothing(Eto)
        Eto = alloc(spec, TT, tshape)
    else
        (eltype(Eto) === TT && size(Eto) == tshape) || error(
            "the time-domain buffer handed to this transform is a $(typeof(Eto)) of "*
            "size $(size(Eto)), but it needs a $TT array of size $tshape.")
        fill!(Eto, zero(TT))
    end
    Pto = similar(Eto)
    Eωo, Eto, Pto
end

#= The oversampled real-space noise field, shared by both Cartesian transforms: ω -> t on
   the oversampled grid and k⊥ -> r⊥ in the same multi-axis inverse transform, exactly the
   transform the field itself takes on every step. The noise is a state-unit quantity, so
   it carries the same 1/Eref the state does. =#
function freenoise(spec, ::Type{TT}, noise_field, grid, ωshape, tshape, IFT,
                   scaling) where {TT}
    isnothing(noise_field) && return nothing, nothing
    size(noise_field) == ωshape || error(
        "the noise field is $(size(noise_field)), but this transform's frequency-domain "*
        "shape is $ωshape. Build it with the state's shape, polarisation axis included.")
    CT = Complex{realtype(spec)}
    Eωo_noise = alloc(spec, CT, ωshape_oversampled(grid, ωshape))
    Et_noise = alloc(spec, TT, tshape)
    nf = isunity(scaling) ? noise_field : noise_field ./ scaling.Eref
    to_time!(Et_noise, todevice(spec, nf), Eωo_noise, IFT)
    Et_noise, alloc(spec, TT, tshape)
end

ωshape_oversampled(grid, ωshape) = (length(grid.ωo), ωshape[2:end]...)

"""
    TransFree

Transform E(ω) -> Pₙₗ(ω) for 3-D free-space propagation on a
[`Grid.FreeGrid`](@ref Luna.Grid.FreeGrid).

Parametric in its buffer array type, like [`TransRadial`](@ref): every buffer and grid
mirror is allocated with [`Luna.alloc`](@ref)/[`Luna.todevice`](@ref) from the run's
[`Luna.DeviceSpec`](@ref), so the same code runs on the host, in reduced precision and on
a device. The `(t, x, y)` <-> `(ω, kx, ky)` transform is a single plan over the region
`(1, 3, 4)` on every backend.

# Fields
- `gv`: mirror of the grid vectors the kernels broadcast against.
- `prefac`: the z-independent part of the frequency-domain normalisation,
  `ωwin·(-iω)·Pref`, precombined on the host.
- `Eωo`, `Pωo`: the oversampled frequency-domain buffer, held under both names because
  they are the **same array** (see [`freebuffers`](@ref)).
- `Et_noise`: precomputed time-domain noise on the oversampled real-space grid
  `(nto, npol, nx, ny)` for the modified shot-noise model, or `nothing`.
- `Et_nl`: preallocated buffer for the combined field + noise, passed to `Et_to_Pt!`. The
  propagating field (`Eto`) is never modified.
- `scaling`: the [`Luna.UnitScaling`](@ref) the state and the polarisation are in.
"""
struct TransFree{TT, ωT, FTT, IFTT, nT, rT, gT, gvT, xygT, dT, iT, pT, eT, nlT}
    FT::FTT # joint (t, x, y) -> (ω, kx, ky) transform, region (1, 3, 4)
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    normfun::nT # callable returning the normalisation factor
    resp::rT # nonlinear responses (tuple of callables)
    grid::gT # host grid, for metadata and for anything not in a kernel
    gv::gvT # mirror of the grid vectors the kernels broadcast against
    xygrid::xygT # transverse grid
    densityfun::dT # callable which returns density
    Pto::TT # buffer for oversampled time-domain NL polarisation
    Eto::TT # buffer for oversampled time-domain field
    Eωo::ωT # buffer for oversampled frequency-domain field
    Pωo::ωT # === Eωo: the same buffer under the name the polarisation pass uses
    idcs::iT # CartesianIndices for Et_to_Pt! to iterate over
    prefac::pT # ωwin*(-im*ω)*Pref: the z-independent normalisation factor
    Et_noise::eT # time-domain noise for modified shot-noise model, or nothing
    Et_nl::nlT # buffer for field+noise passed to Et_to_Pt!, or nothing
    scaling::UnitScaling # units the state and the polarisation are expressed in
end

function show(io::IO, t::TransFree)
    grid = "grid type: $(typeof(t.grid))"
    samples = "time grid size: $(length(t.grid.t)) / $(length(t.grid.to))"
    resp = "responses: "*join([string(typeof(ri)) for ri in t.resp], "\n    ")
    y = "y grid: $(minimum(t.xygrid.y)) to $(maximum(t.xygrid.y)), N=$(length(t.xygrid.y))"
    x = "x grid: $(minimum(t.xygrid.x)) to $(maximum(t.xygrid.x)), N=$(length(t.xygrid.x))"
    out = join(["TransFree", grid, samples, y, x, resp], "\n  ")
    print(io, out)
end

"""
    TransFree(TT, grid, xygrid, FT, responses, densityfun, normfun, pol=false; kwargs...)
    TransFree(grid, xygrid, FT, responses, densityfun, normfun, pol=false; kwargs...)

Construct a `TransFree` to calculate the reciprocal-domain nonlinear polarisation for 3-D
free-space propagation. `TT` is the time-domain element type (`Float64`/`Float32` on a
`Grid.RealGrid`, complex on a `Grid.EnvGrid`); the form without it picks it from the grid.
`FT` is the oversampled `(t, x, y) -> (ω, kx, ky)` plan.

# Keyword arguments
- `noise_field=nothing`: optional `(nω, npol, nx, ny)` frequency/k-space noise field for
  the modified shot-noise model. When provided, it is converted to the oversampled
  real-space time domain and stored as `Et_noise`. Generate with
  [`Fields.generate_noise_field`](@ref Luna.Fields.generate_noise_field).
- `spec=HostSpec()`: the [`Luna.DeviceSpec`](@ref) the buffers, mirrors and responses live
  on.
- `scaling=UNIT_SCALING`: the [`Luna.UnitScaling`](@ref) the state and the nonlinear
  polarisation are expressed in. The responses are converted to it with
  [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale), the noise field is divided by `Eref`
  and `Pref` is folded into `prefac`.
- `Eto=nothing`: an oversampled time-domain block of exactly the transform's shape and
  element type, **taken over** as its `Eto` rather than allocated ([`freebuffers`](@ref)).
  `Luna.setup` passes the array it planned `FT` against, so that block is not left to the
  garbage collector; the caller must not use it for anything else afterwards.
"""
function TransFree(TT, grid, xygrid::Grid.FreeGrid, FT, responses, densityfun, normfun,
                   pol=false; noise_field=nothing, spec=HostSpec(), scaling=UNIT_SCALING,
                   Eto=nothing)
    np = pol ? 2 : 1
    Nx, Ny = length(xygrid.x), length(xygrid.y)
    ωshape = (length(grid.ω), np, Nx, Ny)
    tshape = (length(grid.to), np, Nx, Ny)
    IFT = Utils.plan_ift(FT)
    Eωo, Eto, Pto = freebuffers(spec, TT, ωshape_oversampled(grid, ωshape), tshape,
                                Eto)
    gv = gridvectors(grid, spec)
    prefac = fsprefac(grid, spec, scaling)
    Et_noise, Et_nl = freenoise(spec, TT, noise_field, grid, ωshape, tshape, IFT, scaling)
    responses = Nonlinear.rescale_responses(Tuple(responses), spec, scaling, Eto)
    check_norm(normfun, spec, scaling)
    resparrays = Nonlinear.resident_arrays_all(responses)
    assert_resident(spec, Eωo, Eto, Pto, prefac, gv.ω, gv.ωwin, gv.twin, gv.towin,
                    gv.sidx, Et_noise, Et_nl, resparrays...)
    TransFree(FT, IFT, normfun, responses, grid, gv, xygrid, densityfun, Pto, Eto,
              Eωo, Eωo, CartesianIndices((Nx, Ny)), prefac, Et_noise, Et_nl, scaling)
end

function TransFree(grid::Grid.RealGrid, args...; spec=HostSpec(), kwargs...)
    TransFree(realtype(spec), grid, args...; spec, kwargs...)
end

function TransFree(grid::Grid.EnvGrid, args...; spec=HostSpec(), kwargs...)
    TransFree(Complex{realtype(spec)}, grid, args...; spec, kwargs...)
end

"""
    (t::TransFree)(nl, Eωk, z)

Calculate the reciprocal-domain (ω-kx-ky-space) nonlinear response due to the field `Eωk`
and place the result in `nl`.
"""
(t::TransFree)(nl, Eωk, z) = freetransform!(t, nl, Eωk, z)

"""
    TransFree2D

Transform E(ω) -> Pₙₗ(ω) for 2-D (x-z) free-space propagation on a
[`Grid.Free2DGrid`](@ref Luna.Grid.Free2DGrid). The 2-D counterpart of
[`TransFree`](@ref), with the same fields and the same `Pωo === Eωo` aliasing; its
transform region is `(1, 3)`.
"""
struct TransFree2D{TT, ωT, FTT, IFTT, nT, rT, gT, gvT, xgT, dT, iT, pT, eT, nlT}
    FT::FTT # joint (t, x) -> (ω, kx) transform, region (1, 3)
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    normfun::nT # callable returning the normalisation factor
    resp::rT # nonlinear responses (tuple of callables)
    grid::gT # host grid, for metadata and for anything not in a kernel
    gv::gvT # mirror of the grid vectors the kernels broadcast against
    xgrid::xgT # transverse grid
    densityfun::dT # callable which returns density
    Pto::TT # buffer for oversampled time-domain NL polarisation
    Eto::TT # buffer for oversampled time-domain field
    Eωo::ωT # buffer for oversampled frequency-domain field
    Pωo::ωT # === Eωo: the same buffer under the name the polarisation pass uses
    idcs::iT # CartesianIndices for Et_to_Pt! to iterate over
    prefac::pT # ωwin*(-im*ω)*Pref: the z-independent normalisation factor
    Et_noise::eT # time-domain noise for modified shot-noise model, or nothing
    Et_nl::nlT # buffer for field+noise passed to Et_to_Pt!, or nothing
    scaling::UnitScaling # units the state and the polarisation are expressed in
end

function show(io::IO, t::TransFree2D)
    grid = "grid type: $(typeof(t.grid))"
    samples = "time grid size: $(length(t.grid.t)) / $(length(t.grid.to))"
    resp = "responses: "*join([string(typeof(ri)) for ri in t.resp], "\n    ")
    x = "x grid: $(minimum(t.xgrid.x)) to $(maximum(t.xgrid.x)), N=$(length(t.xgrid.x))"
    out = join(["TransFree2D", grid, samples, x, resp], "\n  ")
    print(io, out)
end

"""
    TransFree2D(TT, grid, xgrid, FT, responses, densityfun, normfun, pol=false; kwargs...)
    TransFree2D(grid, xgrid, FT, responses, densityfun, normfun, pol=false; kwargs...)

Construct a `TransFree2D` to calculate the reciprocal-domain nonlinear polarisation for
2-D free-space propagation. Arguments and keyword arguments are those of
[`TransFree`](@ref), with `xgrid` a [`Grid.Free2DGrid`](@ref Luna.Grid.Free2DGrid), `FT`
the oversampled `(t, x) -> (ω, kx)` plan and `noise_field` of shape `(nω, npol, nx)`.
"""
function TransFree2D(TT, grid, xgrid::Grid.Free2DGrid, FT, responses, densityfun, normfun,
                     pol=false; noise_field=nothing, spec=HostSpec(),
                     scaling=UNIT_SCALING, Eto=nothing)
    np = pol ? 2 : 1
    Nx = length(xgrid.x)
    ωshape = (length(grid.ω), np, Nx)
    tshape = (length(grid.to), np, Nx)
    IFT = Utils.plan_ift(FT)
    Eωo, Eto, Pto = freebuffers(spec, TT, ωshape_oversampled(grid, ωshape), tshape,
                                Eto)
    gv = gridvectors(grid, spec)
    prefac = fsprefac(grid, spec, scaling)
    Et_noise, Et_nl = freenoise(spec, TT, noise_field, grid, ωshape, tshape, IFT, scaling)
    responses = Nonlinear.rescale_responses(Tuple(responses), spec, scaling, Eto)
    check_norm(normfun, spec, scaling)
    resparrays = Nonlinear.resident_arrays_all(responses)
    assert_resident(spec, Eωo, Eto, Pto, prefac, gv.ω, gv.ωwin, gv.twin, gv.towin,
                    gv.sidx, Et_noise, Et_nl, resparrays...)
    TransFree2D(FT, IFT, normfun, responses, grid, gv, xgrid, densityfun, Pto, Eto,
                Eωo, Eωo, CartesianIndices((Nx,)), prefac, Et_noise, Et_nl, scaling)
end

function TransFree2D(grid::Grid.RealGrid, args...; spec=HostSpec(), kwargs...)
    TransFree2D(realtype(spec), grid, args...; spec, kwargs...)
end

function TransFree2D(grid::Grid.EnvGrid, args...; spec=HostSpec(), kwargs...)
    TransFree2D(Complex{realtype(spec)}, grid, args...; spec, kwargs...)
end

"""
    (t::TransFree2D)(nl, Eωk, z)

Calculate the reciprocal-domain (ω-kx-space) nonlinear response due to the field `Eωk`
and place the result in `nl`.
"""
(t::TransFree2D)(nl, Eωk, z) = freetransform!(t, nl, Eωk, z)

"""
    freetransform!(t, nl, Eωk, z)

One right-hand side of a Cartesian free-space transform: the joint inverse transform to
`(t, r⊥)`, the nonlinear responses, the temporal apodisation, the joint forward transform
and the frequency-domain normalisation. Written once for [`TransFree`](@ref) and
[`TransFree2D`](@ref), which differ only in how many transverse axes they have.
"""
function freetransform!(t, nl, Eωk, z)
    to_time!(t.Eto, Eωk, t.Eωo, t.IFT) # (ω, k⊥) -> (t, r⊥)
    #= Modified shot-noise: field+noise goes into a separate buffer, so the propagating
       field (`Eto`) is never contaminated. With no noise field `Eto` is passed straight
       through, with no copy. =#
    if !isnothing(t.Et_noise)
        @. t.Et_nl = t.Eto + t.Et_noise
        Et_to_Pt!(t.Pto, t.Et_nl, t.resp, t.densityfun(z), t.idcs; scaling=t.scaling)
    else
        Et_to_Pt!(t.Pto, t.Eto, t.resp, t.densityfun(z), t.idcs; scaling=t.scaling)
    end
    @. t.Pto *= t.gv.towin # apodisation
    to_freq!(nl, t.Pωo, t.Pto, t.FT) # (t, r⊥) -> (ω, k⊥)
    fsnorm!(nl, t.prefac, t.normfun(z))
end

#= reflength! on a transform forwards to its normalisation; defined here, after every
   transform type exists. Transforms which are not free-space have nothing to taper. =#
reflength!(t::Union{TransRadial, TransFree, TransFree2D}, ℓ; kwargs...) = reflength!(t.normfun, ℓ; kwargs...)
reflength!(t, ℓ; kwargs...) = nothing


#=================================================#
#=========  TABULATION OF THE TRANSFORM  =========#
#=================================================#

"""
    tabulate(transform, z0, z1, tol, proto)

The transform with every z-dependent host quantity it evaluates inside the right-hand side
replaced by a table over `[z0, z1]` built to relative tolerance `tol`, on the array type
and precision of the propagating field `proto`.

[`Luna.run`](@ref) calls this when `tabulate_linop=true`. The generic method returns the
transform unchanged: only the mode-averaged transform has such quantities (the propagation
constant `β(z)` and the effective area `Aeff(z)`), and only it is device-capable so far.
A transform which is already z-independent gets a two-node table, which costs nothing and
keeps one code path.

The same generic method covers a caller-supplied normalisation (`Luna.setup`'s `norm!`
keyword), which is passed through unchanged rather than inspected: only Luna's own
normalisations are known to have an `aeff` to tabulate. It takes the `aeff` keyword the
known ones take so that the call site does not have to know which it has.
"""
tabulate(t, z0, z1, tol, proto; aeff=nothing) = t

function tabulate(t::TransModeAvg, z0, z1, tol, proto)
    atab = _aefftab(t.aeff, nothing, z0, z1, tol)
    #= `Luna.setup` passes the same `aeff` callable to the transform and to the
       normalisation, so one table serves both and `Modes.Aeff` is called once per node
       rather than twice. A caller who supplied two different ones, or a normalisation
       which is not one of Luna's, gets `nothing` and tabulates its own. =#
    shared = t.aeff === _normaeff(t.norm!) ? atab : nothing
    TransModeAvg(t.Pto, t.Eto, t.Eωo, t.Pωo, t.FT, t.IFT, t.resp, t.grid, t.gv,
                 t.densityfun, tabulate(t.norm!, z0, z1, tol, proto; aeff=shared), atab,
                 t.Et_noise, t.Et_nl, t.scaling)
end

#= The `Aeff` callable of a normalisation, or `nothing` for one Luna did not build: a
   caller-supplied `norm!` (`Luna.setup`'s keyword, `NonlinearRHS.check_norm`) need not
   have the field at all, and reading it would abort the run before dispatch could help. =#
_normaeff(n::Union{NormModeAvg, NormModeAvgGNLSE}) = n.aeff
_normaeff(n) = nothing

"""
    _aefftab(f, shared, z0, z1, tol)

The `Aeff` table to use: `shared` if one has already been built for this transform,
otherwise a new one over `[z0, z1]`.

A table which already covers the span is kept as it is. One which does not is rebuilt from
its own source callable rather than from itself: `Luna.prop_capillary` tabulates `Aeff` over
`[0, flength]` so that the statistics hold a table (see [`Luna.run`](@ref)'s
`tabulate_linop`), and the propagation needs it up to one step past the end of the fibre.
Rebuilding through the interpolant instead would place nodes at every kink of the
interpolant it was reading.
"""
_aefftab(f, ::Nothing, z0, z1, tol) = LinearOps.TabulatedScalar(f, z0, z1; tol)
_aefftab(f, tab, z0, z1, tol) = tab
_aefftab(f::LinearOps.TabulatedScalar, ::Nothing, z0, z1, tol) =
    (f.z[1] <= z0 && f.z[end] >= z1) ? f :
        LinearOps.TabulatedScalar(f.src, z0, z1; tol)

#= With `constβ` the propagation constant is already folded into `pre` and there is no
   `βfun!` to tabulate; only the effective area is left. =#
tabulate(n::NormModeAvg{vT, mT, Nothing}, z0, z1, tol, proto; aeff=nothing) where {vT, mT} =
    NormModeAvg(n.pre, n.mask, n.β, n.βfun!, _aefftab(n.aeff, aeff, z0, z1, tol), n.scaling)

function tabulate(n::NormModeAvg, z0, z1, tol, proto; aeff=nothing)
    β = LinearOps.TabulatedVector(n.βfun!, proto, length(n.mask), z0, z1; tol)
    NormModeAvg(n.pre, n.mask, β, n.βfun!, _aefftab(n.aeff, aeff, z0, z1, tol), n.scaling)
end

tabulate(n::NormModeAvgGNLSE, z0, z1, tol, proto; aeff=nothing) =
    NormModeAvgGNLSE(n.pre, n.mask, _aefftab(n.aeff, aeff, z0, z1, tol), n.scaling)

end
