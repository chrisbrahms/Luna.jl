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
import Luna.LinearOps
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
    ModifiedNoise(Eω_noise, phase, grid)

The noise field of the modified shot-noise model (Chen & Wise, arXiv:2410.20567) together
with the propagator which advances it in `z`.

`Eω_noise` is the one-photon-per-mode spectral draw from
[`Fields.generate_noise_field`](@ref Luna.Fields.generate_noise_field), made once before
propagation and never re-randomised. `phase` is an
[`AbstractUnitaryPhase`](@ref Luna.LinearOps.AbstractUnitaryPhase), created with
[`LinearOps.unitary_phase`](@ref Luna.LinearOps.unitary_phase), which supplies the
accumulated phase `Φ(ω, z) = ∫₀ᶻ imag(linop) dz'`.

The noise is advanced by the **unitary** part of the linear operator only:

- **dispersion — yes.** It is unitary, and it is the phase evolution of the medium's own
  vacuum modes. Without it the model is not related to the traditional (`:input`) model by
  the identity `A_s = A_s' + A_noise` on which it rests, and Kerr/FWM magnitudes are wrong by
  up to ~6× per frequency (much more in conservative propagation).
- **loss — no.** Fluctuation–dissipation: loss replenishes vacuum rather than removing it, so
  the local noise level stays at one photon per mode. This is also what allows the model to
  seed spontaneous generation in absorptive media.
- **gain — no.** The source term already injects fresh vacuum at every `z`, which the
  homogeneous propagator amplifies over the remaining length — that integral *is* the ASE.
  Amplifying the noise as well double-counts it.

Because the field is drawn once and only ever propagated deterministically, results remain
step-size independent, which is what the "constant noise field" requirement of the paper
actually protects.

The draw is multiplied by `grid.ωwin` on construction, matching the windowing the traditional
input shot noise receives from the propagation loop.

In free space the draw is restricted to the **propagating** components: the one-photon
vacuum level is that of the radiating modes, and an evanescent transverse wavenumber
(`k⊥ > k(ω)`, where the linear operator is a real decay rather than a phase) is not one.
Left in, those components would keep their full amplitude for the whole propagation (the
unitary phase does not damp them, while the field there decays) and feed the nonlinear
polarisation with a noise floor that does not exist. The free-space transforms zero them on
construction ([`mask_evanescent!`](@ref)), with `k(ω)` at `z = 0`.
"""
struct ModifiedNoise{aT, pT, bT}
    Eω0::aT # z = 0 spectral draw (windowed)
    phase::pT # accumulated unitary phase of the linear operator
    φ::bT # buffer for Φ(ω, z)
    Eω::aT # buffer holding the propagated spectral noise
    lastz::Base.RefValue{Float64}
end

function ModifiedNoise(Eω_noise, phase, grid)
    Eω0 = Eω_noise .* grid.ωwin
    ModifiedNoise(Eω0, phase, zeros(Float64, size(Eω0)), similar(Eω0), Ref(NaN))
end

"""
    noise_at!(noise::ModifiedNoise, z)

Fill and return the propagated spectral noise field at `z`. The result aliases an internal
buffer, so copy it if it needs to outlive the next call.
"""
function noise_at!(n::ModifiedNoise, z)
    if z != n.lastz[]
        n.phase(n.φ, z)
        @. n.Eω = n.Eω0 * cis(n.φ)
        n.lastz[] = z
    end
    n.Eω
end

"""
    noise_at(noise::ModifiedNoise, z)

As [`noise_at!`](@ref) but returns a fresh array.
"""
noise_at(n::ModifiedNoise, z) = copy(noise_at!(n, z))

"""
    retarget_noise(noise, spec, scaling)

The [`ModifiedNoise`](@ref) a transform's kernels advance on `spec`, in `scaling`'s units:
the draw divided by `scaling.Eref` (the noise is a state-unit quantity, like the field) and
converted to `spec`'s complex type and array type, with the phase tables converted alike
(see `retarget_phase`). On the default path -- host arrays, `Float64`, unit scaling -- it is
`noise` itself, so that path is unchanged. `nothing` stays `nothing`.

The transform keeps the caller's object too, in physical units: that is what is saved to the
output, read by [`Processing.withnoise`](@ref Luna.Processing.withnoise) and, on a resumed
cached run, overwritten with the recorded draw (see [`sync_noise!`](@ref)).

In `Float32` the accumulated phase `Φ(z)` is formed from the absolute `z`, so it carries an
absolute error of about `6e-8·|Φ|` rad. The field's propagator does not: a constant operator
is applied over a step's `Δz` only, and a tabulated one stores `Φ` with its secant removed
(`LinearOps.TabulatedLinop`), whereas the noise's phase table holds the absolute `Φ`; the
same secant subtraction would remove this error, if it ever matters. How large `|Φ|` gets
depends on the operator. In a waveguide's co-moving frame it is the dispersion residual, tens to thousands of rad, so ~1e-4 rad. In free space it is
`(β_z − β₁(ω − ω₀) − β_ref)·z`, which near cutoff (`β_z → 0`) approaches `k(ω)·z`: ~4e6 rad
over 0.5 m at 800 nm, i.e. an error of order 0.3 rad on the high-k⊥ components (on axis it is
again only the dispersion residual). And the one-photon noise is itself below the field's
`Float32` resolution near the pulse peak.
"""
retarget_noise(::Nothing, spec, scaling) = nothing

function retarget_noise(n::ModifiedNoise, spec, scaling)
    (arraytype(spec) === Array && realtype(spec) === Float64 && isunity(scaling)) && return n
    T = realtype(spec)
    Eω0 = todevice(spec, Complex{T}.(n.Eω0 ./ scaling.Eref))
    ModifiedNoise(Eω0, retarget_phase(n.phase, spec), alloc(spec, T, size(Eω0)),
                  similar(Eω0), Ref(NaN))
end

"The phase of a `ModifiedNoise` on `spec`'s array type and real type."
retarget_phase(p::LinearOps.ConstUnitaryPhase, spec) =
    LinearOps.ConstUnitaryPhase(todevice(spec, realtype(spec).(p.dφ)))

retarget_phase(p::LinearOps.TabulatedUnitaryPhase, spec) =
    LinearOps.TabulatedUnitaryPhase(p.z, todevice(spec, realtype(spec).(p.Φ)),
                                    todevice(spec, realtype(spec).(p.dΦ)), p.atol, p.err)

"The arrays a `ModifiedNoise` holds, for `assert_resident`; none for `nothing`."
noise_arrays(::Nothing) = ()
noise_arrays(n::ModifiedNoise) = (n.Eω0, n.φ, n.Eω, phase_arrays(n.phase)...)
phase_arrays(p::LinearOps.ConstUnitaryPhase) = (p.dφ,)
phase_arrays(p::LinearOps.TabulatedUnitaryPhase) = (p.Φ, p.dΦ)

"""
    sync_noise!(transform)

After the caller's `transform.noise` has been changed in place (a resumed cached run adopts
the draw recorded in the output), bring the transform's kernel copy `dnoise` up to date.
Nothing to do when the two are the same object, or for a transform without noise.
"""
function sync_noise!(t)
    hasfield(typeof(t), :noise) || return t
    n = t.noise
    isnothing(n) && return t
    n.lastz[] = NaN
    #= a free-space transform's draw is restricted to the propagating components
       (`mask_evanescent!`, a no-op elsewhere), including one adopted from an output =#
    hasfield(typeof(t), :normfun) && mask_evanescent!(n, t.normfun)
    hasfield(typeof(t), :dnoise) || return t
    d = t.dnoise
    d === n && return t
    # host -> device copy of the draw in the transform's units (copyto! handles the transfer)
    copyto!(d.Eω0, eltype(d.Eω0).(n.Eω0 ./ t.scaling.Eref))
    d.lastz[] = NaN
    t
end

"""
    TransModal

Transform E(ω) -> Pₙₗ(ω) for multimode propagation via spatial integration.

# Fields
- `noise`: [`ModifiedNoise`](@ref) for the modified shot-noise model, or `nothing`. It is
  advanced to the current `z` once per call, in `reset!`, and then projected to real space at
  each integration point and combined with the field in a separate buffer (`Er_nl`) for
  nonlinear evaluation. The propagating field (`Er`) is never modified.
- `Er_noise`: preallocated buffer for the real-space time-domain noise, same shape as `Er`.
- `Er_nl`: preallocated buffer for the combined field + noise, passed to `Et_to_Pt!`.
"""
mutable struct TransModal{tsT, lT, TT, FTT, IFTT, rT, gT, dT, ddT, nT, nsT, enT, enlT}
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
    noise::nsT # ModifiedNoise for the modified shot-noise model, or nothing
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
    TransModal(grid, ts, FT, resp, densityfun, norm!; rtol=1e-3, atol=0.0, mfcn=300, full=false, noise=nothing)

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
- `noise=nothing` : optional [`ModifiedNoise`](@ref) for the modified shot-noise model, whose
  `(nω, nmodes)` draw holds independent one-photon-per-mode noise in each mode column.
"""
function TransModal(tT, grid, ts::Modes.ToSpace, FT, resp, densityfun, norm!;
                    rtol=1e-3, atol=0.0, mfcn=512, full=false, noise=nothing)
    Emω = Array{ComplexF64,2}(undef, length(grid.ω), ts.nmodes)
    Erω = Array{ComplexF64,2}(undef, length(grid.ω), ts.npol)
    Erωo = Array{ComplexF64,2}(undef, length(grid.ωo), ts.npol)
    Er = Array{tT,2}(undef, length(grid.to), ts.npol)
    Pr = Array{tT,2}(undef, length(grid.to), ts.npol)
    Prω = Array{ComplexF64,2}(undef, length(grid.ω), ts.npol)
    Prωo = Array{ComplexF64,2}(undef, length(grid.ωo), ts.npol)
    Prmω = Array{ComplexF64,2}(undef, length(grid.ω), ts.nmodes)
    IFT = Utils.plan_ift(FT)
    # For the modified shot-noise model, allocate a buffer for the modal noise field at the
    # current z and one for the real-space time-domain noise. The noise is projected to space
    # at each integration point in Erω_to_Prω!, so we hold it in the modal domain.
    if !isnothing(noise)
        Er_noise = Array{tT,2}(undef, length(grid.to), ts.npol)
        Er_nl = Array{tT,2}(undef, length(grid.to), ts.npol)
    else
        Er_noise = nothing
        Er_nl = nothing
    end
    TransModal(ts, full, Modes.dimlimits(ts.ms[1]), Emω, Erω, Erωo, Er, Pr, Prω, Prωo, Prmω,
               FT, IFT, resp, grid, densityfun, densityfun(0.0), norm!, 0, 0.0, rtol, atol, mfcn,
               similar(Prmω), noise, Er_noise, Er_nl)
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
    # Advance the modified shot-noise field to this z. Done once per RHS evaluation rather
    # than at every cubature point, which is why this path costs essentially nothing.
    isnothing(t.noise) || noise_at!(t.noise, z)
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
    # Modified shot-noise model: project the noise modes (already advanced to t.z by reset!)
    # to real space at this spatial point, convert to the oversampled time domain, and
    # combine with the field in a separate buffer (Er_nl) so the propagating field (Er) is
    # never contaminated.
    if !isnothing(t.noise)
        Modes.to_space!(t.Erω, t.noise.Eω, x, t.ts, z=t.z)
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
- `noise`: [`ModifiedNoise`](@ref) for the modified shot-noise model (Chen & Wise,
  arXiv:2410.20567), or `nothing` for the traditional model.
- `Eωo_noise`/`Et_noise`: buffers holding the noise field at the current `z`, on the
  oversampled frequency and time grids.
- `Et_nl`: preallocated buffer for the combined field + noise. When `noise` is present,
  `Et_nl = Eto + Et_noise` is computed at each step and passed to `Et_to_Pt!`. The
  propagating field (`Eto`) is never modified: the noise feels only the unitary part of the
  linear operator, and never the loss.
"""
struct TransModeAvg{TT, ωT, FTT, IFTT, rT, gT, gvT, dT, nT, aT, nsT, dnT, eoT, eT, nlT}
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
    noise::nsT # ModifiedNoise for the modified shot-noise model (host, physical units), or nothing
    dnoise::dnT # the copy the kernels advance: `noise` itself on the host in Float64, else on `spec`
    Eωo_noise::eoT # buffer for the oversampled spectral noise at the current z, or nothing
    Et_noise::eT # buffer for the time-domain noise at the current z, or nothing
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
- `noise=nothing`: optional [`ModifiedNoise`](@ref) for the modified shot-noise model. When
  provided, the noise field is advanced to the current `z`, converted to the oversampled time
  grid, and injected into the nonlinear operator at every propagation step. The transform
  keeps the caller's object as `noise` (host, physical units: what is saved and what
  [`Processing.withnoise`](@ref Luna.Processing.withnoise) needs) and advances
  `dnoise = retarget_noise(noise, spec, scaling)`, which is the same object on the default host
  `Float64` path and otherwise a copy on `spec` in `scaling`'s units.
- `spec=HostSpec()`: the [`Luna.DeviceSpec`](@ref) the buffers, mirrors and responses live
  on.
- `scaling=UNIT_SCALING`: the [`Luna.UnitScaling`](@ref) the state and the nonlinear
  polarisation are expressed in. The responses are converted to it with
  [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale) and the noise draw is divided by
  `Eref`.
"""
function TransModeAvg(TT, grid, FT, IFT, resp, densityfun, norm!, aeff;
                      noise=nothing, spec=HostSpec(), scaling=UNIT_SCALING)
    nωo = length(grid.ωo)
    nto = length(grid.to)
    CT = Complex{realtype(spec)}
    Eωo = alloc(spec, CT, (nωo,))
    Eto = alloc(spec, TT, (nto,))
    Pto = similar(Eto)
    Pωo = similar(Eωo)
    gv = gridvectors(grid, spec)
    dnoise = retarget_noise(noise, spec, scaling)
    if !isnothing(noise)
        Eωo_noise = alloc(spec, CT, (nωo,))
        Et_noise = alloc(spec, TT, (nto,))
        Et_nl = alloc(spec, TT, (nto,))
    else
        Eωo_noise = nothing
        Et_noise = nothing
        Et_nl = nothing
    end
    #= The four-argument form: `Eto` is the prototype of the block the responses are
       called with, so one which owns buffers can allocate them here, in the run's array
       type, in time for the residency assertion below. =#
    resp = map(r -> Nonlinear.rescale(r, spec, scaling, Eto), Tuple(resp))
    #= Every mirror the transform holds and every array its responses carry, not only the
       ones this transform's own kernels touch: the assertion is what catches a future
       mistake, so it has to cover everything. =#
    resparrays = reduce((a, r) -> (a..., Nonlinear.resident_arrays(r)...), resp; init=())
    assert_resident(spec, Eωo, Eto, Pto, Pωo, gv.ω, gv.ωwin, gv.twin, gv.towin, gv.sidx,
                    Eωo_noise, Et_noise, Et_nl, noise_arrays(dnoise)..., resparrays...)
    TransModeAvg(Pto, Eto, Eωo, Pωo, FT, IFT, resp, grid, gv, densityfun, norm!, aeff,
                 noise, dnoise, Eωo_noise, Et_noise, Et_nl, scaling)
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
    # Modified shot-noise model: advance the noise to this z, then compute field+noise in a
    # separate buffer (Et_nl) so that the propagating field (Eto) is never contaminated. The
    # noise is scaled by the same normalisation factor (nlscale × √Aeff) so it enters in
    # physical units.
    if !isnothing(t.noise)
        to_time!(t.Et_noise, noise_at!(t.dnoise, z), t.Eωo_noise, t.IFT)
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
pressure gradient needs) `βfun!` is called on every evaluation and its result uploaded;
that is the interim arrangement until `gpu/23` tabulates it.
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

function (n::NormModeAvg)(nl, z)
    n.βfun!(n.β.host, z)
    β = upload!(n.β)
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

# Fields
- `noise`: [`ModifiedNoise`](@ref) for the modified shot-noise model, or `nothing`.
- `Eωo_noise`/`Et_noise`: buffers holding the noise field at the current `z`, on the
  oversampled frequency/k-space and time/real-space grids `(nωo, np, nr)`/`(nto, np, nr)`.
- `Et_nl`: preallocated buffer for the combined field + noise, passed to `Et_to_Pt!`. The
  propagating field (`Eto`) is never modified.
"""
struct TransRadial{TT, RGT, FTT, IFTT, nT, rT, gT, dT, iT, nsT, eoT, eT, nlT}
    rgrid::RGT # transverse grid (Grid.RadialGrid: space to k-space)
    FT::FTT # Fourier transform (time to frequency)
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    normfun::nT # Function which returns normalisation factor
    resp::rT # nonlinear responses (tuple of callables)
    grid::gT # time grid
    densityfun::dT # callable which returns density
    Pto_r::Array{TT, 3} # Buffer array for NL polarisation on oversampled time grid
    Pto_k::Array{TT, 3} # Buffer array for NL polarisation on oversampled time grid
    Eto_r::Array{TT, 3} # Buffer array for field on oversampled time grid
    Eto_k::Array{TT, 3} # Buffer array for field on oversampled time grid
    Eωo::Array{ComplexF64, 3} # Buffer array for field on oversampled frequency grid
    Pωo::Array{ComplexF64, 3} # Buffer array for NL polarisation on oversampled frequency grid
    idcs::iT # CartesianIndices for Et_to_Pt! to iterate over
    Tfwd::Matrix{TT} # forward Hankel transform matrix
    Tbwd::Matrix{TT} # backward Hankel transform matrix
    noise::nsT # ModifiedNoise for the modified shot-noise model, or nothing
    Eωo_noise::eoT # buffer for the oversampled spectral noise at the current z, or nothing
    Et_noise::eT # buffer for the time-domain noise at the current z, or nothing
    Et_nl::nlT # buffer for field+noise passed to Et_to_Pt!, or nothing
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
    TransRadial(TT, grid, rgrid, FT, responses, densityfun, normfun; noise=nothing)

Construct a `TransRadial` to calculate the reciprocal-domain nonlinear polarisation.
`rgrid` is a [`Grid.RadialGrid`](@ref Luna.Grid.RadialGrid).

# Keyword arguments
- `noise=nothing`: optional [`ModifiedNoise`](@ref) for the modified shot-noise model, whose
  draw is an `(nω, np, nk)` frequency/k-space noise field, the same size as the `Eω` returned
  by [`Luna.setup`](@ref) (`np` polarisation components). At each step it is advanced to the
  current `z` and converted to the real-space time domain `(nto, np, nr)` via inverse FFT and
  the backward Hankel matrix, like the field. [`Fields.generate_noise_field`](@ref
  Luna.Fields.generate_noise_field) makes only `(nω[, nmodes])` draws, so a free-space draw
  is built by the caller.
"""
function TransRadial(TT, grid, rgrid::Grid.RadialGrid, FT, responses, densityfun, normfun,
                     pol=false; noise=nothing)
    noise = mask_evanescent!(noise, normfun) # the vacuum noise of the propagating components
    np = pol ? 2 : 1
    N = rgrid.N
    IFT = Utils.plan_ift(FT)
    Eωo = zeros(ComplexF64, (length(grid.ωo), np, N))
    Eto_r = zeros(TT, (length(grid.to), np, N))
    Pto_r = similar(Eto_r)
    Eto_k = similar(Eto_r)
    Pto_k = similar(Eto_r)
    Pωo = similar(Eωo)
    idcs = CartesianIndices(size(Pto_r)[3:end])
    #= Our own copies of the grid's transform matrices in the type we multiply: a GEMM
       needs both operands in the same element type (and on devices it is required). =#
    Tfwd = convert(Matrix{TT}, rgrid.Tfwd)
    Tbwd = convert(Matrix{TT}, rgrid.Tbwd)
    # buffers for the noise at the current z: (ω, k) -> (t, k) -> (t, r) each step, through
    # the same Tbwd as the field
    if !isnothing(noise)
        Eωo_noise = zeros(ComplexF64, (length(grid.ωo), np, N))
        Et_noise = zeros(TT, (length(grid.to), np, N))
        Et_nl = zeros(TT, (length(grid.to), np, N))
    else
        Eωo_noise = nothing
        Et_noise = nothing
        Et_nl = nothing
    end
    TransRadial(rgrid, FT, IFT, normfun, responses, grid, densityfun, Pto_r, Pto_k, Eto_r, Eto_k, Eωo, Pωo, idcs,
                Tfwd, Tbwd, noise, Eωo_noise, Et_noise, Et_nl)
end

# accept a Hankel.QDHT as before, converting it (with a deprecation warning)
function TransRadial(TT::Type, grid, q::Grid.HankelTransform, args...; kwargs...)
    TransRadial(TT, grid, Grid.RadialGrid(q), args...; kwargs...)
end

function TransRadial(grid::Grid.RealGrid, args...; kwargs...)
    TransRadial(Float64, grid, args...; kwargs...)
end

function TransRadial(grid::Grid.EnvGrid, args...; kwargs...)
    TransRadial(ComplexF64, grid, args...; kwargs...)
end

"""
    (t::TransRadial)(nl, Eω, z)

Calculate the reciprocal-domain (ω-k-space) nonlinear response due to the field `Eω` and
place the result in `nl`
"""
function (t::TransRadial)(nl, Eω, z)
    to_time!(t.Eto_k, Eω, t.Eωo, t.IFT) # transform ω -> t
    # transform Eto k -> r
    # iterate over polarisation directions (either 1:2 or just 1)
    for ip in axes(t.Eto_k, 2)
        mul!(view(t.Eto_r, :, ip, :), view(t.Eto_k, :, ip, :), t.Tbwd)
    end
    # Modified shot-noise: advance the noise to this z, take it ω -> t (into Et_nl, used
    # as scratch for the k-space noise) and k -> r (into Et_noise), then compute
    # field+noise in Et_nl so the propagating field (Eto_r) is never contaminated.
    if !isnothing(t.noise)
        to_time!(t.Et_nl, noise_at!(t.noise, z), t.Eωo_noise, t.IFT)
        for ip in axes(t.Et_nl, 2)
            mul!(view(t.Et_noise, :, ip, :), view(t.Et_nl, :, ip, :), t.Tbwd)
        end
        @. t.Et_nl = t.Eto_r + t.Et_noise
        Et_to_Pt!(t.Pto_r, t.Et_nl, t.resp, t.densityfun(z), t.idcs)
    else
        Et_to_Pt!(t.Pto_r, t.Eto_r, t.resp, t.densityfun(z), t.idcs)
    end
    @. t.Pto_r *= t.grid.towin # apodisation
    # transform Pto r -> k
    for ip in axes(t.Pto_k, 2)
        mul!(view(t.Pto_k, :, ip, :), view(t.Pto_r, :, ip, :), t.Tfwd)
    end
    to_freq!(nl, t.Pωo, t.Pto_k, t.FT) # transform t -> ω
    nl .*= t.grid.ωwin .* (-im.*t.grid.ω)./(2 .* t.normfun(z))
end

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
- `kperp2`, `kidcs`: squared transverse wavevector and the indices of the k axes
- `out`: the normalisation array (`ComplexF64`, since ``\\beta_z`` is complex below cutoff)
- `ℓ`, `κmax`, `kwin`: taper parameters (see above)
- `constant`: if `true`, `out` is computed once and reused (the index does not depend on `z`)
"""
mutable struct FreeSpaceNorm{gT, sT, nT, kT, iT, oT, wT}
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
end

npol(nfun::Tuple, grid) = 2 # crystal optics: (nfunx, nfuny)
npol(nfun, grid) = length(nfun(grid.ω[findfirst(grid.sidx)]; z=0)) # 1 if single index, 2 if nx, ny

function FreeSpaceNorm(grid, spacegrid, nfun; constant)
    kperp2, kidcs = transverse_k2(spacegrid)
    np = npol(nfun, grid)
    out = zeros(ComplexF64, (length(grid.ω), np, size(kidcs)...))
    kwin = ones(Float64, size(kidcs))
    FreeSpaceNorm(grid, spacegrid, nfun, kperp2, kidcs, out, 0.0, Inf, kwin, constant, false)
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

# isotropic: nfun(ω; z) -> n or (nx, ny), the same k⊥ for every polarisation
fillnorm!(nf::FreeSpaceNorm, z) =
    fillβsq!(nf.out, nf, z, (βsq, ω, wk) -> normfactor(nf, βsq, ω, wk), 1)

#= Fill `out` (shaped like `nf.out`, `(nω, npol, k⊥...)`) with `f(βsq, ω, wk)` of the
   longitudinal wavenumber squared `βsq = k(ω)² - k⊥²` at `z`, and with `outside` at ω = 0
   and outside the simulation band. Shared by the normalisation and by the propagating
   mask of the modified shot noise, so that the two agree on which components are
   evanescent. =#
function fillβsq!(out, nf::FreeSpaceNorm, z, f, outside)
    ω = nf.grid.ω
    for ii in nf.kidcs
        for iω in eachindex(ω)
            if ω[iω] == 0 || !nf.grid.sidx[iω]
                out[iω, :, ii] .= outside
                continue
            end
            for (ip, n) in enumerate(nf.nfun(ω[iω]; z))
                βsq = (real(n)*ω[iω]/PhysData.c)^2 - nf.kperp2[ii]
                out[iω, ip, ii] = f(βsq, ω[iω], nf.kwin[ii])
            end
        end
    end
    out
end

#= crystal optics: nfunx(λ, δθ; z) depends on the internal angle, which depends on kx only,
   so the angle is found once per (ω, kx) and reused along ky. For Free2DGrid the k axes are
   (Nkx,) and the trailing ky index below is the (allowed) singleton 1. =#
function fillβsq!(out, nf::FreeSpaceNorm{<:Any, <:Any, <:Tuple}, z, f, outside)
    nfunx, nfuny = nf.nfun
    ω = nf.grid.ω
    kx = nf.spacegrid.kx
    for iω in eachindex(ω)
        if ω[iω] == 0 || !nf.grid.sidx[iω]
            out[iω, :, nf.kidcs] .= outside
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
                out[iω, 1, ix, iy] = f(ksq_xpol - kperp2, ω[iω], wk)
                out[iω, 2, ix, iy] = f(ksq_ypol - kperp2, ω[iω], wk)
            end
        end
    end
    out
end

"""
    propagating_mask(nf::FreeSpaceNorm, z)

`true` at the propagating components of a free-space field, `(nω, npol, k⊥...)`, i.e.
where `β_z² = k(ω)² − k⊥² ≥ 0` at `z` (with the same `k(ω)` the normalisation uses, crystal
optics included), `false` at the evanescent ones. `true` at `ω = 0` and outside the
simulation band, where the windows remove everything anyway.
"""
propagating_mask(nf::FreeSpaceNorm, z) =
    fillβsq!(Array{Bool}(undef, size(nf.out)), nf, z, (βsq, ω, wk) -> βsq >= 0, true)

"""
    mask_evanescent!(noise, normfun)

Restrict the modified shot-noise draw of a free-space propagation to its propagating
components: zero `noise.Eω0` wherever [`propagating_mask`](@ref)`(normfun, 0)` is `false`
(see [`ModifiedNoise`](@ref) for why). Called by the free-space transforms on construction,
so the caller's `ModifiedNoise` is changed in place and the draw recorded in the output
(`"noise_field"`) is the masked one. A no-op without noise or for a normalisation which is
not a [`FreeSpaceNorm`](@ref), where the transverse wavenumbers are not known. With a
z-dependent index the mask is that of `z = 0`.
"""
mask_evanescent!(noise, normfun) = noise

function mask_evanescent!(n::ModifiedNoise, nf::FreeSpaceNorm)
    prop = propagating_mask(nf, 0.0)
    size(prop) == size(n.Eω0) || error(
        "the noise draw is $(size(n.Eω0)), but the free-space field is $(size(prop)). "*
        "Build it with the state's shape, polarisation axis included.")
    n.Eω0 .*= prop
    n.lastz[] = NaN
    n
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
norm_radial(grid, rg::Grid.RadialGrid, nfun) = FreeSpaceNorm(grid, rg, nfun; constant=false)
norm_radial(grid, q::Grid.HankelTransform, nfun) = norm_radial(grid, Grid.RadialGrid(q), nfun)
norm_free(grid, xygrid::Grid.FreeGrid, nfun) = FreeSpaceNorm(grid, xygrid, nfun; constant=false)
norm_free2D(grid, xgrid::Grid.Free2DGrid, nfun) = FreeSpaceNorm(grid, xgrid, nfun; constant=false)

"""
    const_norm_radial(grid, q, nfun)
    const_norm_free(grid, xygrid, nfun)
    const_norm_free2D(grid, xgrid, nfun)

Make the normalisation factor ([`FreeSpaceNorm`](@ref)) for a `z`-independent refractive
index, computed once and reused. `nfun(λ)` takes wavelength; for crystal optics pass
`(nfunx, nfuny)` with `nfunx(λ, δθ)` and `nfuny(λ)`.
"""
const_norm_radial(grid, rg::Grid.RadialGrid, nfun) = FreeSpaceNorm(grid, rg, _zfun(nfun); constant=true)
const_norm_radial(grid, q::Grid.HankelTransform, nfun) = const_norm_radial(grid, Grid.RadialGrid(q), nfun)
const_norm_free(grid, xygrid::Grid.FreeGrid, nfun) = FreeSpaceNorm(grid, xygrid, _zfun(nfun); constant=true)
const_norm_free2D(grid, xgrid::Grid.Free2DGrid, nfun) = FreeSpaceNorm(grid, xgrid, _zfun(nfun); constant=true)

# wrap a z-independent index function in the (ω; z) / (λ, δθ; z), (λ; z) forms
_zfun(nfun) = (ω; z) -> nfun(wlfreq(ω))
function _zfun(nfuns::Tuple)
    nfunx, nfuny = nfuns
    ((λ, δθ; z) -> nfunx(λ, δθ), (λ; z) -> nfuny(λ))
end

"""
    TransFree

Transform E(ω) -> Pₙₗ(ω) for 3D free-space propagation.

# Fields
- `noise`: [`ModifiedNoise`](@ref) for the modified shot-noise model, or `nothing`.
- `Eωo_noise`/`Et_noise`: buffers holding the noise field at the current `z`, on the
  oversampled frequency/k-space and time/real-space grids `(nωo, np, nx, ny)`/
  `(nto, np, nx, ny)`, the shapes of the field's buffers.
- `Et_nl`: preallocated buffer for the combined field + noise, passed to `Et_to_Pt!`. The
  propagating field (`Eto`) is never modified.
"""
mutable struct TransFree{TT, FTT, IFTT, nT, rT, gT, xygT, dT, iT, nsT, eoT, eT, nlT}
    FT::FTT # 3D Fourier transform (space to k-space and time to frequency)
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    normfun::nT # Function which returns normalisation factor
    resp::rT # nonlinear responses (tuple of callables)
    grid::gT # time grid
    xygrid::xygT
    densityfun::dT # callable which returns density
    Pto::Array{TT, 4} # buffer for oversampled time-domain NL polarisation
    Eto::Array{TT, 4} # buffer for oversampled time-domain field
    Eωo::Array{ComplexF64, 4} # buffer for oversampled frequency-domain field
    Pωo::Array{ComplexF64, 4} # buffer for oversampled frequency-domain NL polarisation
    scale::Float64 # scale factor to be applied during oversampling
    idcs::iT # iterating over these slices Eto/Pto into Vectors, one at each position
    noise::nsT # ModifiedNoise for the modified shot-noise model, or nothing
    Eωo_noise::eoT # buffer for the oversampled spectral noise at the current z, or nothing
    Et_noise::eT # buffer for the time-domain noise at the current z, or nothing
    Et_nl::nlT # buffer for field+noise passed to Et_to_Pt!, or nothing
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
    TransFree(TT, scale, grid, xygrid, FT, responses, densityfun, normfun, pol=false; noise=nothing)

Construct a `TransFree` to calculate the reciprocal-domain nonlinear polarisation for 3D
free-space propagation.

# Keyword arguments
- `noise=nothing`: optional [`ModifiedNoise`](@ref) for the modified shot-noise model, whose
  draw is an `(nω, np, nx, ny)` frequency/k-space noise field, the same size as the `Eω`
  returned by [`Luna.setup`](@ref). At each step it is advanced to the current `z` and
  converted to the real-space oversampled time domain `(nto, np, nx, ny)` like the field.
  [`Fields.generate_noise_field`](@ref Luna.Fields.generate_noise_field) makes only
  `(nω[, nmodes])` draws, so a free-space draw is built by the caller.
"""
function TransFree(TT, scale, grid, xygrid, FT, responses, densityfun, normfun, pol=false;
                   noise=nothing)
    noise = mask_evanescent!(noise, normfun) # the vacuum noise of the propagating components
    Ny = length(xygrid.y)
    Nx = length(xygrid.x)
    Eωo = zeros(ComplexF64, (length(grid.ωo), pol ? 2 : 1, Nx, Ny))
    Eto = zeros(TT, (length(grid.to), pol ? 2 : 1, Nx, Ny))
    Pto = similar(Eto)
    Pωo = similar(Eωo)
    idcs = CartesianIndices((Nx, Ny))
    # buffers for the noise at the current z, shaped like the field's (ω/t, pol, x, y)
    if !isnothing(noise)
        Eωo_noise = similar(Eωo)
        Et_noise = similar(Eto)
        Et_nl = similar(Eto)
    else
        Eωo_noise = nothing
        Et_noise = nothing
        Et_nl = nothing
    end
    TransFree(FT, Utils.plan_ift(FT), normfun, responses, grid, xygrid, densityfun,
              Pto, Eto, Eωo, Pωo, scale, idcs, noise, Eωo_noise, Et_noise, Et_nl)
end

function TransFree(grid::Grid.RealGrid, args...; kwargs...)
    N = length(grid.ω)
    No = length(grid.ωo)
    scale = (No-1)/(N-1)
    TransFree(Float64, scale, grid, args...; kwargs...)
end

function TransFree(grid::Grid.EnvGrid, args...; kwargs...)
    N = length(grid.ω)
    No = length(grid.ωo)
    scale = No/N
    TransFree(ComplexF64, scale, grid, args...; kwargs...)
end

"""
    (t::TransFree)(nl, Eω, z)

Calculate the reciprocal-domain (ω-kx-ky-space) nonlinear response due to the field `Eω`
and place the result in `nl`.
"""
function (t::TransFree)(nl, Eωk, z)
    to_time!(t.Eto, Eωk, t.Eωo, t.IFT) # transform (ω, kx, ky) -> (t, x, y)
    # Modified shot-noise: advance the noise to this z ((ω,kx,ky) → (t,x,y)) and compute
    # field+noise in a separate buffer (Et_nl) so the propagating field (Eto) is never
    # contaminated.
    if !isnothing(t.noise)
        to_time!(t.Et_noise, noise_at!(t.noise, z), t.Eωo_noise, t.IFT)
        @. t.Et_nl = t.Eto + t.Et_noise
        Et_to_Pt!(t.Pto, t.Et_nl, t.resp, t.densityfun(z), t.idcs)
    else
        Et_to_Pt!(t.Pto, t.Eto, t.resp, t.densityfun(z), t.idcs)
    end
    @. t.Pto *= t.grid.towin # apodisation
    to_freq!(nl, t.Pωo, t.Pto, t.FT) # transform (t, x, y) -> (ω, kx, ky)
    nl .*= t.grid.ωwin .* (-im.*t.grid.ω)./(2 .* t.normfun(z))
end

mutable struct TransFree2D{TT, FTT, IFTT, nT, rT, gT, xgT, dT, iT}
    FT::FTT # 2D Fourier transform (space to k-space and time to frequency)
    IFT::IFTT # explicit inverse of FT (see Utils.plan_ift)
    normfun::nT # Function which returns normalisation factor
    resp::rT # nonlinear responses (tuple of callables)
    grid::gT # time grid
    xgrid::xgT
    densityfun::dT # callable which returns density
    Pto::Array{TT, 3} # buffer for oversampled time-domain NL polarisation
    Eto::Array{TT, 3} # buffer for oversampled time-domain field
    Eωo::Array{ComplexF64, 3} # buffer for oversampled frequency-domain field
    Pωo::Array{ComplexF64, 3} # buffer for oversampled frequency-domain NL polarisation
    scale::Float64 # scale factor to be applied during oversampling
    idcs::iT # iterating over these slices Eto/Pto into Vectors, one at each position
end

function show(io::IO, t::TransFree2D)
    grid = "grid type: $(typeof(t.grid))"
    samples = "time grid size: $(length(t.grid.t)) / $(length(t.grid.to))"
    resp = "responses: "*join([string(typeof(ri)) for ri in t.resp], "\n    ")
    x = "x grid: $(minimum(t.xgrid.x)) to $(maximum(t.xgrid.x)), N=$(length(t.xgrid.x))"
    out = join(["TransFree2D", grid, samples, x, resp], "\n  ")
    print(io, out)
end

function TransFree2D(TT, scale, grid, xgrid, FT, responses, densityfun, normfun, pol=false)
    Nx = length(xgrid.x)
    Eωo = zeros(ComplexF64, (length(grid.ωo), pol ? 2 : 1, Nx))
    Eto = zeros(TT, (length(grid.to), pol ? 2 : 1, Nx))
    Pto = similar(Eto)
    Pωo = similar(Eωo)
    idcs = CartesianIndices(size(Pto)[3:end])
    TransFree2D(FT, Utils.plan_ift(FT), normfun, responses, grid, xgrid, densityfun,
              Pto, Eto, Eωo, Pωo, scale, idcs)
end

"""
    TransFree2D(grid, xygrid, FT, responses, densityfun, normfun)

Construct a `TransFree2D` to calculate the reciprocal-domain nonlinear polarisation.

# Arguments
- `grid::AbstractGrid` : the grid used in the simulation
- `xgrid` : the spatial grid (instances of [`Grid.FreeGrid`](@ref))
- `FT::FFTW.Plan` : the 2D (t-x) Fourier transform for the oversampled time grid
- `responses` : `Tuple` of response functions
- `densityfun` : callable which returns the gas density as a function of `z`
- `normfun` : normalisation factor as fctn of `z`, can be created via [`norm_free`](@ref)
"""
function TransFree2D(grid::Grid.RealGrid, args...)
    N = length(grid.ω)
    No = length(grid.ωo)
    scale = (No-1)/(N-1)
    TransFree2D(Float64, scale, grid, args...)
end

function TransFree2D(grid::Grid.EnvGrid, args...)
    N = length(grid.ω)
    No = length(grid.ωo)
    scale = No/N
    TransFree2D(ComplexF64, scale, grid, args...)
end

"""
    (t::TransFree2D)(nl, Eω, z)

Calculate the reciprocal-domain (ω-kx-space) nonlinear response due to the field `Eω`
and place the result in `nl`.
"""
function (t::TransFree2D)(nl, Eωk, z)
    # TODO: this can probably be combined with the case for TransFree
    to_time!(t.Eto, Eωk, t.Eωo, t.IFT) # transform (ω, kx) -> (t, x)
    Et_to_Pt!(t.Pto, t.Eto, t.resp, t.densityfun(z), t.idcs) # add up responses
    @. t.Pto *= t.grid.towin # apodisation
    to_freq!(nl, t.Pωo, t.Pto, t.FT) # transform (t, x) -> (ω, kx)
    nl .*= t.grid.ωwin .* (-im.*t.grid.ω)./(2 .* t.normfun(z))
end

#= reflength! on a transform forwards to its normalisation; defined here, after every
   transform type exists. Transforms which are not free-space have nothing to taper. =#
reflength!(t::Union{TransRadial, TransFree, TransFree2D}, ℓ; kwargs...) = reflength!(t.normfun, ℓ; kwargs...)
reflength!(t, ℓ; kwargs...) = nothing

end
