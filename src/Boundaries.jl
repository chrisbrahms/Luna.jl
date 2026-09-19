"""
    Boundaries

Absorbing boundaries with **rate semantics**.

Luna's spectral and temporal absorbing boundaries used to be applied by multiplying the
solution by a fixed profile `W` after every accepted step. That makes the cumulative
absorption `W^N`, where `N` is the number of steps the adaptive controller happens to take,
so the answer depended on `rtol`. Tightening the tolerance does not converge on the taper
the profile describes: at the step counts of a typical run `W^N` is not a taper at all but a
brick wall at the edge of the flat region, and that brick wall is the limit it approaches.
It is why the old scheme eroded genuine spectral wings. An absorbing boundary must instead
have an absorption *coefficient per unit propagation distance* `α`, so that
traversing a distance `L` attenuates by a fixed factor however that path is subdivided.

`α` here is a **power** absorption coefficient, as everywhere else in Luna (`Modes.α`,
`LinearOps`): power falls as `exp(-αL)`, the field as `exp(-αL/2)`, and a linear operator
carries `-α/2`. This module turns the historical taper profiles into such coefficients,

    α(x) = -2 log(W(x)) / ℓ

where `ℓ` is a reference length: propagating through `ℓ` attenuates the *field* by exactly
the historical profile once, since the profiles were applied to the field. `ℓ` defaults to
`zmax/N`, i.e. "the historical window is applied `N` times over the whole propagation".

The two rates are then applied by different mechanisms, chosen by what the existing
machinery can do *exactly*:

- `α_ω` is diagonal in ω, so it is folded into the linear operator ([`addloss`](@ref)) and
  the interaction-picture propagator applies `exp(-α_ω Δz/2)` exactly, for whatever
  sub-interval the stepper uses, at no extra cost, and consistently with dense output.
- `α_t` is diagonal in `t`, so it cannot ride the propagator. It is applied as an
  exponential split-step factor `exp(-α_t Δz/2)` in `Luna.run`'s `stepfun`. Those factors
  telescope exactly to `exp(-α_t L/2)` regardless of the step layout, so the *total*
  absorption is still a rate. Interleaving the two parts this way would only be exact if
  they commuted, which they do not; the resulting error is first order in the step size and
  is confined to the collar, where the field is being destroyed anyway.

The *hard* band limit — the 0/1 part of the spectral window — is not applied here at all.
`NonlinearRHS` zeroes the nonlinear polarisation outside `grid.sidx` and the linear operator
is zero there, so once the field has been band-limited at the start of `Luna.run` nothing
can put anything back.

# Free space

A free-space simulation has two more finite axes, and one more decaying process, all of
which end up in the same place: a negative real part of the linear operator over the
transverse wavevector `k⊥`.

- **Evanescent decay.** Components with `k⊥ > k(ω)` do not propagate; the exact one-way
  solution is `exp(-κz)` with `κ = sqrt(k⊥² - k²)`, which `LinearOps` writes into the
  operator uncapped (see [`LinearOps.βz`](@ref)). At the low-frequency end of a grid (where
  `k → 0`) `κ` reaches `k⊥,max`, i.e. `1e5–1e6 /m` for micron transverse sampling. That is
  physics, not a boundary, so it is applied in every `boundary` mode. What it needs from
  this module is the same treatment as the spectral absorber: the decay rate of the
  operator is clamped at `MAX_αℓ/(2ℓ)` ([`clampdecay`](@ref)) and the nonlinear source in
  those channels is tapered by the matching factor `exp(-κℓ)`
  ([`NonlinearRHS.FreeSpaceNorm`](@ref), set up through [`NonlinearRHS.reflength!`](@ref)),
  so that the interaction picture never amplifies by more than `exp(MAX_αℓ/2)`. Components
  with `κℓ ≫ 1` thereby lose their source, which is the correct limit — their driven
  amplitude `S/κ` vanishes — and components with `κℓ ≲ 1` are integrated exactly. In
  `:rate` mode `ℓ` is the absorber reference length; in `:none` and `:legacy` it is
  `max_dz`, since the only requirement is `Δz ≤ ℓ`. Before this the operator capped `κ` at
  an arbitrary `200 /m` and the source prefactor below cutoff was an arbitrary `1.0`, so
  free-space `:legacy` results are *not* bit-identical to historical ones.
- **k-space window** (`:rate` only). Where `k⊥,max < k(ω)` — the whole optical band on a
  typical grid — nothing attenuates content piling up at the edge of the transverse k
  grid, exactly as at the edge of the frequency window. A Planck taper over the outer
  `kcollar` fraction of `|k⊥|` ([`kprofile`](@ref)) becomes a rate ([`rate`](@ref)), rides
  the propagator through [`addloss_k`](@ref) and tapers the source through the same
  `FreeSpaceNorm`. Where the evanescent region reaches the grid edge this is redundant
  with `κ`, and harmless.
- **Transverse collar** (`:rate` only). The analogue of the temporal collar: a beam
  reaching the transverse edge wraps around (FFT grids) or reflects off the `E(R) = 0`
  wall of the QDHT. A rate over the outer part of the transverse grid ([`rprofile`](@ref))
  is applied as a split-step factor per accepted step ([`RadialCollar`](@ref),
  [`CartesianCollar`](@ref)), with the same first-order, collar-confined
  non-commutation error as the temporal collar. As for the temporal collar, what it removes
  is measured and reported once if it becomes noticeable.

# Step size

`Luna.run` caps `max_dz` (and `init_dz`) at the reference length `ℓ`. This is not about
accuracy, it is about not confusing the step-size controller.

RK45 works in the interaction picture, and `fbar!` back-propagates the right-hand side with
`exp(-linop Δz)` — which, for a linop that decays, means it *amplifies* by `exp(+α Δz/2)`.
The elements it amplifies most are at the band edge, where the absorber is deepest. Those
elements are not independent of the rest, because `RK45.weaknorm` sums the error estimate
over every element and divides by the global norm of the solution, so an amplified band-edge
element can come to dominate the error estimate and drive the controller into repeated
rejections.

What saves it is that the nonlinear drive at the band edge is itself proportional to
`ωwin = exp(-α ℓ/2)`. The amplified drive therefore scales as `exp(α(Δz - ℓ)/2)`, which
stays below 1 for any `α` as long as `Δz ≤ ℓ`. Hence the cap, which costs nothing in
practice: it requires at least `boundary_N` steps, and real runs take far more.

See also `Luna.run`'s `boundary` keyword argument.
"""
module Boundaries

import ..Maths
import ..Grid
import ..NonlinearRHS
import Hankel
import LinearAlgebra: mul!, ldiv!
import Logging
import Printf: @sprintf

"Default number of applications of the historical window profile over the propagation."
const DEFAULT_N = 20

"Default minimum temporal collar width, as a fraction of the full time window."
const DEFAULT_TCOLLAR = 0.05

"Default width of the k-space absorber collar, as a fraction of the largest `|k⊥|`."
const DEFAULT_KCOLLAR = 0.1

"Default width of the transverse absorber collar of a radial grid, as a fraction of `R`."
const DEFAULT_RCOLLAR = 0.1

"Fraction of the pulse the temporal absorber may remove before `Luna.run` warns about it."
const DEFAULT_WARNFRAC = 1e-3

"""
Cap on `α*ℓ`, the depth of the absorber over one reference length.

Two jobs. It keeps `α` finite at all — the tapers reach exactly zero at the band edge, and
`Inf*0` is `NaN` — and, given the step-size cap described in the module docstring, it keeps
`exp(α Δz/2)` finite in the interaction picture. At 60 the field is attenuated by at most
`exp(-30) ≈ 1e-13` over one reference length, far below anything physically meaningful.

The same number caps the evanescent decay of a free-space operator ([`clampdecay`](@ref)):
a field decay rate of `MAX_αℓ/(2ℓ)`, so that the clamped components still fall by
`exp(-30)` per reference length.
"""
const MAX_αℓ = 60.0

"""
    reflength(zmax, N, ℓ)

The absorber reference length in metres. `ℓ` wins if given; otherwise `zmax/N`, `zmax`
being the propagation length passed to [`Luna.run`](@ref).

Both are validated here so that every call site is protected. A non-positive length is not
merely useless, it is silently destructive: `ℓ == 0` makes [`rate`](@ref) evaluate `0/0` in
the interior and fills the grid with NaN, and `ℓ < 0` inverts the clamp bounds so that
`clamp` returns the (negative) upper bound everywhere — an absorber which amplifies, over
the whole grid rather than just the collar. To switch the absorbers off, use
`boundary=:none`.
"""
function reflength(zmax, N, ℓ)
    if isnothing(ℓ)
        (N > 0 && isfinite(N)) || error(
            "boundary_N must be finite and positive, got $N")
        zmax/N
    else
        (ℓ > 0 && isfinite(ℓ)) || error(
            "boundary_length must be finite and positive, got $ℓ")
        float(ℓ)
    end
end

"""
    rate(W, ℓ)

Convert an apodisation profile `W ∈ [0, 1]` (1 in the interior, falling to 0 at the
boundary) into a power absorption coefficient in 1/m, such that propagating through `ℓ`
attenuates the *field* by `W` — that is, `exp(-αℓ/2) == W`, since the historical profiles
multiplied the field. Exactly zero where `W == 1`, and clamped at `MAX_αℓ/ℓ` so that
`W == 0` gives a large but finite coefficient.
"""
rate(W, ℓ) = clamp.(.-2 .* log.(W) ./ ℓ, 0.0, MAX_αℓ/ℓ)

"""
    spectral_rate(grid, zmax; N=DEFAULT_N, ℓ=nothing)

Power absorption coefficient in 1/m across the frequency grid, derived from
`grid.ωwin`. Exactly zero outside `grid.sidx`, where the linear operator is zero by
construction and there is nothing left to absorb. `zmax` is the propagation length, used
only to derive the reference length from `N` (see [`reflength`](@ref)).
"""
function spectral_rate(grid, zmax; N=DEFAULT_N, ℓ=nothing)
    α = rate(grid.ωwin, reflength(zmax, N, ℓ))
    α[.!grid.sidx] .= 0
    α
end

"""
    tcollarwidth(grid, collar)

Width (in seconds) of the temporal absorber collar at each end of the time window.

`grid.twin`'s own collar is only the *slack* between the requested `trange` and the
realised power-of-two window, so it can be exactly zero — for a `trange` that already lands
on a power of two, `grid.twin` is identically 1 except for a single sample, i.e. those
grids have no temporal absorber at all. The collar used here is therefore the wider of that
slack and `collar` times the full time window.
"""
function tcollarwidth(grid, collar)
    t = grid.t
    tmin, tmax = extrema(t)
    i1 = findfirst(isequal(1.0), grid.twin)
    i2 = findlast(isequal(1.0), grid.twin)
    slack = (isnothing(i1) || isnothing(i2)) ? 0.0 : min(t[i1] - tmin, tmax - t[i2])
    max(slack, collar*(tmax - tmin))
end

"""
    tprofile(grid; collar=DEFAULT_TCOLLAR)

The temporal absorber profile: the pointwise minimum of `grid.twin` and a Planck taper with
a guaranteed non-degenerate collar (see [`tcollarwidth`](@ref)). Taking the minimum means
the absorber is never weaker than the historical profile, and equals it exactly whenever
the natural collar is already wide enough.

Deliberately not stored on the grid: `Grid.to_dict`/`from_dict` serialise every grid field
and `from_dict` requires them all, so adding one would make existing output files
unreadable.
"""
function tprofile(grid; collar=DEFAULT_TCOLLAR)
    t = grid.t
    tmin, tmax = extrema(t)
    w = tcollarwidth(grid, collar)
    min.(grid.twin, Maths.planck_taper(t, tmin, tmin + w, tmax - w, tmax))
end

"""
    temporal_rate(grid, zmax; N=DEFAULT_N, ℓ=nothing, collar=DEFAULT_TCOLLAR)

Power absorption coefficient in 1/m across the time grid. `zmax` is the propagation
length, used only to derive the reference length from `N` (see [`reflength`](@ref)).
"""
temporal_rate(grid, zmax; N=DEFAULT_N, ℓ=nothing, collar=DEFAULT_TCOLLAR) =
    rate(tprofile(grid; collar), reflength(zmax, N, ℓ))

"""
    addloss(linop, α)

Add the power absorption coefficient `α` (a vector over ω) to a linear operator as
`-α/2`, which is how loss enters a linop everywhere in Luna (see `LinearOps`). Returns a
linop of the same kind, so that `RK45.make_prop!` dispatch and `Luna.linoptype` are
unchanged. Handles both forms Luna uses: a materialised `AbstractArray` and a closure
`linop!(out, z)`.

ω is axis 1 for every linop shape — `(Nω,)` mode-averaged, `(Nω, nmodes)` multimode,
`(Nω, q.N)` radial, `(Nω, Nky, Nkx)` free-space — so one broadcast covers all of them.

The array method **allocates a copy** rather than subtracting in place:
`Interface.prop_capillary_args` is documented as being for repeated simulations in an
identical fibre, i.e. the same linop is deliberately reused across `Luna.run` calls, and an
in-place subtraction would compound the absorber on every reuse.
"""
addloss(linop::AbstractArray, α) = linop .- α./2

function addloss(linop!, α)
    function linop_absorbing!(out, z)
        linop!(out, z)
        out .-= α./2
        out
    end
end

# ---------------------------------------------------------------------------- free space

"""
    spacegrid(transform)

The transverse grid of a free-space transform (`Hankel.QDHT`, `Grid.FreeGrid` or
`Grid.Free2DGrid`), or `nothing` for a modal transform. The free-space parts of the
boundaries are set up only when this is not `nothing`.
"""
spacegrid(t::NonlinearRHS.TransRadial) = t.QDHT
spacegrid(t::NonlinearRHS.TransFree) = t.xygrid
spacegrid(t::NonlinearRHS.TransFree2D) = t.xgrid
spacegrid(t) = nothing

# Planck taper over |k|: 1 up to (1-collar)*K and 0 at K, the largest |k| on the axis.
function ktaper(k, collar)
    K = maximum(abs, k)
    Maths.planck_taper(abs.(k), -K, -(1 - collar)*K, (1 - collar)*K, K)
end

"""
    kprofile(spacegrid, collar)

The k-space absorber profile over the transverse k axes of `spacegrid`: 1 in the interior,
falling to 0 over the outer `collar` fraction of the largest wavevector the grid holds
(`q.K` for a `QDHT`, the Nyquist wavevector for FFT grids). For a `FreeGrid` it is the
product of the tapers along `kx` and `ky`, matching the square edge of the grid. Sized
like the k axes, `(Nk,)` or `(Nkx, Nky)`, so that it broadcasts over a linop's trailing
axes ([`addloss_k`](@ref)) and matches [`NonlinearRHS.FreeSpaceNorm`](@ref)'s `kwin`.
`collar = 0` gives a profile of ones, i.e. no k-space absorber.
"""
kprofile(q::Hankel.QDHT, collar) = Maths.planck_taper(
    q.k, -q.K, -(1 - collar)*q.K, (1 - collar)*q.K, q.K)
kprofile(sg::Grid.Free2DGrid, collar) = ktaper(sg.kx, collar)
kprofile(sg::Grid.FreeGrid, collar) = ktaper(sg.kx, collar) .* ktaper(sg.ky, collar)'

"""
    rprofile(spacegrid, collar)

The transverse absorber profile in real space. For a `QDHT` it is a Planck taper from
`(1 - collar)R` to the aperture `R`, where the transform imposes `E(R) = 0`. The Cartesian
grids already carry the window they were built with (`window_factor` in
[`Grid.FreeGrid`](@ref) and [`Grid.Free2DGrid`](@ref) extends the box to make room for
it), so for those `collar` is ignored and that window is returned.
"""
rprofile(q::Hankel.QDHT, collar) = Maths.planck_taper(
    q.r, -q.R, -(1 - collar)*q.R, (1 - collar)*q.R, q.R)
rprofile(sg::Grid.Free2DGrid, collar) = copy(sg.xwin)
rprofile(sg::Grid.FreeGrid, collar) = dropdims(sg.xywin; dims=1)

"""
    addloss_k(linop, α)

Like [`addloss`](@ref), for a power absorption coefficient `α` defined over the transverse
k axes of a free-space operator: the operator has shape `(Nω, Npol, Nk...)` and `α` has
shape `Nk...`, so it is broadcast along ω and polarisation.
"""
addloss_k(linop::AbstractArray, α) = linop .- reshape(α, 1, 1, size(α)...)./2

function addloss_k(linop!, α)
    αr = reshape(α, 1, 1, size(α)...)
    function linop_kabsorbing!(out, z)
        linop!(out, z)
        out .-= αr./2
        out
    end
end

"""
    clampdecay(linop, ratemax)

Limit the decay of a free-space linear operator — its negative real part, which is the
evanescent `-κ` ([`LinearOps.βz`](@ref)) — to the *field* rate `ratemax` in 1/m, keeping the
imaginary part. Returns a new operator of the same kind as [`addloss`](@ref) does.

`κ` is exact and unbounded in `LinearOps`, and the interaction picture amplifies the
nonlinear source in a channel by `exp(κ Δz)`, which overflows for the `κ` of any fine
transverse grid. Clamping alone would inflate the driven amplitude of those channels from
`S/κ` to `S/ratemax`; it is only correct together with the matching taper of the source
([`NonlinearRHS.reflength!`](@ref)), which removes the source from exactly the channels the
clamp touches. See "Free space" in the module docstring.
"""
clampdecay(linop::AbstractArray, ratemax) = @. complex(max(real(linop), -ratemax), imag(linop))

function clampdecay(linop!, ratemax)
    function linop_clamped!(out, z)
        linop!(out, z)
        @. out = complex(max(real(out), -ratemax), imag(out))
        out
    end
end

"""
    evanescent(linop, transform, ℓ; kwin=nothing)

Make the evanescent channels of a free-space operator safe for the interaction-picture
stepper taking steps up to `ℓ`: clamp the decay of `linop` at `MAX_αℓ/(2ℓ)`
([`clampdecay`](@ref)) and taper the nonlinear source of `transform` over `ℓ` to match
([`NonlinearRHS.reflength!`](@ref)), including the k-space window profile `kwin` if given.
Returns the new operator. A modal `transform` has no evanescent channels and `linop` is
returned unchanged.
"""
function evanescent(linop, transform, ℓ; kwin=nothing)
    isnothing(spacegrid(transform)) && return linop
    ratemax = MAX_αℓ/(2ℓ)
    NonlinearRHS.reflength!(transform, ℓ; κmax=ratemax, kwin)
    clampdecay(linop, ratemax)
end

"""
    RadialCollar(q, αr, Eω)

Transverse absorbing boundary for radially symmetric propagation: the power rate `αr` over
`q.r`, applied as `exp(-αr Δz/2)` per accepted step. The QDHT imposes `E(R) = 0`, a hard
wall which reflects whatever reaches the aperture, and the collar absorbs it first. It is
applied to `Eω` directly (the collar is diagonal in ω) with one inverse and one forward
Hankel transform along the last axis, into the buffer `buf` sized like `Eω`.
"""
struct RadialCollar{qT, bT}
    q::qT
    αr::Vector{Float64}
    ridcs::Vector{Int} # only the collar is ever ≠ 1
    weight::Vector{Float64} # radial integration weights, to measure what is removed
    buf::bT
    removed::Base.RefValue{Float64}
    reference::Base.RefValue{Float64}
    warned::Base.RefValue{Bool}
end

RadialCollar(q, αr, Eω) = RadialCollar(q, αr, findall(>(0), αr), copy(q.scaleR), similar(Eω),
                                       Ref(0.0), Ref(0.0), Ref(false))

# applied before the temporal collar, in (ω, k⊥) space
function apply_kspace!(c::RadialCollar, Eω, Δz)
    ldiv!(c.buf, c.q, Eω) # (ω, pol, k) -> (ω, pol, r)
    d = ndims(c.buf)
    if c.reference[] == 0
        c.reference[] = sum(i -> c.weight[i]*sum(abs2, selectdim(c.buf, d, i)), axes(c.buf, d))
    end
    removed = 0.0
    for i in c.ridcs
        fac = exp(-c.αr[i]*Δz/2)
        s = selectdim(c.buf, d, i)
        removed += c.weight[i]*sum(abs2, s)*(1 - fac^2)
        s .*= fac
    end
    c.removed[] += removed
    mul!(Eω, c.q, c.buf)
    nothing
end
apply_kspace!(c, Eω, Δz) = nothing

"""
    CartesianCollar(αxy)

Transverse absorbing boundary for the Cartesian free-space grids: the power rate `αxy` over
the spatial axes `(Nx,)` or `(Nx, Ny)`, applied as `exp(-αxy Δz/2)` per accepted step. The
Fourier transform of those grids is joint in `(t, x[, y])`, so the collar is applied in
the same real-space pass as the temporal collar, at no extra transform cost.
"""
struct CartesianCollar{N}
    αxy::Array{Float64, N}
    idcs::Vector{CartesianIndex{N}} # only the collar is ever ≠ 1
    removed::Base.RefValue{Float64}
    reference::Base.RefValue{Float64}
    warned::Base.RefValue{Bool}
end

CartesianCollar(αxy) = CartesianCollar(
    αxy, vec(collect(CartesianIndices(αxy)))[vec(αxy .> 0)], Ref(0.0), Ref(0.0), Ref(false))

# applied after the temporal collar, in (t, x[, y]) space
function apply_realspace!(c::CartesianCollar, Et, Δz)
    c.reference[] == 0 && (c.reference[] = sum(abs2, Et))
    removed = 0.0
    for J in c.idcs
        fac = exp(-c.αxy[J]*Δz/2)
        s = view(Et, :, :, J)
        removed += sum(abs2, s)*(1 - fac^2)
        s .*= fac
    end
    c.removed[] += removed
    nothing
end
apply_realspace!(c, Et, Δz) = nothing

"""
    spatialcollar(spacegrid, αr, grid, Et)

The transverse absorber functor for `spacegrid`, given the power rate `αr` over its real
space. `grid` and `Et` size the buffer the radial collar needs.
"""
spatialcollar(q::Hankel.QDHT, αr, grid, Et) = RadialCollar(
    q, αr, zeros(ComplexF64, (length(grid.ω), size(Et)[2:end]...)))
spatialcollar(sg::Union{Grid.Free2DGrid, Grid.FreeGrid}, αr, grid, Et) = CartesianCollar(αr)

# --------------------------------------------------------------------------- application

#= The three ways a boundary can be applied per accepted step. Each is a functor rather than
   a closure so that `setup` has one return type per mode, the state each needs is named,
   and a new kind of boundary (the spatial windows a free-space simulation wants) is another
   struct rather than another branch inside `Luna.run`. Each ends by calling `output`,
   because that is what `RK45.solve` expects of a `stepfun`. =#

"""
    RateAbsorber(αt, Et, FT, output, z0; spatial=nothing)

Applies the temporal absorber as `exp(-αt Δz/2)` over the distance actually travelled.
Successive factors telescope to `exp(-αt L/2)` however the solver subdivides the
propagation, which is the whole point of rate semantics. The spectral absorber is not here:
it rides the propagator, having been folded into the linear operator by [`addloss`](@ref).

In free space `spatial` is the transverse collar ([`RadialCollar`](@ref) or
[`CartesianCollar`](@ref)), applied in the same step; the k-space absorber rides the
propagator like the spectral one.
"""
struct RateAbsorber{tT, fT, oT, sT}
    αt::Vector{Float64}
    tidcs::Vector{Int} # only the collar is ever ≠ 1, so only those need touching each step
    tfac::Vector{Float64}
    Et::tT
    FT::fT
    output::oT
    spatial::sT
    zprev::Base.RefValue{Float64}
    removed::Base.RefValue{Float64} # running total of |E|² taken out of the collar
    reference::Base.RefValue{Float64} # |E|² over the whole window, at the first step
    warned::Base.RefValue{Bool}
    warnfrac::Float64
end

RateAbsorber(αt, Et, FT, output, z0; warnfrac=DEFAULT_WARNFRAC, spatial=nothing) = RateAbsorber(
    αt, findall(>(0), αt), ones(Float64, length(αt)), Et, FT, output, spatial, Ref(float(z0)),
    Ref(0.0), Ref(0.0), Ref(false), warnfrac)

function (b::RateAbsorber)(Eω, z, dz, interpolant)
    Δz = z - b.zprev[]
    b.zprev[] = z
    if Δz > 0
        apply_kspace!(b.spatial, Eω, Δz) # radial collar: transforms to r and back itself
        @inbounds for i in b.tidcs
            b.tfac[i] = exp(-b.αt[i]*Δz/2) # αt is a power coefficient, tfac hits the field
        end
        #= An inverse real FFT overwrites its input, so between these two lines Eω holds
           whatever FFTW left there and must not be read; the forward transform refills it. =#
        ldiv!(b.Et, b.FT, Eω)
        b.reference[] == 0 && (b.reference[] = sum(abs2, b.Et))
        #= Apply the collar and measure what it took out, in one pass. Only the collar can
           change, so this also does less work than multiplying the whole array. The trailing
           index covers every field shape: none for mode-averaged, modes, or transverse. =#
        removed = 0.0
        @inbounds for J in CartesianIndices(size(b.Et)[2:end]), i in b.tidcs
            before = abs2(b.Et[i, J])
            b.Et[i, J] *= b.tfac[i]
            removed += before - abs2(b.Et[i, J])
        end
        b.removed[] += removed
        apply_realspace!(b.spatial, b.Et, Δz) # Cartesian collar: Et is already (t, x[, y])
        mul!(Eω, b.FT, b.Et)
        warn_maybe(b, z)
        warn_maybe(b.spatial, z, b.warnfrac)
    end
    b.output(Eω, z, dz, interpolant)
end

"""
    warn_maybe(b::RateAbsorber, z)

Warn, once, if the temporal absorber has eaten a noticeable fraction of the pulse.

This is a measurement of what the boundary actually removed, not a prediction from the
dispersion. A prediction cannot work well here: at setup time the only thing available is
the group delay across the simulation band, and on a broadband grid that is dominated by
the band edges — where a capillary has an enormous group delay and nothing but shot noise
to carry it. Such a check fires on almost every realistic run, and a warning which is
usually wrong is one people learn to ignore.

What the user needs to know is whether light they care about is leaving the time window,
and the absorber is the thing that finds out. Measuring costs nothing extra, since the
collar is being multiplied anyway, and it fires when it becomes true rather than
speculating beforehand.
"""
function warn_maybe(b::RateAbsorber, z)
    (b.warned[] || b.reference[] == 0) && return nothing
    frac = b.removed[]/b.reference[]
    frac > b.warnfrac || return nothing
    b.warned[] = true
    Logging.@warn(@sprintf(
        "Temporal absorbing boundary has removed %.2g%% of the pulse by z = %.3g m. Light \
         is reaching the edge of the time window and being absorbed; if that is not \
         intended, widen `trange`. (Reported once.)", 100frac, z))
    nothing
end

# the same, for the transverse collar of a free-space simulation
function warn_maybe(c::Union{RadialCollar, CartesianCollar}, z, warnfrac)
    (c.warned[] || c.reference[] == 0) && return nothing
    frac = c.removed[]/c.reference[]
    frac > warnfrac || return nothing
    c.warned[] = true
    Logging.@warn(@sprintf(
        "Transverse absorbing boundary has removed %.2g%% of the beam by z = %.3g m. Light \
         is reaching the edge of the transverse grid and being absorbed; if that is not \
         intended, enlarge the grid. (Reported once.)", 100frac, z))
    nothing
end
warn_maybe(c::Nothing, z, warnfrac) = nothing

"""
    LegacyAbsorber(grid, Et, FT, output)

The historical scheme: multiply the solution by the fixed profiles once per accepted step.
Kept only so that results from before rate semantics can be reproduced exactly.
"""
struct LegacyAbsorber{gT, tT, fT, oT}
    grid::gT
    Et::tT
    FT::fT
    output::oT
end

function (b::LegacyAbsorber)(Eω, z, dz, interpolant)
    Eω .*= b.grid.ωwin
    ldiv!(b.Et, b.FT, Eω) # destroys Eω, see RateAbsorber
    b.Et .*= b.grid.twin
    mul!(Eω, b.FT, b.Et)
    b.output(Eω, z, dz, interpolant)
end

"""
    NoAbsorber(output)

No boundary at all. The field is still band-limited once at the start of `Luna.run` and the
nonlinear polarisation is still band-limited in `NonlinearRHS`, but nothing stops energy
wrapping around the time window or piling up at the edge of the frequency window.
"""
struct NoAbsorber{oT}
    output::oT
end

(b::NoAbsorber)(Eω, z, dz, interpolant) = b.output(Eω, z, dz, interpolant)

"Report the absorbing-boundary configuration."
function log_setup(grid, zmax, ℓ, collar, sg=nothing, kcollar=0, rcollar=0)
    w = tcollarwidth(grid, collar)
    trange = maximum(grid.t) - minimum(grid.t)
    #= The clamp in `rate` is not worth reporting: it bites only where the profile is
       already below exp(-30), and even there the attenuation over the whole propagation is
       exp(-30 zmax/ℓ). What is worth reporting is the strength the user actually chose, so
       quote the attenuation over the propagation at the half-way point of a taper. =#
    Logging.@info(@sprintf(
        "Absorbing boundaries: rate-based, reference length %.3g m (%.3g applications of \
         the window profile over %.3g m; a 50%% point of the taper attenuates by %.1e over \
         the propagation). Temporal collar %.3g fs, %.1f%% of the time window.",
        ℓ, zmax/ℓ, zmax, 0.5^(zmax/ℓ), w*1e15, 100*w/trange))
    isnothing(sg) && return nothing
    rdesc = sg isa Hankel.QDHT ? @sprintf("%.1f%% of the aperture", 100rcollar) :
                                 "the grid's own window"
    Logging.@info(@sprintf(
        "Free-space boundaries: k-space collar %.1f%% of the largest k⊥, transverse \
         collar %s. Evanescent decay clamped at %.3g /m and its source tapered over %.3g m.",
        100kcollar, rdesc, MAX_αℓ/(2ℓ), ℓ))
end

"""
    setup(boundary, grid, transform, linop, Et, FT, output, z0, zmax, max_dz, init_dz;
          N=DEFAULT_N, ℓ=nothing, collar=DEFAULT_TCOLLAR,
          kcollar=DEFAULT_KCOLLAR, rcollar=DEFAULT_RCOLLAR, warnfrac=DEFAULT_WARNFRAC)

Everything `Luna.run` needs in order to apply absorbing boundaries, as a named tuple
`(; stepfun, linop, max_dz, init_dz, ℓ)`. `ℓ` is the reference length actually used, or
`nothing` for the modes which do not have one.

Three of those are returned because setting up an absorber genuinely changes them, and it is
clearer to hand them back than to mutate them from inside a branch:

- `linop` gains the spectral absorber, which the interaction-picture propagator then applies
  exactly (see [`addloss`](@ref)). In free space it also gains the k-space absorber
  ([`addloss_k`](@ref)) and, in *every* mode, the clamp of its evanescent decay
  ([`evanescent`](@ref)), the nonlinear source of `transform` being tapered to match.
- `max_dz` is capped at the reference length `ℓ`, and `init_dz` with it, so that the
  step-size controller is not thrown by the amplified band-edge elements the interaction
  picture produces. See "Step size" in the module docstring.

`transform` is used only to find the transverse grid ([`spacegrid`](@ref)) and to taper its
normalisation; for a modal transform it is untouched. `kcollar` and `rcollar` are the
k-space and transverse collar widths ([`kprofile`](@ref), [`rprofile`](@ref)). `zmax` is
the propagation length, which sets the reference length when `N` is used rather than `ℓ`.
"""
function setup(boundary, grid, transform, linop, Et, FT, output, z0, zmax, max_dz, init_dz;
               N=DEFAULT_N, ℓ=nothing, collar=DEFAULT_TCOLLAR,
               kcollar=DEFAULT_KCOLLAR, rcollar=DEFAULT_RCOLLAR, warnfrac=DEFAULT_WARNFRAC)
    ℓabs = nothing
    sg = spacegrid(transform)
    if boundary === :rate
        ℓabs = reflength(zmax, N, ℓ)
        if max_dz > ℓabs
            Logging.@info(@sprintf(
                "Reducing max_dz from %.3g m to the absorber reference length %.3g m.",
                max_dz, ℓabs))
            max_dz = ℓabs
        end
        init_dz = min(init_dz, max_dz)
        αt = temporal_rate(grid, zmax; N, ℓ, collar)
        spatial = nothing
        if isnothing(sg)
            linop = addloss(linop, spectral_rate(grid, zmax; N, ℓ))
        else
            #= Order matters: the clamp must see only the physical decay, not the absorbers
               added after it. The source taper carries the k-window at its clamped depth,
               never zero, so the division in FreeSpaceNorm stays finite. =#
            Wk = kprofile(sg, kcollar)
            linop = evanescent(linop, transform, ℓabs; kwin=max.(Wk, exp(-MAX_αℓ/2)))
            linop = addloss(linop, spectral_rate(grid, zmax; N, ℓ))
            linop = addloss_k(linop, rate(Wk, ℓabs))
            spatial = spatialcollar(sg, rate(rprofile(sg, rcollar), ℓabs), grid, Et)
        end
        log_setup(grid, zmax, ℓabs, collar, sg, kcollar, rcollar)
        stepfun = RateAbsorber(αt, Et, FT, output, z0; warnfrac, spatial)
    elseif boundary === :legacy || boundary === :none
        if boundary === :legacy
            Logging.@warn(
                "boundary=:legacy applies the absorbing boundaries once per accepted step, " *
                "so the absorption depends on the step count and the result depends on rtol.")
            stepfun = LegacyAbsorber(grid, Et, FT, output)
        else
            stepfun = NoAbsorber(output)
        end
        #= No absorbers, but the evanescent channels still need their clamp and taper, and
           the only requirement on the reference length is that no step exceeds it. =#
        linop = evanescent(linop, transform, min(max_dz, zmax))
    else
        error("boundary must be :rate, :legacy or :none, not $boundary")
    end
    (; stepfun, linop, max_dz, init_dz, ℓ=ℓabs)
end

end
