module Stats
import Luna
import Luna: Maths, Grid, Modes, Utils, settings, PhysData, Fields, Processing, Ionisation,
             Boundaries
import Luna.PhysData: wlfreq, c, ε_0
import Luna.NonlinearRHS: TransModal, TransModalFixed, TransModeAvg, TransRadial,
                          TransFree, TransFree2D, Erω_to_Prω!,
                          integral_error!, has_error_estimate
import Luna.Nonlinear: PlasmaCumtrapz
import Luna.Capillary: MarcatiliMode
import FFTW
import LinearAlgebra: mul!
import Printf: @sprintf
import Logging
import Logging: @warn
import Base.Broadcast
import Luna.RK45: _zipreduce

#=================================================#
#==========  THE STATISTICS PROTOCOL  ============#
#=================================================#

#= Every statistic is a callable with the signature `(d, Eω, Et, z, dz)` which mutates the
   dictionary `d`. The ones Luna ships are callable structs rather than closures so that
   they can carry the three traits below; a user statistic is usually a closure, for which
   the fallbacks apply.

   Each struct has two branches: the host branch, which is the code Luna has always run and
   is kept arithmetically identical to it, and the device branch, one or more reductions
   and broadcasts over the state as the stepper holds it. Which one runs is decided by
   `Utils.isdevice`, not by precision, so a `Float32` run on the host is also unchanged.

   On a device the state is *scaled* (`e = E/E_ref`, `Luna.UnitScaling`), because nothing
   has unscaled it: `Luna.ScaledOutput` only unscales into a host buffer, which is the copy
   these statistics exist to avoid. The device branches therefore apply `E_ref` themselves,
   always to the scalar result of a reduction and in `Float64`, never to the field. See
   `StatsContext`. =#

"""
    device_capable(f) -> Bool

Whether the statistics function `f` can be evaluated on the propagating state as the
stepper holds it: a device array, in the scaled units of [`Luna.UnitScaling`](@ref).

`false` for anything Luna does not know about, in particular a user-written closure, which
is why a statistics set containing one makes `Luna.ScaledOutput` copy the state to the host
on every step whose statistics fire.

The trait is structural: it does not depend on the array type a statistic was prepared
for, only on whether its algorithm has a device form. `fwhm_r` and
`mode_reconstruction_error` are the two statistics in the default sets which do not.

Which path a statistic actually *takes* is a separate, fixed decision -- see
[`collect_stats`](@ref) and `Stats._onstate`. For a whole set,
[`device_capable`](@ref)`(::StatsCollector)` reports the path, because that is what
`Luna.ScaledOutput` has to know.
"""
device_capable(f) = false

"""
    needs_time(f) -> Bool

Whether the statistics function `f` reads the time-domain field `Et`.
[`collect_stats`](@ref) does the inverse transform which produces `Et` only when at least
one of its statistics does. `true` for anything Luna does not know about, so a user
statistic always gets a valid `Et`.
"""
needs_time(f) = true

"""
    statlabel(f) -> String

A short name for the statistics function `f`, used in the one-time warning
`Luna.ScaledOutput` emits when a statistic forces a device-to-host copy. The type's name
for a callable struct; for an anonymous closure this is the gensym Julia gave it, which is
the best available.
"""
statlabel(f) = string(nameof(typeof(f)))

"""
    prepare(f, ctx::StatsContext)

The statistics function `f` with its buffers allocated for the run described by `ctx` (the
state's array type and precision, and the unit scaling it is expressed in). Called by
[`collect_stats`](@ref) on every statistic it is given; the fallback returns `f`
unchanged, which is what a user closure needs.
"""
prepare(f, ctx) = f

#= Which branch a statistic takes is decided once, at `prepare`, and carried in its
   `ondevice` field -- not read off the array it is handed. The two can disagree (review
   round 1, finding 1: an `HDF5Output` with a resume cache used to hand a host copy to a
   set built for the device, which threw on Metal and silently applied `E_ref^2` twice on
   JLArrays), and a statistic which reads the array type cannot tell a legitimate call
   from a mistake. `Utils.isdevice` is a type-level trait, so the check costs nothing once
   the method is specialised. =#
@inline function _onstate(f, x)
    f.ondevice === Utils.isdevice(x) || _branchmismatch(f, x)
    f.ondevice
end

@noinline _branchmismatch(f, x) = error(
    "the statistic $(statlabel(f)) was prepared for a $(f.ondevice ? "device" : "host") "*
    "state but was called with a $(typeof(x)). The array a statistics set is called with "*
    "is fixed for the propagation: `Stats.collect_stats` builds the set for it and "*
    "`Luna.stats_device_capable` tells `Luna.ScaledOutput` which one it is.")

#= `mapreduce` over a lazy `Broadcasted` rather than over `a .* b`, which materialises a
   field-sized temporary before reducing (GPU_PLAN.md section 11, the amendment which added
   `RK45._zipreduce` for the stepper norms). `RK45._zipreduce` is the whole-array form;
   this is the same thing reducing along the frequency axis only, for a state with more
   than one column. =#
@inline _zipreduce1(f, op, init, arrs...) = mapreduce(
    identity, op, Broadcast.instantiate(Broadcast.broadcasted(f, arrs...));
    init=init, dims=1)

#= The same over an arbitrary set of dimensions, which is what a state with transverse
   axes needs: the frequency axis and the transverse axes are reduced together and the
   polarisation axis is kept. =#
@inline _zipreducedims(f, op, init, dims, arrs...) = mapreduce(
    identity, op, Broadcast.instantiate(Broadcast.broadcasted(f, arrs...));
    init=init, dims=dims)

@inline _wabs2(w, e) = w*abs2(e)
@inline _w2abs2(w1, w2, e) = w1*w2*abs2(e)
@inline _wmul(w, x) = w*x
@inline _w2mul(w1, w2, x) = w1*w2*x

"""
    StatsContext(grid, Eω, Eref=1.0)

What a statistic has to know about the run it is being built for:

- `grid`, the frequency/time grid;
- `proto`, a prototype of the propagating state (its array type, element type and shape);
- `spec`, the [`Luna.DeviceSpec`](@ref) derived from `proto`, which is what mirrors and
  buffers are allocated with;
- `Eref`, the field unit the state is expressed in ([`Luna.UnitScaling`](@ref)). `1.0` for
  every `Float64` run. **It applies to the device branches only**: on the host the state
  has already been unscaled by `Luna.ScaledOutput` before the statistics see it;
  [`collect_stats`](@ref) passes `1.0` when it builds a set for the host.
- `ondevice`, whether the set being built will be called with the device state. Each
  statistic stores it and branches on it, rather than on the array type of what it is
  handed (see `Stats._onstate`).

[`collect_stats`](@ref) builds one and passes it to [`prepare`](@ref).
"""
struct StatsContext{G, A, S}
    grid::G
    proto::A
    spec::S
    Eref::Float64
    ondevice::Bool
end

StatsContext(grid, Eω, Eref=1.0, ondevice=Utils.isdevice(Eω)) =
    StatsContext(grid, Eω, _specof(Eω), Float64(Eref), ondevice)

#= The `DeviceSpec` a state array belongs to. `Luna.todevice`/`Luna.alloc` are written
   against a spec rather than an array, and the statistics are the only place which has
   the array but not the spec. `.wrapper` strips the element type and dimensionality from
   the array type (`MtlArray{ComplexF32, 1, ...}` -> `MtlArray`), which is what
   `DeviceSpec` holds. =#
_specof(x::AbstractArray) = Luna.DeviceSpec(_arraytype(x), real(eltype(x)))
_arraytype(::Array) = Array
_arraytype(x::AbstractArray) = Base.typename(typeof(x)).wrapper

"The number of time-domain samples the analytic field has on `grid` (see `plan_analytic`)."
_ntime(grid::Grid.RealGrid) = (length(grid.ω) - 1)*2
_ntime(grid::Grid.EnvGrid) = length(grid.ω)

"The shape of the time-domain buffer for a state shaped like `ctx.proto`."
_timedims(ctx::StatsContext) = (_ntime(ctx.grid), size(ctx.proto)[2:end]...)

#=================================================#
#===========  SPECTRAL-DOMAIN WEIGHTS  ===========#
#=================================================#

#= `Fields.energyfuncs(grid)[2]` integrates |Eω|^2 over the frequency axis with
   `NumericalIntegration`'s `SimpsonEven` (a `RealGrid`) or a plain `sum` (an `EnvGrid`).
   `SimpsonEven` is the alternative extended Simpson rule: every sample with weight 1
   except the first four and last four, which carry 17/48, 59/48, 43/48 and 49/48. As a
   weight vector that is one fused reduction on any array type, where the library's
   version is a scalar loop over a host `Vector`.

   The weights are only ever used on a device; the host branch calls `energyfun_ω` itself,
   so the default CPU path is bit-for-bit what it was. =#
function _spectral_weights(grid::Grid.RealGrid)
    nω = length(grid.ω)
    nω >= 8 || return nothing # fewer samples than the rule has distinct weights
    w = ones(nω)
    w[1] = w[nω] = 17/48
    w[2] = w[nω-1] = 59/48
    w[3] = w[nω-2] = 43/48
    w[4] = w[nω-3] = 49/48
    prefac = 2π/(grid.ω[end]^2) * (grid.ω[2] - grid.ω[1])
    (w, prefac)
end

function _spectral_weights(grid::Grid.EnvGrid)
    nω = length(grid.ω)
    δω = grid.ω[2] - grid.ω[1]
    Δω = nω*δω
    (ones(nω), 2π*δω/(Δω^2))
end

_spectral_weights(grid) = nothing

#= A deterministic probe spectrum, so that the check below neither touches the global RNG
   nor depends on it. =#
function _probespectrum(nω)
    n = 1:nω
    @. exp(-((n - nω/3)/(nω/7))^2) * cis(0.37*n) + 1e-3*cis(1.1*n)
end

"A deterministic probe state of shape `shape`, from the same samples as `_probespectrum`."
_probestate(shape) = reshape(_probespectrum(prod(shape)), shape)

"""
    EnergyWeights

The device form of an energy functional: a scalar `prefac` and one weight vector per
weighted axis, already reshaped to broadcast against the state, such that

    prefac * sum(w[1] .* w[2] .* ... .* abs2.(Eω); dims=dims)

is the functional applied to each column of the state. `dims` is `nothing` for a state
with no column axis, where the reduction is over the whole array and the result is a
scalar. Built and checked by [`_devenergy`](@ref).
"""
struct EnergyWeights{W, D}
    w::W
    dims::D
    prefac::Float64
end

#= Which kernel the reduction folds over depends only on how many weight vectors there
   are, which is a property of the tuple's type. =#
_ekernel(::Tuple{Any}) = _wabs2
_ekernel(::Tuple{Any, Any}) = _w2abs2

"""
    _energy_weights(grid, spacegrid)

`(weights, axes, prefactor, dims)` for the spectral energy functional of `grid` over
`spacegrid`: the host weight vectors, the axis of the state each one belongs to, the
scalar prefactor, and the axes the reduction runs over. `nothing` where there is no such
form.

`spacegrid` is `nothing` for a mode-averaged or multimode state, whose functional
integrates over frequency alone; a [`Grid.RadialGrid`](@ref Luna.Grid.RadialGrid) adds
the quadrature weights of the reciprocal radial axis
([`Grid.integrate_k`](@ref Luna.Grid.integrate_k)); the Cartesian transverse grids add
nothing but a prefactor, because their functional is a plain sum over every axis --
including the frequency axis, which is why they do not take this grid's Simpson weights.
"""
function _energy_weights(grid, ::Nothing)
    ws = _spectral_weights(grid)
    isnothing(ws) && return nothing
    w, prefac = ws
    ((w,), (1,), prefac, (1,))
end

function _energy_weights(grid, rg::Grid.RadialGrid)
    ws = _spectral_weights(grid)
    isnothing(ws) && return nothing
    w, prefacω = ws
    ((w, rg.wk), (1, 3), 2π*c*ε_0/2 * prefacω, (1, 3))
end

function _energy_weights(grid, sg::Grid.Free2DGrid)
    ws = _spectral_weights(grid)
    isnothing(ws) && return nothing
    _, prefacω = ws
    ((ones(length(grid.ω)),), (1,),
     c*ε_0/2 * prefacω * _kfactor(sg.kx), (1, 3))
end

function _energy_weights(grid, sg::Grid.FreeGrid)
    ws = _spectral_weights(grid)
    isnothing(ws) && return nothing
    _, prefacω = ws
    ((ones(length(grid.ω)),), (1,),
     c*ε_0/2 * prefacω * _kfactor(sg.kx) * _kfactor(sg.ky), (1, 3, 4))
end

_energy_weights(grid, spacegrid) = nothing

"The `2πδk/Δk²` each transverse axis of a Cartesian free-space energy functional carries."
function _kfactor(k)
    δk = k[2] - k[1]
    2π*δk/(length(k)*δk)^2
end

"""
    _devenergy(ctx, spacegrid, energyfun_ω)

The [`EnergyWeights`](@ref) which reproduce `energyfun_ω` on a state shaped like
`ctx.proto`, or `nothing` if there are none for this grid or if `energyfun_ω` is not the
spectral energy functional of `ctx.grid` over `spacegrid`.

The agreement is *checked* on a probe state rather than assumed, because [`energy`](@ref)
and [`energy_window`](@ref) take the functional as an argument: a caller who passes
something else (their own functional, or the one for a different transverse grid) gets a
host-only statistic rather than a silently different quantity.
"""
function _devenergy(ctx::StatsContext, spacegrid, energyfun_ω, window=nothing)
    ews = _energy_weights(ctx.grid, spacegrid)
    isnothing(ews) && return nothing
    wh, waxes, prefac, rdims = ews
    nd = ndims(ctx.proto)
    all(<=(nd), waxes) && all(<=(nd), rdims) || return nothing
    #= One column of the state -- the array `energyfun_ω` is written against -- padded
       with a singleton polarisation axis, so that the weights can be broadcast against
       it exactly as they will be against the state. =#
    cshape = (size(ctx.proto, 1), size(ctx.proto)[3:end]...)
    pshape = (size(ctx.proto, 1), 1, size(ctx.proto)[3:end]...)
    Ep = _probestate(pshape)
    ref = try
        energyfun_ω(reshape(Ep, cshape))
    catch
        return nothing
    end
    (ref isa Real && isfinite(ref) && ref > 0) || return nothing
    whr = map((v, ax) -> reshape(v, _axisshape(length(v), ax, nd)), wh, waxes)
    got = prefac*sum(_ekernel(whr).(whr..., Ep))
    abs(got - ref) <= 1e-10*ref || return nothing
    #= A spectral window multiplies the field, so it enters the reduction squared and
       folds into the frequency weights: one vector and one reduction, and no windowed
       copy of the field. It does not affect the check above, which is the identity
       between the unwindowed functional and the weights. =#
    wh = isnothing(window) ? wh : ((wh[1] .* abs2.(window)), Base.tail(wh)...)
    wd = map((v, ax) -> reshape(Luna.todevice(ctx.spec, v),
                                _axisshape(length(v), ax, nd)), wh, waxes)
    EnergyWeights(wd, nd == 1 ? nothing : rdims, prefac)
end

"The shape which puts a vector of length `n` on axis `ax` of an `nd`-dimensional array."
_axisshape(n, ax, nd) = ntuple(i -> i == ax ? n : 1, nd)

#=================================================#
#===============  THE STATISTICS  ================#
#=================================================#

struct CentreOfMass{V, W}
    ω::V        # the frequency axis, as `Maths.moment` takes it (host)
    ωd::W       # the same on the state's array type and precision
    ondevice::Bool
end

"""
    ω0(grid)

Create stats function to calculate the centre of mass (first moment) of the spectral power
density.
"""
ω0(grid) = CentreOfMass(grid.ω, grid.ω, false)

prepare(f::CentreOfMass, ctx::StatsContext) =
    CentreOfMass(f.ω, Luna.todevice(ctx.spec, f.ω), ctx.ondevice)

device_capable(::CentreOfMass) = true
needs_time(::CentreOfMass) = false

function (f::CentreOfMass)(d, Eω, Et, z, dz)
    if _onstate(f, Eω)
        #= The two sums `Maths.moment` takes, each as one reduction over a lazy
           `Broadcasted`: nothing field-sized is materialised. The ratio does not depend
           on the unit scaling. =#
        RT = real(eltype(Eω))
        if ndims(Eω) > 1
            num = Array(_zipreduce1(_wabs2, +, zero(RT), f.ωd, Eω))
            den = Array(sum(abs2, Eω; dims=1))
            d["ω0"] = _dropfreq(Float64.(num) ./ Float64.(den))
        else
            num = _zipreduce(_wabs2, +, zero(RT), f.ωd, Eω)
            d["ω0"] = Float64(num)/Float64(sum(abs2, Eω))
        end
    else
        d["ω0"] = _dropfreq(Maths.moment(f.ω, abs2.(Eω); dim=1))
    end
end

#= What a reduction along the frequency axis leaves behind depends on the rank of what
   was reduced, and only on that: a mode-averaged state (nω,) leaves a scalar, a modal or
   an on-axis state (nω, ncols) one value per column, and a state with transverse axes
   keeps them. This used to be `Stats.squeeze`, with one method for each of the first two
   cases and a `MethodError` for anything else -- which is why the default statistics did
   not run on a radial or free-space state at all, on the host as much as on a device. =#
_dropfreq(x::AbstractArray{<:Any, 1}) = x[1]
_dropfreq(x::AbstractArray) = dropdims(x; dims=1)

#= The `i`th column of a state: everything at index `i` of the second axis, which is the
   mode axis of a modal state and the polarisation axis of a radial or free-space one.
   `copy` rather than a view, because the energy functionals reshape what they are given.
   For a state of rank 2 this is `Eω[:, i]`, which is what the statistics below did
   before they had to cope with transverse axes as well. =#
_column(Eω, i) = copy(selectdim(Eω, 2, i))

#= A statistic which is one number per column: a scalar for a state with no column axis,
   otherwise a vector over it, with everything past the second axis integrated or reduced
   by `f` itself. =#
function _bycolumn(f, Eω)
    ndims(Eω) == 1 && return f(Eω)
    [f(_column(Eω, i)) for i in axes(Eω, 2)]
end

struct SpectralEnergy{F, S, W}
    energyfun_ω::F
    sg::S           # the transverse grid the functional integrates over, or `nothing`
    ew::W           # the device form of the functional, or `nothing`
    Eref2::Float64  # E_ref^2: the energy is quadratic in the field
    key::String
    ondevice::Bool
end

"""
    energy(grid, energyfun_ω)
    energy(grid, spacegrid, energyfun_ω)

Create stats function to calculate the total energy.

The second form is for a radial or free-space state, whose energy functional
(`Fields.energyfuncs(grid, spacegrid)[2]`) integrates over the transverse axes as well;
`spacegrid` is the transverse grid the state is held on and is what tells the device
branch which axes to reduce over.

On a device the integral is the same rule as `energyfun_ω`'s written as one weighted
reduction over a lazy `Broadcasted` (see [`EnergyWeights`](@ref) and `Stats._devenergy`);
the scalar prefactor and the square of the unit scaling are applied to the result, in
`Float64`.
"""
energy(grid, energyfun_ω) =
    SpectralEnergy(energyfun_ω, nothing, nothing, 1.0, "energy", false)

energy(grid, spacegrid, energyfun_ω) =
    SpectralEnergy(energyfun_ω, spacegrid, nothing, 1.0, "energy", false)

function prepare(f::SpectralEnergy, ctx::StatsContext)
    ew = _devenergy(ctx, f.sg, f.energyfun_ω)
    isnothing(ew) && return SpectralEnergy(f.energyfun_ω, f.sg, nothing, 1.0, f.key, false)
    SpectralEnergy(f.energyfun_ω, f.sg, ew, ctx.Eref^2, f.key, ctx.ondevice)
end

device_capable(f::SpectralEnergy) = !isnothing(f.ew)
needs_time(::SpectralEnergy) = false
statlabel(f::SpectralEnergy) = f.key

function (f::SpectralEnergy)(d, Eω, Et, z, dz)
    if _onstate(f, Eω)
        d[f.key] = _weightedenergy(f.ew, f.Eref2, Eω)
    else
        d[f.key] = _bycolumn(f.energyfun_ω, Eω)
    end
end

#= One reduction over a lazy `Broadcasted`, so nothing field-sized is materialised. The
   prefactors are applied on the host, in Float64, to the result: in a scaled Float32 run
   the field is deliberately of order one and the physical units (up to 1e20 apart) would
   not survive being put back on it. =#
function _weightedenergy(ew, Eref2, Eω)
    isnothing(ew) && error(
        "this energy statistic was built from an energy functional with no device form, "*
        "so it cannot be evaluated on a $(typeof(Eω)).")
    RT = real(eltype(Eω))
    fac = ew.prefac*Eref2
    isnothing(ew.dims) &&
        return fac*Float64(_zipreduce(_ekernel(ew.w), +, zero(RT), ew.w..., Eω))
    #= `vec`: the reduction leaves the axes it reduced as singletons, and the result is
       one number per column, in the order the column axis has them. =#
    s = Array(_zipreducedims(_ekernel(ew.w), +, zero(RT), ew.dims, ew.w..., Eω))
    [fac*Float64(x) for x in vec(s)]
end

"""
    energy_λ(grid, energyfun_ω, λlims; label)
    energy_λ(grid, spacegrid, energyfun_ω, λlims; label)

Create stats function to calculate the energy in a wavelength region given by `λlims`.
If `label` is omitted, the stats dataset is named by the wavelength limits. The second
form is for a radial or free-space state, as for [`energy`](@ref).
"""
function energy_λ(grid, energyfun_ω, λlims; label=nothing, winwidth=0)
    window, label = _λwindow(grid, λlims, label, winwidth)
    energy_window(grid, energyfun_ω, window; label=label)
end

@doc (@doc energy_λ)
function energy_λ(grid, spacegrid, energyfun_ω, λlims; label=nothing, winwidth=0)
    window, label = _λwindow(grid, λlims, label, winwidth)
    energy_window(grid, spacegrid, energyfun_ω, window; label=label)
end

function _λwindow(grid, λlims, label, winwidth)
    λlims = collect(λlims)
    ωmin, ωmax = extrema(wlfreq.(λlims))
    window = Maths.planck_taper(grid.ω, ωmin-winwidth, ωmin, ωmax, ωmax+winwidth)
    if isnothing(label)
        λnm = 1e9.*λlims
        label = @sprintf("%.2fnm_%.2fnm", minimum(λnm), maximum(λnm))
    end
    (window, label)
end

struct SpectralEnergyWindow{F, S, V, W}
    energyfun_ω::F
    sg::S           # the transverse grid the functional integrates over, or `nothing`
    window::V       # the host window, as given
    ew::W           # the device form with the window folded in, or `nothing`
    Eref2::Float64
    key::String
    ondevice::Bool
end

"""
    energy_window(grid, energyfun_ω, window; label)
    energy_window(grid, spacegrid, energyfun_ω, window; label)

Create stats function to calculate the energy filtered by a `window`. The stats dataset will
be named `energy_[label]`. The second form is for a radial or free-space state, as for
[`energy`](@ref).
"""
energy_window(grid, energyfun_ω, window::AbstractVector{<:Real}; label) =
    SpectralEnergyWindow(energyfun_ω, nothing, window, nothing, 1.0, "energy_$label",
                         false)

energy_window(grid, spacegrid, energyfun_ω, window::AbstractVector{<:Real}; label) =
    SpectralEnergyWindow(energyfun_ω, spacegrid, window, nothing, 1.0, "energy_$label",
                         false)

function prepare(f::SpectralEnergyWindow, ctx::StatsContext)
    ew = _devenergy(ctx, f.sg, f.energyfun_ω, f.window)
    isnothing(ew) && return SpectralEnergyWindow(f.energyfun_ω, f.sg, f.window, nothing,
                                                 1.0, f.key, false)
    SpectralEnergyWindow(f.energyfun_ω, f.sg, f.window, ew, ctx.Eref^2, f.key,
                         ctx.ondevice)
end

device_capable(f::SpectralEnergyWindow) = !isnothing(f.ew)
needs_time(::SpectralEnergyWindow) = false
statlabel(f::SpectralEnergyWindow) = f.key

function (f::SpectralEnergyWindow)(d, Eω, Et, z, dz)
    if _onstate(f, Eω)
        d[f.key] = _weightedenergy(f.ew, f.Eref2, Eω)
    else
        d[f.key] = _bycolumn(Eωi -> f.energyfun_ω(Eωi.*f.window), Eω)
    end
end

struct PeakPower{V}
    t::V
    Eref2::Float64
    ondevice::Bool
end

"""
    peakpower(grid)

Create stats function to calculate the peak power.
"""
peakpower(grid) = PeakPower(grid.t, 1.0, false)

prepare(f::PeakPower, ctx::StatsContext) = PeakPower(f.t, ctx.Eref^2, ctx.ondevice)

device_capable(::PeakPower) = true

function (f::PeakPower)(d, Eω, Et, z, dz)
    if _onstate(f, Et)
        if ndims(Et) > 1
            # `maximum(abs2, Et; dims=1)`, not `maximum(abs2.(Et), dims=1)`: the second
            # materialises the whole block first
            pp = Array(maximum(abs2, Et; dims=1))
            d["peakpower"] = [f.Eref2*Float64(pp[1, i]) for i in axes(pp, 2)]
            d["peakpower_allmodes"] = f.Eref2*Float64(maximum(sum(abs2, Et; dims=2)))
        else
            d["peakpower"] = f.Eref2*Float64(maximum(abs2, Et))
        end
    elseif ndims(Et) > 1
        d["peakpower"] = dropdims(maximum(abs2.(Et), dims=1), dims=1)
        d["peakpower_allmodes"] = maximum(eachindex(f.t)) do ii
            sum(abs2, Et[ii, :])
        end
    else
        d["peakpower"] = maximum(abs2, Et)
    end
end

"""
    peakpower(grid, Eω, window; label)

Create stats function to calculate the peak power within a frequency range defined by the
window function `window`. `window` must have the same length as `grid.ω`. The stats
dataset is labeled as `peakpower_[label]`.

Host-only: it plans its own inverse transform and holds its own buffers, which are built
for the `Eω` given here.
"""
function peakpower(grid, Eω, window::Vector{<:Real}; label)
    Etbuf, analytic! = plan_analytic(grid, Eω) # output buffer and function for inverse FT
    Eωbuf = similar(Eω) # buffer for Eω with window applied
    Pt = zeros((length(grid.t), size(Eωbuf, 2)))
    key = "peakpower_$label"
    function addstat!(d, Eω, Et, z, dz)
        Eωbuf .= Eω .* window
        analytic!(Etbuf, Eωbuf)
        if ndims(Etbuf) > 1
            Pt .= abs2.(Etbuf)
            d[key] = dropdims(maximum(Pt, dims=1), dims=1)
            d[key*"_allmodes"] = maximum(eachindex(grid.t)) do ii
                sum(Pt[ii, :]; dims=2)
            end
        else
            d[key] = maximum(abs2, Etbuf)
        end
    end
end

"""
    peakpower(grid, Eω, λlims; label=nothing)

Create stats function to calculate the peak power within a frequency range defined by the
wavelength limits `λlims`. If `label` is given, the stats dataset is labeled as
`peakpower_[label]`, otherwise `label` is created automatically from `λlims`.
"""
function peakpower(grid, Eω, λlims::NTuple{2, <:Real}; label=nothing, winwidth=:auto)
    window = Processing.ωwindow_λ(grid.ω, λlims; winwidth=winwidth)
    if isnothing(label)
        λnm = 1e9.*λlims
        label = @sprintf("%.2fnm_%.2fnm", minimum(λnm), maximum(λnm))
    end
    peakpower(grid, Eω, window; label=label)
end


struct PeakIntensityAeff{A}
    aeff::A
    Eref2::Float64
    ondevice::Bool
end

"""
    peakintensity(grid, aeff)

Create stats function to calculate the mode-averaged peak intensity given the effective area
`aeff(z)`.
"""
peakintensity(grid, aeff) = PeakIntensityAeff(aeff, 1.0, false)

prepare(f::PeakIntensityAeff, ctx::StatsContext) =
    PeakIntensityAeff(f.aeff, ctx.Eref^2, ctx.ondevice)

device_capable(::PeakIntensityAeff) = true

function (f::PeakIntensityAeff)(d, Eω, Et, z, dz)
    if _onstate(f, Et)
        d["peakintensity"] = f.Eref2*Float64(maximum(abs2, Et))/f.aeff(z)
    else
        d["peakintensity"] = maximum(abs2, Et)/f.aeff(z)
    end
end

struct PeakIntensityModes{T, B}
    tospace::T
    Et0::B          # host buffer for the projected field
    Eref2::Float64
    ondevice::Bool
end

"""
    peakintensity(grid, modes; components=:y)

Create stats function to calculate the peak intensity for several modes.

Device-capable for a single mode, where projecting onto the transverse plane is a scalar
factor per polarisation component and comes out of the reduction.
"""
function peakintensity(grid, modes::Modes.ModeCollection; components=:y)
    tospace = Modes.ToSpace(modes, components=components)
    Et0 = zeros(ComplexF64, (length(grid.t), tospace.npol))
    PeakIntensityModes(tospace, Et0, 1.0, false)
end

prepare(f::PeakIntensityModes, ctx::StatsContext) =
    PeakIntensityModes(f.tospace, f.Et0, ctx.Eref^2,
                       ctx.ondevice && device_capable(f))

device_capable(f::PeakIntensityModes) = f.tospace.nmodes == 1

function (f::PeakIntensityModes)(d, Eω, Et, z, dz)
    if _onstate(f, Et)
        d["peakintensity"] = c*ε_0/2 * f.Eref2 * _onaxisfactor(f.tospace, z) *
                             Float64(maximum(abs2, Et))
    else
        tospace = f.tospace
        Modes.to_space!(f.Et0, Et, (0, 0), tospace; z=z)
        if tospace.npol > 1
            d["peakintensity"] = c*ε_0/2 * maximum(axes(f.Et0, 1)) do ii
                sum(abs2, f.Et0[ii, :])
            end
        else
            d["peakintensity"] = c*ε_0/2 * maximum(abs2, f.Et0)
        end
    end
end

#= The single mode's on-axis field, summed in quadrature over the polarisation components
   kept. `Modes.to_space!` multiplies the field by that mode field sample by sample; with
   one mode the modulus squared of the result is a constant times `abs2(Et)` and the
   constant comes out of the maximum. Evaluated per step on the host because the mode's
   transverse profile can depend on z; it is a handful of scalar operations. =#
function _onaxisfactor(ts::Modes.ToSpace, z)
    ts.nmodes == 1 || error(
        "the on-axis peak intensity of $(ts.nmodes) modes does not reduce to a scalar "*
        "factor and has no device form.")
    E = Modes.Exy(ts.ms[1], (0, 0), z=z)[ts.indices]
    sum(abs2, E)
end

struct FWHMt{V, B, H}
    t::V
    Pd::B       # |Et|^2 on the state's array type
    Ph::H       # the same on the host
    ondevice::Bool
end

"""
    fwhm_t(grid)

Create stats function to calculate the temporal FWHM (pulse duration) for mode average.

On a device only `|E|^2` is copied to the host; the FWHM itself is then found by the same
root-finding on the same samples as on the host.
"""
fwhm_t(grid) = FWHMt(grid.t, nothing, nothing, false)

function prepare(f::FWHMt, ctx::StatsContext)
    ctx.ondevice || return FWHMt(f.t, nothing, nothing, false)
    Pd = Luna.alloc(ctx.spec, Luna.realtype(ctx.spec), _timedims(ctx))
    FWHMt(f.t, Pd, Array{Luna.realtype(ctx.spec)}(undef, size(Pd)), true)
end

device_capable(::FWHMt) = true

function (f::FWHMt)(d, Eω, Et, z, dz)
    Pt = if _onstate(f, Et)
        f.Pd .= abs2.(Et)
        copyto!(f.Ph, f.Pd)
        f.Ph
    else
        abs2.(Et)
    end
    if ndims(Pt) > 1
        Ptsum = dropdims(sum(Pt; dims=2); dims=2)
        d["fwhm_t_min"] = [Maths.fwhm(f.t, Pt[:, i], method=:linear)
                          for i=1:size(Pt, 2)]
        d["fwhm_t_max"] = [Maths.fwhm(f.t, Pt[:, i], method=:linear, minmax=:max)
                          for i=1:size(Pt, 2)]
        d["fwhm_t_min_allmodes"] = Maths.fwhm(f.t, Ptsum, method=:linear)
        d["fwhm_t_max_allmodes"] = Maths.fwhm(f.t, Ptsum, method=:linear, minmax=:max)
    else
        d["fwhm_t_min"] = Maths.fwhm(f.t, Pt, method=:linear, minmax=:min)
        d["fwhm_t_max"] = Maths.fwhm(f.t, Pt, method=:linear, minmax=:max)
    end
end

struct FWHMr{T, B}
    tospace::T
    Eω0::B
end

"""
    fwhm_r(grid, modes; components=:y)

Create stats function to calculate the radial FWHM (aka beam size) in a modal propagation.

Host-only: the beam size is found by root-finding on a function of the radius, each
evaluation of which projects the whole spectrum onto one transverse point. There is no
reduction over the state to move to a device, and the modal transform which produces such
a state is host-only in any case.
"""
function fwhm_r(grid, modes; components=:y)
    tospace = Modes.ToSpace(modes, components=components)
    FWHMr(tospace, zeros(ComplexF64, (length(grid.ω), tospace.npol)))
end

device_capable(::FWHMr) = false
needs_time(::FWHMr) = false

function (f::FWHMr)(d, Eω, Et, z, dz)
    function g(r)
        Modes.to_space!(f.Eω0, Eω, (r, 0), f.tospace; z=z)
        sum(abs2, f.Eω0)
    end
    d["fwhm_r"] = 2*Maths.hwhm(g)
end

struct ElectronDensityAeff{R, RD, D, A, V, B, W}
    ratefunc::R     # the rate as given: host, physical units
    ratedev::RD     # the same in the run's precision and array type
    dfun::D
    aeff::A
    oversampling::Int
    t::V            # the time axis of the grid, host (what `oversample` is given)
    δt::Float64     # the step of the *oversampled* axis
    frac::Vector{Float64}   # host buffer: the rate, then the ionisation fraction
    rate::B         # device buffer for the rate
    w::W            # trapezoid weights on the state's array type
    Eref::Float64
    ondevice::Bool
end

"""
    electrondensity(grid, ionrate, dfun, aeff; oversampling=1)

Create stats function to calculate the maximum electron density in mode average.

If oversampling > 1, the field is oversampled before the calculation
!!! warning
    Oversampling can lead to a significant performance hit

Device-capable when `oversampling == 1` and the rate has a device kernel
([`Ionisation.device_capable`](@ref Luna.Ionisation.device_capable)). The device branch
evaluates the rate straight off the state, with the intensity conversion and the unit
scaling folded into the field reference the kernel multiplies each sample by, and takes
the ionisation integral as one weighted reduction: only the end point of the cumulative
integral is ever read, and that end point is the trapezoid rule over the whole window.
"""
function electrondensity(grid::Grid.RealGrid, ionrate!, dfun, aeff; oversampling=1)
    to, _ = Maths.oversample(grid.t, complex(grid.t), factor=oversampling)
    ElectronDensityAeff(ionrate!, ionrate!, dfun, aeff, oversampling, grid.t,
                        to[2]-to[1], similar(to), nothing, nothing, 1.0, false)
end

function prepare(f::ElectronDensityAeff, ctx::StatsContext)
    (ctx.ondevice && device_capable(f)) || return f
    n = length(f.frac)
    w = ones(n)
    w[1] = w[n] = 0.5
    ElectronDensityAeff(f.ratefunc, Ionisation.device_rate(f.ratefunc, ctx.spec),
                        f.dfun, f.aeff, f.oversampling, f.t, f.δt, f.frac,
                        Luna.alloc(ctx.spec, Luna.realtype(ctx.spec), (n,)),
                        Luna.todevice(ctx.spec, w), ctx.Eref, true)
end

device_capable(f::ElectronDensityAeff) =
    f.oversampling == 1 && Ionisation.device_capable(f.ratefunc)

function (f::ElectronDensityAeff)(d, Eω, Et, z, dz)
    if _onstate(f, Et)
        #= `Eref/sqrt(...)` is what the kernel multiplies each sample by: it takes the
           scaled state to the physical field and the physical field to the one which
           drives the rate. The host branch divides the field itself instead. =#
        Erate = f.Eref/sqrt(ε_0*c*f.aeff(z)/2)
        rk = Ionisation.ratekernel(f.ratedev, Luna.scalar(f.rate, Erate))
        f.rate .= rk.(real.(Et))
        ratemax = Float64(maximum(f.rate))
        # one reduction over a lazy `Broadcasted`, so no temporary the size of `rate`
        intg = f.δt*Float64(_zipreduce(_wmul, +, zero(eltype(f.rate)), f.w, f.rate))
        d["electrondensity"] = (1 - exp(-intg))*f.dfun(z)
        d["peak_ionisation_rate"] = ratemax
    else
        # note: oversampling returns its arguments without any work done if factor==1
        _, Eto = Maths.oversample(f.t, Et, factor=f.oversampling)
        @. Eto /= sqrt(ε_0*c*f.aeff(z)/2)
        ratemax = _ionfrac!(f.frac, f.ratefunc, real(Eto), f.δt)
        d["electrondensity"] = f.frac[end]*f.dfun(z)
        d["peak_ionisation_rate"] = ratemax
    end
end

#= `out` holds the time-dependent ionisation fraction on return; the maximum ionisation
   rate is the return value. The serial `Maths.cumtrapz!` is deliberate: this is the host
   branch, and it is the arithmetic Luna has always done here. =#
function _ionfrac!(out, ionrate!, Et, δt)
    ionrate!(out, Et)
    ratemax = maximum(out)
    Maths.cumtrapz!(out, δt) # in-place cumulative integration
    @. out = 1 - exp(-out)
    return ratemax
end

struct ElectronDensityModes{R, D, T, V, B}
    ratefunc::R
    dfun::D
    tospace::T
    oversampling::Int
    t::V            # the time axis of the grid, host (what `oversample` is given)
    δt::Float64     # the step of the *oversampled* axis
    frac::Vector{Float64}
    Et0::B
end

"""
    electrondensity(grid, ionrate, dfun, modes; oversampling=1)

Create stats function to calculate the maximum electron density for multimode simulations.

If oversampling > 1, the field is oversampled before the calculation
!!! warning
    Oversampling can lead to a significant performance hit

Host-only: the field has to be projected onto the transverse plane first, and the modal
transform which produces the state is host-only in any case.
"""
function electrondensity(grid::Grid.RealGrid, ionrate!, dfun,
                         modes::Modes.ModeCollection;
                         components=:y, oversampling=1)
    to, _ = Maths.oversample(grid.t, complex(grid.t), factor=oversampling)
    tospace = Modes.ToSpace(modes, components=components)
    ElectronDensityModes(ionrate!, dfun, tospace, oversampling, grid.t, to[2]-to[1],
                         similar(to), zeros(ComplexF64, (length(to), tospace.npol)))
end

device_capable(::ElectronDensityModes) = false

function (f::ElectronDensityModes)(d, Eω, Et, z, dz)
    # note: oversampling returns its arguments without any work done if factor==1
    _, Eto = Maths.oversample(f.t, Et, factor=f.oversampling)
    Modes.to_space!(f.Et0, Eto, (0, 0), f.tospace; z=z)
    if f.tospace.npol > 1
        ratemax = _ionfrac!(f.frac, f.ratefunc,
                            hypot.(real(f.Et0[:, 1]), real(f.Et0[:, 2])), f.δt)
    else
        ratemax = _ionfrac!(f.frac, f.ratefunc, real(f.Et0[:, 1]), f.δt)
    end
    d["electrondensity"] = f.frac[end]*f.dfun(z)
    d["peak_ionisation_rate"] = ratemax
end

struct ModeReconstructionError{T, A, B}
    t::T
    Prω_recon::A
    difference::A
    nl::B
end

"""
    mode_reconstruction_error(t::TransModal)

Create a stats function to calculate and collect the mode reconstruction error in the
induced polarisation on axis at every step.

Host-only: it evaluates the modal transform itself, which is host code.
"""
function mode_reconstruction_error(t::TransModal)
    Prω_recon = similar(t.Prω)
    ModeReconstructionError(t, Prω_recon, similar(Prω_recon), similar(t.Emω))
end

device_capable(::ModeReconstructionError) = false
needs_time(::ModeReconstructionError) = false

function (f::ModeReconstructionError)(d, Eω, Et, z, dz)
    t = f.t
    t(f.nl, Eω, z)
    x = (0.0, 0.0) # on-axis coordinate
    # reconstruct
    Modes.to_space!(f.Prω_recon, f.nl, x, t.ts, z=z)
    # in going to modes and back we've picked up two factors of the mode normalisation
    f.Prω_recon .*= 1/2*sqrt(PhysData.ε_0/PhysData.μ_0)
    Erω_to_Prω!(t, x)
    f.difference .= f.Prω_recon .- t.Prω
    d["mode_reconstruction_error"] =
        sqrt(sum(abs2, f.difference))/sqrt(sum(abs2, f.Prω_recon))
    d["transverse_points"] = float(t.ncalls) # convert to Float64 to enable NaN padding
    d["transverse_integral_error_abs"] = sqrt(sum(abs2, t.err)/length(t.err))
    d["transverse_integral_error_rel"] =
        d["transverse_integral_error_abs"]/sqrt(sum(abs2, f.nl)/length(f.nl))
end

struct TransverseIntegralError{TT, A, B}
    t::TT           # the TransModalFixed
    nl::A           # the modal polarisation, where the transform's buffers live
    Eh::B           # host staging buffer, in the transform's element type
    Ein::A          # the state in the transform's units and array type (=== Eh on a host
                    # transform, so that nothing is copied there)
    Eref::Float64   # the unit scaling the transform works in
    invEref::Float64
    haserr::Bool
    points::Float64
end

"""
    transverse_integral_error(t::TransModalFixed)

Create a stats function which records the embedded error estimate of the fixed transverse
quadrature rule, and the number of nodes that rule uses. It is what
[`mode_reconstruction_error`](@ref) is for the adaptive transform, and it records the same
three datasets as that one's error part:

- `transverse_points`: the node count of the rule, which is fixed for the propagation;
- `transverse_integral_error_abs`: the root-mean-square over `(ω, mode)` of
  `P_coarse - P_fine`, the difference between the embedded coarse rule and the rule the
  propagation uses
  ([`NonlinearRHS.integral_error!`](@ref Luna.NonlinearRHS.integral_error!));
- `transverse_integral_error_rel`: the same divided by the root-mean-square of the
  polarisation itself.

Both error datasets are `NaN` when the rule has no embedded coarse rule
([`NonlinearRHS.has_error_estimate`](@ref Luna.NonlinearRHS.has_error_estimate)), i.e.
without `modal_kronrod=true` on a radial or Cartesian rule. `transverse_points` is
recorded either way.

The transform is evaluated once per call on the accepted state, as
[`mode_reconstruction_error`](@ref) does, so the estimate belongs to the state which was
saved rather than to whichever stage the stepper last called.

Host-only: like every other statistic in the multimode default set. The transform itself
may be on a device, in which case the state is staged into the transform's units and array
type first and the error is scaled back to physical units afterwards.
"""
function transverse_integral_error(t::TransModalFixed)
    nl = similar(t.err)
    Eh = Array{eltype(t.err)}(undef, size(t.err))
    Ein = Utils.isdevice(t.err) ? similar(t.err) : Eh
    Eref = Float64(Luna.runscaling(t).Eref)
    TransverseIntegralError(t, nl, Eh, Ein, Eref, inv(Eref), has_error_estimate(t),
                            float(t.ncalls))
end

device_capable(::TransverseIntegralError) = false
needs_time(::TransverseIntegralError) = false

function (f::TransverseIntegralError)(d, Eω, Et, z, dz)
    #= The state arrives on the host in physical units (the set is host-only, so
       `Luna.ScaledOutput` has unscaled it); the transform works in its own units and on
       its own array type. `invEref` is exactly 1 for every Float64 run. =#
    @. f.Eh = Eω * f.invEref
    f.Ein === f.Eh || copyto!(f.Ein, f.Eh)
    f.t(f.nl, f.Ein, z)
    d["transverse_points"] = f.points
    if f.haserr
        err = integral_error!(f.t)
        rms = sqrt(_meanabs2(err))
        d["transverse_integral_error_abs"] = f.Eref*rms
        d["transverse_integral_error_rel"] = rms/sqrt(_meanabs2(f.nl))
    else
        d["transverse_integral_error_abs"] = NaN
        d["transverse_integral_error_rel"] = NaN
    end
end

"The mean of `abs2` over a whole array, as one reduction, in `Float64`."
_meanabs2(x) = Float64(sum(abs2, x))/length(x)

#=================================================#
#======  RADIAL AND FREE-SPACE STATISTICS  =======#
#=================================================#

#= A radial or free-space state lives in transverse *reciprocal* space: `(nω, npol, nk)`
   on a `Grid.RadialGrid`, `(nω, npol, nkx)` on a `Grid.Free2DGrid` and
   `(nω, npol, nkx, nky)` on a `Grid.FreeGrid`. Two things follow.

   The total energy is a functional of the whole state, because Parseval's theorem is what
   makes the transverse integral available there at all
   (`Fields.energyfuncs(grid, spacegrid)`); that is the `spacegrid` method of `energy`
   above.

   Everything which is a property of the pulse rather than of the beam -- its peak
   intensity, its duration, the centre of mass of its spectrum -- is taken **on the
   propagation axis**, because the transverse average of a duration is not a duration. The
   field on axis is a fixed linear combination of the transverse samples, so it is one
   weighted reduction over the transverse axes, after which the statistics above apply
   unchanged to an `(nω[, npol])` field. `OnAxis` is that wrapper. =#

struct OnAxisProjection{W, D, S, B, F}
    w::W            # one weight array per transverse axis, shaped to broadcast
    dims::D         # the transverse axes, reduced over
    shape::S        # the shape of the on-axis spectral field
    Et0::B          # buffer for the on-axis analytic time-domain field
    analytic!::F
end

"""
    _axisweights(x)

The weights which take a field held as its DFT along a Cartesian transverse axis to its
sample at `x = 0`.

`FFTW.ifft` is `E[n] = (1/N) Σ_m Ek[m] exp(2πi(n-1)(m-1)/N)`, which is the transform that
takes a free-space state back to real space, so the on-axis sample is that sum with `n`
the index at which `x` is zero. [`Grid.FreeGrid`](@ref Luna.Grid.FreeGrid) and
[`Grid.Free2DGrid`](@ref Luna.Grid.Free2DGrid) centre their transverse axes, so `n - 1` is
`N/2` and the weights are exactly `±1/N` -- `cispi` of an integer is exactly `±1`, so
nothing rounds. An axis which does not contain `x = 0` gives genuinely complex weights,
which is the correct interpolation to the axis rather than an approximation of it.
"""
function _axisweights(x)
    N = length(x)
    n0 = argmin(abs.(x))
    [cispi(2*(n0 - 1)*(m - 1)/N)/N for m in 1:N]
end

"""
    _onaxis_vectors(spacegrid)

`(weights, axes)`: the host weight vector for each transverse axis of `spacegrid` and the
axis of the state it belongs to. On a radial grid this is
[`Grid.onaxis`](@ref Luna.Grid.onaxis), the integral over the reciprocal radial axis; on a
Cartesian grid it is [`_axisweights`](@ref) per axis.
"""
function _onaxis_vectors(rg::Grid.RadialGrid)
    rg.order == 0 || throw(DomainError(
        rg.order, "the on-axis field is only defined for a 0-order RadialGrid"))
    ((rg.wk,), (3,))
end

_onaxis_vectors(sg::Grid.Free2DGrid) = ((_axisweights(sg.x),), (3,))
_onaxis_vectors(sg::Grid.FreeGrid) = ((_axisweights(sg.x), _axisweights(sg.y)), (3, 4))

function _projection(spacegrid, ctx::StatsContext)
    vecs, axes_ = _onaxis_vectors(spacegrid)
    nd = ndims(ctx.proto)
    w = map((v, ax) -> reshape(Luna.todevice(ctx.spec, v),
                               _axisshape(length(v), ax, nd)), vecs, axes_)
    nω = size(ctx.proto, 1)
    npol = size(ctx.proto, 2)
    shape = npol == 1 ? (nω,) : (nω, npol)
    proto0 = Luna.alloc(ctx.spec, eltype(ctx.proto), shape)
    Et0, analytic! = plan_analytic(ctx.grid, proto0)
    OnAxisProjection(w, axes_, shape, Et0, analytic!), proto0
end

#= One reduction over a lazy `Broadcasted`: nothing the size of the state is
   materialised, and what is allocated per call is the `(nω, npol)` result, which is
   smaller than the state by the number of transverse points. =#
function _project(p::OnAxisProjection, Eω)
    red = _zipreducedims(_pkernel(p.w), +, zero(eltype(Eω)), p.dims, p.w..., Eω)
    reshape(red, p.shape)
end

_pkernel(::Tuple{Any}) = _wmul
_pkernel(::Tuple{Any, Any}) = _w2mul

struct OnAxis{S, P, F}
    spacegrid::S
    proj::P         # `nothing` until `prepare`
    f::F
end

"""
    onaxis(f, spacegrid)

The statistics function `f` evaluated on the field on the propagation axis, for a state
held in the transverse reciprocal space of `spacegrid`.

`f` is handed the on-axis spectral field, of shape `(nω)` for one polarisation component
and `(nω, npol)` for more, and its analytic time-domain field -- exactly the shapes the
mode-averaged and multimode statistics are written for, so `ω0(grid)`, `fwhm_t(grid)`,
`peakintensity(grid)` and `electrondensity(grid, ...)` all work through it.

The projection is one weighted reduction over the transverse axes and one inverse
transform of the result, which is `(nω, npol)` rather than the whole state. It is done
per wrapper, so a set with three on-axis statistics reduces the state three times; the
reductions are a small part of the step which produced it, and sharing the result between
statistics would mean caching it against `z`, which a repeated `z` would silently get
wrong.
"""
onaxis(f, spacegrid) = OnAxis(spacegrid, nothing, f)

function prepare(o::OnAxis, ctx::StatsContext)
    proj, proto0 = _projection(o.spacegrid, ctx)
    ctx0 = StatsContext(ctx.grid, proto0, ctx.Eref, ctx.ondevice)
    OnAxis(o.spacegrid, proj, prepare(o.f, ctx0))
end

device_capable(o::OnAxis) = device_capable(o.f)
statlabel(o::OnAxis) = statlabel(o.f)*" (on axis)"

#= The projection does its own, much smaller, inverse transform, so an on-axis statistic
   never reads the collector's `Et` -- which for a free-space state is the largest array
   in the whole statistics set. =#
needs_time(::OnAxis) = false

function (o::OnAxis)(d, Eω, Et, z, dz)
    Eω0 = _project(o.proj, Eω)
    needs_time(o.f) && o.proj.analytic!(o.proj.Et0, Eω0)
    o.f(d, Eω0, o.proj.Et0, z, dz)
end

struct PeakIntensityField
    Eref2::Float64
    ondevice::Bool
end

"""
    peakintensity(grid)

Create stats function to calculate the peak intensity of a field which is already an
electric field in V/m -- a radial or free-space state projected onto the propagation axis
([`onaxis`](@ref)) -- rather than a mode-averaged state, whose intensity needs an
effective area. Polarisation components are summed in quadrature.
"""
peakintensity(grid) = PeakIntensityField(1.0, false)

prepare(f::PeakIntensityField, ctx::StatsContext) =
    PeakIntensityField(ctx.Eref^2, ctx.ondevice)

device_capable(::PeakIntensityField) = true

function (f::PeakIntensityField)(d, Eω, Et, z, dz)
    fac = _onstate(f, Et) ? f.Eref2 : 1.0
    d["peakintensity"] = c*ε_0/2 * fac * Float64(_peakabs2(Et))
end

#= Both branches are the same two reductions -- there is no older host expression for
   this statistic to be a transcription of -- and differ only in the unit scaling the
   device branch has to put back. =#
_peakabs2(Et::AbstractVector) = maximum(abs2, Et)
_peakabs2(Et) = size(Et, 2) == 1 ? maximum(abs2, Et) : maximum(sum(abs2, Et; dims=2))

struct ElectronDensityField{R, RD, D, V, B, W}
    ratefunc::R     # the rate as given: host, physical units
    ratedev::RD     # the same in the run's precision and array type
    dfun::D
    oversampling::Int
    t::V            # the time axis of the grid, host
    δt::Float64     # the step of the *oversampled* axis
    frac::Vector{Float64}   # host buffer: the rate, then the ionisation fraction
    rate::B         # device buffer for the rate
    w::W            # trapezoid weights on the state's array type
    Eref::Float64
    ondevice::Bool
end

"""
    electrondensity(grid, ionrate, dfun; oversampling=1)

Create stats function to calculate the maximum electron density from a field which is
already an electric field in V/m -- a radial or free-space state projected onto the
propagation axis ([`onaxis`](@ref)). Two polarisation components drive the rate through
their quadrature sum, as in the multimode method.

Device-capable on the same terms as the mode-averaged method: `oversampling == 1` and a
rate with a device kernel. The device branch takes the ionisation integral as one
weighted reduction, since only the end point of the cumulative integral is ever read and
that end point is the trapezoid rule over the whole window.
"""
function electrondensity(grid::Grid.RealGrid, ionrate!, dfun; oversampling=1)
    to, _ = Maths.oversample(grid.t, complex(grid.t), factor=oversampling)
    ElectronDensityField(ionrate!, ionrate!, dfun, oversampling, grid.t, to[2]-to[1],
                         similar(to), nothing, nothing, 1.0, false)
end

function prepare(f::ElectronDensityField, ctx::StatsContext)
    (ctx.ondevice && device_capable(f)) || return f
    n = length(f.frac)
    w = ones(n)
    w[1] = w[n] = 0.5
    ElectronDensityField(f.ratefunc, Ionisation.device_rate(f.ratefunc, ctx.spec),
                         f.dfun, f.oversampling, f.t, f.δt, f.frac,
                         Luna.alloc(ctx.spec, Luna.realtype(ctx.spec), (n,)),
                         Luna.todevice(ctx.spec, w), ctx.Eref, true)
end

device_capable(f::ElectronDensityField) =
    f.oversampling == 1 && Ionisation.device_capable(f.ratefunc)

function (f::ElectronDensityField)(d, Eω, Et, z, dz)
    if _onstate(f, Et)
        #= The field already drives the rate; the only conversion is the unit scaling,
           which is folded into the reference the kernel multiplies each sample by. =#
        rk = Ionisation.ratekernel(f.ratedev, Luna.scalar(f.rate, f.Eref))
        _ratebroadcast!(f.rate, rk, Et)
        ratemax = Float64(maximum(f.rate))
        intg = f.δt*Float64(_zipreduce(_wmul, +, zero(eltype(f.rate)), f.w, f.rate))
        d["electrondensity"] = (1 - exp(-intg))*f.dfun(z)
        d["peak_ionisation_rate"] = ratemax
    else
        # note: oversampling returns its arguments without any work done if factor==1
        _, Eto = Maths.oversample(f.t, Et, factor=f.oversampling)
        ratemax = _ionfrac!(f.frac, f.ratefunc, _drivefield(Eto), f.δt)
        d["electrondensity"] = f.frac[end]*f.dfun(z)
        d["peak_ionisation_rate"] = ratemax
    end
end

"The real field which drives the ionisation rate: the quadrature sum over polarisation."
_drivefield(Et::AbstractVector) = real.(Et)
_drivefield(Et) = size(Et, 2) == 1 ? real.(view(Et, :, 1)) :
                  hypot.(real.(view(Et, :, 1)), real.(view(Et, :, 2)))

_ratebroadcast!(rate, rk, Et::AbstractVector) = (rate .= rk.(real.(Et)))

function _ratebroadcast!(rate, rk, Et)
    if size(Et, 2) == 1
        rate .= rk.(real.(view(Et, :, 1)))
    else
        rate .= rk.(hypot.(real.(view(Et, :, 1)), real.(view(Et, :, 2))))
    end
end

struct BeamProfile{S, B, X, A, M, W}
    sg::S           # the transverse grid
    Er::B           # host buffer for the transverse real-space field
    xfrm!::X        # the inverse transverse transform, host
    axis::A         # the transverse coordinate the width is measured on
    mask::M         # which samples are inside the absorber collar, or `nothing`
    measure::W      # transverse integration weights, or `nothing` for a uniform grid
end

"""
    beam_profile(grid, spacegrid; collar=Boundaries.DEFAULT_RCOLLAR)

Create stats function to record the transverse size of the beam and the fraction of its
energy which has reached the absorbing collar at the edge of the transverse grid.

It records, from the transverse fluence profile `Σ_ω,pol |E(ω, pol, r)|²` in **real**
space:

- `fwhm_r` on a [`Grid.RadialGrid`](@ref Luna.Grid.RadialGrid), from the profile mirrored
  about the axis ([`Grid.rsymmetric`](@ref Luna.Grid.rsymmetric)), with the on-axis sample
  supplied by [`Grid.onaxis`](@ref Luna.Grid.onaxis);
- `fwhm_x` on a [`Grid.Free2DGrid`](@ref Luna.Grid.Free2DGrid), and `fwhm_x`/`fwhm_y` on a
  [`Grid.FreeGrid`](@ref Luna.Grid.FreeGrid), from the cuts through the axis;
- `collar_energy_fraction`, the fraction of the transverse integral of the profile which
  lies where the transverse absorber ([`Boundaries.rprofile`](@ref
  Luna.Boundaries.rprofile)) is active. `collar` is the collar width as a fraction of the
  aperture and should match `Luna.run`'s `rcollar`; it is ignored on the Cartesian grids,
  whose absorber profile is the grid's own window. `collar=nothing` skips the dataset.

Host-only, and the one expensive statistic in the free-space set: it applies the inverse
transverse transform to the whole state, which is a `N×N` matrix product on a radial grid
and an FFT over the transverse axes on a Cartesian one. That is a fraction of the two
such transforms each of the six right-hand-side evaluations of a step already does, but
it is not a reduction, and it needs one more host buffer the size of the state.
"""
function beam_profile(grid, spacegrid; collar=Boundaries.DEFAULT_RCOLLAR)
    BeamProfile(spacegrid, nothing, nothing, _widthaxis(spacegrid),
                _collarmask(spacegrid, collar), _measure(spacegrid))
end

device_capable(::BeamProfile) = false
needs_time(::BeamProfile) = false

_widthaxis(rg::Grid.RadialGrid) = Grid.rsymmetric(rg)
_widthaxis(sg::Grid.Free2DGrid) = sg.x
_widthaxis(sg::Grid.FreeGrid) = (sg.x, sg.y)

_measure(rg::Grid.RadialGrid) = rg.wr
_measure(::Grid.Free2DGrid) = nothing
_measure(::Grid.FreeGrid) = nothing

_collarmask(sg, ::Nothing) = nothing
_collarmask(sg, collar) = Boundaries.rprofile(sg, collar) .< 1

function prepare(f::BeamProfile, ctx::StatsContext)
    CT = eltype(ctx.proto)
    Er = Array{CT}(undef, size(ctx.proto))
    BeamProfile(f.sg, Er, _transverse_inverse(f.sg, Er), f.axis, f.mask, f.measure)
end

#= The inverse transverse transform, on the host. A radial grid's is the backward Hankel
   matrix applied along the last axis, held in the state's element type so that `mul!`
   reaches a BLAS `gemm` rather than the generic fallback a real matrix against a complex
   block would take; a Cartesian grid's is an inverse FFT over the transverse axes. =#
function _transverse_inverse(rg::Grid.RadialGrid, Er)
    Tb = convert(Matrix{eltype(Er)}, rg.Tbwd)
    (out, Eω) -> Grid.radial_matmul!(out, Eω, Tb)
end

_transverse_inverse(sg::Grid.Free2DGrid, Er) = _planned_ifft(Er, (3,))
_transverse_inverse(sg::Grid.FreeGrid, Er) = _planned_ifft(Er, (3, 4))

function _planned_ifft(buf, dims)
    Utils.loadFFTwisdom()
    iFT = FFTW.plan_ifft(buf, dims, flags=settings["fftw_flag"])
    Utils.saveFFTwisdom()
    (out, Eω) -> mul!(out, iFT, Eω)
end

function (f::BeamProfile)(d, Eω, Et, z, dz)
    f.xfrm!(f.Er, Eω)
    P = dropdims(sum(abs2, f.Er; dims=(1, 2)); dims=(1, 2))
    _beamwidth!(d, f, P, Eω)
    isnothing(f.mask) && return nothing
    tot = _transverse_sum(f.measure, P)
    d["collar_energy_fraction"] = _transverse_sum(f.measure, P, f.mask)/tot
    nothing
end

_transverse_sum(::Nothing, P) = sum(P)
_transverse_sum(w, P) = sum(w .* P)
_transverse_sum(::Nothing, P, mask) = sum(P[mask])
_transverse_sum(w, P, mask) = sum((w .* P)[mask])

#= The radial profile is sampled at `rg.r`, which does not include the axis, so the
   on-axis sample is taken from the state itself -- one more weighted reduction, the same
   one `Grid.symmetric` makes when it mirrors a real-space field. =#
function _beamwidth!(d, f::BeamProfile{<:Grid.RadialGrid}, P, Eω)
    E0 = Grid.onaxis(f.sg, Eω; dim=ndims(Eω))
    Psym = vcat(reverse(P), sum(abs2, E0), P)
    d["fwhm_r"] = Maths.fwhm(f.axis, Psym; method=:linear, minmax=:max)
    nothing
end

function _beamwidth!(d, f::BeamProfile{<:Grid.Free2DGrid}, P, Eω)
    d["fwhm_x"] = Maths.fwhm(f.axis, P; method=:linear, minmax=:max)
    nothing
end

function _beamwidth!(d, f::BeamProfile{<:Grid.FreeGrid}, P, Eω)
    x, y = f.axis
    ix, iy = argmin(abs.(x)), argmin(abs.(y))
    d["fwhm_x"] = Maths.fwhm(x, P[:, iy]; method=:linear, minmax=:max)
    d["fwhm_y"] = Maths.fwhm(y, P[ix, :]; method=:linear, minmax=:max)
    nothing
end

#= The statistics which read nothing but `z`: they never touch the state, so they run
   wherever it lives. =#

struct Density{D}
    dfun::D
end

"""
    density(dfun)

Create stats function to capture the gas density as defined by `dfun(z)`
"""
density(dfun) = Density(dfun)
device_capable(::Density) = true
needs_time(::Density) = false
(f::Density)(d, Eω, Et, z, dz) = d["density"] = f.dfun(z)

struct Pressure{D, G}
    dfun::D
    gas::G
end

"""
    pressure(dfun, gas)

Create stats function to capture the pressure. Like [`density`](@ref) but converts to
pressure. A `Tuple` of gases (a mixture) records one dataset per component.
"""
pressure(dfun, gas) = Pressure(dfun, gas)
device_capable(::Pressure) = true
needs_time(::Pressure) = false

(f::Pressure)(d, Eω, Et, z, dz) = d["pressure"] = PhysData.pressure(f.gas, f.dfun(z))

function (f::Pressure{D, <:Tuple})(d, Eω, Et, z, dz) where {D}
    dens = f.dfun(z)
    for (di, gi) in zip(dens, f.gas)
        d["pressure_$gi"] = PhysData.pressure(gi, di)
    end
end

struct CoreRadius{A}
    a::A
end

"""
    core_radius(a)

Create stats function to capture core radius as defined by `a` (either a `Number` or a
callable `a(z)`)
"""
core_radius(a) = CoreRadius(a)
device_capable(::CoreRadius) = true
needs_time(::CoreRadius) = false

(f::CoreRadius{<:Number})(d, Eω, Et, z, dz) = d["core_radius"] = f.a
(f::CoreRadius)(d, Eω, Et, z, dz) = d["core_radius"] = f.a(z)

mutable struct ZDW{M, L}
    mode_s::M
    λ00::L      # the last ZDW found, as the next root-finding guess
end

"""
    zdw(mode; λmin, λmax)

Create stats function to capture the zero-dispersion wavelength (ZDW) of one mode or of
each of several modes. The previous step's ZDW is the starting guess for the next.

!!! warning
    Since [`Modes.zdw`](@ref) is based on root-finding of a derivative, this can be slow!
"""
function zdw(mode::Modes.AbstractMode; λmin=100e-9, λmax=3000e-9)
    λ00 = Modes.zdw(mode; λmin=λmin, λmax=λmax, z=0)
    ZDW(mode, ismissing(λ00) ? λmin : λ00)
end

function zdw(modes; λmin=100e-9, λmax=3000e-9)
    λ00 = zeros(length(modes))
    for (ii, mode) in enumerate(modes)
        tmp = Modes.zdw(mode; λmin=λmin, λmax=λmax, z=0)
        λ00[ii] = ismissing(tmp) ? λmin : tmp
    end
    ZDW(modes, λ00)
end

device_capable(::ZDW) = true
needs_time(::ZDW) = false

function (f::ZDW{<:Modes.AbstractMode})(d, Eω, Et, z, dz)
    d["zdw"] = missnan(Modes.zdw(f.mode_s, f.λ00; z=z))
    f.λ00 = d["zdw"]
end

function (f::ZDW)(d, Eω, Et, z, dz)
    d["zdw"] = [missnan(Modes.zdw(f.mode_s[ii], f.λ00[ii]; z=z))
                for ii in eachindex(f.mode_s)]
    f.λ00 = d["zdw"]
end

"A ZDW which does not depend on `z`, which is what a constant linear operator implies."
struct ConstZDW{T}
    zdw::T
end

device_capable(::ConstZDW) = true
needs_time(::ConstZDW) = false
(f::ConstZDW)(d, Eω, Et, z, dz) = d["zdw"] = f.zdw

"""
    UserStat(f, label)

A statistics function the user supplied through `Stats.default`'s `userfuns`, wrapped so
that it has a name. It has no device form, so it is what makes the whole set run on the
host; the wrapper exists only so that the one-time warning can say `userfuns[1]` instead
of the gensym an anonymous closure's type carries.
"""
struct UserStat{F}
    f::F
    label::String
end

(u::UserStat)(d, Eω, Et, z, dz) = u.f(d, Eω, Et, z, dz)
statlabel(u::UserStat) = u.label

# convert missing to NaN
missnan(x) = ismissing(x) ? NaN : x

function zdz!(d, Eω, Et, z, dz)
    d["z"] = z
    d["dz"] = dz
end

device_capable(::typeof(zdz!)) = true
needs_time(::typeof(zdz!)) = false

#=================================================#
#===========  THE ANALYTIC TRANSFORM  ============#
#=================================================#

"""
    plan_analytic(grid, Eω)

Plan a transform from the frequency-domain field `Eω` to the analytic time-domain field.

Returns both a buffer for the analytic field and a closure to do the transform.

The buffers and the plan are of `Eω`'s own array type and precision, so this is one
inverse FFT on whatever the state lives on. On a device the result is in the state's own
(scaled) units; the statistics which read it apply `E_ref` to their scalar results.
"""
function plan_analytic(grid::Grid.EnvGrid, Eω)
    Eta = similar(Eω)
    iFT = _plan_ifft(Eω, copy(Eω))
    function analytic!(Eta, Eω)
        mul!(Eta, iFT, Eω) # for envelope fields, we only need to do the inverse transform
    end
    return Eta, analytic!
end

function plan_analytic(grid::Grid.RealGrid, Eω)
    s = collect(size(Eω))
    s[1] = (length(grid.ω) - 1)*2 # e.g. for 4097 rFFT samples, we need 8192 FFT samples
    Eta = similar(Eω, Tuple(s))
    Eωa = zero(Eta)
    iFT = _plan_ifft(Eω, Eωa)
    #= The rFFT has a sample at +fs/2 and the FFT does not, so the copy stops one short of
       the end of the rFFT axis; the rest of `Eωa` is zero and is never written. One
       masked broadcast rather than the scalar loop this used to be. =#
    n = size(Eω, 1) - 1
    hi = ntuple(_ -> Colon(), ndims(Eω) - 1)
    two = Luna.scalar(Eωa, 2.0)
    function analytic!(Eta, Eω)
        # copy across to the FFT-sampled buffer, doubling for the analytic signal
        view(Eωa, 1:n, hi...) .= two .* view(Eω, 1:n, hi...)
        mul!(Eta, iFT, Eωa) # now do the inverse transform
    end
    return Eta, analytic!
end

#= The inverse complex transform, planned on `buf` (which is of the state's array type).
   On the host this is `FFTW.plan_ifft` as it always was, with the wisdom logic and Luna's
   planning flags; on a device it is the generic `AbstractFFTs` planner, which is what
   `Utils.plan_ft`/`Utils.plan_ift` choose between. `inv` of a forward complex plan is the
   same `ScaledPlan` that `plan_ifft` returns, so the host arithmetic is unchanged. =#
function _plan_ifft(Eω, buf)
    Utils.isdevice(Eω) && return Utils.plan_ift(Utils.plan_ft(buf, 1))
    Utils.loadFFTwisdom()
    iFT = FFTW.plan_ifft(buf, 1, flags=settings["fftw_flag"])
    Utils.saveFFTwisdom()
    iFT
end

#=================================================#
#==============  THE COLLECTOR  ==================#
#=================================================#

struct StatsCollector{F, B, A}
    funcs::F
    Et::B
    analytic!::A
    needtime::Bool
    devicecapable::Bool
    hostlist::Vector{String}
end

"""
    STATS_DEVICE_MINLEN

The number of elements a state has to have before [`collect_stats`](@ref) evaluates the
statistics on the device rather than on a host copy of it. Below it the copy is cheaper.

The cost the device path adds is a fixed number of device-to-host round trips -- one per
statistic that ends in a scalar read, six for the default set -- and a round trip does not
depend on the size of the state. The cost it removes is one transfer of the state, which
does. Measured on an M1 Pro through Metal (`benchmark/stats.jl`), a round trip is ~400 µs
and the default set costs 2.7-3.6 ms on the device for 1025 to 16385 elements against
0.34-1.31 ms for the copy plus the host branches, and it is still the slower of the two on
a 128-column, 131200-element state (18.8 ms against 16.1 ms). The transfer reaches the
~2.7 ms the round trips cost at a few million elements, which is what this threshold is.

The rule is **the size of the state and nothing else**: an earlier version also took the
device path for any state with more than one column, which the column sweep in
`benchmark/stats.jl` does not support -- `fwhm_t` copies the time-domain intensity to the
host on either path and its per-column root-finding is host work either way, so extra
columns alone do not make the device path pay. `stats_device=:device` forces it regardless.

Re-measured on `gpu/int-E`, once the radial and free-space transforms could produce a
multi-column device state (M1 Pro, Metal 1.11, the device-capable part of the default set,
one call):

| state | elements | host + copy | device |
| --- | ---: | ---: | ---: |
| radial, 256 radial points | 33024 | 2.03 ms | 333 ms |
| radial, 1024 radial points | 132096 | 7.65 ms | 1.44 s |
| 3-D, 64 x 64 | 524288 | 662 ms | 615 ms |
| 3-D, 128 x 128 | 2097152 | 5.23 s | 4.88 s |

The threshold stays where it was: the device path does not pay at any of these sizes. On a
radial state it is two orders of magnitude worse, and the whole of that is [`fwhm_t`](@ref)
-- its device branch reduces the *scaled* field, so the columns far off axis
underflow to exactly zero in `Float32` and the host root-finding which follows is far
slower on them than on the small but non-zero numbers the host path gives it. On a
many-column state both paths are dominated by that same per-column host root-finding, so
`stats_period` (or `Output.nostats`) is the lever there rather than this switch.
"""
const STATS_DEVICE_MINLEN = 1 << 22

"Whether the device path is worth taking for a state this size; see `STATS_DEVICE_MINLEN`."
_devicepays(x) = length(x) >= STATS_DEVICE_MINLEN

"""
    collect_stats(grid, Eω, funcs...; Eref=1.0, stats_device=:auto)

Create a callable which collects statistics from the individual functions in `funcs`.

Each function given will be called with the arguments `(d, Eω, Et, z, dz)`, where
- d -> dictionary to store statistics values. each `func` should **mutate** this
- Eω -> frequency-domain field
- Et -> analytic time-domain field
- z -> current propagation distance
- dz -> current stepsize

`Eω` is a prototype of the propagating state: its array type, element type and shape are
what the buffers and the inverse transform are built for, and every statistic is
[`prepare`](@ref)d for it. `Eref` is the field unit the state is expressed in
([`Luna.UnitScaling`](@ref)), which is `1.0` for every `Float64` run and is applied by the
device branches only -- on the host the state is in physical units by the time the
statistics see it. [`default`](@ref) takes it from the transform.

The inverse transform which produces `Et` is planned, allocated and applied only when at
least one of `funcs` reads it ([`needs_time`](@ref)); `Et` is `nothing` otherwise.

# Which array the set is built for

Everything -- buffers, mirrors, the inverse plan, and each statistic's `ondevice` flag --
is built for one array, decided here and reported by
[`Luna.stats_device_capable`](@ref), which is what `Luna.ScaledOutput` uses to decide
whether to copy the state to the host. The device state is used when

- it is a device array, **and**
- every one of `funcs` has a device form ([`device_capable`](@ref)), **and**
- `stats_device` allows it.

`stats_device` is `:auto` (the default), `:device` or `:host`. Under `:auto` the device
state is used only when it has at least [`STATS_DEVICE_MINLEN`](@ref) elements; that
docstring has the measurement behind the threshold. A mode-averaged state is far below it,
and there the device path costs more than the copy it avoids: every statistic which ends
in a device-to-host transfer costs the same round trip whatever the size of the state
(~400 µs on an M1 Pro through Metal), and the default set makes six of them, against one
transfer of a few tens of kilobytes for the host path. `:device` overrides the size test
(the capability test still applies); `:host` builds for the host whatever the state is.

Otherwise everything is built for a host copy in physical units, with `Eref = 1`, because
that is what `Luna.ScaledOutput` will hand it.
"""
function collect_stats(grid, Eω, funcs...; Eref=1.0, stats_device=:auto)
    # make sure z and dz are recorded
    if !(zdz! in funcs)
        funcs = (funcs..., zdz!)
    end
    stats_device in (:auto, :device, :host) || throw(ArgumentError(
        "stats_device must be :auto, :device or :host, got :$stats_device"))
    ondev = Utils.isdevice(Eω) && stats_device !== :host &&
            (stats_device === :device || _devicepays(Eω))
    #= A set built for the host is called with a host copy in physical units, so it gets
       the host prototype and `Eref = 1`. Building it for the device instead would leave
       its buffers and its plan on the device while the field is not, which JLArrays
       tolerates silently and real hardware does not (review round 1, finding 1). =#
    ctx = StatsContext(grid, ondev ? Eω : Luna.tohost(Eω), ondev ? Eref : 1.0, ondev)
    prepped = map(f -> prepare(f, ctx), funcs)
    capable = all(device_capable, prepped)
    if ondev && !capable
        ondev = false
        ctx = StatsContext(grid, Luna.tohost(Eω), 1.0, false)
        prepped = map(f -> prepare(f, ctx), funcs)
    end
    hostlist = String[statlabel(f) for f in prepped if !device_capable(f)]
    _logstatspath(Eω, ondev, hostlist, stats_device)
    needtime = any(needs_time, prepped)
    #= The analytic transform of the whole state is the largest thing a statistics set
       allocates -- four times the state on a `RealGrid`, which for a 3-D free-space run
       is the whole memory budget -- and the free-space set does not read it at all: its
       on-axis statistics transform the projected `(nω, npol)` field instead. So it is
       planned only when something asks for it. =#
    Et, analytic! = needtime ? plan_analytic(grid, ctx.proto) : (nothing, nothing)
    StatsCollector(prepped, Et, analytic!, needtime, ondev, hostlist)
end

#= One line per propagation, and only when there is a choice to report: on a host run
   there is nothing to say. `Luna.ScaledOutput`'s one-time warning covers only the case
   where a statistic has no device form at all, so the other two reasons for the host
   path are said here. =#
function _logstatspath(Eω, ondev, hostlist, stats_device)
    Utils.isdevice(Eω) || return nothing
    if ondev
        Logging.@info("Per-step statistics run on the device.")
    elseif !isempty(hostlist)
        Logging.@info(
            "Per-step statistics run on the host: "*join(hostlist, ", ")*" "*
            (length(hostlist) == 1 ? "has" : "have")*" no device form, so the field is "*
            "copied from the device on every step the statistics fire.")
    elseif stats_device === :host
        Logging.@info("Per-step statistics run on the host (stats_device=:host).")
    else
        Logging.@info(
            "Per-step statistics run on the host: this state has $(length(Eω)) elements, "*
            "below Stats.STATS_DEVICE_MINLEN, where copying it down costs less than the "*
            "device reductions. Pass stats_device=:device to override.")
    end
    nothing
end

function (c::StatsCollector)(Eω, z, dz)
    d = Dict{String, Any}()
    c.needtime && c.analytic!(c.Et, Eω)
    for func in c.funcs
        func(d, Eω, c.Et, z, dz)
    end
    return d
end

"""
    device_capable(c::StatsCollector) -> Bool

For a whole set this is the *path* it was built for, not only what its members are capable
of: `false` when every statistic has a device form but [`collect_stats`](@ref) chose the
host anyway (a single small column, or `stats_device=:host`). `Luna.ScaledOutput` reads it
to decide whether to copy the state down, so it has to answer "which array will this set
be called with".
"""
device_capable(c::StatsCollector) = c.devicecapable
statlabel(::StatsCollector) = "statistics"

"""
    host_statistics(f) -> Vector{String}

The names of the statistics in `f` which have no device form, so cannot be evaluated on
the propagating state where it lives. Empty for a statistics set which runs entirely on a
device. `Luna.ScaledOutput` names them in the warning it emits when they force a
device-to-host copy of the state on every step.
"""
host_statistics(f) = device_capable(f) ? String[] : String[statlabel(f)]
host_statistics(c::StatsCollector) = c.hostlist

# The hooks `Luna.ScaledOutput` reaches the two traits above through; see `Device.jl`.
Luna.stats_device_capable(c::StatsCollector) = device_capable(c)
Luna.stats_host_list(c::StatsCollector) = host_statistics(c)

#=================================================#
#=============  THE DEFAULT SETS  ================#
#=================================================#

"""
    default(grid, Eω, mode_s, linop, transform; kwargs...)

The default statistics set for a mode-averaged (`mode_s::Modes.AbstractMode`) or
multimode (`mode_s::Modes.ModeCollection`) propagation.

`Eω` is the propagating state itself, not a host copy of it: the statistics are built for
its array type and precision, and the unit scaling they have to undo is taken from
`transform`.

`stats_device` (`:auto`, `:device` or `:host`) is passed to [`collect_stats`](@ref), which
documents what it chooses and why. `prop_capillary` does not take it directly; reach it
through `stats_kwargs`.

For a multimode propagation, `mode_error=true` (the default) adds the transverse-integral
diagnostic of whichever transform is in use: [`mode_reconstruction_error`](@ref) for the
adaptive cubature ([`NonlinearRHS.TransModal`](@ref Luna.NonlinearRHS.TransModal)) and
[`transverse_integral_error`](@ref) for the fixed quadrature rule
([`NonlinearRHS.TransModalFixed`](@ref Luna.NonlinearRHS.TransModalFixed)).
"""
function default(grid, Eω, mode::Modes.AbstractMode, linop, transform;
                 windows=nothing, gas=nothing, onaxis=false, userfuns=Any[],
                 stats_device=:auto)
    _, energyfunω = Fields.energyfuncs(grid)
    funs = [ω0(grid), energy(grid, energyfunω), peakpower(grid),
            fwhm_t(grid), zdw_linop(mode, linop),
            density(transform.densityfun)]
    if !isnothing(gas)
        push!(funs, pressure(transform.densityfun, gas))
    end
    if onaxis
        push!(funs, peakintensity(grid, (mode,)))
    else
        push!(funs, peakintensity(grid, transform.aeff))
    end
    for resp in transform.resp
        if resp isa PlasmaCumtrapz
            if onaxis
                push!(funs, electrondensity(grid, resp.ratefunc, transform.densityfun, (mode,)))
            else
                push!(funs, electrondensity(grid, resp.ratefunc, transform.densityfun, transform.aeff))
            end
        end
    end
    if !isnothing(windows)
        for win in windows
            push!(funs, energy_λ(grid, energyfunω, win))
        end
    end
    _adduserfuns!(funs, userfuns)
    collect_stats(grid, Eω, funs...;
                  Eref=Luna.runscaling(transform).Eref, stats_device)
end

@doc (@doc default)
function default(grid, Eω, modes::Modes.ModeCollection, linop, transform;
                 windows=nothing, gas=nothing, mode_error=true, userfuns=Any[],
                 stats_device=:auto)
    _, energyfunω = Fields.energyfuncs(grid)
    pol = transform.ts.indices == 1:2 ? :xy : transform.ts.indices == 1 ? :x : :y
    funs = [ω0(grid), energy(grid, energyfunω), peakpower(grid),
            peakintensity(grid, modes, components=pol), fwhm_t(grid),
            zdw_linop(modes, linop), density(transform.densityfun),
            fwhm_r(grid, modes; components=pol)]
    if !isnothing(gas)
        push!(funs, pressure(transform.densityfun, gas))
    end
    if mode_error
        push!(funs, _mode_error_stat(transform))
    end
    for resp in transform.resp
        if resp isa PlasmaCumtrapz
            ed = electrondensity(grid, resp.ratefunc, transform.densityfun, modes,
                                 components=pol)
            push!(funs, ed)
        end
    end
    if !isnothing(windows)
        for win in windows
            push!(funs, energy_λ(grid, energyfunω, win))
        end
    end
    _adduserfuns!(funs, userfuns)
    collect_stats(grid, Eω, funs...;
                  Eref=Luna.runscaling(transform).Eref, stats_device)
end

"""
    default(grid, Eω, transform, linop; kwargs...)

The default statistics set for a radial or free-space propagation, i.e. for a
[`NonlinearRHS.TransRadial`](@ref Luna.NonlinearRHS.TransRadial),
[`NonlinearRHS.TransFree`](@ref Luna.NonlinearRHS.TransFree) or
[`NonlinearRHS.TransFree2D`](@ref Luna.NonlinearRHS.TransFree2D).

It records

- `energy`, the total energy, from the transverse energy functional
  `Fields.energyfuncs(grid, spacegrid)[2]` -- one number per polarisation component;
- `peakintensity`, `ω0`, `fwhm_t_min` and `fwhm_t_max` **on the propagation axis** (see
  [`onaxis`](@ref)), because those are properties of the pulse and a transverse average of
  them is not one;
- `fwhm_r` (or `fwhm_x`/`fwhm_y` on a Cartesian transverse grid) and
  `collar_energy_fraction`, the beam size and how much of the beam has reached the
  absorbing collar (see [`beam_profile`](@ref));
- `density`, and `pressure` when `gas` is given;
- `electrondensity` and `peak_ionisation_rate` on axis, when the responses contain a
  `Nonlinear.PlasmaCumtrapz`;
- the energy in each wavelength region in `windows`;
- `z` and `dz`.

`Eω` is the propagating state itself, not a host copy of it, as for the modal methods.
`linop` is accepted for symmetry with those and is not used: a free-space linear operator
has no mode and so no zero-dispersion wavelength to record.

# Keyword arguments
- `windows=nothing`: wavelength regions to record the energy in, as for the modal methods
- `gas=nothing`: record the pressure of this gas (or `Tuple` of gases) as well as the
  density
- `userfuns=Any[]`: extra statistics functions, which are handed the whole state
- `collar=Boundaries.DEFAULT_RCOLLAR`: the width of the transverse absorber collar the
  energy fraction is measured in; should match `Luna.run`'s `rcollar`. `nothing` skips
  that dataset.
- `beam_profile=true`: record the beam size and the collar fraction at all. This is the
  one statistic in the set which has no device form, so it is also what decides whether
  the set can be evaluated on a device state; `false` makes the rest of the set
  device-capable.
- `stats_device=:auto`: as for the modal methods, see [`collect_stats`](@ref)
"""
function default(grid, Eω, transform::Union{TransRadial, TransFree, TransFree2D}, linop;
                 windows=nothing, gas=nothing, userfuns=Any[],
                 collar=Boundaries.DEFAULT_RCOLLAR, beam_profile=true,
                 stats_device=:auto)
    sg = _spacegrid(transform)
    _, energyfunω = Fields.energyfuncs(grid, sg)
    funs = Any[energy(grid, sg, energyfunω),
               onaxis(ω0(grid), sg),
               onaxis(peakintensity(grid), sg),
               onaxis(fwhm_t(grid), sg),
               density(transform.densityfun)]
    # `Stats.` because the keyword argument of the same name shadows the function here
    if beam_profile
        push!(funs, Stats.beam_profile(grid, sg; collar))
    end
    if !isnothing(gas)
        push!(funs, pressure(transform.densityfun, gas))
    end
    for resp in transform.resp
        if resp isa PlasmaCumtrapz
            push!(funs, onaxis(electrondensity(grid, resp.ratefunc,
                                               transform.densityfun), sg))
        end
    end
    if !isnothing(windows)
        for win in windows
            push!(funs, energy_λ(grid, sg, energyfunω, win))
        end
    end
    _adduserfuns!(funs, userfuns)
    collect_stats(grid, Eω, funs...;
                  Eref=Luna.runscaling(transform).Eref, stats_device)
end

"The transverse grid a radial or free-space transform holds."
_spacegrid(t::TransRadial) = t.rgrid
_spacegrid(t::TransFree) = t.xygrid
_spacegrid(t::TransFree2D) = t.xgrid

#= What `mode_error=true` means depends on which transverse integral the transform uses.
   The adaptive one can reconstruct the polarisation at a single transverse point and
   carries the cubature's own error estimate; the fixed rule has neither, but it has an
   embedded coarse rule, which is the same kind of diagnostic. Both record
   `transverse_points` and `transverse_integral_error_abs`/`_rel`; only the adaptive one
   records `mode_reconstruction_error`. =#
_mode_error_stat(transform::TransModal) = mode_reconstruction_error(transform)
_mode_error_stat(transform::TransModalFixed) = transverse_integral_error(transform)

#= Each user statistic is wrapped in a `UserStat` carrying `userfuns[i]` as its name, so
   that the one-time warning about the host copy names the thing the user wrote rather
   than the gensym an anonymous closure's type carries. The duplicate check compares the
   functions as given, not the wrappers. =#
function _adduserfuns!(funs, userfuns)
    seen = Any[]
    for (idx, uf) in enumerate(userfuns)
        if (uf in funs) || any(u -> u === uf, seen)
            @warn("userfun $idx is already present in the default set and will be ignored")
        else
            push!(seen, uf)
            push!(funs, UserStat(uf, "userfuns[$idx]"))
        end
    end
    funs
end

# For constant linop, ZDW is also constant
zdw_linop(mode::Modes.AbstractMode, linop::AbstractArray) =
    ConstZDW(missnan(Modes.zdw(mode)))

zdw_linop(modes, linop::AbstractArray) =
    ConstZDW([missnan(Modes.zdw(mode)) for mode in modes])

zdw_linop(mode_s, linop) = zdw(mode_s)

end
