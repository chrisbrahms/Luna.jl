module Stats
import Luna
import Luna: Maths, Grid, Modes, Utils, settings, PhysData, Fields, Processing, Ionisation
import Luna.PhysData: wlfreq, c, ε_0
import Luna.NonlinearRHS: TransModal, TransModeAvg, Erω_to_Prω!
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

@inline _wabs2(w, e) = w*abs2(e)
@inline _wmul(w, x) = w*x

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

"""
    _devenergy(grid, energyfun_ω)

`(weights, prefactor)` such that `prefactor*sum(weights .* abs2.(Eω))` is
`energyfun_ω(Eω)`, or `nothing` if there is no such pair for this grid or if
`energyfun_ω` is not the grid's own spectral energy functional.

The agreement is *checked* on a probe spectrum rather than assumed, because
[`energy`](@ref) and [`energy_window`](@ref) take the functional as an argument: a caller
who passes something else (`Fields.energyfuncs(grid, spacegrid)[2]`, or their own) gets a
host-only statistic rather than a silently different quantity.
"""
function _devenergy(grid, energyfun_ω)
    ws = _spectral_weights(grid)
    isnothing(ws) && return nothing
    w, prefac = ws
    Ep = _probespectrum(length(w))
    ref = try
        energyfun_ω(Ep)
    catch
        return nothing
    end
    (ref isa Real && isfinite(ref) && ref > 0) || return nothing
    abs(prefac*sum(w .* abs2.(Ep)) - ref) <= 1e-10*ref || return nothing
    (w, prefac)
end

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
            d["ω0"] = squeeze(Float64.(num) ./ Float64.(den))
        else
            num = _zipreduce(_wabs2, +, zero(RT), f.ωd, Eω)
            d["ω0"] = Float64(num)/Float64(sum(abs2, Eω))
        end
    else
        d["ω0"] = squeeze(Maths.moment(f.ω, abs2.(Eω); dim=1))
    end
end

squeeze(ω0::Array{T, 1}) where T = ω0[1]
squeeze(ω0::Array{T, 2}) where T = ω0[1, :]

struct SpectralEnergy{F, W}
    energyfun_ω::F
    w::W            # spectral weights on the state's array type, or `nothing`
    prefac::Float64
    Eref2::Float64  # E_ref^2: the energy is quadratic in the field
    key::String
    ondevice::Bool
end

"""
    energy(grid, energyfun_ω)

Create stats function to calculate the total energy.

On a device the integral is the same rule as `energyfun_ω`'s written as one weighted
reduction over a lazy `Broadcasted` (see `Stats._devenergy`); the scalar prefactor and the
square of the unit scaling are applied to the result, in `Float64`.
"""
energy(grid, energyfun_ω) = SpectralEnergy(energyfun_ω, nothing, 1.0, 1.0, "energy", false)

function prepare(f::SpectralEnergy, ctx::StatsContext)
    we = _devenergy(ctx.grid, f.energyfun_ω)
    isnothing(we) && return SpectralEnergy(f.energyfun_ω, nothing, 1.0, 1.0, f.key, false)
    w, prefac = we
    SpectralEnergy(f.energyfun_ω, Luna.todevice(ctx.spec, w), prefac, ctx.Eref^2, f.key,
                   ctx.ondevice)
end

device_capable(f::SpectralEnergy) = !isnothing(f.w)
needs_time(::SpectralEnergy) = false
statlabel(f::SpectralEnergy) = f.key

function (f::SpectralEnergy)(d, Eω, Et, z, dz)
    if _onstate(f, Eω)
        d[f.key] = _weightedenergy(f.w, f.prefac*f.Eref2, Eω)
    elseif ndims(Eω) > 1
        d[f.key] = [f.energyfun_ω(Eω[:, i]) for i=1:size(Eω, 2)]
    else
        d[f.key] = f.energyfun_ω(Eω)
    end
end

#= One reduction along the frequency axis over a lazy `Broadcasted`, so nothing
   field-sized is materialised. The prefactors are applied on the host, in Float64, to the
   result: in a scaled Float32 run the field is deliberately of order one and the physical
   units (up to 1e20 apart) would not survive being put back on it. =#
function _weightedenergy(w, fac, Eω)
    isnothing(w) && error(
        "this energy statistic was built from an energy functional with no device form, "*
        "so it cannot be evaluated on a $(typeof(Eω)).")
    RT = real(eltype(Eω))
    if ndims(Eω) > 1
        s = Array(_zipreduce1(_wabs2, +, zero(RT), w, Eω))
        return [fac*Float64(s[1, i]) for i in axes(s, 2)]
    end
    fac*Float64(_zipreduce(_wabs2, +, zero(RT), w, Eω))
end

"""
    energy_λ(grid, energyfun_ω, λlims; label)

Create stats function to calculate the energy in a wavelength region given by `λlims`.
If `label` is omitted, the stats dataset is named by the wavelength limits.
"""
function energy_λ(grid, energyfun_ω, λlims; label=nothing, winwidth=0)
    λlims = collect(λlims)
    ωmin, ωmax = extrema(wlfreq.(λlims))
    window = Maths.planck_taper(grid.ω, ωmin-winwidth, ωmin, ωmax, ωmax+winwidth)
    if isnothing(label)
        λnm = 1e9.*λlims
        label = @sprintf("%.2fnm_%.2fnm", minimum(λnm), maximum(λnm))
    end
    energy_window(grid, energyfun_ω, window; label=label)
end

struct SpectralEnergyWindow{F, V, W}
    energyfun_ω::F
    window::V       # the host window, as given
    w::W            # spectral weights times window^2, on the state's array type
    prefac::Float64
    Eref2::Float64
    key::String
    ondevice::Bool
end

"""
    energy_window(grid, energyfun_ω, window; label)

Create stats function to calculate the energy filtered by a `window`. The stats dataset will
be named `energy_[label]`.
"""
energy_window(grid, energyfun_ω, window::AbstractVector{<:Real}; label) =
    SpectralEnergyWindow(energyfun_ω, window, nothing, 1.0, 1.0, "energy_$label", false)

function prepare(f::SpectralEnergyWindow, ctx::StatsContext)
    we = _devenergy(ctx.grid, f.energyfun_ω)
    isnothing(we) && return SpectralEnergyWindow(f.energyfun_ω, f.window, nothing,
                                                 1.0, 1.0, f.key, false)
    w, prefac = we
    #= The window multiplies the field, so it enters the reduction squared and folds into
       the weights: one vector and one reduction, no windowed copy of the field. =#
    SpectralEnergyWindow(f.energyfun_ω, f.window,
                         Luna.todevice(ctx.spec, w .* abs2.(f.window)),
                         prefac, ctx.Eref^2, f.key, ctx.ondevice)
end

device_capable(f::SpectralEnergyWindow) = !isnothing(f.w)
needs_time(::SpectralEnergyWindow) = false
statlabel(f::SpectralEnergyWindow) = f.key

function (f::SpectralEnergyWindow)(d, Eω, Et, z, dz)
    if _onstate(f, Eω)
        d[f.key] = _weightedenergy(f.w, f.prefac*f.Eref2, Eω)
    elseif ndims(Eω) > 1
        d[f.key] = [f.energyfun_ω(Eω[:, i].*f.window) for i=1:size(Eω, 2)]
    else
        d[f.key] = f.energyfun_ω(Eω.*f.window)
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

The number of elements a *single-column* state has to have before [`collect_stats`](@ref)
evaluates the statistics on the device rather than on a host copy of it. Below it the copy
is cheaper.

The cost the device path adds is a fixed number of device-to-host round trips -- one per
statistic that ends in a scalar read, six for the default set -- and a round trip does not
depend on the grid size. The cost it removes is one transfer of the state, which does.
Measured on an M1 Pro through Metal (`benchmark/stats.jl`), a round trip is ~400 µs and
the default set costs 2.7-3.6 ms on the device for 1025 to 16385 elements, against
0.34-1.31 ms for the copy plus the host branches. The transfer only reaches the ~2.7 ms
the round trips cost at a few million elements, which is what this threshold is; a
single-column state that large is not something Luna produces, so in practice a
mode-averaged run always takes the host path. `stats_device=:device` overrides it.

A state with **more than one column** takes the device path whatever its length. That is a
forward-looking rule, not one this branch can measure: no transform here produces a
multi-column device state (the radial, free-space and multimode transforms are host-only
until `gpu/20`-`gpu/22`). On a synthetic multi-column state the device path is still the
slower of the two at 16 and 128 columns, because `fwhm_t` copies the time-domain intensity
to the host on either path and its per-column root-finding is host work either way. See
`PR_24-stats-device.md`; it is worth re-measuring when a device-capable multi-column
transform exists.
"""
const STATS_DEVICE_MINLEN = 1 << 22

"The number of columns (transverse points, modes) of a state."
_ncols(x::AbstractArray) = ndims(x) < 2 ? 1 : prod(size(x)[2:end])

"Whether the device path is worth taking for a state of this shape; see `STATS_DEVICE_MINLEN`."
_devicepays(x) = _ncols(x) > 1 || length(x) >= STATS_DEVICE_MINLEN

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

The inverse transform which produces `Et` is *applied* once per call, and only when at
least one of `funcs` reads it ([`needs_time`](@ref)); it is planned and its buffers are
allocated either way.

# Which array the set is built for

Everything -- buffers, mirrors, the inverse plan, and each statistic's `ondevice` flag --
is built for one array, decided here and reported by
[`Luna.stats_device_capable`](@ref), which is what `Luna.ScaledOutput` uses to decide
whether to copy the state to the host. The device state is used when

- it is a device array, **and**
- every one of `funcs` has a device form ([`device_capable`](@ref)), **and**
- `stats_device` allows it.

`stats_device` is `:auto` (the default), `:device` or `:host`. Under `:auto` the device
state is used only when it has more than one column or at least
[`STATS_DEVICE_MINLEN`](@ref) elements; that docstring has the measurement behind the
rule. A single-column mode-averaged state is below both, and there the device path costs
more than the copy it avoids: every statistic which ends in a device-to-host transfer
costs the same round trip whatever the grid size (~400 µs on an M1 Pro through Metal), and
the default set makes six of them, against one transfer of a few tens of kilobytes for the
host path. `:device` overrides the shape test (the capability test still applies); `:host`
builds for the host whatever the state is.

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
    Et, analytic! = plan_analytic(grid, ctx.proto)
    StatsCollector(prepped, Et, analytic!, any(needs_time, prepped), ondev, hostlist)
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
            "Per-step statistics run on the host: this state is a single column of "*
            "$(length(Eω)) elements, below Stats.STATS_DEVICE_MINLEN, where copying it "*
            "down costs less than the device reductions. Pass stats_device=:device to "*
            "override.")
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
        push!(funs, mode_reconstruction_error(transform))
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
