module Interface
using Luna
import Luna.PhysData: wlfreq, roomtemp
import Luna: Grid, Modes, Output, Fields, Boundaries, NonlinearRHS
import Random: AbstractRNG, GLOBAL_RNG
import Logging: @info, @debug

module Pulses

import Luna: Fields, Output, Processing, Capillary

export AbstractPulse, CustomPulse, GaussPulse, SechPulse, DataPulse, LunaPulse

abstract type AbstractPulse end

struct CustomPulse{fT<:Fields.TimeField} <: AbstractPulse
    mode::Symbol
    polarisation
    field::fT
end

"""
    CustomPulse(;λ0, energy=nothing, power=nothing, ϕ=Float64[],
                mode=:lowest, polarisation=:linear, propagator=nothing)

A custom pulse defined by a function for use with `prop_capillary`, with either energy or
peak power specified.

# Keyword arguments
- `λ0::Number`: the central wavelength
- `Itshape::function`: a function `I(t)`` which defines the intensity/power envelope of the
                       pulse as a function of time `t`. Note that the normalisation of this
                       envelope is irrelevant as it will be re-scaled by `energy` or `power`.
- `energy::Number`: the pulse energy.
- `power::Number`: the pulse peak power (**after** applying any spectral phases).
- `ϕ::Vector{Number}`: spectral phases (CEP, group delay, GDD, TOD, ...).
- `mode::Symbol`: Mode in which this input should be coupled. Can be `:lowest` for the
                  lowest-order mode in the simulation, or a mode designation
                  (e.g. `:HE11`, `:HE12`, `:TM01`, etc.). Defaults to `:lowest`.
- `polarisation`: Can be `:linear`, `:x`, `:y`, `:circular`, or an ellipticity number -1 ≤ ε ≤ 1,
                  where ε=-1 corresponds to left-hand circular, ε=1 to right-hand circular,
                  and ε=0 to linear polarisation.
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
"""
function CustomPulse(;mode=:lowest, polarisation=:linear, propagator=nothing, kwargs...)
    CustomPulse(mode, polarisation,
                Fields.PropagatedField(propagator, Fields.PulseField(;kwargs...)))
end

struct GaussPulse{fT<:Fields.TimeField} <: AbstractPulse
    mode::Symbol
    polarisation
    field::fT
end

"""
    GaussPulse(;λ0, τfwhm, energy=nothing, power=nothing, ϕ=Float64[], m=1,
               mode=:lowest, polarisation=:linear, propagator=nothing)

A (super)Gaussian pulse for use with `prop_capillary`, with either energy or peak power
specified.

# Keyword arguments
- `λ0::Number`: the central wavelength.
- `τfwhm::Number`: the pulse duration (power/intensity FWHM).
- `energy::Number`: the pulse energy.
- `power::Number`: the pulse peak power (**after** applying any spectral phases).
- `ϕ::Vector{Number}`: spectral phases (CEP, group delay, GDD, TOD, ...).
- `m::Int`: super-Gaussian parameter (the power in the Gaussian exponent is 2m).
            Defaults to 1.
- `mode::Symbol`: Mode in which this input should be coupled. Can be `:lowest` for the
                  lowest-order mode in the simulation, or a mode designation
                  (e.g. `:HE11`, `:HE12`, `:TM01`, etc.). Defaults to `:lowest`.
- `polarisation`: Can be `:linear`, `:x`, `:y`, `:circular`, or an ellipticity number -1 ≤ ε ≤ 1,
                  where ε=-1 corresponds to left-hand circular, ε=1 to right-hand circular,
                  and ε=0 to linear polarisation.
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
"""
function GaussPulse(;mode=:lowest, polarisation=:linear, propagator=nothing, kwargs...)
    GaussPulse(mode, polarisation,
               Fields.PropagatedField(propagator, Fields.GaussField(;kwargs...)))
end

struct SechPulse{fT<:Fields.TimeField} <: AbstractPulse
    mode::Symbol
    polarisation
    field::fT
end

"""
    SechPulse(;λ0, τfwhm=nothing, τw=nothing, energy=nothing, power=nothing, ϕ=Float64[],
               mode=:lowest, polarisation=:linear, propagator=nothing)

A sech²(τ/τw) pulse for use with `prop_capillary`, with either `energy` or peak `power`
specified, and duration given either as `τfwhm` or `τw`.

# Keyword arguments
- `λ0::Number`: the central wavelength.
- `τfwhm::Number`: the pulse duration (power/intensity FWHM).
- `τw::Number`: "natural" pulse duration of a sech²(τ/τw) pulse.
- `energy::Number`: the pulse energy.
- `power::Number`: the pulse peak power (**after** applying any spectral phases).
- `ϕ::Vector{Number}`: spectral phases (CEP, group delay, GDD, TOD, ...)
- `mode::Symbol`: Mode in which this input should be coupled. Can be `:lowest` for the
                  lowest-order mode in the simulation, or a mode designation
                  (e.g. `:HE11`, `:HE12`, `:TM01`, etc.). Defaults to `:lowest`.
- `polarisation`: Can be `:linear`, `:x`, `:y`, `:circular`, or an ellipticity number -1 ≤ ε ≤ 1,
                  where ε=-1 corresponds to left-hand circular, ε=1 to right-hand circular,
                  and ε=0 to linear polarisation.
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
"""
function SechPulse(;mode=:lowest, polarisation=:linear, propagator=nothing, kwargs...)
    SechPulse(mode, polarisation,
              Fields.PropagatedField(propagator, Fields.SechField(;kwargs...)))
end

struct DataPulse{fT<:Fields.TimeField} <: AbstractPulse
    mode::Symbol
    polarisation
    field::fT
end

#TODO add peak power to DataPulses
"""
    DataPulse(ω, Iω, ϕω; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)
    DataPulse(ω, Eω; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)
    DataPulse(fpath; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)

A custom pulse defined by tabulated data to be used with `prop_capillary`.

# Data input options
- `ω, Iω, ϕω`: arrays of angular frequency `ω` (units rad/s), spectral energy density `Iω`
               and spectral phase `ϕω`. `ϕω` should be unwrapped.
- `ω, Eω`: arrays of angular frequency `ω` (units rad/s) and the complex frequency-domain
           field `Eω`.
- `fpath`: a string containing the path to a file which contains 3 columns:
    Column 1: frequency (units of Hertz)
    Column 2: spectral energy density
    Column 3: spectral phase (unwrapped)

# Keyword arguments
- `energy::Number`: the pulse energy
- `λ0::Number`: the central wavelength (optional; defaults to the centre of mass of the
                given spectral energy density).
- `ϕ::Vector{Number}`: spectral phases (CEP, group delay, GDD, TOD, ...) to be applied to the
                       pulse (in addition to any phase already present in the data).
- `mode::Symbol`: Mode in which this input should be coupled. Can be `:lowest` for the
                  lowest-order mode in the simulation, or a mode designation
                  (e.g. `:HE11`, `:HE12`, `:TM01`, etc.). Defaults to `:lowest`.
- `polarisation`: Can be `:linear`, `:x`, `:y`, `:circular`, or an ellipticity number -1 ≤ ε ≤ 1,
                  where ε=-1 corresponds to left-hand circular, ε=1 to right-hand circular,
                  and ε=0 to linear polarisation.
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
"""
function DataPulse(ω::AbstractVector, Iω, ϕω;
                   mode=:lowest, polarisation=:linear, propagator=nothing, kwargs...)
    DataPulse(mode, polarisation,
              Fields.PropagatedField(propagator, Fields.DataField(ω, Iω, ϕω; kwargs...)))
end

function DataPulse(ω, Eω;
                   mode=:lowest, polarisation=:linear, propagator=nothing, kwargs...)
    DataPulse(mode, polarisation,
              Fields.PropagatedField(propagator, Fields.DataField(ω, Eω; kwargs...)))
end

function DataPulse(fpath;
                   mode=:lowest, polarisation=:linear, propagator=nothing, kwargs...)
    DataPulse(mode, polarisation,
              Fields.PropagatedField(propagator, Fields.DataField(fpath; kwargs...)))
end

"""
    LunaPulse(output; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)

A pulse defined to be used with `prop_capillary` which comes from a previous `Luna`
propagation simulation.

For multi-mode simulations, only the lowest-order modes is transferred.

# Arguments
- `output::AbstractOutput`: output from a previous `Luna` simulation.

# Keyword arguments
- `energy::Number`: the pulse energy. When transferring multi-mode simulations this defines the **total** energy.
- `scale_energy`: if given instead of `energy`, scale the field from `output` by this number. Defaults to 1, so giving `energy` is **not** required. For multi-mode simulations, this can also be a `Vector` with the same number of elements as the number of modes, in which case the energy of each mode is scaled by the corresponding number.
- `λ0::Number`: the central wavelength (optional; defaults to the centre of mass of the
                given spectral energy density).
- `ϕ::Vector{Number}`: spectral phases (CEP, group delay, GDD, TOD, ...) to be applied to the
                       pulse (in addition to any phase already present in the data).
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
"""
function LunaPulse(o::Output.AbstractOutput; energy=nothing, scale_energy=nothing, kwargs...)
    ω = o["grid"]["ω"]
    t = o["grid"]["t"]
    τ = length(t) * (t[2] - t[1])/2 # middle of old time window
    Eω = o["Eω"]
    if ndims(Eω) == 2
        # mode-averaged
        Eωm = Eω[:, end]
        eout = Processing.energy(o)[end]
        e = make_energies(energy, scale_energy, eout)
        return DataPulse(ω, Eωm .* exp.(1im .* ω .* τ); energy=e, kwargs...)
    elseif ndims(Eω) == 3
        # multi-mode
        modes = Processing.makemodes(o; warn_dispersion=false)
        symbols = makesymbol.(modes)
        eout = Processing.energy(o)[:, end]
        es = make_energies(energy, scale_energy, eout)
        return [DataPulse(ω, Eω[:, ii, end] .* exp.(1im .* ω .* τ); mode=symbols[ii], energy=es[ii], kwargs...) for ii in eachindex(modes)]
    end
end

makesymbol(mode::Capillary.MarcatiliMode) = Symbol("$(mode.kind)$(mode.n)$(mode.m)")

make_energies(energy::Number, scale_energy::Nothing, eout) = eout ./ sum(eout) .* energy
make_energies(energy::Nothing, scale_energy, eout) = eout .* scale_energy
make_energies(energy::Nothing, scale_energy::Nothing, eout) = eout

struct GaussBeamPulse{pT, NmT} <: AbstractPulse
    waist::Float64
    timepulse::pT
    polarisation
    Nmodes::NmT
end

"""
    GaussBeamPulse(waist, timepulse, Nmodes=:all)

A pulse whose shape in time is defined by the `timepulse::AbstractPulse`, and whose modal content is calculated by considering the overlap of an ideal Gaussian laser beam with 1/e² radius `waist` with the modes of the waveguide. `Nmodes` determines how many of the available modes to couple to. By default (`Nmodes=:all`) all modes are taken into account, but this can lead to numerical inaccuracies.
"""
function GaussBeamPulse(waist, timepulse, Nmodes=:all)
    GaussBeamPulse(waist, timepulse, timepulse.polarisation, Nmodes)
end

end


"""
    prop_capillary(radius, flength, gas, pressure; λ0, λlims, trange, kwargs...)

Simulate pulse propagation in a hollow fibre using the capillary model.

# Mandatory arguments
- `radius`: Core radius of the fibre. Can be a `Number` for constant radius, or a function
    `a(z)` which returns the `z`-dependent radius.
- `flength::Number`: Length of the fibre.
- `gas::Symbol`: Filling gas species.
- `pressure`: Gas pressure. Can be a `Number` for constant pressure, a 2-`Tuple` of `Number`s
    for a simple pressure gradient, or a `Tuple` of `(Z, P)` where `Z` and `P`
    contain `z` positions and the pressures at those positions.
- `λ0`: (keyword argument) the reference wavelength for the simulation. For simple
    single-pulse inputs, this is also the central wavelength of the input pulse.
- `λlims::Tuple{<:Number, <:Number}`: The wavelength limits for the simulation grid.
- `trange::Number`: The total width of the time grid. To make the number of samples a
    power of 2, the actual grid used is usually bigger.

# Grid options
- `envelope::Bool`: Whether to use envelope fields for the simulation. Defaults to `false`.
    By default, envelope simulations ignore third-harmonic generation.
    Plasma has not yet been implemented for envelope fields.
- `δt::Number`: Time step on the fine grid used for the nonlinear interaction. By default,
    this is determined by the wavelength grid. If `δt` is given **and smaller** than the
    required value, it is used instead.

# Input pulse options
A single pulse in the lowest-order mode can be specified by the keyword arguments below.
More complex inputs can be defined by a single `AbstractPulse` or a `Vector{AbstractPulse}`.
In this case, all keyword arguments except for `λ0` are ignored.

- `λ0`: Central wavelength
- `τfwhm`: The pulse duration as defined by the full width at half maximum.
- `τw`: The "natural" pulse duration. Only available if pulseshape is `sech`.
- `ϕ`: Spectral phases to be applied to the transform-limited pulse. Elements are
    the usual polynomial phases ϕ₀ (CEP), ϕ₁ (group delay), ϕ₂ (GDD), ϕ₃ (TOD), etc.
- `energy`: Pulse energy.
- `power`: Peak power **after any spectral phases are added**.
- `pulseshape`: Shape of the transform-limited pulse. Can be `:gauss` for a Gaussian pulse
    or `:sech` for a sech² pulse.
- `polarisation`: Polarisation of the input pulse. Can be `:linear` (default), `:x`, `:y`,
    `:circular`, or an ellipticity number -1 ≤ ε ≤ 1, where ε=-1 corresponds to left-hand circular,
    ε=1 to right-hand circular, and ε=0 to linear polarisation. The major axis for
    elliptical polarisation is always the y-axis.
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
- `shotnoise`: Whether and how to include quantum noise. Can be one of:
    - `true` (default) -- same as `:modified`.
    - `false` -- disable all noise.
    - `:modified` -- use the modified shot-noise model of Chen & Wise
      (arXiv:2410.20567), where a constant noise field enters the nonlinear operator
      at every step but is excluded from dispersion. This prevents artificial FWM
      phase-matching and elevated noise floor artefacts.
    - `:input` -- use traditional one-photon-per-mode shot noise added to the input
      field at `z = 0`.
    See the [Noise model](@ref) documentation for details.
- `rng`: Random number generator for noise field generation. Defaults to `GLOBAL_RNG`.
    Pass a seeded RNG (e.g. `MersenneTwister(seed)`) for reproducible noise realisations,
    or different seeds for ensemble/shot-to-shot statistics.

# Modes options
- `modes`: Defines which modes are included in the propagation. Can be any of:
    - a single mode signifier (default: :HE11), which leads to mode-averaged propagation
        (as long as all inputs are linearly polarised).
    - a `Dict` mode signifier with keys `:kind`, `:n`, and `:m`, e.g. `Dict(:kind=>:HE, :n=>1, :m=>1)`
    - a `Tuple` of mode signifiers (`Symbol`s or `Dict`s), which leads to multi-mode propagation in those modes.
    - a `Number` `N` of modes, which simply creates the first `N` `HE` modes.
    Note that when elliptical or circular polarisation is included, each mode is present
    twice in the output, once for `x` and once for `y` polarisation.
- `modal_integral::Symbol`: how the transverse integral of the nonlinear polarisation is
    evaluated in a multimode simulation. `:adaptive` (the default) uses an adaptive
    cubature rule, which chooses its own transverse points to reach
    `radial_integral_rtol` and runs on the host in double precision. `:fixed` uses a fixed
    Gauss quadrature rule of `nr` (and, for the full 2-D integral, `nθ`) nodes, which
    costs the same on every step, is the multimode transform which runs on a device or in
    reduced precision, and is a different discretisation of the same integral -- it agrees
    with `:adaptive` to the accuracy of the rule rather than to rounding. Ignored for
    mode-averaged propagation, which has no transverse integral. See
    [`NonlinearRHS.TransModalFixed`](@ref Luna.NonlinearRHS.TransModalFixed).
- `radial_integral_rtol::Number`: relative tolerance of the adaptive transverse integral
    (`modal_integral=:adaptive` only).
- `modal_nr::Int`, `modal_nθ::Int`, `modal_kronrod::Bool`: the quadrature rule of
    `modal_integral=:fixed`: the number of nodes along r (or x) and θ (or y), and whether
    to use a Gauss-Kronrod rule in r so that the rule carries an embedded error estimate.
    `modal_nr` has to resolve the transverse structure of the highest mode, and
    `modal_nθ` has to be at least `4h+1` for modes of azimuthal order up to `h` (so at
    least 5 for an HE₁ₘ set). They are `nr`, `nθ` and `kronrod` on [`Luna.setup`](@ref),
    where the context leaves no room for confusion.
- `model::Symbol`: Can be `:full`, which includes the full complex refractive index of the cladding
    in the effective index of the mode, or `:reduced`, which uses the simpler model more
    commonly seen in the literature. See `Luna.Capillary` for more details.
    Defaults to `:full`.
- `loss::Bool`: Whether to include propagation loss. Defaults to `true`.
- `temperature::Number`: Temperature of the gas in Kelvin. Defaults to room temperature.

# Nonlinear interaction options
- `kerr`: Whether to include the Kerr effect. Defaults to `true`.
- `raman`: Whether to include the Raman effect. Defaults to `false`.
- `plasma`: Can be one of
    - `:ADK` -- include plasma using the ADK ionisation rate.
    - `:PPT` -- include plasma using the PPT ionisation rate.
    - `true` (default) -- same as `:PPT`.
    - `false` -- ignore plasma.
    Note that plasma is only available for full-field simulations.
- `PPT_options::Dict{Symbol, Any}`: when using the PPT ionisation rate for the
    plasma nonlinearity, this allows for fine-tuning of the options in calculating
    the ionisation. See [`IonRatePPTAccel`](@ref Ionisation.IonRatePPTAccel) for possible
    keyword arguments.
- `preionfrac::Float64`: fraction of the gas that is pre-ionised before the pulse. Defaults to `0.0`.
    Note that this is a very simplistic model of pre-ionisation and should be used with
    caution.
- `thg::Bool`: Whether to include third-harmonic generation. Defaults to `true` for
    full-field simulations and to `false` for envelope simulations.
If `raman` is `true`, then the following options apply:
    - `rotation::Bool = true`: whether to include the rotational Raman contribution
    - `vibration::Bool = true`: whether to include the vibrational Raman contribution

# Output options
- `stats_kwargs::Dict{Symbol, Any}`: a dictionary of keyword arguments to `Stats.default`
- `saveN::Integer`: Number of points along z at which to save the field.
- `filepath`: If `nothing` (default), create a `MemoryOutput` to store the simulation results
    only in the working memory. If not `nothing`, should be a file path as a `String`,
    and the results are saved in a file at this location. If `scan` is passed, `filepath`
    determines the output **directory** for the scan instead.
- `scan`: A `Scan` instance defining a parameter scan. If `scan` is given`, a
    `Output.ScanHDF5Output` is used to automatically name and populate output files of
    the scan. `scanidx` must also be given.
- `scanidx`: Current scan index within a scan being run. Only used when `scan` is passed.
- `filename`: Can be used to to overwrite the scan name when running a parameter scan.
    The running `scanidx` will be appended to this filename. Ignored if no `scan` is given.
- `status_period::Number`: Interval (in seconds) between printed status updates.
- `boundary::Symbol=:rate`: How the absorbing boundaries at the edges of the frequency and
    time windows are applied. `:rate` treats them as an absorption rate per unit distance,
    so the total absorption depends only on the propagation distance and not on how many
    steps the solver took. `:legacy` reproduces the historical behaviour, in which the
    windows were applied once per accepted step and the result therefore depended on
    `rtol`. `:none` disables them. See [`Luna.run`](@ref).
- `boundary_N::Real`: Absorber strength, expressed as the number of times the historical
    window profile is applied over the whole propagation length.
- `boundary_length`: Absorber reference length in metres, overriding `boundary_N`.
- `tcollar::Real`: Minimum width of the temporal absorber collar, as a fraction of the time
    window.
- `linop_integral::Symbol=:auto`: how the integral `Φ(z) = ∫linop dz'` of a tapered
    or pressure-graded capillary's z-dependent linear operator is obtained. The stepper
    propagates the linear part by `exp(Φ(t2) − Φ(t1))`, which is exact.
    `:tabulated` builds a table of `Φ` at setup, along with the propagation constant and
    effective area, so that the propagation does no host work per stage -- which is what
    makes such a run go entirely on a device; `:quadrature` integrates the operator over
    each step instead, with no table and no setup pass but around fifteen host
    evaluations of the operator per stage; `:auto` (the default) tabulates unless the
    table would not fit in a byte budget, in which case it falls back to the quadrature and
    says so. See [`Luna.run`](@ref). A mode-averaged capillary always tabulates.

    A uniform fibre has a constant operator, which is already exact in the propagator, and
    is not affected by this keyword at all. `prop_gnlse` does not accept the keyword: its
    operator is always constant.
- `linop_tol::Real`: the tolerance `Φ` is computed to, in radians. See [`Luna.run`](@ref).
- `tabulate_linop`: **deprecated** and ignored; see [`Luna.run`](@ref).
- `device`: where to run: `:cpu`, `:auto`, `:metal`, `:cuda` or a [`Luna.DeviceSpec`](@ref).
    `nothing` (the default) means "not specified": it becomes `Luna.device_request()`,
    i.e. `Luna.settings["device"]` as the user set it (`:cpu` if nothing was set and
    nothing loaded, `:auto` once a GPU package has been `using`d), only when the
    propagation has a device path -- mode-averaged (`modes` a single mode), or multimode
    with `modal_integral=:fixed` -- *and* every nonlinear response it
    was built with is device-capable (`Nonlinear.device_capable`; true for the Kerr
    responses including the no-THG one, the plasma, Raman and χ⁽²⁾ responses, false for
    anything user-written). Otherwise it stays on
    the CPU, whatever `Luna.settings["device"]` says, exactly as before this keyword
    existed -- loading a GPU package must never turn a working default call into an
    error. An *explicit* `device` or `precision` request which cannot be honoured
    (multimode propagation with `modal_integral=:adaptive`, a response with no device
    kernel, and [`prop_gnlse`](@ref)) errors naming the fix (`device=:cpu` or
    `Luna.set_device(:cpu)`), rather than being silently narrowed to the CPU, run on the
    host through `Nonlinear.HostResponse` at every step, or failing with an unrelated
    `MethodError`. See the "Running on a GPU" page (`docs/src/gpu.md`).
- `precision`: `Float32` to run in reduced precision (on the CPU or on a device),
    `Float64` for double, `nothing` (default) for whatever `device` resolves to
    normally (`Float64` on the CPU, `Float32` on Metal). A `Float32` run is scaled (see
    [`Luna.UnitScaling`](@ref)); the output is unscaled automatically and saved as
    `ComplexF32`. `precision` alone does not select the CPU: if a GPU package is loaded
    and `device` is left at its default, `precision=Float32` runs on the GPU (whatever
    `Luna.settings["device"]` resolves to) in that precision, not on the CPU. Pass
    `device=:cpu` (or call `Luna.set_device(:cpu)`) to force the host.
- `stats_period::Real=1`: collect the default statistics less often than every accepted
    step (see [`Output.PeriodicStats`](@ref) and [`Output.maybe_periodic`](@ref)). An
    integer (the default, `1`) collects every `stats_period`-th accepted step; a
    non-integer value collects every time the propagation distance has advanced by at
    least `stats_period` metres. On a device this is still the lever it always was: a
    mode-averaged state is a single column, where computing the default statistics from a
    host copy of the field costs less than the device reductions, so the field is copied
    down on every step they fire on (see [`Stats.collect_stats`](@ref
    Luna.Stats.collect_stats) and `Stats.STATS_DEVICE_MINLEN`; pass
    `stats_kwargs=Dict(:stats_device => :device)` to override the choice). A statistics
    function added through `stats_kwargs[:userfuns]` has no device form at all; `Luna.run`
    warns once per propagation, naming it, when one forces the copy.
"""
function prop_capillary(args...; status_period=5, kwargs...)
    Eω, grid, linop, transform, FT, output = prop_capillary_args(args...; kwargs...)
    #= args[2] is `flength`, the second positional argument of prop_capillary_args: the
       grid no longer carries the propagation length, so it is passed to Luna.run here.
       `makeoutput` builds the save grid from the named `flength`, so if the positional
       layout ever changed the consistency check in Luna.run would raise rather than
       silently propagate the wrong distance. =#
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=args[2], status_period,
             boundary_kwargs(kwargs)...)
    output
end

#= Pick the options which belong to Luna.run -- the absorbing boundaries and the
   tabulation of the linear operator -- out of the user's keyword arguments so they can be
   forwarded to it. They are declared on the *_args functions rather than here so that
   there is one set of defaults and so that saveargs records them; whatever the user did not
   pass simply falls through to Luna.run's own defaults. =#
boundary_kwargs(kwargs) = NamedTuple(
    k => v for (k, v) in pairs(kwargs)
    if k in (:boundary, :boundary_N, :boundary_length, :tcollar,
             :linop_integral, :linop_tol, :tabulate_linop))

#= Error, naming the fix, when an *explicit* `device`/`precision` request cannot be
   honoured well because `resp` (mode-averaged only; multimode/radial go through
   `_cpu_only!` instead) contains a response with no device kernel.
   Does nothing for the CPU/Float64 default, whatever `resp` contains.

   Such a response is not impossible on a device: `Nonlinear.rescale` wraps it in a
   `Nonlinear.HostResponse`, which copies the whole block to the host and back at every
   right-hand side. That is the low-level interface's hackability fallback, and it is far
   slower than the CPU run the caller almost certainly wanted; the simple interface
   therefore refuses rather than silently produces it. =#
function _check_responses_device_capable!(device, precision, resp)
    spec = Luna.withprecision(Luna.resolve_device(device), precision)
    (Luna.arraytype(spec) === Array && Luna.realtype(spec) === Float64) && return nothing
    #= `nameof` rather than the full type: a plasma response's type parameters run to
       several hundred characters and would bury the fix past a wrapped paragraph. =#
    bad = unique(string.(nameof.(typeof.(
        Iterators.filter(!Nonlinear.device_capable, resp)))))
    isempty(bad) && return nothing
    error("this propagation includes a response with no device kernel "*
          "($(join(bad, ", "))), so it cannot run on a device or in reduced precision "*
          "without falling back to the host at every step. Pass device=:cpu (or call "*
          "Luna.set_device(:cpu)) and leave `precision` unset to run on the CPU "*
          "instead, or remove the response for a device-capable run. The low-level "*
          "interface will run it through "*
          "`Nonlinear.HostResponse` if that is really what you want.")
end

"""
    prop_capillary_args(radius, flength, gas, pressure; λ0, λlims, trange, kwargs...)

Prepare to simulate pulse propagation in a hollow fibre using the capillary model. This
function takes the same arguments as `prop_capillary` but instead or running the
simulation and returning the output, it returns the required arguments for `Luna.run`,
which is useful for repeated simulations in an indentical fibre with different initial
conditions.

The propagation length is not among them: run the propagation with

```julia
Eω, grid, linop, transform, FT, output = prop_capillary_args(args...; kwargs...)
Luna.run(Eω, grid, linop, transform, FT, output; zmax=flength)
```
"""
function prop_capillary_args(radius, flength, gas, pressure;
                        λlims, trange, envelope=false, thg=nothing, δt=1,
                        λ0, τfwhm=nothing, τw=nothing, ϕ=Float64[],
                        power=nothing, energy=nothing,
                        pulseshape=:gauss, polarisation=:linear, propagator=nothing,
                        pulses=nothing,
                        shotnoise=true,
                        rng=GLOBAL_RNG,
                        modes=:HE11, model=:full, loss=true,
                        radial_integral_rtol=1e-3, modal_integral=:adaptive,
                        modal_nr=NonlinearRHS.FIXED_NR,
                        modal_nθ=NonlinearRHS.FIXED_Nθ, modal_kronrod=false,
                        raman=nothing, kerr=true, plasma=nothing,
                        stats_kwargs=Dict{Symbol, Any}(),
                        PPT_options=Dict{Symbol, Any}(), preionfrac=0.0,
                        rotation=true, vibration=true, temperature=roomtemp,
                        saveN=201, filepath=nothing,
                        scan=nothing, scanidx=nothing, filename=nothing,
                        boundary=:rate, boundary_N=Boundaries.DEFAULT_N,
                        boundary_length=nothing, tcollar=Boundaries.DEFAULT_TCOLLAR,
                        linop_integral=:auto,
                        linop_tol=LinearOps.DEFAULT_LINOP_TOL,
                        tabulate_linop=nothing,
                        device=nothing, precision=nothing, stats_period=1)

    #= Validated here rather than only inside `Luna.run`, so that a misspelled symbol fails
       before the grid, the FFT plans, the input field and the statistics are built.
       `Luna.run` resolves it again, which is idempotent. =#
    linop_integral = Luna._linop_integral(linop_integral, tabulate_linop)

    # do we have energy in the orthogonal polarisation states, or just the fundamental?
    # if so, we need to treat double the number of modes
    both_modes = needpol(polarisation, pulses)
    @info "Orthogonal polarisation modes are "* (both_modes ? "required." : "not required.")
    #= need to treat vector fields if:
        a) we have both polarisation states in the field AND/OR
        b) the modes themselves contain x and y polarisation components
    =#
    pol = both_modes || needpol_modes(modes)
    @info "Vector fields are "* (pol ? "required." : "not required.")

    plasma = isnothing(plasma) ? !envelope : plasma
    thg = isnothing(thg) ? !envelope : thg

    grid = makegrid(λ0, λlims, trange, envelope, thg, δt)
    mode_s = makemode_s(
        modes, flength, radius, gas, pressure, temperature, model, loss, both_modes)
    check_orth(mode_s)
    density = makedensity(flength, gas, pressure, temperature)
    resp = makeresponse(grid, gas, raman, kerr, plasma, thg, pol, rotation, vibration,
                        PPT_options, preionfrac, temperature)
    inputs = makeinputs(mode_s, λ0, pulses, τfwhm, τw, ϕ,
                        power, energy, pulseshape, polarisation, propagator)
    inputs, noise_field = makenoise(grid, mode_s, inputs, shotnoise, rng)
    #= `device=nothing` means "not specified". It resolves to `Luna.device_request()`,
       i.e. `Luna.settings["device"]` as the user set it -- so an untouched call follows a
       loaded GPU package exactly as the low-level interface does -- only when the
       propagation has a device path (`hasdevicepath` below: mode-averaged, or multimode
       with `modal_integral=:fixed`) *and* every response it was built with is
       device-capable (`Nonlinear.device_capable`, true for the Kerr, plasma, Raman and
       χ⁽²⁾ responses). Otherwise it resolves to the CPU regardless of
       `Luna.settings["device"]`, exactly as gpu/10 hardcoded, so that loading a GPU
       package does not turn a silent, working default call -- multimode on the adaptive
       transverse integral, or a user-written response -- into an error: only an
       *explicit* `device`/`precision` request reaches `Luna.setup_modal`'s check (the
       adaptive transverse integral) or `_check_responses_device_capable!`'s (a response
       with no device kernel).

       An explicit `precision` counts as an explicit request even with `device` left
       unspecified: `precision=Float32` with a response which has no device kernel is a
       scaled run in which that response falls back to the host, which is not what the
       caller asked for, so it is refused with the same message. =#
    #= Which geometries have a device path at all: mode-averaged, and multimode on the
       fixed quadrature rule. The adaptive transverse integral is driven by `Cubature`,
       which is host scalar code returning `Vector{Float64}`. =#
    hasdevicepath = (mode_s isa Modes.AbstractMode) || (modal_integral === :fixed)
    if hasdevicepath && !(isnothing(device) && isnothing(precision))
        _check_responses_device_capable!(something(device, Luna.HostSpec()), precision,
                                         resp)
    end
    devicereq = if !isnothing(device)
        device
    elseif hasdevicepath && all(Nonlinear.device_capable, resp)
        Luna.device_request()
    else
        Luna.HostSpec()
    end
    #= `aefftol`/`aeffspan`: with `linop_integral=:tabulated` and a z-dependent operator,
       the effective area is tabulated here, before `Stats.default` closes over it
       (below), rather than only inside `Luna.run`. `Luna.run` rebinds its own `transform`
       and leaves this one alone, so a statistics function built from `transform.aeff`
       would otherwise keep calling `Modes.Aeff` -- memoised on `(mode, z)`, so once per
       accepted step forever -- for a tapered fibre. The span is the fibre; the
       propagation needs `Aeff` up to one step past the end, and `NonlinearRHS.tabulate`
       rebuilds a wider table for that from this one's source callable.

       A uniform fibre is left alone: its `Aeff` is a constant, `Luna.run` tabulates
       nothing for it, and reading the constant off a two-node table would be
       `(1-s)f + sf` rather than `f` -- a rounding difference for no gain. =#
    tabvalues = (linop_integral !== :quadrature) &&
                (const_linop(radius, pressure) === Val(false))
    linop, Eω, transform, FT = setup(grid, mode_s, density, resp, inputs, pol,
                                     radial_integral_rtol, const_linop(radius, pressure);
                                     noise_field, thg, device=devicereq, precision,
                                     aefftol=(tabvalues ? linop_tol : nothing),
                                     aeffspan=(0.0, float(flength)),
                                     modal_integral, nr=modal_nr, nθ=modal_nθ,
                                     kronrod=modal_kronrod)
    #= The state itself, not a host copy of it: `Stats.default` builds its buffers and
       plans its inverse transform for the array type and precision `Eω` has, and takes
       the unit scaling from `transform`, so that the default statistics run on the
       device with the state left where it is. `_statskwargs` turns off the mode
       reconstruction error for `modal_integral=:fixed`, which has no single-point
       machinery to compute it from. =#
    stats = Stats.default(grid, Eω, mode_s, linop, transform;
                          gas=gas, _statskwargs(transform, stats_kwargs)...)
    stats = Output.maybe_periodic(stats, stats_period)
    output = makeoutput(flength, saveN, stats, filepath, scan, scanidx, filename)

    saveargs(output; radius, flength, gas, pressure, λlims, trange, envelope, thg, δt,
        λ0, τfwhm, τw, ϕ, power, energy, pulseshape, polarisation, propagator, pulses,
        shotnoise, modes, model, loss, raman, kerr, plasma, PPT_options,
        modal_integral, modal_nr, modal_nθ, modal_kronrod,
        temperature, saveN, filepath, filename,
        boundary, boundary_N, boundary_length, tcollar, linop_integral, linop_tol,
        device, precision, stats_period)

    return Eω, grid, linop, transform, FT, output
end

check_orth(mode::Modes.AbstractMode) = nothing
function check_orth(modes)
    if length(modes) > 1
        if !Modes.orthonormal(modes)
            ms = join(modes, "\n")
            error("The selected modes do not form an orthonormal set:\n$ms")
        end
    end
end

function saveargs(output; kwargs...)
    d = Dict{String, String}()
    for (k, v) in kwargs
        d[string(k)] = string(v)
    end
    output(d; group="prop_capillary_args")
end

function needpol(pol)
    if pol == :linear
        return false
    elseif pol in (:circular, :x, :y)
        return true
    else
        error("Polarisation must be :linear, :circular, :x/:y, or an ellipticity, not $pol")
    end
end

needpol(pol::Number) = true
needpol(pulse::Pulses.AbstractPulse) = needpol(pulse.polarisation)
needpol(pulses::Vector{<:Pulses.AbstractPulse}) = any(needpol, pulses)

needpol(pol, pulses::Nothing) = needpol(pol)
needpol(pol, pulse::Pulses.AbstractPulse) = needpol(pulse)
needpol(pol, pulses) = any(needpol, pulses)

needpol_modes(mode::Symbol) = false # mode average
needpol_modes(modes::Number) = false # only HE1m modes

function needpol_modes(modes::Tuple)
    any(modes) do mode
        md = parse_mode(mode)
        md[:kind] ≠ :HE || md[:n] > 1
    end
end


const_linop(radius::Number, pressure::Number) = Val(true)
const_linop(radius, pressure) = Val(false)

function makegrid(λ0, λlims, trange, envelope, thg, δt)
    if envelope
        isnothing(thg) && (thg = false)
        Grid.EnvGrid(λ0, λlims, trange; δt, thg)
    else
        Grid.RealGrid(λ0, λlims, trange; δt)
    end
end

makegrid(λ0::Tuple, args...) = makegrid(λ0[1], args...)

function parse_mode(mode)
    ms = String(mode)
    kind_string = ms[1:2]
    if length(ms) > 4
        throw(DomainError(mode, "Ambiguous mode designation $mode. Pass modes as `Dict`s to disambiguate, e.g. Dict(:kind => :HE, :n => 1, :m => 12)."))
    else
        nstring = ms[3]
        mstring = ms[4]
    end
    Dict(:kind => Symbol(kind_string), :n => parse(Int, nstring), :m => parse(Int, mstring))
end

parse_mode(mode::Dict) = mode

function makemodes_pol(both, args...; kwargs...)
    if both
        if kwargs[:kind] == :HE
            return [Capillary.MarcatiliMode(args...; ϕ=0.0, kwargs...),
                    Capillary.MarcatiliMode(args...; ϕ=π/(2*kwargs[:n]), kwargs...)]
        else # TE/TM: there is only one mode
            return [Capillary.MarcatiliMode(args...; ϕ=0.0, kwargs...)]
        end
    else
        Capillary.MarcatiliMode(args...; kwargs...)
    end
end

function makemode_s(mode::Union{Symbol, Dict}, flength, radius, gas, pressure::Number, temperature, model, loss, both)
    makemodes_pol(both, radius, gas, pressure; T=temperature, model, loss, parse_mode(mode)...)
end

function makemode_s(mode::Union{Symbol, Dict}, flength, radius, gas, pressure::Tuple{<:Number, <:Number},
                    temperature, model, loss, both)
    coren, _ = Capillary.gradient(gas, flength, pressure..., T=temperature)
    makemodes_pol(both, radius, coren; model, loss, parse_mode(mode)...)
end

function makemode_s(mode::Union{Symbol, Dict}, flength, radius, gas, pressure, temperature, model, loss, both)
    Z, P = pressure
    coren, _ = Capillary.gradient(gas, Z, P, T=temperature)
    makemodes_pol(both, radius, coren; model, loss, parse_mode(mode)...)
end

function makemode_s(modes::Int, args...)
    _flatten([makemode_s(Dict(:kind => :HE, :n => 1, :m => m), args...) for m=1:modes])
end

function makemode_s(modes::Tuple, args...)
    _flatten([makemode_s(m, args...) for m in modes])
end

# Iterators.flatten recursively flattens arrays of arrays, but can't handle scalars
_flatten(modes::Vector{<:AbstractArray}) = collect(Iterators.flatten(modes))
_flatten(mode) = mode

function makedensity(flength, gas, pressure::Number, temperature)
    ρ0 = PhysData.density(gas, pressure, temperature)
    z -> ρ0
end

function makedensity(flength, gas, pressure::Tuple{<:Number, <:Number}, temperature)
    _, density = Capillary.gradient(gas, flength, pressure..., T=temperature)
    density
end

function makedensity(flength, gas, pressure, temperature)
    _, density = Capillary.gradient(gas, pressure..., T=temperature)
    density
end

function makeresponse(grid::Grid.RealGrid, gas, raman, kerr, plasma, thg, pol,
                      rotation, vibration, PPT_options, preionfrac, temperature)
    out = Any[]
    if kerr
        if thg
            push!(out, Nonlinear.Kerr_field(PhysData.γ3_gas(gas)))
        else
            push!(out, Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas), length(grid.to)))
        end
    end
    makeplasma!(out, grid, gas, plasma, pol, PPT_options, preionfrac)
    if isnothing(raman)
        raman = gas in (:N2, :H2, :D2, :N2O, :CH4, :SF6)
    end
    if raman
        @info("Including the Raman response (due to molecular gas choice).")
        rr = Raman.raman_response(grid.to, gas;
            rotation, vibration, temp=temperature)
        if thg
            push!(out, Nonlinear.RamanPolarField(grid.to, rr))
        else
            push!(out, Nonlinear.RamanPolarField(grid.to, rr, thg=false))
        end
    end
    Tuple(out)
end

function makeplasma!(out, grid, gas, plasma::Bool, pol,
                     PPT_options, preionfrac)
    # simple true/false => default to PPT for atoms, ADK for molecules
    if ~plasma
        return
    end
    if gas in (:H2, :D2, :N2O, :CH4, :SF6)
        @info("Using ADK ionisation rate (due to molecular gas choice).")
        model = :ADK
    else
        @info("Using PPT ionisation rate.")
        model = :PPT
    end
    makeplasma!(out, grid, gas, model, pol, PPT_options, preionfrac)
end

function makeplasma!(out, grid, gas, plasma::Symbol, pol,
                     PPT_options, preionfrac)
    ionpot = PhysData.ionisation_potential(gas)
    if plasma == :ADK
        ionrate = Ionisation.IonRateADK(gas)
    elseif plasma == :PPT
        ionrate = Ionisation.IonRatePPTCached(gas, grid.referenceλ;
                                                    PPT_options...)
    else
        throw(DomainError(plasma, "Unknown ionisation rate $plasma."))
    end
    Et = pol ? Array{Float64}(undef, length(grid.to), 2) : grid.to
    push!(out, Nonlinear.PlasmaCumtrapz(grid.to, Et, ionrate, ionpot; preionfrac))
end

function makeresponse(grid::Grid.EnvGrid, gas, raman, kerr, plasma, thg, pol,
                      rotation, vibration, PPT_options, preionfrac, temperature)
    plasma && error("Plasma response for envelope fields has not been implemented yet.")
    isnothing(thg) && (thg = false)
    out = Any[]
    if kerr
        if thg
            ω0 = wlfreq(grid.referenceλ)
            r = Nonlinear.Kerr_env_thg(PhysData.γ3_gas(gas), ω0, grid.to)
            push!(out, r)
        else
            push!(out, Nonlinear.Kerr_env(PhysData.γ3_gas(gas)))
        end
    end
    if isnothing(raman)
        raman = gas in (:N2, :H2, :D2, :N2O, :CH4, :SF6)
    end
    if raman
        @info("Including the Raman response (due to molecular gas choice).")
        rr = Raman.raman_response(grid.to, gas;
            rotation, vibration, temp=temperature)
        push!(out, Nonlinear.RamanPolarEnv(grid.to, rr))
    end
    Tuple(out)
end

getAeff(mode::Modes.AbstractMode) = Modes.Aeff(mode)
getAeff(modes) = Modes.Aeff(modes[1])

function makeinputs(mode_s, λ0, pulses::Nothing, τfwhm, τw, ϕ, power, energy,
                    pulseshape, polarisation, propagator)
    if pulseshape == :gauss
        return makeinputs(mode_s, λ0, Pulses.GaussPulse(;λ0, τfwhm, power=power, energy=energy,
                          polarisation, ϕ, propagator))
    elseif pulseshape == :sech
        return makeinputs(mode_s, λ0, Pulses.SechPulse(;λ0, τfwhm, τw, power=power, energy=energy,
                          polarisation, ϕ, propagator))
    else
        error("Valid pulse shapes are :gauss and :sech")
    end
end

function makeinputs(mode_s, λ0, pulses, args...)
    makeinputs(mode_s, λ0, pulses)
end

function findmode(mode_s, pulse)
    if pulse.mode == :lowest
        if pulse.polarisation == :linear
            return [1]
        else
            return [1, 2]
        end
    else
        md = parse_mode(pulse.mode)
        return _findmode(mode_s, md)
    end
end

function _findmode(mode_s::AbstractArray, md)
    return findall(mode_s) do m
        (m.kind == md[:kind]) && (m.n == md[:n]) && (m.m == md[:m])
    end
end


function makeinputs(mode_s, λ0, pulse::Pulses.GaussBeamPulse)
    k = 2π/λ0
    gauss = Fields.normalised_gauss_beam(k, pulse.waist)
    ovlps = [Modes.overlap(mi, gauss) for mi in selectmodes(mode_s, pulse.Nmodes, pulse.polarisation)]
    fields = Any[]
    if pulse.polarisation == :linear
        for (modeidx, ovlp) in enumerate(ovlps)
            energyfac = abs2(ovlp)
            phase = -angle(ovlp)
            sf = scalefield(pulse.timepulse.field, energyfac, phase)
            push!(fields, (mode=modeidx, fields=(sf,)))
        end
    else
        fy, fx = ellfields(pulse.timepulse)
        for (idx, ovlp) in enumerate(ovlps[1:2:end])
            energyfac = abs2(ovlp)
            phase = -angle(ovlp)
            sfy = scalefield(fy, energyfac, phase)
            sfx = scalefield(fx, energyfac, phase)
            push!(fields, (mode=2idx-1, fields=(sfy,)))
            push!(fields, (mode=2idx, fields=(sfx,)))
        end
    end
    Tuple(fields)
end

function selectmodes(mode_s, Nmodes, pol)
    if pol == :linear
        mode_s[1:Nmodes]
    else
        mode_s[1:2Nmodes]
    end
end

selectmodes(mode_s, Nmodes::Symbol, pol) = mode_s

function scalefield(f::Fields.PulseField, fac, phase)
    Fields.PulseField(f.λ0, nmult(f.energy, fac), nmult(f.power, fac), addphase(f.ϕ, phase), f.Itshape)
end

function scalefield(f::Fields.DataField, fac, phase)
    Fields.DataField(f.ω, f.Iω, f.ϕω, nmult(f.energy, fac), addphase(f.ϕ, phase), f.λ0)
end

function scalefield(f::Fields.PropagatedField, fac, phase)
    Fields.PropagatedField(f.propagator!, scalefield(f.field, fac, phase))
end

function addphase(ϕ, phase)
    if phase == 0
        return copy(ϕ)
    end
    if length(ϕ) == 0
        return [phase]
    else
        out = copy(ϕ)
        out[1] += phase
        return out
    end
end

_findmode(mode_s, md) = _findmode([mode_s], md)

function makeinputs(mode_s, λ0, pulse::Pulses.AbstractPulse)
    idcs = findmode(mode_s, pulse)
    (length(idcs) > 0) || error("Mode $(pulse.mode) not found in mode list: $mode_s")
    if pulse.polarisation == :linear || pulse.polarisation == :x
        ((mode=idcs[1], fields=(pulse.field,)),)
    elseif pulse.polarisation == :y
        ((mode=idcs[2], fields=(pulse.field,)),)
    else
        (length(idcs) == 2) || error("Modes not set up for circular/elliptical polarisation")
        f1, f2 = ellfields(pulse)
        ((mode=idcs[1], fields=(f1,)), (mode=idcs[2], fields=(f2,)))
    end
end

function makeinputs(mode_s, λ0, pulses::AbstractVector)
    i = Tuple(collect(Iterators.flatten([makeinputs(mode_s, λ0, pii) for pii in pulses])))
    @debug join(string.(i), "\n")
    return i
end

ellphase(ϕ, pol::Symbol) = ellphase(ϕ, 1.0)
ellphase(ϕ, ε) = addphase(ϕ, π/2 * sign(ε))

ellfac(pol::Symbol) = (1/2, 1/2) # circular
function ellfac(ε::Number)
    (-1 <= ε <= 1) || throw(DomainError(ε, "Ellipticity must be between -1 and 1."))
    (1-ε^2/(1+ε^2), ε^2/(1+ε^2))
end
# sqrt(px/py) = ε => px = ε^2*py; px+py = 1 => px = ε^2*(1-px) => px = ε^2/(1+ε^2)

nmult(x::Nothing, fac) = x
nmult(x, fac) = x*fac

function ellfields(pulse::Union{Pulses.CustomPulse, Pulses.GaussPulse, Pulses.SechPulse})
    f = pulse.field
    py, px = ellfac(pulse.polarisation)
    f1 = Fields.PulseField(f.λ0, nmult(f.energy, py), nmult(f.power, py), f.ϕ, f.Itshape)
    f2 = Fields.PulseField(f.λ0, nmult(f.energy, px), nmult(f.power, px),
                           ellphase(f.ϕ, pulse.polarisation), f.Itshape)
    f1, f2
end

ellfields(pulse::Pulses.DataPulse) = ellfields(pulse, pulse.field)

function ellfields(pulse::Pulses.DataPulse, pf::Fields.PropagatedField)
    f = pf.field
    py, px = ellfac(pulse.polarisation)
    f1 = Fields.DataField(f.ω, f.Iω, f.ϕω, nmult(f.energy, py), f.ϕ, f.λ0)
    f2 = Fields.DataField(f.ω, f.Iω, f.ϕω, nmult(f.energy, px),
                          ellphase(f.ϕ, pulse.polarisation), f.λ0)
    Fields.PropagatedField(pf.propagator!, f1), Fields.PropagatedField(pf.propagator!, f2)
end

function ellfields(pulse::Pulses.DataPulse, pf)
    f = pf
    py, px = ellfac(pulse.polarisation)
    f1 = Fields.DataField(f.ω, f.Iω, f.ϕω, nmult(f.energy, py), f.ϕ, f.λ0)
    f2 = Fields.DataField(f.ω, f.Iω, f.ϕω, nmult(f.energy, px),
                          ellphase(f.ϕ, pulse.polarisation), f.λ0)
    f1, f2
end

function makenoise(grid, mode_s, inputs, shotnoise::Bool, rng)
    shotnoise ? makenoise(grid, mode_s, inputs, :modified, rng) : (inputs, nothing)
end

function makenoise(grid, mode_s, inputs, shotnoise::Symbol, rng)
    if shotnoise == :modified
        nm = mode_s isa AbstractArray ? length(mode_s) : 1
        noise_field = Fields.generate_noise_field(grid; rng, nmodes=nm)
        @info("Modified shot-noise model enabled. Traditional input shot noise is " *
              "disabled (noise enters through nonlinear operator instead).")
        return (inputs, noise_field)
    elseif shotnoise == :input
        inputs = _add_input_shotnoise(inputs, mode_s, rng)
        return (inputs, nothing)
    else
        throw(DomainError(shotnoise, "Unknown shotnoise=$shotnoise. Use true, false, :modified, or :input."))
    end
end

function _add_input_shotnoise(inputs, mode::Modes.AbstractMode, rng)
    (inputs..., (mode=1, fields=(Fields.ShotNoise(rng),)))
end

function _add_input_shotnoise(inputs, modes, rng)
    (inputs..., [(mode=ii, fields=(Fields.ShotNoise(rng),)) for ii in eachindex(modes)]...)
end

#=
For envelope grids the `thg` flag must also be passed to the linear operator, which then
uses a reference frame transparent to the carrier-mixing THG response (see
`LinearOps.getω0`). For real grids the linops take no `thg` choice (`thg` only selects
the response), so no keyword is passed.
=#
linopkw(grid::Grid.RealGrid, thg) = NamedTuple()
linopkw(grid::Grid.EnvGrid, thg) = (; thg)

#= The effective-area callable `Luna.setup` is given: `Modes.Aeff` directly, or a table
   over `aeffspan` when `aefftol` is a tolerance rather than `nothing`. Tabulating here
   rather than only in `Luna.run` is what stops the memoised `Modes.Aeff` `Dict` growing
   once per accepted step for a tapered fibre when the statistics are on -- see the call
   site in `prop_capillary_args`. =#
function makeaeff(mode, aefftol::Nothing, aeffspan)
    z -> Modes.Aeff(mode, z=z)
end

function makeaeff(mode, aefftol, aeffspan)
    #= `quiet`: this table spans the fibre and the statistics of the last accepted step are
       recorded a fraction of a step past the end of it, by design. The table the
       propagation uses is rebuilt by `NonlinearRHS._aefftab` over the whole of what the
       stepper can reach, and that one is not quiet. =#
    LinearOps.TabulatedScalar(z -> Modes.Aeff(mode, z=z), aeffspan...; tol=aefftol,
                              quiet=true)
end

#= `rtol` and the `modal_integral`/quadrature keywords describe the transverse integral
   of a multimode propagation, which mode-averaged propagation does not have; they are
   accepted and ignored here, as `rtol` always was. =#
function setup(grid, mode::Modes.AbstractMode, density, responses, inputs, pol, rtol,
               c::Val{true}; noise_field=nothing, thg=LinearOps.thg_default(grid),
               device=Luna.device_request(), precision=nothing,
               aefftol=nothing, aeffspan=(0.0, 1.0), modal_integral=:adaptive,
               nr=nothing, nθ=nothing, kronrod=nothing)
    @info("Using mode-averaged propagation.")
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, mode, grid.referenceλ;
                                                    linopkw(grid, thg)...)

    #= The operator is constant, so `βfun!` is too: the normalisation folds it in once at
       setup instead of calling it on every right-hand side. `device`/`precision` are the
       caller's request (`prop_capillary`'s own keywords, default `Luna.device_request()`
       so an untouched call reproduces today's behaviour); `Luna.setup` resolves them. =#
    Eω, transform, FT = Luna.setup(grid, density, responses, inputs,
                                   βfun!, makeaeff(mode, aefftol, aeffspan);
                                   noise_field, constβ=true, device, precision)
    linop, Eω, transform, FT
end

function setup(grid, mode::Modes.AbstractMode, density, responses, inputs, pol, rtol,
               c::Val{false}; noise_field=nothing, thg=LinearOps.thg_default(grid),
               device=Luna.device_request(), precision=nothing,
               aefftol=nothing, aeffspan=(0.0, 1.0), modal_integral=:adaptive,
               nr=nothing, nθ=nothing, kronrod=nothing)
    @info("Using mode-averaged propagation.")
    linop, βfun! = LinearOps.make_linop(grid, mode, grid.referenceλ;
                                        linopkw(grid, thg)...)

    # `device`/`precision`: see the constant-operator branch above
    Eω, transform, FT = Luna.setup(grid, density, responses, inputs,
                                   βfun!, makeaeff(mode, aefftol, aeffspan);
                                   noise_field, device, precision)
    linop, Eω, transform, FT
end

needfull(modes) = !all(modes) do mode
    (mode.kind == :HE) && (mode.n == 1)
end

#= `prop_gnlse` is the one thing the simple interface builds which is not device- or
   reduced-precision-capable (GPU_PLAN.md Group E): its `Luna.setup` call takes no
   `device`/`precision` keyword at all. This accepts and validates them instead of
   erroring with an unhelpful "no keyword argument device", so a `device`/`precision`
   request which does not resolve to the CPU in `Float64` -- the only thing that path can
   produce -- gets a message naming the actual limitation. Multimode propagation is
   checked by `Luna.setup` itself, which knows whether the adaptive or the fixed
   transverse integral was asked for. `prop_capillary`/`prop_gnlse` never build a
   free-space run, so the device paths of `TransRadial` (gpu/20), `TransFree` and
   `TransFree2D` (gpu/21) are reached through the low-level interface only. =#
function _cpu_only!(device, precision, what)
    spec = Luna.withprecision(Luna.resolve_device(device), precision)
    (Luna.arraytype(spec) === Array && Luna.realtype(spec) === Float64) || error(
        "$what does not run on a device or in reduced precision yet (GPU_PLAN.md Group "*
        "E); got $(spec). Use prop_capillary for a device-capable mode-averaged or "*
        "multimode (modal_integral=:fixed) run, or device=Luna.HostSpec() (the default) "*
        "here.")
    nothing
end

#= The `device`/`precision` request is passed straight to `Luna.setup`, which refuses
   anything but the host in `Float64` for `modal_integral=:adaptive` and names
   `modal_integral=:fixed` as the fix. `aefftol`/`aeffspan` are accepted and ignored:
   multimode propagation has no single effective area, and `norm_modal` does not use
   one. =#
function setup(grid, modes, density, responses, inputs, pol, rtol, c::Val{true};
               noise_field=nothing, thg=LinearOps.thg_default(grid),
               device=Luna.HostSpec(), precision=nothing, modal_integral=:adaptive,
               nr=NonlinearRHS.FIXED_NR, nθ=NonlinearRHS.FIXED_Nθ, kronrod=false,
               aefftol=nothing, aeffspan=(0.0, 1.0))
    nf = needfull(modes)
    @info(nf ? "Using full 2-D modal integral." : "Using radial modal integral.")
    linop = LinearOps.make_const_linop(grid, modes, grid.referenceλ; linopkw(grid, thg)...)
    Eω, transform, FT = Luna.setup(grid, density, responses, inputs, modes,
                                   pol ? :xy : :y; full=nf, rtol, noise_field,
                                   modal_integral, nr, nθ, kronrod, device, precision)
    linop, Eω, transform, FT
end

function setup(grid, modes, density, responses, inputs, pol, rtol, c::Val{false};
               noise_field=nothing, thg=LinearOps.thg_default(grid),
               device=Luna.HostSpec(), precision=nothing, modal_integral=:adaptive,
               nr=NonlinearRHS.FIXED_NR, nθ=NonlinearRHS.FIXED_Nθ, kronrod=false,
               aefftol=nothing, aeffspan=(0.0, 1.0))
    nf = needfull(modes)
    @info(nf ? "Using full 2-D modal integral." : "Using radial modal integral.")
    linop = LinearOps.make_linop(grid, modes, grid.referenceλ; linopkw(grid, thg)...)
    Eω, transform, FT = Luna.setup(grid, density, responses, inputs, modes,
                                   pol ? :xy : :y; full=nf, rtol, noise_field,
                                   modal_integral, nr, nθ, kronrod, device, precision)
    linop, Eω, transform, FT
end

"""
    _statskwargs(transform, stats_kwargs)

The keyword arguments for `Stats.default`. `Stats.mode_reconstruction_error` is written
against [`NonlinearRHS.TransModal`](@ref Luna.NonlinearRHS.TransModal): it re-evaluates
the transform at one transverse point and compares the result with the modal
reconstruction, which needs the adaptive transform's single-point machinery, and it
records the cubature's own error estimate. The fixed quadrature rule has neither. Its own
embedded error estimate
([`NonlinearRHS.integral_error!`](@ref Luna.NonlinearRHS.integral_error!)) becomes a
statistic in a later branch of the GPU work; until then a `modal_integral=:fixed` run
collects the other default statistics and not this one. An explicit `mode_error` in
`stats_kwargs` is left alone.
"""
_statskwargs(transform, stats_kwargs) = stats_kwargs

function _statskwargs(transform::NonlinearRHS.TransModalFixed, stats_kwargs)
    if haskey(stats_kwargs, :mode_error)
        stats_kwargs[:mode_error] && error(
            "stats_kwargs[:mode_error] = true, but the mode reconstruction error "*
            "statistic is defined only for the adaptive transverse integral "*
            "(modal_integral=:adaptive): it re-evaluates the transform at a single "*
            "transverse point and records the cubature's own error estimate, neither "*
            "of which the fixed quadrature rule has. Leave `mode_error` out (it is off "*
            "by default with modal_integral=:fixed) or pass modal_integral=:adaptive.")
        return stats_kwargs
    end
    _mergekw(stats_kwargs)
end

#= `stats_kwargs` is splatted everywhere else, so a `NamedTuple` works as well as the
   documented `Dict`; `merge` of a `NamedTuple` with a `Dict` is a `MethodError`. The
   default goes first in both, so that anything in `stats_kwargs` wins -- there is
   nothing there to win, by the `haskey` above, but that is the order a reader expects. =#
_mergekw(kw::AbstractDict) = merge(Dict{Symbol, Any}(:mode_error => false), kw)
_mergekw(kw) = merge((; mode_error=false), NamedTuple(pairs(kw)))

function makeoutput(flength, saveN, stats, filepath::Nothing, scan::Nothing, scanidx, filename)
    Output.MemoryOutput(0, flength, saveN, stats)
end

function makeoutput(flength, saveN, stats, filepath, scan::Nothing, scanidx, filename)
    Output.HDF5Output(filepath, 0, flength, saveN, stats)
end

function makeoutput(flength, saveN, stats, filepath, scan, scanidx, filename)
    isnothing(scanidx) && error("scanidx must be passed along with scan.")
    Output.ScanHDF5Output(scan, scanidx, 0, flength, saveN, stats;
                          fdir=filepath, fname=filename)
end

"""
    prop_gnlse(γ, flength, βs; λ0, λlims, trange, kwargs...)

Simulate pulse propagation using the GNLSE.

# Mandatory arguments
- `γ::Number`: The nonlinear coefficient.
- `flength::Number`: Length of the fibre.
- `βs`: The Taylor expansion of the propagation constant about `λ0`.
- `λ0`: (keyword argument) the reference wavelength for the simulation. For simple
    single-pulse inputs, this is also the central wavelength of the input pulse.
- `λlims::Tuple{<:Number, <:Number}`: The wavelength limits for the simulation grid.
- `trange::Number`: The total width of the time grid. To make the number of samples a
    power of 2, the actual grid used is usually bigger.

# Grid options
- `δt::Number`: Time step on the fine grid used for the nonlinear interaction. By default,
    this is determined by the wavelength grid. If `δt` is given **and smaller** than the
    required value, it is used instead.

# Input pulse options
A single pulse can be specified by the keyword arguments below.
More complex inputs can be defined by a single `AbstractPulse` or a `Vector{AbstractPulse}`.
In this case, all keyword arguments except for `λ0` are ignored.
Note that the current GNLSE model is single mode only.

- `λ0`: Central wavelength
- `τfwhm`: The pulse duration as defined by the full width at half maximum.
- `τw`: The "natural" pulse duration. Only available if pulseshape is `sech`.
- `ϕ`: Spectral phases to be applied to the transform-limited pulse. Elements are
    the usual polynomial phases ϕ₀ (CEP), ϕ₁ (group delay), ϕ₂ (GDD), ϕ₃ (TOD), etc.
- `energy`: Pulse energy.
- `power`: Peak power **after any spectral phases are added**.
- `pulseshape`: Shape of the transform-limited pulse. Can be `:gauss` for a Gaussian pulse
    or `:sech` for a sech² pulse.
- `polarisation`: Polarisation of the input pulse. Can be `:linear` (default), `:circular`,
    or an ellipticity number -1 ≤ ε ≤ 1, where ε=-1 corresponds to left-hand circular,
    ε=1 to right-hand circular, and ε=0 to linear polarisation. The major axis for
    elliptical polarisation is always the y-axis.
- `propagator`: A function `propagator!(Eω, grid)` which **mutates** its first argument to
                apply an arbitrary propagation to the pulse before the simulation starts.
- `shotnoise`: Whether and how to include quantum noise. Can be one of:
    - `true` (default) -- same as `:modified`.
    - `false` -- disable all noise.
    - `:modified` -- use the modified shot-noise model of Chen & Wise
      (arXiv:2410.20567), where a constant noise field enters the nonlinear operator
      at every step but is excluded from dispersion. This prevents artificial FWM
      phase-matching and elevated noise floor artefacts.
    - `:input` -- use traditional one-photon-per-mode shot noise added to the input
      field at `z = 0`.
    See the [Noise model](@ref) documentation for details.
- `rng`: Random number generator for noise field generation. Defaults to `GLOBAL_RNG`.
    Pass a seeded RNG (e.g. `MersenneTwister(seed)`) for reproducible noise realisations,
    or different seeds for ensemble/shot-to-shot statistics.

# GNLSE options
- `shock::Bool`: Whether to include the shock derivative term. Default is `true`.
- `raman::Bool`: Whether to include the Raman effect. Defaults to `true`.
- `ramanmodel`; which Raman model to use, defaults to `:sdo` which uses a simple
   damped oscillator model, defined `τ1` and `τ2` (which default to values commonly
   used for silica). `ramanmodel` can also be set to `:SiO2` which uses the more
   advanced model of Hollenbeck and Cantrell.
- `loss`: the power loss [dB/m]. Defaults to 0.
- `fr`: fractional Raman contribution to `γ`. Defaults to `fr = 0.18`.
- `τ1`: the Raman oscillator period.
- `τ2`: the Raman damping time.

# Output options
- `saveN::Integer`: Number of points along z at which to save the field.
- `filepath`: If `nothing` (default), create a `MemoryOutput` to store the simulation results
    only in the working memory. If not `nothing`, should be a file path as a `String`,
    and the results are saved in a file at this location. If `scan` is passed, `filepath`
    determines the output **directory** for the scan instead.
- `scan`: A `Scan` instance defining a parameter scan. If `scan` is given`, a
    `Output.ScanHDF5Output` is used to automatically name and populate output files of
    the scan. `scanidx` must also be given.
- `scanidx`: Current scan index within a scan being run. Only used when `scan` is passed.
- `filename`: Can be used to to overwrite the scan name when running a parameter scan.
    The running `scanidx` will be appended to this filename. Ignored if no `scan` is given.
- `status_period::Number`: Interval (in seconds) between printed status updates.
- `boundary::Symbol=:rate`: How the absorbing boundaries at the edges of the frequency and
    time windows are applied. `:rate` treats them as an absorption rate per unit distance,
    so the total absorption depends only on the propagation distance and not on how many
    steps the solver took. `:legacy` reproduces the historical behaviour, in which the
    windows were applied once per accepted step and the result therefore depended on
    `rtol`. `:none` disables them. See [`Luna.run`](@ref).
- `boundary_N::Real`: Absorber strength, expressed as the number of times the historical
    window profile is applied over the whole propagation length.
- `boundary_length`: Absorber reference length in metres, overriding `boundary_N`.
- `tcollar::Real`: Minimum width of the temporal absorber collar, as a fraction of the time
    window.
- `device`, `precision`: accepted for symmetry with [`prop_capillary`](@ref), but
    `prop_gnlse` is not device- or reduced-precision-capable (its normalisation is built
    before the unit scaling is known); anything other than the default `Luna.HostSpec()`
    in `Float64` errors.
- `stats_period::Real=1`: collect the default statistics less often than every accepted
    step: an integer (default `1`) every `stats_period`-th accepted step, a non-integer
    value every `stats_period` metres of propagation (see [`Output.PeriodicStats`](@ref)
    and [`Output.maybe_periodic`](@ref)).
"""
function prop_gnlse(args...; status_period=5, kwargs...)
    Eω, grid, linop, transform, FT, output = prop_gnlse_args(args...; kwargs...)
    #= args[2] is `flength`, the second positional argument of prop_gnlse_args: the grid no
       longer carries the propagation length, so it is passed to Luna.run here. `makeoutput`
       builds the save grid from the named `flength`, so if the positional layout ever
       changed the consistency check in Luna.run would raise rather than silently propagate
       the wrong distance. =#
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=args[2], status_period,
             boundary_kwargs(kwargs)...)
    output
end

"""
    prop_gnlse_args(γ, flength, βs; λ0, λlims, trange, kwargs...)

Prepare to simulate pulse propagation using the GNLSE. This
function takes the same arguments as `prop_gnlse` but instead or running the
simulation and returning the output, it returns the required arguments for `Luna.run`,
which is useful for repeated simulations in an indentical fibre with different initial
conditions.

The propagation length is not among them: run the propagation with

```julia
Eω, grid, linop, transform, FT, output = prop_gnlse_args(args...; kwargs...)
Luna.run(Eω, grid, linop, transform, FT, output; zmax=flength)
```
"""
function prop_gnlse_args(γ, flength, βs; λ0, λlims, trange,
                        δt=1, τfwhm=nothing, τw=nothing, ϕ=Float64[],
                        power=nothing, energy=nothing,
                        pulseshape=:gauss, propagator=nothing,
                        pulses=nothing,
                        shotnoise=true, shock=true,
                        rng=GLOBAL_RNG,
                        loss=0.0, raman=true, fr=0.18,
                        ramanmodel=:sdo, τ1=12.2e-15, τ2=32e-15,
                        saveN=201, filepath=nothing,
                        scan=nothing, scanidx=nothing, filename=nothing,
                        boundary=:rate, boundary_N=Boundaries.DEFAULT_N,
                        boundary_length=nothing, tcollar=Boundaries.DEFAULT_TCOLLAR,
                        device=Luna.HostSpec(), precision=nothing, stats_period=1)
    #= `device` defaults to the CPU outright (not `Luna.device_request()`): prop_gnlse is
       never device-capable, so it must keep giving the CPU answer whatever
       `Luna.settings["device"]` says, exactly as gpu/10 hardcoded, and only refuse when
       the caller explicitly asks for something else (`_cpu_only!`, just below).

       Unlike `prop_capillary`, `prop_gnlse` builds its own normalisation
       (`norm_mode_average_gnlse`) before the unit scaling is known (`Luna.setup`
       derives it from the peak of the input field, once the transform is built), so it
       cannot yet be handed a `spec`/`scaling` matching a real device or a reduced
       precision the way `Luna.setup`'s own default normalisation is. Refuse rather than
       build a normalisation silently wrong by a factor of `Pref`/`Eref`; the GNLSE device
       path is not in this branch's exit criteria (GPU_PLAN.md §6 Group B). =#
    _cpu_only!(device, precision, "prop_gnlse")
    envelope = true
    thg = false
    polarisation=:linear
    grid = makegrid(λ0, λlims, trange, envelope, thg, δt)
    mode_s = SimpleFibre.SimpleMode(PhysData.wlfreq(λ0), βs; loss)
    aeff = z -> 1.0
    density = z -> 1.0
    linop, βfun!, β1, αfun = LinearOps.make_const_linop(grid, mode_s, λ0)
    k0 = 2π/λ0
    n2 = γ/k0*aeff(0.0)
    # factor of 4/3 below compensates for the factor of 3/4 in Nonlinear.jl, as
    # n2 and γ are usually defined for the envelope case already
    χ3 = 4/3 * (1 - fr) * n2 * (PhysData.ε_0*PhysData.c)
    resp = Any[Nonlinear.Kerr_env(χ3)]
    if raman
        # factor of 2 here compensates for factor 1/2 in Nonlinear.jl as fr is
        # defined for the envelope case already
        χ3R = 2 * fr * n2 * (PhysData.ε_0*PhysData.c)
        if ramanmodel == :SiO2
            push!(resp, Nonlinear.RamanPolarEnv(grid.to, Raman.raman_response(grid.to, :SiO2,
                                                                              χ3R * PhysData.ε_0)))
        elseif ramanmodel == :sdo
            if isnothing(τ1) || isnothing(τ2)
                error("for :sdo ramanmodel you must specify τ1 and τ2")
            end
            push!(resp, Nonlinear.RamanPolarEnv(grid.to,
                Raman.CombinedRamanResponse(grid.to,
                    [Raman.RamanRespNormedSingleDampedOscillator(χ3R * PhysData.ε_0, 1/τ1, τ2)])))
        else
            error("unrecognised value for ramanmodel")
        end
    end
    resp = Tuple(resp)

    inputs = makeinputs(mode_s, λ0, pulses, τfwhm, τw, ϕ,
                        power, energy, pulseshape, polarisation, propagator)
    inputs, noise_field = makenoise(grid, mode_s, inputs, shotnoise, rng)

    norm! = NonlinearRHS.norm_mode_average_gnlse(grid, aeff; shock)
    #= `device=Luna.HostSpec()`: `_cpu_only!` above has already refused anything else, so
       this is not a silent narrowing of the caller's request. =#
    Eω, transform, FT = Luna.setup(grid, density, resp, inputs, βfun!, aeff;
                                   norm!, noise_field, device=Luna.HostSpec())
    stats = Stats.default(grid, Eω, mode_s, linop, transform)
    stats = Output.maybe_periodic(stats, stats_period)
    output = makeoutput(flength, saveN, stats, filepath, scan, scanidx, filename)

    saveargs(output; γ, flength, βs, λlims, trange, envelope, thg, δt,
        λ0, τfwhm, τw, ϕ, power, energy, pulseshape, polarisation, propagator, pulses,
        shotnoise, shock, loss, raman, ramanmodel, fr, τ1, τ2,
        saveN, filepath, filename,
        boundary, boundary_N, boundary_length, tcollar,
        device, precision, stats_period)

    return Eω, grid, linop, transform, FT, output
end

end
