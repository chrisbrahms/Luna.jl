module Luna
import FFTW
import Hankel
import Logging
import LinearAlgebra: mul!, ldiv!
Logging.disable_logging(Logging.BelowMinLevel)

"""
    Luna.settings

Dictionary of global settings for `Luna`.
"""
settings = Dict{String, Any}("fftw_flag" => FFTW.PATIENT,
                             "fftw_threads" => 0,
                             "fftw_wisdom" => true)

"""
    set_fftw_mode(mode)

Set FFTW planning mode for all FFTW transform planning in `Luna`.

Possible values for `mode` are `:estimate`, `:measure`, `:patient`, and `:exhaustive`.
The initial value upon loading `Luna` is `:patient`

# Examples
```jldoctest
julia> Luna.set_fftw_mode(:patient)
0x00000020
```
"""
function set_fftw_mode(mode)
    s = uppercase(string(mode))
    flag = getfield(FFTW, Symbol(s))
    settings["fftw_flag"] = flag
end

"""
    set_fftw_threads(nthr)

Set number of threads to be used by FFTW. If set to `0`, the number of threads used by
FFTW is determined automatically (see `Utils.FFTWthreads`).
"""
function set_fftw_threads(nthr=0)
    settings["fftw_threads"] = nthr
    FFTW.set_num_threads(Utils.FFTWthreads())
end

"""
    set_fftw_wisdom(enabled::Bool)

Enable (`true`, the default) or disable (`false`) the on-disk FFTW wisdom cache.

When enabled, `Utils.loadFFTwisdom()` imports accumulated FFTW wisdom from a file in
`Utils.cachedir()` before planning and `Utils.saveFFTwisdom()` writes it back
afterwards, so that expensive planning modes (`:measure`, `:patient`, `:exhaustive`, see
[`set_fftw_mode`](@ref)) only have to be paid for once per transform shape. When disabled,
both functions do nothing, so the plan FFTW produces depends only on the planning mode and
the transform shape. The FFTW thread count is re-asserted either way.

Disabling wisdom is needed to make runs reproducible: the wisdom file is shared by every
process using the same Julia depot, so wisdom written by an unrelated `:patient` run (for
example Luna's own precompilation) silently changes the plan an `:estimate` run gets, and
with it the order of the floating-point operations. Turning wisdom off also calls
`FFTW.forget_wisdom()` once, so that wisdom already imported into the running process does
not leak into plans made later.

Turning wisdom back on does not re-import the file; the next call to
`Utils.loadFFTwisdom()` (which `Luna.setup` makes) does that.

Like [`set_fftw_mode`](@ref) and [`set_fftw_threads`](@ref), this changes `settings` in the
calling process only: `Scans` workers load `Luna` fresh and start from the defaults, so a
scan that has to be reproducible needs `@everywhere Luna.set_fftw_wisdom(false)` — and the
same for the other two — after the workers exist.
"""
function set_fftw_wisdom(enabled::Bool)
    settings["fftw_wisdom"] = enabled
    enabled || FFTW.forget_wisdom()
    enabled
end

function __init__()
    set_fftw_threads()
end

#= Device.jl is included after Output.jl (rather than right after Utils.jl, where its
   core -- DeviceSpec, alloc, todevice -- would also work) because it defines
   `ScaledOutput`, the device/precision boundary `Luna.run` wraps an output in, and that
   needs `Output.MemoryOutput`, `Output.HDF5Output`, `Output.willsave` and
   `Output.check_cache` to already exist. Nothing between here and there needs Device.jl's
   own definitions. =#
include("Utils.jl")
include("Scans.jl")
include("Output.jl")
include("Device.jl")
include("Maths.jl")
include("PhysData.jl")
include("Grid.jl")
include("Modes.jl")
include("Fields.jl")
include("RK45.jl")
include("LinearOps.jl")
include("Capillary.jl")
include("Antiresonant.jl")
include("RectModes.jl")
include("StepIndexFibre.jl")
include("SimpleFibre.jl")
include("Nonlinear.jl")
include("Ionisation.jl")
include("NonlinearRHS.jl")
include("Boundaries.jl")
include("Processing.jl")
include("Stats.jl")
include("Polarisation.jl")
include("Tools.jl")
include("Plotting.jl")
include("Raman.jl")
include("SFA.jl")
include("Interface.jl")

prop_capillary = Interface.prop_capillary
prop_gnlse = Interface.prop_gnlse
Pulses = Interface.Pulses

Scan = Scans.Scan
runscan = Scans.runscan
makefilename = Scans.makefilename
addvariable! = Scans.addvariable!

export Utils, Scans, Output, Maths, PhysData, Grid, RK45, Modes, Capillary, RectModes,
       Nonlinear, Ionisation, NonlinearRHS, LinearOps, Boundaries, Stats, Polarisation,
       Tools, Plotting, Raman, Antiresonant, Fields, Processing, Interface, SFA,
       prop_capillary, prop_gnlse, Pulses, Scan, runscan, makefilename, addvariable!,
       StepIndexFibre, SimpleFibre

# for a tuple of TimeFields we assume all inputs are for mode 1
function doinput_sm(grid, inputs::Tuple{Vararg{T} where T <: Fields.TimeField}, FT)
    out = fill(0.0 + 0.0im, length(grid.ω))
    for field in inputs
        out .+= field(grid, FT)
    end
    return out
end

# for a single Fields.TimeField we assume a single input for mode 1
function doinput_sm(grid, inputs::Fields.TimeField, FT)
    doinput_sm(grid, (inputs,), FT)
end

function doinput_sm(grid, inputs::Tuple{Vararg{T} where T <: NamedTuple{<:Any, <:Tuple{Vararg{Any}}}}, FT)
    if any([i.mode ≠ 1 for i in inputs])
        error("For mode-averaged propagation, all inputs must be in 1st mode.")
    end
    inputs_flat = Tuple(Iterators.flatten([i.fields for i in inputs]))
    doinput_sm(grid, inputs_flat, FT)
end

"""
    setup(grid, densityfun, responses, inputs, βfun!, aeff; kwargs...)

Set up a mode-averaged propagation: plan the transforms, build the initial
frequency-domain field from `inputs`, and return `(Eω, transform, FT)`.

# Keyword arguments
- `norm!`: the normalisation of the nonlinear polarisation. Built with
    [`NonlinearRHS.norm_mode_average`](@ref Luna.NonlinearRHS.norm_mode_average) for the
    chosen device and precision if not given; a normalisation supplied here must have
    been built for the same ones.
- `noise_field=nothing`: frequency-domain noise field for the modified shot-noise model.
- `constβ=false`: declare that `βfun!` does not depend on `z`, which lets the
    normalisation fold it in once instead of calling it on every right-hand side. True
    whenever the linear operator is constant.
- `device`: where to run, as `:cpu`, `:auto`, `:metal`, `:cuda` or a
    [`Luna.DeviceSpec`](@ref). The default is `Luna.settings["device"]` as the user set
    it (`:cpu` if the key is absent), which a loaded GPU package sets to `:auto`; it is
    resolved here, so that the log line can say what was asked for as well as what was
    chosen. See [`Luna.set_device`](@ref).
- `precision=nothing`: `Float32` to run in single precision on the chosen device,
    `Float64` for double, `nothing` for whatever the device's own spec says (`Float64` on
    the CPU, `Float32` on Metal). A `Float32` run is scaled (see
    [`Luna.UnitScaling`](@ref)) and its state is stored and saved in `ComplexF32`.
"""
function setup(grid::Grid.RealGrid, densityfun, responses, inputs, βfun!, aeff; kwargs...)
    setup_mode_average(grid, densityfun, responses, inputs, βfun!, aeff; kwargs...)
end

@doc (@doc setup)
function setup(grid::Grid.EnvGrid, densityfun, responses, inputs, βfun!, aeff; kwargs...)
    setup_mode_average(grid, densityfun, responses, inputs, βfun!, aeff; kwargs...)
end

"The time-domain element type for a grid at real precision `T`."
timetype(::Grid.RealGrid, ::Type{T}) where {T} = T
timetype(::Grid.EnvGrid, ::Type{T}) where {T} = Complex{T}

"""
    runscaling(transform)

The [`UnitScaling`](@ref) `transform` was built with, or `UNIT_SCALING` (the
identity) for a transform which does not carry one. `NonlinearRHS.TransModeAvg`,
`NonlinearRHS.TransRadial`, `NonlinearRHS.TransFree`, `NonlinearRHS.TransFree2D` and
`NonlinearRHS.TransModalFixed` do; `NonlinearRHS.TransModal` (the adaptive transverse
integral) cannot be device- or reduced-precision-capable at all -- its cubature driver is
host scalar code returning `Vector{Float64}` -- and always runs at `E_ref = 1`. Used
by [`run`](@ref) to decide whether the output needs [`ScaledOutput`](@ref).
"""
runscaling(transform) = UNIT_SCALING
runscaling(transform::NonlinearRHS.TransModeAvg) = transform.scaling
runscaling(transform::NonlinearRHS.TransRadial) = transform.scaling
runscaling(transform::NonlinearRHS.TransFree) = transform.scaling
runscaling(transform::NonlinearRHS.TransFree2D) = transform.scaling
runscaling(transform::NonlinearRHS.TransModalFixed) = transform.scaling

function setup_mode_average(grid, densityfun, responses, inputs, βfun!, aeff;
                            norm! = nothing, noise_field=nothing, constβ=false,
                            device=device_request(), precision=nothing)
    #= The *unresolved* request is the default and the resolution happens here, so that
       `log_device` can tell `:auto` which found no GPU from a plain `:cpu` -- which is
       what a `Scans` worker that only did `using Luna` sees. =#
    spec = withprecision(resolve_device(device), precision)
    T = realtype(spec)
    log_device(spec, device)
    Logging.@info("Setting up and planning FFTs...")
    flush(stderr)
    Utils.loadFFTwisdom()
    #= The input fields are built on the host in Float64 (`Fields` uses host FFTs and
       scalar code), so the transform they need is planned on the host whatever the run
       uses. On the default CPU path it is also the transform `setup` returns. =#
    xh = Array{timetype(grid, Float64)}(undef, length(grid.t))
    FTh = Utils.plan_ft(xh, 1)
    Eωh = doinput_sm(grid, inputs, FTh)
    scaling = unitscaling(T, () -> FTh \ Eωh, PhysData.ε_0)
    xo = alloc(spec, timetype(grid, T), (length(grid.to),))
    FTo = Utils.plan_ft(xo, 1)
    IFTo = Utils.plan_ift(FTo)
    FT = (arraytype(spec) === Array && T === Float64) ? FTh :
         Utils.plan_ft(alloc(spec, timetype(grid, T), (length(grid.t),)), 1)
    Utils.plan_ift(FT) # create inverse FT plans now, so wisdom is saved
    Utils.plan_ift(FTh)
    if isnothing(norm!)
        norm! = NonlinearRHS.norm_mode_average(grid, βfun!, aeff; spec, scaling, constβ)
    else
        NonlinearRHS.check_norm(norm!, spec, scaling)
    end
    transform = NonlinearRHS.TransModeAvg(grid, FTo, IFTo, responses, densityfun, norm!,
                                          aeff; noise_field, spec, scaling)
    Eω = todevice(spec, isunity(scaling) ? Eωh : Eωh ./ scaling.Eref)
    Utils.saveFFTwisdom()
    Logging.@info("Setup finished.")
    flush(stderr)
    Eω, transform, FT
end

# for a tuple of NamedTuple's with tuple fields we assume all is well
function doinput_mm!(Eω, grid, inputs::Tuple{Vararg{T} where T <: NamedTuple{<:Any, <:Tuple{Vararg{Any}}}}, FT)
    for input in inputs
        out = @view Eω[:, input.mode]
        for field in input.fields
            out .+= field(grid, FT)
        end
    end
end

# for a tuple of TimeFields we assume all inputs are for mode 1
function doinput_mm!(Eω, grid, inputs::Tuple{Vararg{T} where T <: Fields.TimeField}, FT)
    doinput_mm!(Eω, grid, ((mode=1, fields=inputs),), FT)
end

# for a single Fields.TimeField we assume a single input for mode 1
function doinput_mm!(Eω, grid, inputs::Fields.TimeField, FT)
    doinput_mm!(Eω, grid, ((mode=1, fields=(inputs,)),), FT)
end

"""
    setup(grid, densityfun, responses, inputs, modes, components; kwargs...)

Set up a multimode (modal) propagation: plan the transforms, build the initial
frequency-domain field from `inputs`, and return `(Eω, transform, FT)`. `modes` is a
collection of [`Modes.AbstractMode`](@ref Luna.Modes.AbstractMode)s and `components` is
`:x`, `:y` or `:xy`.

# Keyword arguments
- `modal_integral=:adaptive`: how the transverse integral of the nonlinear polarisation is
    evaluated. `:adaptive` builds a
    [`NonlinearRHS.TransModal`](@ref Luna.NonlinearRHS.TransModal), which drives an
    adaptive cubature rule on the host; `:fixed` builds a
    [`NonlinearRHS.TransModalFixed`](@ref Luna.NonlinearRHS.TransModalFixed), which uses a
    fixed quadrature rule and is the multimode transform which runs on a device or in
    reduced precision. `:fixed` is a different discretisation of the same integral, so it
    agrees with `:adaptive` to the accuracy of the quadrature rather than to rounding.
- `full=false`: use the full 2-D transverse integral rather than the radial one.
- `norm!`: the normalisation of the nonlinear polarisation, built with
    [`NonlinearRHS.norm_modal`](@ref Luna.NonlinearRHS.norm_modal) for the chosen device
    and precision if not given.
- `rtol=1e-3`, `atol=0.0`, `mfcn=512`: cubature tolerances and evaluation limit
    (`:adaptive` only).
- `maxbatch=$(NonlinearRHS.MODAL_MAXBATCH)`: the largest number of transverse points
    evaluated in one block (`:adaptive` only); see
    [`NonlinearRHS.MODAL_MAXBATCH`](@ref Luna.NonlinearRHS.MODAL_MAXBATCH).
- `nr=$(NonlinearRHS.FIXED_NR)`, `nθ=$(NonlinearRHS.FIXED_Nθ)`, `kronrod=false`,
    `zconstant=nothing`: the quadrature rule (`:fixed` only); see
    [`NonlinearRHS.TransModalFixed`](@ref Luna.NonlinearRHS.TransModalFixed).
- `noise_field=nothing`: `(nω, nmodes)` noise field for the modified shot-noise model.
- `device`, `precision`: where to run and in what precision, as for the mode-averaged
    `setup` above. A device or a `Float32` run needs `modal_integral=:fixed`.

!!! note "Statistics with `modal_integral=:fixed`"
    `Stats.default`'s `mode_error=true` records the diagnostic the transform has: the
    mode reconstruction error and the cubature's error estimate for
    [`NonlinearRHS.TransModal`](@ref Luna.NonlinearRHS.TransModal), and the embedded
    Gauss--Kronrod estimate
    ([`Stats.transverse_integral_error`](@ref Luna.Stats.transverse_integral_error)) for
    a [`NonlinearRHS.TransModalFixed`](@ref Luna.NonlinearRHS.TransModalFixed). The
    latter is `NaN` unless the rule was built with `kronrod=true`.
"""
function setup(grid::Grid.RealGrid, densityfun, responses, inputs,
               modes::Modes.ModeCollection, components; kwargs...)
    setup_modal(grid, densityfun, responses, inputs, modes, components; kwargs...)
end

@doc (@doc setup)
function setup(grid::Grid.EnvGrid, densityfun, responses, inputs,
               modes::Modes.ModeCollection, components; kwargs...)
    setup_modal(grid, densityfun, responses, inputs, modes, components; kwargs...)
end

function setup_modal(grid, densityfun, responses, inputs, modes, components;
                     modal_integral=:adaptive, full=false, norm! = nothing,
                     rtol=1e-3, atol=0.0, mfcn=512,
                     maxbatch=NonlinearRHS.MODAL_MAXBATCH,
                     nr=NonlinearRHS.FIXED_NR, nθ=NonlinearRHS.FIXED_Nθ, kronrod=false,
                     zconstant=nothing, noise_field=nothing,
                     device=device_request(), precision=nothing)
    modal_integral in (:adaptive, :fixed) || error(
        "modal_integral must be :adaptive or :fixed, got $(repr(modal_integral))")
    spec = withprecision(resolve_device(device), precision)
    T = realtype(spec)
    #= The adaptive driver is `Cubature`, which is host scalar code and returns the
       integral and its error estimate as `Vector{Float64}`. Refused here, where the
       alternative can be named, rather than in the transform constructor. =#
    if modal_integral === :adaptive && !(arraytype(spec) === Array && T === Float64)
        error("the adaptive transverse integral (modal_integral=:adaptive, the default) "*
              "runs on the host in Float64: its cubature driver is host scalar code and "*
              "returns Vector{Float64}. It cannot run on $(spec). Pass "*
              "modal_integral=:fixed for the fixed-quadrature transform, which runs on a "*
              "device and in reduced precision, or device=:cpu to stay on the host.")
    end
    log_device(spec, device)
    Logging.@info("Setting up and planning FFTs...")
    flush(stderr)
    ts = Modes.ToSpace(modes, components=components)
    Utils.loadFFTwisdom()
    nmodes = length(modes)
    TTh = timetype(grid, Float64)
    #= The input fields are built on the host in Float64 (`Fields` uses host FFTs and
       scalar code), so the transform they need is planned on the host whatever the run
       uses. =#
    FTt = Utils.plan_ft(Array{TTh}(undef, length(grid.t)), 1)
    Eωh = zeros(ComplexF64, length(grid.ω), nmodes)
    doinput_mm!(Eωh, grid, inputs, FTt)
    #= The unit scaling needs the time-domain input field, which a Float64 run never asks
       for -- hence the thunk, which is also where the host transform it needs is
       planned. =#
    scaling = unitscaling(T, () -> Utils.plan_ft(
        Array{TTh}(undef, length(grid.t), nmodes), 1) \ Eωh, PhysData.ε_0)
    if isnothing(norm!)
        norm! = NonlinearRHS.norm_modal(grid; spec, scaling)
    else
        NonlinearRHS.check_norm(norm!, spec, scaling)
    end
    TT = timetype(grid, T)
    transform = if modal_integral === :adaptive
        NonlinearRHS.TransModal(TT, grid, ts, responses, densityfun, norm!;
                                rtol, atol, mfcn, full, noise_field, maxbatch)
    else
        NonlinearRHS.TransModalFixed(TT, grid, ts, responses, densityfun, norm!;
                                     full, nr, nθ, kronrod, noise_field, zconstant,
                                     spec, scaling)
    end
    # the transform of the state itself, which `Luna.run` gives to the absorbers
    FT = Utils.plan_ft(alloc(spec, TT, (length(grid.t), nmodes)), 1)
    Utils.plan_ift(FT) # create the inverse plan now, so the wisdom is saved
    Eω = todevice(spec, isunity(scaling) ? Eωh : Eωh ./ scaling.Eref)
    Utils.saveFFTwisdom()
    Logging.@info("Setup finished.")
    flush(stderr)
    Eω, transform, FT
end

function doinputs_fs!(Eωk, grid, spacegrid::Grid.TransverseGrid, FT,
                   inputs::Tuple{Vararg{T} where T <: Fields.SpatioTemporalField})
    for field in inputs
        Eωki = field(grid, spacegrid, FT)
        if size(Eωk, 2) == 2
            Eωk .+= Eωki
        else
            # take y-polarisation if only 1 polarisation specified
            Eωk .+= Eωki[:, [2], :, :] # use array index [2] to preserve dimensionality
        end
    end
end

function doinputs_fs!(Eωk, grid, spacegrid::Grid.TransverseGrid, FT,
                   inputs::Fields.SpatioTemporalField)
    doinputs_fs!(Eωk, grid, spacegrid, FT, (inputs,))
end

#= Radial simulations used to be set up with a Hankel.QDHT. Convert, so that scripts
   written against the old interface keep working; Grid.RadialGrid warns once.

   These take the same six concrete positional arguments as the RadialGrid methods below
   rather than `args...`: a `setup(grid::TimeGrid, q::QDHT, args...)` shim is ambiguous with
   the six-argument mode-averaged `setup(grid::RealGrid, densityfun, responses, inputs,
   βfun!, aeff)`, which is exactly the call a legacy radial script makes. =#
function setup(grid::Grid.RealGrid, q::Grid.HankelTransform,
               densityfun, normfun, responses, inputs; kwargs...)
    setup(grid, Grid.RadialGrid(q), densityfun, normfun, responses, inputs; kwargs...)
end

function setup(grid::Grid.EnvGrid, q::Grid.HankelTransform,
               densityfun, normfun, responses, inputs; kwargs...)
    setup(grid, Grid.RadialGrid(q), densityfun, normfun, responses, inputs; kwargs...)
end

"""
    setup(grid, rg::Grid.RadialGrid, densityfun, normfun, responses, inputs; kwargs...)

Set up a radially symmetric free-space propagation: plan the transforms, build the initial
`(ω, polarisation, k⊥)` field from `inputs`, and return `(Eωk, transform, FT)`.

# Keyword arguments
- `noise_field=nothing`: `(nω, npol, nk)` frequency/k-space noise field for the modified
  shot-noise model.
- `device`, `precision`: where and in what precision to run; see the mode-averaged
  [`setup`](@ref) for the meaning. `normfun` is built by the caller, before the device is
  known, so it is retargeted here with
  [`NonlinearRHS.retarget`](@ref Luna.NonlinearRHS.retarget).

The input fields are built on the host in `Float64` (`Fields` is host scalar code), so the
transform they need is planned on the host whatever the run uses; the returned `FT` is the
plan on the *state's* array type, which is what the absorbing boundaries apply.

!!! note "Statistics"
    [`Stats.default`](@ref Luna.Stats.default)`(grid, Eωk, transform, linop)` builds the
    default statistics set for the returned transform: the total energy, the peak
    intensity, `ω0` and the duration on the propagation axis, the beam size, the energy
    in the absorbing collar, and the electron density with plasma.
"""
function setup(grid::Grid.RealGrid, rg::Grid.RadialGrid,
               densityfun, normfun, responses, inputs; kwargs...)
    setup_radial(Float64, grid, rg, densityfun, normfun, responses, inputs; kwargs...)
end

@doc (@doc setup)
function setup(grid::Grid.EnvGrid, rg::Grid.RadialGrid,
               densityfun, normfun, responses, inputs; kwargs...)
    setup_radial(ComplexF64, grid, rg, densityfun, normfun, responses, inputs; kwargs...)
end

function setup_radial(::Type{TH}, grid, rg::Grid.RadialGrid,
                      densityfun, normfun, responses, inputs;
                      noise_field=nothing, device=device_request(),
                      precision=nothing) where {TH}
    spec = withprecision(resolve_device(device), precision)
    T = realtype(spec)
    log_device(spec, device)
    Logging.@info("Setting up and planning FFTs...")
    flush(stderr)
    Utils.loadFFTwisdom()
    #= The normalisation is built by the caller (it is a positional argument), so it does
       not know the device; move it before anything asks it for its shape. =#
    normfun = NonlinearRHS.retarget(normfun, spec)
    np = size(normfun(0), 2) # number of polarisation directions (1 or 2)
    tshape = (length(grid.t), np, rg.N)
    ωshape = (length(grid.ω), np, rg.N)
    # host plans and buffers: the input fields are built in Float64 on the host
    FTh = Utils.plan_ft(zeros(TH, tshape), 1)
    Eωk = Grid.to_kspace(rg, zeros(ComplexF64, ωshape))
    # plan FFT for xy polarisation for field creation
    FT_xy = Utils.plan_ft(zeros(TH, (length(grid.t), 2, rg.N)), 1)
    doinputs_fs!(Eωk, grid, rg, FT_xy, inputs)
    #= The unit scaling needs the peak of the physical time-domain field, which for a
       radial run is the input taken back to (t, pol, r). Only evaluated for Float32. =#
    #= `copy` because a real inverse FFTW plan overwrites its input, and `Eωk` is the
       state. Only evaluated at all for a Float32 run. =#
    scaling = unitscaling(T, () -> Grid.to_rspace(rg, FTh \ copy(Eωk)), PhysData.ε_0)
    xo = alloc(spec, timetype(grid, T), (length(grid.to), np, rg.N))
    FTo = Utils.plan_ft(xo, 1)
    FT = (arraytype(spec) === Array && T === Float64) ? FTh :
         Utils.plan_ft(alloc(spec, timetype(grid, T), tshape), 1)
    Utils.plan_ift(FT) # create inverse FT plans now, so wisdom is saved
    Utils.plan_ift(FTh)
    Utils.plan_ift(FTo)
    transform = NonlinearRHS.TransRadial(
        grid, rg, FTo, responses, densityfun, normfun, np > 1;
        noise_field, spec, scaling)
    Eωk = todevice(spec, isunity(scaling) ? Eωk : Eωk ./ scaling.Eref)
    Utils.saveFFTwisdom()
    Logging.@info("Setup finished.")
    flush(stderr)
    Eωk, transform, FT
end

"""
    setup(grid, xygrid::Grid.FreeGrid, densityfun, normfun, responses, inputs; kwargs...)
    setup(grid, xgrid::Grid.Free2DGrid, densityfun, normfun, responses, inputs; kwargs...)

Set up a 3-D or 2-D Cartesian free-space propagation: plan the transforms, build the
initial `(ω, polarisation, k⊥...)` field from `inputs`, and return `(Eωk, transform, FT)`.

The joint time-and-space transform is one multi-axis plan — region `(1, 3, 4)` in 3-D and
`(1, 3)` in 2-D, the polarisation axis skipped — on every backend.

# Keyword arguments
- `noise_field=nothing`: `(nω, npol, nk...)` frequency/k-space noise field for the
  modified shot-noise model.
- `device`, `precision`: where and in what precision to run; see the mode-averaged
  [`setup`](@ref) for the meaning. `normfun` is built by the caller, before the device is
  known, so it is retargeted here with
  [`NonlinearRHS.retarget`](@ref Luna.NonlinearRHS.retarget).

The input fields are built on the host in `Float64` (`Fields` is host scalar code), so the
transform they need is planned on the host whatever the run uses; the returned `FT` is the
plan on the *state's* array type, which is what the absorbing boundaries apply.

!!! note "Statistics"
    [`Stats.default`](@ref Luna.Stats.default)`(grid, Eωk, transform, linop)` builds the
    default statistics set for the returned transform: the total energy, the peak
    intensity, `ω0` and the duration on the propagation axis, the beam size, the energy
    in the absorbing collar, and the electron density with plasma.
"""
function setup(grid::Grid.RealGrid, xygrid::Grid.FreeGrid,
               densityfun, normfun, responses, inputs; kwargs...)
    setup_free(Float64, grid, xygrid, densityfun, normfun, responses, inputs; kwargs...)
end

@doc (@doc setup)
function setup(grid::Grid.EnvGrid, xygrid::Grid.FreeGrid,
               densityfun, normfun, responses, inputs; kwargs...)
    setup_free(ComplexF64, grid, xygrid, densityfun, normfun, responses, inputs; kwargs...)
end

@doc (@doc setup)
function setup(grid::Grid.RealGrid, xgrid::Grid.Free2DGrid,
               densityfun, normfun, responses, inputs; kwargs...)
    setup_free(Float64, grid, xgrid, densityfun, normfun, responses, inputs; kwargs...)
end

@doc (@doc setup)
function setup(grid::Grid.EnvGrid, xgrid::Grid.Free2DGrid,
               densityfun, normfun, responses, inputs; kwargs...)
    setup_free(ComplexF64, grid, xgrid, densityfun, normfun, responses, inputs; kwargs...)
end

#= The transverse shape and FFT region of a Cartesian free-space grid: `(Nx, Ny)` and
   `(1, 3, 4)` in 3-D, `(Nx,)` and `(1, 3)` in 2-D. Everything else in `setup_free` is
   written once for both. =#
freeshape(xygrid::Grid.FreeGrid) = (length(xygrid.x), length(xygrid.y))
freeshape(xgrid::Grid.Free2DGrid) = (length(xgrid.x),)
freeregion(::Grid.FreeGrid) = (1, 3, 4)
freeregion(::Grid.Free2DGrid) = (1, 3)

freetransform(grid, spacegrid::Grid.FreeGrid, args...; kwargs...) =
    NonlinearRHS.TransFree(grid, spacegrid, args...; kwargs...)
freetransform(grid, spacegrid::Grid.Free2DGrid, args...; kwargs...) =
    NonlinearRHS.TransFree2D(grid, spacegrid, args...; kwargs...)

function setup_free(::Type{TH}, grid, spacegrid, densityfun, normfun, responses, inputs;
                    noise_field=nothing, device=device_request(),
                    precision=nothing) where {TH}
    spec = withprecision(resolve_device(device), precision)
    T = realtype(spec)
    log_device(spec, device)
    Logging.@info("Setting up and planning FFTs...")
    flush(stderr)
    Utils.loadFFTwisdom()
    #= The normalisation is built by the caller (it is a positional argument), so it does
       not know the device; move it before anything asks it for its shape. =#
    normfun = NonlinearRHS.retarget(normfun, spec)
    np = size(normfun(0), 2) # number of polarisation directions (1 or 2)
    xyshape = freeshape(spacegrid)
    region = freeregion(spacegrid)
    tshape = (length(grid.t), np, xyshape...)
    ωshape = (length(grid.ω), np, xyshape...)
    # host plans and buffers: the input fields are built in Float64 on the host
    FTh = Utils.plan_ft(Array{TH}(undef, tshape), region)
    Eωk = zeros(ComplexF64, ωshape)
    # plan the transform for xy polarisation for field creation
    FT_xy = Utils.plan_ft(Array{TH}(undef, (length(grid.t), 2, xyshape...)), region)
    doinputs_fs!(Eωk, grid, spacegrid, FT_xy, inputs)
    #= The unit scaling needs the peak of the physical time-domain field, which the joint
       inverse transform gives directly. `copy` because a real inverse FFTW plan
       overwrites its input, and `Eωk` is the state. Only evaluated for Float32. =#
    scaling = unitscaling(T, () -> FTh \ copy(Eωk), PhysData.ε_0)
    xo = alloc(spec, timetype(grid, T), (length(grid.to), np, xyshape...))
    FTo = Utils.plan_ft(xo, region)
    FT = (arraytype(spec) === Array && T === Float64) ? FTh :
         Utils.plan_ft(alloc(spec, timetype(grid, T), tshape), region)
    Utils.plan_ift(FT) # create inverse FT plans now, so wisdom is saved
    Utils.plan_ift(FTh)
    Utils.plan_ift(FTo)
    #= `xo` is handed to the transform rather than left to the garbage collector: it is
       an oversampled block of exactly the shape and type the transform's own `Eto` has,
       and it has served its purpose once `FTo` is planned. =#
    transform = freetransform(grid, spacegrid, FTo, responses, densityfun, normfun, np > 1;
                              noise_field, spec, scaling, Eto=xo)
    Eωk = todevice(spec, isunity(scaling) ? Eωk : Eωk ./ scaling.Eref)
    Utils.saveFFTwisdom()
    Logging.@info("Setup finished.")
    flush(stderr)
    Eωk, transform, FT
end

linoptype(l::AbstractArray) = "constant"
linoptype(l::LinearOps.TabulatedLinop) = "tabulated"
linoptype(l) = "variable"

gridtype(g::Grid.RealGrid) = "field-resolved"
gridtype(g::Grid.EnvGrid) = "envelope"
gridtype(g) = "unknown"

simtype(g, t, l) = Dict("field" => gridtype(g),
                        "transform" => string(t),
                        "linop" => linoptype(l))

function save_modeinfo_maybe(output, t::NonlinearRHS.AbstractTransModal)
    pol = t.ts.indices == 1:2 ? "xy" : t.ts.indices == 1 ? "x" : "y"
    modeinfos = unnest([Modes.modeinfo(m) for m in t.ts.ms])
    output(modeinfos; group="modes")
    output("polarisation", pol)
end

function unnest(dicts)
    out = Dict{String, Any}()
    for k in keys(dicts[1]) # assuming all dicts have the same keys
        out[sym2string(k)] = [sym2string(di[k]) for di in dicts]
    end
    out
end

sym2string(sym::Symbol) = string(sym)
sym2string(other) = other

save_modeinfo_maybe(output, t) = nothing

"""
    run(Eω, grid, linop, transform, FT, output; zmax, kwargs...)

Run the propagation over a distance `zmax`.

# Arguments
- `zmax::Real`: the propagation length in metres. Required: the grid no longer carries it.
    It is converted to `Float64` and written to `output` as a top-level entry.

    If `output` is an output object whose save condition is a fixed grid
    (`Output.GridCondition`, which is what `MemoryOutput(zmin, zmax, saveN)` and
    `HDF5Output(path, zmin, zmax, saveN)` build), `zmax` must agree with the end of that
    save grid. That check, and the guard which stops the `zmax` entry being written twice
    when an `HDF5Output` propagation is resumed, both need the output object itself. If
    `output` is a closure or another wrapper around one -- the way several outputs are
    driven at once -- neither applies: the save grid is not checked, and a resumed run warns
    that the file already has the entry and overwrites it with the same value.
- `max_dz::Real=zmax/2`: the largest step the solver may take.

# Absorbing boundaries
- `boundary::Symbol=:rate`: how the absorbing boundaries at the edges of the frequency and
    time windows are applied.
    - `:rate` applies them as an absorption *coefficient per unit distance*: the spectral
      absorber is folded into `linop` (so the propagator applies it exactly) and the
      temporal absorber is applied to the field as `exp(-α_t Δz/2)`, `α` being a power
      coefficient as it is elsewhere in Luna. The total absorption over a distance depends
      only on that distance, not on how many steps the solver took.
    - `:legacy` reproduces the historical behaviour — multiplying the solution by
      `grid.ωwin` and `grid.twin` after every accepted step. This makes the cumulative
      absorption depend on the step count and hence on `rtol`, and tightening the tolerance
      approaches a hard truncation at the edge of the flat region rather than the taper the
      profiles describe. Provided only to reproduce older results.
    - `:none` disables the absorbers. The nonlinear polarisation is still band-limited in
      `NonlinearRHS` and the field is still band-limited once at the start, but nothing
      stops energy wrapping around the time window or piling up at the edge of the
      frequency window.
- `boundary_N::Real=$(Boundaries.DEFAULT_N)`: absorber strength, expressed as the number of
    times the historical window profile is applied over the whole propagation length. The
    reference length is `ℓ = zmax/boundary_N`.
- `boundary_length=nothing`: reference length `ℓ` in metres, overriding `boundary_N`.
- `tcollar::Real=$(Boundaries.DEFAULT_TCOLLAR)`: width of the temporal absorber collar as a
    fraction of the time window. Only used if it is wider than the collar `grid.twin`
    already has (which can be nearly zero, depending on how `trange` rounds up).
- `kcollar::Real=$(Boundaries.DEFAULT_KCOLLAR)`: free space only; width of the k-space
    absorber collar as a fraction of the largest transverse wavevector on the grid. `0`
    disables it.
- `rcollar::Real=$(Boundaries.DEFAULT_RCOLLAR)`: free space only; width of the transverse
    absorber collar of a radial ([`Grid.RadialGrid`](@ref)) grid as a fraction of its
    aperture `R`. The Cartesian
    grids use the window they were built with (`window_factor`) instead.

In free space the evanescent part of `linop` is also made safe for the stepper in every
`boundary` mode: its decay is clamped and the nonlinear source is tapered to match, over
the absorber reference length (`:rate`) or `max_dz` (`:none`, `:legacy`). See "Free
space" in [`Luna.Boundaries`](@ref) and [`NonlinearRHS.FreeSpaceNorm`](@ref).

See [`Luna.Boundaries`](@ref) for the rationale.

# Tabulating a z-dependent linear operator
- `tabulate_linop::Bool=false`: tabulate the z-dependent quantities of the propagation --
    the integrated linear operator `Φ(z) = ∫ linop dz'`
    ([`LinearOps.TabulatedLinop`](@ref Luna.LinearOps.TabulatedLinop)), and the
    propagation constant `β(z)` and effective area `Aeff(z)` of a mode-averaged transform
    -- on adaptively placed z nodes at setup, instead of evaluating them on the host at
    every stage of every step. A tapered or pressure-graded run then does no host work
    inside the right-hand side, which is what a device run needs; on the CPU it replaces
    the per-stage `Modes.neff` loop with an interpolation.

    It **changes the discretisation**: the propagator becomes `exp(Φ(t2) − Φ(t1))`, the
    exact interaction-picture propagator of the linear part over the step, where the
    default is `exp(linop(t2)·(t2 − t1))`, a one-point rule. Both converge to the same
    solution as the step shrinks and the difference is largest where the operator varies
    fastest within a step -- the entrance of a `p₀ = 0` pressure gradient. It is off by
    default for that reason, and a run which uses it is not comparable element by element
    with one which does not.

    A constant operator is unaffected: it is already exact in the propagator and is not
    tabulated. The tables are built for the propagation: the `transform` object the caller
    passed in is not modified, so a statistics function built from `transform.aeff` before
    the run keeps calling the untabulated one (once per accepted step, on the host, where
    the statistics already are).
- `linop_tol::Real=$(LinearOps.DEFAULT_LINOP_TOL)`: the tolerance the nodes are placed to
    satisfy, in radians for the integrated operator (absolute) and relative for `β` and
    `Aeff`. The tables cost `2·length(Eω)·nnodes` numbers for the operator, so a tolerance
    far below the solver's own `rtol` buys nothing and costs memory.
"""
function run(Eω, grid,
             linop, transform, FT, output;
             zmax=nothing, min_dz=0, max_dz=nothing, init_dz=1e-4, z0=0.0,
             rtol=1e-6, atol=1e-10, safety=0.9, norm=RK45.weaknorm,
             status_period=1,
             boundary=:rate, boundary_N=Boundaries.DEFAULT_N, boundary_length=nothing,
             tcollar=Boundaries.DEFAULT_TCOLLAR, kcollar=Boundaries.DEFAULT_KCOLLAR,
             rcollar=Boundaries.DEFAULT_RCOLLAR,
             tabulate_linop=false, linop_tol=LinearOps.DEFAULT_LINOP_TOL)

    isnothing(zmax) && error(
        "Luna.run requires the propagation length as the keyword argument zmax, e.g. "*
        "Luna.run(Eω, grid, linop, transform, FT, output; zmax=flength). The grid no "*
        "longer stores it.")
    #= `grid.zmax` was a Float64 field, so the length reached the absorbers, the stepper and
       the output as a Float64 whatever the caller wrote. Keep that. =#
    zmax = float(zmax)
    isnothing(max_dz) && (max_dz = zmax/2)
    #= The save grid and zmax are two separate statements of the propagation length; check
       here that they cannot drift apart. Only possible when `output` is the output object
       itself -- a closure wrapping one hides the save condition, and then this is skipped. =#
    if hasproperty(output, :save_cond) && output.save_cond isa Output.GridCondition
        zend = output.save_cond.grid[end]
        zend ≈ zmax || error(
            "zmax ($zmax m) does not agree with the end of the output's save grid "*
            "($zend m). The output was created for a different propagation length; pass "*
            "the same length to both.")
    end

    #= Absorbing boundaries used to be host scalar code, so a device run needed
       boundary=:none. `Boundaries.RateAbsorber`/`LegacyAbsorber` are now broadcasts and
       reductions over mirrored arrays (see `Boundaries.jl`), so every `boundary` mode
       works on a device. Every transform except `TransModal` can produce a device `Eω`:
       `TransModeAvg`, `TransRadial`, `TransFree`, `TransFree2D` and `TransModalFixed`
       (`modal_integral=:fixed`) all take `device`/`precision`. `TransModal` cannot take
       them at all (its cubature driver is host scalar code returning `Vector{Float64}`)
       and always builds a host array. =#

    #= Et is the time-domain buffer the absorbers and the transverse collar work on --
       nothing about its *contents* matters here, only its shape and element type, since
       every absorber overwrites it before reading it back. Building it with `similar`
       rather than an actual inverse transform means this costs nothing and needs no
       `\` on the state's plan (which the generic device planners do not define; see
       `Utils.plan_ift`). =#
    Et = similar(Eω, timetype(grid, real(eltype(Eω))), (length(grid.t), size(Eω)[2:end]...))

    #= The unit scaling the state and the polarisation are expressed in (`Luna.jl`'s
       `unitscaling`, GPU_PLAN.md 4.1): the identity for every transform which does not
       carry one. `NonlinearRHS.TransModeAvg`, `TransRadial`, `TransFree2D` and
       `TransFree` all do; `TransModal` does not. Needed here only to decide whether the
       output needs unscaling. =#
    scaling = runscaling(transform)

    #= `Output.jl` stays device-unaware (`ScaledOutput`'s docstring): wrap whenever the
       state lives on a device or the run is scaled (`E_ref != 1`, every `Float32` run,
       host or device). The default CPU path in `Float64` is untouched -- `output` is the
       same object it always was. Must happen after the zmax/save_cond check above,
       which needs the wrapped output's own `save_cond` field, and before `check_cache`
       and `Boundaries.setup`, both of which now see host, physical-unit data whenever
       they read `y` off the wrapped output. =#
    if isdevice(Eω) || !isunity(scaling)
        output = ScaledOutput(output, Eω, scaling.Eref)
    end

    # check_cache does nothing except for HDF5Outputs
    Eωc, zc, dzc = Output.check_cache(output, Eω, z0, init_dz)
    if zc > z0
        Logging.@info("Found cached propagation. Resuming...")
        #= `check_cache` reads the cached field back from the file as a host array in
           physical units (`ScaledOutput` never touches the cache write itself, only
           gates when it happens) -- so it has to be rescaled and put back on the state's
           array type and precision. On the default CPU path this returns the array it
           was given and `isunity(scaling)` skips the division. =#
        Eωr = isunity(scaling) ? Eωc : Eωc ./ scaling.Eref
        Eω, z0, init_dz = upload_like(Eω, Eωr), zc, dzc
    end

    #= NOTE: this must come after check_cache, which can move z0 and init_dz: the temporal
       absorber measures the distance it is applied over from z0. =#
    absorber = Boundaries.setup(boundary, grid, transform, linop, Et, FT, output, z0,
                                zmax, max_dz, init_dz;
                                N=boundary_N, ℓ=boundary_length, collar=tcollar,
                                kcollar, rcollar)
    stepfun = absorber.stepfun
    linop, max_dz, init_dz = absorber.linop, absorber.max_dz, absorber.init_dz

    #= Band-limit the field once, here, rather than after every step. `NonlinearRHS` zeroes
       the nonlinear polarisation outside `grid.sidx` and the linear operator is zero there,
       so nothing can put anything back; all this removes is the numerical dust the input
       transform leaves outside the band. Not for :legacy, which must reproduce the
       historical scheme exactly -- its per-step `ωwin` multiply does the same job from the
       first step anyway.

       One masked broadcast rather than logical indexing, so that it runs on any array
       type. `ifelse` leaves the in-band elements untouched and writes an exact zero
       elsewhere, which is what the indexed assignment did. =#
    if boundary !== :legacy
        sidx = reshape(mask_like(Eω, grid.sidx), :, ntuple(_ -> 1, ndims(Eω) - 1)...)
        z0c = zero(eltype(Eω))
        @. Eω = ifelse(sidx, Eω, z0c)
    end

    #= The linear operator is built on the host, in Float64, and wrapped by
       `Boundaries.setup` -- so the upload has to happen after it. A constant operator is
       uploaded once, here; a z-dependent one is either tabulated just below or evaluated
       into a host buffer per stage by `RK45.make_prop!`. On the default CPU path this
       returns the operator unchanged. =#
    linop = upload_like(Eω, linop)

    #= Tabulation (GPU_PLAN.md section 4.5 layer 2), opt-in. The table has to cover every
       z the stepper can ask about: `RK45.solve` runs `while tn <= tmax`, so the last step
       starts at or before `zmax` and ends up to `max_dz` past it, and that step's stages
       are what the last saved plane is interpolated from. `max_dz` is the absorber's,
       which is the one the stepper will be given. Both tabulations happen after
       `Boundaries.setup`, so the operator includes the absorber and the evanescent clamp.

       A constant operator is already exact in the propagator and is left alone; the
       transform's own z-dependent quantities are tabulated either way, which for a
       constant operator is a two-node table and no change to the arithmetic that matters. =#
    if tabulate_linop
        #= `init_dz` as well as `max_dz`: the first step is taken at `init_dz` before
           `steplims!` has had a chance to clamp it, and `Boundaries.setup` only reduces
           it to `max_dz` for `boundary=:rate`. =#
        ztab = zmax + max(max_dz, init_dz)
        transform = NonlinearRHS.tabulate(transform, z0, ztab, linop_tol, Eω)
        if !(linop isa AbstractArray)
            linop = LinearOps.TabulatedLinop(linop, Eω, z0, ztab; tol=linop_tol)
        end
    end

    output(Grid.to_dict(grid), group="grid")
    #= Written once: on a resumed HDF5 propagation it is already in the file. An output
       which cannot be queried (a closure, say) is written to again and warns. =#
    Output.hasdata(output, "zmax") || output("zmax", zmax)
    st = simtype(grid, transform, linop)
    st["boundary"] = string(boundary)
    boundary === :rate && (st["boundary_length"] = string(absorber.ℓ))
    sg = Boundaries.spacegrid(transform)
    if !isnothing(sg)
        st["kcollar"] = string(kcollar)
        st["rcollar"] = string(rcollar)
        output(Grid.to_dict(sg), group="spacegrid")
    end
    output(st, group="simulation_type")
    save_modeinfo_maybe(output, transform)

    flush(stderr) # flush std error once before starting to show setup steps
    RK45.solve_precon(
        transform, linop, Eω, z0, init_dz, zmax, stepfun=stepfun,
        max_dt=max_dz, min_dt=min_dz,
        rtol=rtol, atol=atol, safety=safety, norm=norm,
        status_period=status_period)
end

# run some code for precompilation
Logging.with_logger(Logging.NullLogger()) do
    prop_capillary(125e-6, 0.3, :He, 1.0; λ0=800e-9, energy=1e-9,
                    τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12, saveN=11,
                    PPT_options=Dict(:cache => false))
    prop_capillary(125e-6, 0.3, :He, (1.0, 0); λ0=800e-9, energy=1e-9,
                    τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12, saveN=11,
                    PPT_options=Dict(:cache => false))
    prop_capillary(125e-6, 0.3, :He, 1.0; λ0=800e-9, energy=1e-9,
                    τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12, saveN=11,
                    modes=4, PPT_options=Dict(:cache => false))
    p = Tools.capillary_params(120e-6, 10e-15, 800e-9, 125e-6, :He, P=1.0)

    # gnlse_sol.jl example but with N=1 and 100th of the fibre length
    γ = 0.1
    β2 = -1e-26
    τ0 = 280e-15
    τfwhm = (2*log(1 + sqrt(2)))*τ0
    fr = 0.18
    P0 = abs(β2)/((1-fr)*γ*τ0^2)
    flength = π*τ0^2/abs(β2)/100
    βs =  [0.0, 0.0, β2]
    λ0 = 835e-9
    λlims = [450e-9, 8000e-9]
    trange = 4e-12
    output = prop_gnlse(γ, flength, βs; λ0, τfwhm, power=P0, pulseshape=:sech, λlims, trange,
                        raman=true, shock=true, fr, shotnoise=true, saveN=11)
end

end # module
