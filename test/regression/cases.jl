#= Regression case matrix.

   Each case is a small propagation which is run twice: once with a fixed step sequence
   (`min_dz == max_dz == init_dz`, which bypasses the step-size controller through
   `RK45.steplims!`) and once with the adaptive controller. The fixed-step runs make a
   difference attributable to the operation that produced it; the adaptive runs are the
   sanity check that the controller still takes the same decisions.

   The cases are deliberately small (each runs in a few seconds) and have `shotnoise=false`
   so that they are deterministic. `saveN=11`.

   This file is `include`d both by the gate (`test/test_regression.jl`), which runs it
   against the working tree, and by `test/regression/run_cases.jl`, which the baseline
   generator copies into a worktree of an older commit. It must therefore only use API that
   exists on the oldest commit it will ever be run against, or go through the compatibility
   shim below.
=#
module RegressionCases

using Luna
import Luna: Boundaries, Capillary, Fields, Grid, Interface, LinearOps, Modes,
             Nonlinear, NonlinearRHS, Output, PhysData, RectModes, Stats
#= Hankel is a direct dependency of Luna, so import it directly rather than as `Luna.Hankel`:
   `src/Luna.jl` keeps an `import Hankel` only so that existing scripts still resolve
   `Luna.Hankel`, and a later branch removing that dead import must not break the gate. =#
import Hankel
import Logging

export CASES, overrides, runcase

# =======================================================================================
# COMPATIBILITY SHIM
#
# `generate.jl` copies this file into a worktree of the baseline commit, so the same
# `cases.jl` has to build the same propagations on `evanescent` (fdf8dbe3) and on the
# Group A branches, which changed two APIs:
#
#   gpu/01-zmax      `Grid.RealGrid`/`Grid.EnvGrid` lose their leading `zmax` argument and
#                    their `zmax` field; `Luna.run` takes `zmax` as a required keyword;
#                    `Boundaries.setup` takes it as a positional argument after `z0`.
#   gpu/02-radialgrid `Grid.RadialGrid` replaces `Hankel.QDHT` as the radial transverse
#                    grid.
#
# The branch is detected from the API itself, not from a version or a commit, so nothing
# here has to be updated when a branch is merged. Every helper reduces to exactly the call
# the pre-Group-A `cases.jl` made when the old API is what is present, so a baseline
# generated with this file from `evanescent` is bit-identical to one generated with the
# original file.
#
# REMOVE THIS BLOCK, and the calls to it, once `evanescent` is no longer a comparison base
# for the regression gate -- i.e. once every baseline in use is from `gpu/int-A` or later.
# =======================================================================================

"`true` when `Grid.RealGrid` still carries a `zmax` field, i.e. before `gpu/01-zmax`."
const GRID_HAS_ZMAX = hasfield(Grid.RealGrid, :zmax)

"`true` when `Grid.RadialGrid` exists, i.e. from `gpu/02-radialgrid` on."
const HAS_RADIALGRID = isdefined(Grid, :RadialGrid)

"""
    makegrid(GT, zmax, referenceλ, λ_lims, trange; kwargs...)

`Grid.RealGrid`/`Grid.EnvGrid` with or without the leading `zmax` argument, whichever the
loaded Luna wants. `GT` is `Grid.RealGrid` or `Grid.EnvGrid`.

Do not pass `δt` as a keyword while `evanescent` is still a baseline commit: there
`RealGrid` takes it as its fifth *positional* argument (`EnvGrid` takes it as a keyword on
both sides), so `δt=...` would reach the deprecated `RealGrid` method as an unsupported
keyword. No case needs it; `thg`, which the BBO envelope case does pass, is a keyword on
both sides and goes through `kwargs...` unchanged.
"""
makegrid(GT, zmax, referenceλ, λ_lims, trange; kwargs...) =
    GRID_HAS_ZMAX ? GT(zmax, referenceλ, λ_lims, trange; kwargs...) :
                    GT(referenceλ, λ_lims, trange; kwargs...)

"""
    runkw(zmax)

The `zmax` keyword for `Luna.run` as a `NamedTuple`: `(; zmax)` from `gpu/01-zmax` on, and
empty before it, where `Luna.run` reads the length off the grid and rejects the keyword.
"""
runkw(zmax) = GRID_HAS_ZMAX ? NamedTuple() : (; zmax=zmax)

"""
    radialgrid(R, N)

The radial transverse grid: `Grid.RadialGrid(R, N)` from `gpu/02-radialgrid` on, and
`Hankel.QDHT(R, N, dim=3)` before it. `dim=3` is the axis Luna's free-space arrays put the
radial coordinate on; a `RadialGrid` always transforms along the last dimension.
"""
radialgrid(R, N) = HAS_RADIALGRID ? Grid.RadialGrid(R, N) : Hankel.QDHT(R, N, dim=3)

#= The cases store `Eω` exactly as `Luna.run` produces it, in k-space, and the gate compares
   it there, so nothing here needs the inverse transform (`Grid.to_rspace` after
   `gpu/02-radialgrid`, `q \\ Eω` before it). =#

"""
    absorber_setup(boundary, grid, transform, linop, Et, FT, output, z0, zmax, max_dz,
                   init_dz; kwargs...)

`Boundaries.setup` with `zmax` in the positional list from `gpu/01-zmax` on, and without it
before, where `Boundaries.setup` reads the length off the grid. Used by `benchmark/run.jl`,
which builds an absorber outside `Luna.run`.
"""
function absorber_setup(boundary, grid, transform, linop, Et, FT, output, z0, zmax,
                        max_dz, init_dz; kwargs...)
    if GRID_HAS_ZMAX
        Boundaries.setup(boundary, grid, transform, linop, Et, FT, output, z0,
                         max_dz, init_dz; kwargs...)
    else
        Boundaries.setup(boundary, grid, transform, linop, Et, FT, output, z0, zmax,
                         max_dz, init_dz; kwargs...)
    end
end

# ============================ end of the compatibility shim ============================

"Number of steps in the fixed-step mode. Must be at least `Boundaries.DEFAULT_N` so that
 the `:rate` absorber does not reduce `max_dz` below the requested step."
const NSTEPS = 20

"Number of saved points along z, for every case."
const SAVEN = 11

"`status_period` for the runs: large enough that nothing is printed."
const STATUS_PERIOD = 1e6

"""
    Case(name, zmax, prepare, runkwargs=NamedTuple())

One regression case.

- `prepare()` builds the propagation and returns `(Eω, grid, linop, transform, FT, output)`,
  i.e. exactly the arguments of `Luna.run`, as a fresh set of arrays every call.
- `runkwargs` are the case-specific keyword arguments for `Luna.run`: the absorbing-boundary
  options, which have to match the ones the case was built with.
- `zmax` is the propagation length, which the step-size overrides are derived from.

`prepare` is separate from the run so that `benchmark/run.jl` can time the pieces of a case
without propagating it.
"""
struct Case
    name::String
    zmax::Float64
    prepare::Function
    runkwargs::NamedTuple
end

Case(name, zmax, prepare) = Case(name, zmax, prepare, NamedTuple())

"""
    overrides(case, mode)

Run overrides for `case` in `mode`, which is either `:fixed` (a fixed step sequence) or
`:adaptive` (the step-size controller, started from a small step).
"""
function overrides(c::Case, mode::Symbol)
    if mode === :fixed
        h = c.zmax/NSTEPS
        (min_dz=h, max_dz=h, init_dz=h)
    elseif mode === :adaptive
        (init_dz=c.zmax/1000,)
    else
        error("Unknown regression mode $mode")
    end
end

"""
    runcase(case, mode; perturb=0.0)

Run `case` in `mode` and return the `Output.MemoryOutput`. `perturb` multiplies the initial
frequency-domain field by `1 + perturb` before the propagation; that is how
`test/regression/sensitivity.jl` measures the one-ulp sensitivity.
"""
function runcase(c::Case, mode::Symbol; perturb=0.0)
    quiet() do
        Eω, grid, linop, transform, FT, output = c.prepare()
        perturb == 0 || (Eω .*= (1 + perturb))
        Luna.run(Eω, grid, linop, transform, FT, output;
                 status_period=STATUS_PERIOD, runkw(c.zmax)...,
                 c.runkwargs..., overrides(c, mode)...)
        output
    end
end

#= Suppress `Info` and below, not everything: a warning raised by a later branch -- the
   absorber's "removed more than `warnfrac` of the pulse", a deprecation from a changed
   signature -- has to be visible in the gate's output, and a `NullLogger` would swallow it. =#
"Run `f()` with `Info`-level logging and below suppressed. Warnings and errors still print."
quiet(f) = Logging.with_logger(f, Logging.ConsoleLogger(stderr, Logging.Warn))

#= The absorbing-boundary options have to be given both to `prop_capillary_args` (which
   records them) and to `Luna.run` (which uses them). `Luna.prop_capillary` does this with
   a private helper; repeating the key list here keeps this file independent of it. =#
const BOUNDARY_KEYS = (:boundary, :boundary_N, :boundary_length, :tcollar)
boundary_kwargs(kw) = NamedTuple(k => v for (k, v) in pairs(kw) if k in BOUNDARY_KEYS)

"""
    capillary_case(name, flength, args...; kwargs...)

A case built from `Interface.prop_capillary_args`. `args` and `kwargs` are what
`prop_capillary` would be called with.
"""
function capillary_case(name, args...; kwargs...)
    flength = args[2]
    prepare = () -> Interface.prop_capillary_args(args...; saveN=SAVEN, shotnoise=false,
                                                  kwargs...)
    Case(name, float(flength), prepare, boundary_kwargs(kwargs))
end

"""
    lowlevel_case(name, zmax, setup)

A case built from the low-level interface. `setup()` must return
`(Eω, grid, linop, transform, FT, output)`.
"""
lowlevel_case(name, zmax, setup) = Case(name, float(zmax), setup)

# ---------------------------------------------------------------------------------------
# Shared parameters for the capillary cases
# ---------------------------------------------------------------------------------------
const A_CAP = 125e-6
const L_CAP = 0.1
const Λ0 = 800e-9
const ΛLIMS = (200e-9, 3e-6)
const TRANGE = 0.5e-12
const ΤFWHM = 10e-15
const NOCACHE = Dict{Symbol, Any}(:cache => false) # keep the shared PPT cache out of it

# ---------------------------------------------------------------------------------------
# Low-level case setups
# ---------------------------------------------------------------------------------------

"Mode-averaged field propagation in a two-component gas mixture (He/Ne), Kerr only."
function setup_mixture()
    a = 13e-6
    gases = (:He, :Ne)
    pressures = (1.0, 1.0)
    L = 0.02
    grid = makegrid(Grid.RealGrid, L, Λ0, ΛLIMS, TRANGE)
    m = Capillary.MarcatiliMode(a, gases, pressures; loss=false)
    aeff(z) = Modes.Aeff(m; z=z)
    dens = [PhysData.density(g, p) for (g, p) in zip(gases, pressures)]
    densityfun = let dens=dens
        z -> dens
    end
    responses = ((Nonlinear.Kerr_field(PhysData.γ3_gas(gases[1])),),
                 (Nonlinear.Kerr_field(PhysData.γ3_gas(gases[2])),))
    inputs = Fields.GaussField(λ0=Λ0, τfwhm=ΤFWHM, energy=1e-6)
    linop, βfun!, β1, αfun = LinearOps.make_const_linop(grid, m, Λ0)
    Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.default(grid, Eω, m, linop, transform; gas=gases)
    output = Output.MemoryOutput(0, L, SAVEN, statsfun)
    Eω, grid, linop, transform, FT, output
end

#= `z` and `dz` are all the free-space cases record. `gpu/25-modal-error-stat` added a
   `Stats.default` for the radial and free-space transforms, so they *could* record the
   whole set; they do not, because every per-case tolerance in `tolerances.jl` comes from a
   one-ulp sensitivity study of the quantities actually compared, and adding a class of
   statistics to a case means re-measuring its tolerance and regenerating every baseline.
   `Stats.collect_stats` with no functions appends `Stats.zdz!`, which is enough for the
   fixed-step check that the step sequence really was imposed, and for the step-count check
   in the adaptive mode.

   `gpu/int-E2` considered adding the default set to the two radial cases only, which
   `gpu/25`'s review suggested, and did not: this file has to run unchanged on the baseline
   commits, and none of the three in use (`b641025b`, `782f55d1`, `fdf8dbe3`) has a
   `Stats.default` for a radial transform. Recording those statistics would make the gate
   impossible to run against any commit before `gpu/25`, which is the comparison the
   project's regression record is made of. It becomes possible once the oldest baseline in
   use is `gpu/int-E2` or later. =#
"`z` and `dz` only: the statistics a free-space geometry can record."
freestats(grid, Eω) = Stats.collect_stats(grid, Eω)

#= Free-space parameters, following test/test_freespace.jl but on smaller grids. =#
const R_FREE = 1.0e-3
const L_FREE = 0.15
const W0_FREE = 200e-6
const GAS_FREE = :Ar
const P_FREE = 1.0
const ΛLIMS_FREE = (400e-9, 2000e-9)
const TRANGE_FREE = 0.1e-12

"""
Pulse energy of `setup_radial_field_raman`. The Kerr-only free-space cases run at 1 pJ,
where a Raman term would be numerical dust; at 50 µJ, with `W0_FREE` and `ΤFWHM`, the Raman
polarisation changes the field at the end of the propagation by 9.3e-03 relative, so the
case measures the Raman response rather than a Kerr-only propagation. Measured against the
same case with the Raman response removed, at 5, 20, 50 and 100 µJ: 8.6e-04, 3.5e-03,
9.3e-03, 2.1e-02, with the adaptive step count (23) reproducible under a one-ulp
perturbation at all four. 50 µJ is 1.6 GW, comfortably below the critical power for
self-focusing in nitrogen.
"""
const RAMAN_FREE_ENERGY = 5e-5

"""
    setup_free(grid, zmax, sg, normfun, responses; gas, pres, energy)

Free-space propagation of a Gaussian beam focusing from `propz = -L_FREE`, on the
transverse grid `sg` (a radial grid or a `Grid.FreeGrid`) with the matching `normfun` and
nonlinear `responses`. `zmax` is the propagation length, which the grid no longer carries.
`gas`, `pres` and `energy` default to the argon/1 bar/1 pJ the Kerr-only cases use.
"""
function setup_free(grid, zmax, sg, normfun, responses;
                    gas=GAS_FREE, pres=P_FREE, energy=1e-12)
    nfunλ = PhysData.ref_index_fun(gas, pres)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    #= `thg` defaults to `true` for a `RealGrid` (required) and `false` for an `EnvGrid`,
       which matches the Kerr responses used below. =#
    linop = LinearOps.make_const_linop(grid, sg, nfun)
    dens0 = PhysData.density(gas, pres)
    densityfun = let dens0=dens0
        z -> dens0
    end
    inputs = Fields.GaussGaussField(;λ0=Λ0, τfwhm=ΤFWHM, energy,
                                     w0=W0_FREE, propz=-L_FREE)
    Eω, transform, FT = Luna.setup(grid, sg, densityfun, normfun, responses, inputs)
    output = Output.MemoryOutput(0, zmax, SAVEN, freestats(grid, Eω))
    Eω, grid, linop, transform, FT, output
end

function setup_radial_field()
    grid = makegrid(Grid.RealGrid, L_FREE, Λ0, ΛLIMS_FREE, TRANGE_FREE)
    q = radialgrid(R_FREE, 32)
    nfunλ = PhysData.ref_index_fun(GAS_FREE, P_FREE)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    setup_free(grid, L_FREE, q, NonlinearRHS.const_norm_radial(grid, q, nfun),
               (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS_FREE)),))
end

function setup_radial_env()
    grid = makegrid(Grid.EnvGrid, L_FREE, Λ0, ΛLIMS_FREE, TRANGE_FREE)
    q = radialgrid(R_FREE, 32)
    nfunλ = PhysData.ref_index_fun(GAS_FREE, P_FREE)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    setup_free(grid, L_FREE, q, NonlinearRHS.const_norm_radial(grid, q, nfun),
               (Nonlinear.Kerr_env(PhysData.γ3_gas(GAS_FREE)),))
end

"""
    setup_radial_field_raman()

Radial Kerr + Raman in nitrogen: the only multi-column Raman case in the matrix. The
batched Raman response owns `(2nt, npol, ncols)` buffers and does one pair of FFTs over the
whole transverse block, which a single-column (mode-averaged or GNLSE) case cannot
exercise. Nitrogen at 1 bar and [`RAMAN_FREE_ENERGY`](@ref), so that the Raman term is a
per-cent-level contribution rather than numerical dust; everything else matches the two
radial Kerr cases.
"""
function setup_radial_field_raman()
    gas, pres = :N2, 1.0
    grid = makegrid(Grid.RealGrid, L_FREE, Λ0, ΛLIMS_FREE, TRANGE_FREE)
    q = radialgrid(R_FREE, 32)
    nfunλ = PhysData.ref_index_fun(gas, pres)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    rr = Raman.raman_response(grid.to, gas)
    setup_free(grid, L_FREE, q, NonlinearRHS.const_norm_radial(grid, q, nfun),
               (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
                Nonlinear.RamanPolarField(grid.to, rr));
               gas, pres, energy=RAMAN_FREE_ENERGY)
end

function setup_free3d_env()
    grid = makegrid(Grid.EnvGrid, L_FREE, Λ0, ΛLIMS_FREE, TRANGE_FREE)
    sg = Grid.FreeGrid(R_FREE, 16, R_FREE, 8)
    nfunλ = PhysData.ref_index_fun(GAS_FREE, P_FREE)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    setup_free(grid, L_FREE, sg, NonlinearRHS.const_norm_free(grid, sg, nfun),
               (Nonlinear.Kerr_env(PhysData.γ3_gas(GAS_FREE)),))
end

#= Type I SHG in BBO on a 2-D Cartesian grid: the χ⁽²⁾ response with two polarisation
   components. Same setup as the "BBO SHG" testset in test/test_freespace.jl, which runs it
   both field-resolved and as an envelope; so do we, because `Chi2Field` and `Chi2Env` are
   separate implementations and GPU_PLAN.md section 6 makes both of them targets. =#
const BBO_THICKNESS = 30e-6
const BBO_θ = deg2rad(29.2) # type I phase-matching angle
const BBO_ϕ = deg2rad(30)
const BBO_λ0 = 800e-9
const BBO_τFWHM = 30e-15
const BBO_W0 = 20e-6
const BBO_ENERGY = 10e-9

"""
    setup_bbo(grid, response)

Type I SHG in BBO on a `Grid.Free2DGrid`, with `response` the χ⁽²⁾ response matching `grid`.
"""
function setup_bbo(grid, response)
    xgrid = Grid.Free2DGrid(4BBO_W0, 2^5)
    nfuns = PhysData.ref_index_fun_xy(:BBO, BBO_θ)
    linop = LinearOps.make_const_linop(grid, xgrid, nfuns)
    normfun = NonlinearRHS.const_norm_free2D(grid, xgrid, nfuns)
    densityfun = z -> 1 # unity density: we're considering a solid
    inputs = Fields.GaussGaussField(;λ0=BBO_λ0, τfwhm=BBO_τFWHM,
                                     energy=BBO_ENERGY/(sqrt(π/2)*BBO_W0), w0=BBO_W0)
    Eω, transform, FT = Luna.setup(grid, xgrid, densityfun, normfun, (response,), inputs)
    output = Output.MemoryOutput(0, BBO_THICKNESS, SAVEN, freestats(grid, Eω))
    Eω, grid, linop, transform, FT, output
end

function setup_bbo_field()
    grid = makegrid(Grid.RealGrid, BBO_THICKNESS, BBO_λ0, (250e-9, 2e-6), 120e-15)
    setup_bbo(grid, Nonlinear.Chi2Field(BBO_θ, BBO_ϕ, PhysData.χ2(:BBO)))
end

function setup_bbo_env()
    grid = makegrid(Grid.EnvGrid, BBO_THICKNESS, BBO_λ0, (250e-9, 2e-6), 120e-15; thg=true)
    setup_bbo(grid, Nonlinear.Chi2Env(BBO_θ, BBO_ϕ, PhysData.χ2(:BBO), grid.ω0, grid.to))
end

#= GNLSE: a second-order (N = 2) soliton, from examples/simple_interface/gnlse_sol.jl.
   `GNLSE_LENGTH = 0.1π τ₀²/|β₂|` is 0.2 soliton periods (z₀ = (π/2) τ₀²/|β₂|), i.e. the
   first compression is in the propagation. =#
const GNLSE_γ = 0.1
const GNLSE_β2 = -1e-26
const GNLSE_τ0 = 280e-15
const GNLSE_FR = 0.18
const GNLSE_LENGTH = 0.1π*GNLSE_τ0^2/abs(GNLSE_β2)

"""
    setup_gnlse(; raman, shock)

An N = 2 sech soliton through `prop_gnlse`. `raman` switches on the `fr`-weighted Raman
response and `shock` the self-steepening term, which are separate code paths in
`Interface.prop_gnlse_args`.
"""
function setup_gnlse(; raman, shock)
    N = 2.0
    P0 = N^2*abs(GNLSE_β2)/((1 - GNLSE_FR)*GNLSE_γ*GNLSE_τ0^2)
    Interface.prop_gnlse_args(
        GNLSE_γ, GNLSE_LENGTH, [0.0, 0.0, GNLSE_β2];
        λ0=835e-9, λlims=(450e-9, 8000e-9), trange=2e-12,
        τfwhm=(2*log(1 + sqrt(2)))*GNLSE_τ0, power=P0, pulseshape=:sech,
        raman, shock, fr=GNLSE_FR, shotnoise=false, saveN=SAVEN)
end

setup_gnlse_sech() = setup_gnlse(raman=false, shock=false)
setup_gnlse_raman_shock() = setup_gnlse(raman=true, shock=true)

"A linear taper of the core radius from `A_CAP` to `3A_CAP/4`."
taper(z) = A_CAP + (0.75A_CAP - A_CAP)*z/L_CAP

# ---------------------------------------------------------------------------------------
# Rectangular multimode case
# ---------------------------------------------------------------------------------------
"""
    A_RECT, B_RECT, L_RECT

Half-widths and length of the rectangular guide of [`setup_rect_modal`](@ref). `a > b`
deliberately: a rectangular guide is the only geometry in Luna with a Cartesian
`Modes.dimlimits`, and `a > b` is the case the Cartesian in-domain test of the adaptive
transverse integral (`NonlinearRHS._points!`) got wrong before `gpu/26-rectmode-fix`. The
aspect ratio 2.5 makes the dropped strip 30 % of the guide's area.
"""
const A_RECT = 100e-6

@doc (@doc A_RECT)
const B_RECT = 40e-6

@doc (@doc A_RECT)
const L_RECT = 0.03

"""
Multimode field propagation in a rectangular guide (`RectModes.RectMode`, argon at 5 bar,
silver cladding) with `a > b`, Kerr only.

The two modes are the `m = 1` and `m = 3` modes of the `x` index -- the coordinate the
guide is wide in, which is the one the in-domain test is applied to -- with the same `y`
index and the same polarisation, so that the Kerr product of the fundamental projects
onto the second mode. At 5 µJ over 3 cm the second mode reaches 7 % of the fundamental's
peak `|Eω|`, and the 20-step fixed run agrees with the adaptive one to 1.7e-3 in the
gate's metric (the maximum difference normalised per mode and per save; 2e-4 normalised
over the whole field), so the case is resolved at the fixed step size. There is no
plasma response: the matrix already has four ionising cases, and this one is here for the
transverse integral.

This is the matrix's only Cartesian transverse domain, and so the only case which reaches
the `:cartesian` branch of `NonlinearRHS._points!`; every other multimode case is polar,
where the test is unchanged.
"""
function setup_rect_modal()
    gas = :Ar
    pres = 5.0
    grid = makegrid(Grid.RealGrid, L_RECT, Λ0, ΛLIMS, TRANGE)
    modes = Tuple(RectModes.RectMode(A_RECT, B_RECT, gas, pres, :Ag; n=1, m=mx, pol=:x)
                  for mx in (1, 3))
    dens0 = PhysData.density(gas, pres)
    densityfun = let dens0=dens0
        z -> dens0
    end
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    inputs = Fields.GaussField(λ0=Λ0, τfwhm=ΤFWHM, energy=5e-6)
    linop = LinearOps.make_const_linop(grid, modes, Λ0)
    Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs, modes, :x;
                                   full=true)
    statsfun = Stats.default(grid, Eω, modes, linop, transform; gas=gas)
    output = Output.MemoryOutput(0, L_RECT, SAVEN, statsfun)
    Eω, grid, linop, transform, FT, output
end

# ---------------------------------------------------------------------------------------
# The case matrix
# ---------------------------------------------------------------------------------------
"""
    CASES

The regression case matrix: a `Vector{Case}`. See `test/regression/README.md`.
"""
const CASES = Case[
    capillary_case("modeavg_field_kerr", A_CAP, L_CAP, :He, 1.0;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   plasma=false),

    #= `thg=false` on a `RealGrid` selects `Nonlinear.Kerr_field_nothg`, which removes the
       third-harmonic term by way of the analytic signal and is what `prop_capillary` uses
       for a field run with `thg=false`. `Kerr_field` is everything else in the matrix. =#
    capillary_case("modeavg_field_nothg", A_CAP, L_CAP, :He, 1.0;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   plasma=false, thg=false),

    #= The four plasma cases are argon at 0.1 bar rather than helium at 1 bar, and at
       hundreds of µJ rather than 800 nJ. Review round 1 of gpu/13-plasma established
       from the baseline files that the original parameters ionise nothing at all --
       `maximum(stats/electrondensity)` was exactly 0.0 in all four cases and both modes,
       so the plasma response contributed nothing and the cases were blind to every
       change to it. These parameters give a peak ionised fraction of 0.03-0.4 %, and
       three of the four are the ones the review measured (gpu-13-plasma-1.md,
       "Recommended replacement parameters"). Everything else about the cases -- core
       radius, length, λ0, λlims, trange, τfwhm -- is unchanged.

       This one is 175 µJ rather than the review's 300 µJ. At 300 µJ the adaptive run
       has no reproducible step sequence: a one-ulp perturbation of the input takes it
       from 92 accepted steps to 98, which the gate reports as a `step count` failure,
       so no tolerance can be measured for its statistics. Measured on this branch
       (adaptive steps unperturbed / perturbed by one ulp, peak ionised fraction):
       300 µJ 92/98, 0.78 %; 250 µJ 77/86, 0.27 %; 200 µJ 48/46, 0.07 %;
       175 µJ 38/38, 0.03 %; 150 µJ 31/30, 0.01 %; 125 µJ 28/28, 0.003 %;
       100 µJ 25/25, 0.001 %. 175 µJ is the most strongly ionising energy whose step
       count is reproducible, and it is reproducible at ±1, ±2 and ±8 ulp and at 1e-14.
       The other three cases are stable at the review's energies. =#
    capillary_case("modeavg_field_plasma", A_CAP, L_CAP, :Ar, 0.1;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=175e-6,
                   plasma=true, PPT_options=NOCACHE),

    capillary_case("modeavg_field_raman", A_CAP, L_CAP, :N2, 0.5;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   raman=true, plasma=false),

    lowlevel_case("modeavg_field_mixture", 0.02, setup_mixture),

    #= ADK rather than PPT. `makeplasma!` takes the model as a `Symbol`, and the default for
       a noble gas is `:PPT`, so nothing else in the matrix reaches `Ionisation.IonRateADK`. =#
    capillary_case("modeavg_field_adk", A_CAP, L_CAP, :Ar, 0.1;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=300e-6,
                   plasma=:ADK),

    #= Elliptically polarised input: two polarisation components, so the vector forms of the
       Kerr and plasma responses and a two-mode `TransModal`. =#
    capillary_case("modeavg_field_vector", A_CAP, L_CAP, :Ar, 0.1;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=150e-6,
                   polarisation=0.5, plasma=true, PPT_options=NOCACHE),

    capillary_case("modeavg_env_kerr", A_CAP, L_CAP, :He, 1.0;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   envelope=true),

    capillary_case("modeavg_env_raman", A_CAP, L_CAP, :N2, 0.5;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   envelope=true, raman=true),

    #= `thg=true` on an `EnvGrid` selects `Nonlinear.Kerr_env_thg`, which carries the
       carrier-phase array and is a different response from `Kerr_env`. =#
    capillary_case("modeavg_env_thg", A_CAP, L_CAP, :He, 1.0;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   envelope=true, thg=true),

    lowlevel_case("gnlse_sech", GNLSE_LENGTH, setup_gnlse_sech),

    lowlevel_case("gnlse_raman_shock", GNLSE_LENGTH, setup_gnlse_raman_shock),

    capillary_case("multimode_field_plasma", A_CAP, L_CAP, :Ar, 0.1;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=150e-6,
                   modes=4, plasma=true, PPT_options=NOCACHE),

    lowlevel_case("rect_modal_field", L_RECT, setup_rect_modal),

    lowlevel_case("radial_field_kerr", L_FREE, setup_radial_field),
    lowlevel_case("radial_env_kerr", L_FREE, setup_radial_env),
    lowlevel_case("radial_field_raman", L_FREE, setup_radial_field_raman),
    lowlevel_case("free3d_env_kerr", L_FREE, setup_free3d_env),
    lowlevel_case("free2d_field_chi2", BBO_THICKNESS, setup_bbo_field),
    lowlevel_case("free2d_env_chi2", BBO_THICKNESS, setup_bbo_env),

    capillary_case("gradient_field_kerr", A_CAP, L_CAP, :He, (2.0, 0.1);
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   plasma=false),

    capillary_case("taper_field_kerr", taper, L_CAP, :He, 1.0;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   plasma=false),

    capillary_case("modeavg_field_legacy", A_CAP, L_CAP, :He, 1.0;
                   λ0=Λ0, λlims=ΛLIMS, trange=TRANGE, τfwhm=ΤFWHM, energy=200e-9,
                   plasma=false, boundary=:legacy),
]

"The two run modes every case is run in."
const MODES = (:fixed, :adaptive)

"Look up a case by name."
function getcase(name::AbstractString)
    idx = findfirst(c -> c.name == name, CASES)
    isnothing(idx) && error("No regression case called $name")
    CASES[idx]
end

end
