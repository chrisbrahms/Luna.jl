#= Whole-propagation cases for the CPU-versus-GPU comparison (`run.jl`).

   Each case is a physically meaningful run, not a fixed number of steps: the time that
   matters to a user is the wall time of the propagation, and that includes the step count
   the adaptive integrator needs, which depends on the physics. The strength of every case
   is set against the critical power for self-focusing (`Tools.Pcr`), so that the
   transverse dynamics are really exercised:

   - `modeavg`: mode-averaged field-resolved Kerr + PPT plasma in an argon-filled capillary
     (the most common Luna run; one transverse column, so a GPU is not expected to win).
   - `modeavg_long`: the same over 20 cm in an 8 ps window, the largest single-column grid.
   - `modal`: 8 HE₁ₘ modes, 125 µm core radius, 1 bar argon, 30 fs at 0.95 P_cr, Kerr +
     plasma. This is the hard case of `benchmark/threads/modal_accuracy.jl`: the higher
     modes are driven by self-focusing, and the adaptive transverse rule grows past 100
     points. `modal_integral=:adaptive` (host only) and `:fixed` with `nr` nodes.
   - `radial`: field-resolved Kerr + PPT plasma in 1 bar argon, a 30 fs pulse at
     0.8 P_cr focused (w0 = 100 µm, focus at the middle of a 10 cm cell), so self-focusing
     and ionisation both act near the focus; `N` radial points.
   - `free3d`: envelope Kerr in 1 bar argon on an `N × N` Cartesian grid, a 30 fs pulse at
     0.8 P_cr focused (w0 = 200 µm, focus at the middle of a 30 cm cell).

   `scale` multiplies the propagation length (1 is the full case; `run.jl` uses a small
   one to compile). Every case returns `(setup_s, run_s, output)`. =#
import Luna: Grid, Fields, LinearOps, Nonlinear, NonlinearRHS, Output, PhysData, Stats,
             Ionisation, Tools, Interface
import Logging

const GAS = :Ar
const PRES = 1.0
const λ0 = 800e-9
const τFWHM = 30e-15

"Critical power for self-focusing of `GAS` at `PRES` and `λ0`, in W."
function pcrit()
    ω = 2π*PhysData.c/λ0
    _, n0, n2 = Tools.getN0n0n2(ω, GAS; P=PRES)
    Tools.Pcr(ω, n0, n2)
end

"Energy of a Gaussian pulse of FWHM `τFWHM` with peak power `frac * pcrit()`."
energy_at(frac) = frac*pcrit()*τFWHM/0.94

"""
Number of statistics points over a case's length, from `STATS_N` (0, the default: every
accepted step, which is Luna's default). The default statistics of the multimode and
free-space runs include members with no device form, so on a GPU they copy the field to
the host every step; `STATS_N > 0` is the `stats_period` a GPU run can use instead.
"""
const STATSN = parse(Int, get(ENV, "STATS_N", "0"))
statskw(L) = STATSN > 0 ? (stats_period=L/STATSN,) : (;)

#= Forwards warnings and errors to the console and picks the step count out of the
   integrator's "Propagation finished in ..., N steps" message. N counts every attempted
   step, rejected ones included, so it is the amount of work; the statistics only see the
   accepted ones, and with `STATS_N > 0` not even all of those. =#
struct StepLogger <: Logging.AbstractLogger
    steps::Base.RefValue{Int}
    inner::Logging.ConsoleLogger
end
StepLogger() = StepLogger(Ref(0), Logging.ConsoleLogger(stderr, Logging.Warn))
Logging.min_enabled_level(::StepLogger) = Logging.Info
Logging.shouldlog(::StepLogger, args...) = true
Logging.catch_exceptions(::StepLogger) = false
function Logging.handle_message(l::StepLogger, level, msg, args...; kw...)
    m = match(r"^Propagation finished in .*, (\d+) steps", string(msg))
    isnothing(m) || (l.steps[] = parse(Int, m[1]))
    level >= Logging.Warn && Logging.handle_message(l.inner, level, msg, args...; kw...)
    nothing
end

"Run `f()` under a [`StepLogger`](@ref); return its result and the step count."
function quiet(f)
    l = StepLogger()
    r = Logging.with_logger(f, l)
    r, l.steps[]
end

function capillary(L, modes; device, precision, kw...)
    (ts, tr, output), n = quiet() do
        ts = @elapsed begin
            Eω, grid, linop, transform, FT, output = Interface.prop_capillary_args(
                125e-6, L, GAS, PRES; λ0, τfwhm=τFWHM, λlims=(150e-9, 4e-6),
                shotnoise=false, modes, device, precision, statskw(L)..., kw...)
        end
        tr = @elapsed Luna.run(Eω, grid, linop, transform, FT, output; zmax=L,
                                status_period=1e9)
        ts, tr, output
    end
    ts, tr, n, output
end

modeavg(; scale=1.0, device, precision) =
    capillary(0.5scale, :HE11; device, precision, energy=150e-6, trange=1e-12)

"The long-window mode-averaged case: an 8 ps window (about 2^16 time points)."
modeavg_long(; scale=1.0, device, precision) =
    capillary(0.2scale, :HE11; device, precision, energy=150e-6, trange=8e-12)

modal(rule, nr; scale=1.0, device, precision) =
    capillary(0.3scale, 8; device, precision, energy=energy_at(0.95), trange=0.5e-12,
              modal_integral=rule, (rule == :fixed ? (modal_nr=nr,) : (;))...)

"""
    freespace(sg, grid, L, w0, frac, responses; device, precision)

Gaussian beam focused to `w0` at `L/2`, with peak power `frac * pcrit()`, propagated over
`L` on the transverse grid `sg` (a `RadialGrid` or `FreeGrid`).
"""
function freespace(sg, grid, L, w0, frac, responses; device, precision)
    (ts, tr, output), n = quiet() do
        nfunλ = PhysData.ref_index_fun(GAS, PRES)
        nfun = (λ; z=0.0) -> nfunλ(λ)
        dens0 = PhysData.density(GAS, PRES)
        densityfun = z -> dens0
        normfun = sg isa Grid.RadialGrid ? NonlinearRHS.const_norm_radial(grid, sg, nfun) :
                                           NonlinearRHS.const_norm_free(grid, sg, nfun)
        inputs = Fields.GaussGaussField(; λ0, τfwhm=τFWHM, energy=energy_at(frac), w0,
                                         propz=-L/2)
        local Eω, transform, FT, linop, output
        ts = @elapsed begin
            linop = LinearOps.make_const_linop(grid, sg, nfun)
            Eω, transform, FT = Luna.setup(grid, sg, densityfun, normfun, responses, inputs;
                                           device, precision)
            statsfun = Stats.default(grid, Eω, linop, transform; gas=GAS)
            STATSN > 0 && (statsfun = Output.PeriodicStats(statsfun, L/STATSN))
            output = Output.MemoryOutput(0, L, 11, statsfun)
        end
        tr = @elapsed Luna.run(Eω, grid, linop, transform, FT, output; zmax=L,
                               status_period=1e9)
        ts, tr, output
    end
    ts, tr, n, output
end

function radial(N; scale=1.0, device, precision)
    L = 0.1scale
    grid = Grid.RealGrid(λ0, (200e-9, 4e-6), 0.25e-12)
    ionrate = Ionisation.IonRatePPTCached(GAS, λ0)
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS)),
                 Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate,
                                          PhysData.ionisation_potential(GAS)))
    freespace(Grid.RadialGrid(1.5e-3, N), grid, L, 100e-6, 0.8, responses;
              device, precision)
end

function free3d(N; scale=1.0, device, precision)
    L = 0.3scale
    grid = Grid.EnvGrid(λ0, (400e-9, 2000e-9), 0.2e-12)
    responses = (Nonlinear.Kerr_env(PhysData.γ3_gas(GAS)),)
    freespace(Grid.FreeGrid(1.2e-3, N), grid, L, 200e-6, 0.8, responses; device, precision)
end
