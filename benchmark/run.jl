#= Per-case CPU benchmarks for the regression case matrix.

   Usage:
       julia --project=benchmark -t 1 benchmark/run.jl [case ...]

   Set the environment up once from the repository root with

       julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

   For every case in `test/regression/cases.jl` this times three things:

   - `rhs`:  one evaluation of the nonlinear right-hand side, `transform(nl, Eω, z)`.
   - `step`: one `RK45.step!` of the preconditioned stepper the propagation uses, which is
             six RHS evaluations plus the linear propagator and the stage combines. A fresh
             stepper is built before every sample, so consecutive samples do the same work.
   - `prop`: a fixed-step propagation of `RegressionCases.NSTEPS` steps through `Luna.run`,
             including setup, the absorbing boundaries and the output.

   `rhs` and `step` are BenchmarkTools minima; `prop` is the minimum of `PROP_SAMPLES`
   whole runs, because a propagation is too slow to sample properly.

   Numbers from this script are recorded in `PR_00-harness.md`. Run it with one Julia
   thread, one FFTW thread and one BLAS thread, as it sets up below, or the numbers are not
   comparable between runs.
=#
using Luna
import Luna: RK45, Boundaries
import LinearAlgebra
import BenchmarkTools: @benchmarkable, run as brun, minimum as bminimum
import Printf: @printf, @sprintf

Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Luna.set_fftw_wisdom(false)

include(joinpath(dirname(@__DIR__), "test", "regression", "cases.jl"))
using .RegressionCases

const ONLY = ARGS

"How many whole propagations to time per case."
const PROP_SAMPLES = 3

"Time budget in seconds for each `rhs` and `step` benchmark."
const BUDGET = 2.0

"""
    prepare(case)

Build `case` and return `(Eω, linop, transform, dz)` with the absorbing boundaries already
folded into `linop`, as `Luna.run` does, and `dz` the fixed step the case is timed at.
"""
function prepare(case)
    RegressionCases.quiet() do
        Eω, grid, linop, transform, FT, output = case.prepare()
        dz = case.zmax/RegressionCases.NSTEPS
        Et = FT \ Eω
        #= `case.runkwargs` carries the same four keys `Interface.boundary_kwargs` forwards
           to `Luna.run`. `:boundary` is positional to `Boundaries.setup` and the other three
           are keywords under different names, so they are mapped rather than splatted. =#
        kw = case.runkwargs
        #= Through the shim in `cases.jl`: `gpu/01-zmax` gave `Boundaries.setup` a `zmax`
           positional argument after `z0`, in place of the `grid.zmax` it used to read. =#
        absorber = RegressionCases.absorber_setup(
            get(kw, :boundary, :rate),
            grid, transform, linop, Et, FT, output, 0.0, case.zmax, dz, dz;
            N=get(kw, :boundary_N, Boundaries.DEFAULT_N),
            ℓ=get(kw, :boundary_length, nothing),
            collar=get(kw, :tcollar, Boundaries.DEFAULT_TCOLLAR))
        (Eω, absorber.linop, transform, min(dz, absorber.max_dz))
    end
end

"""
    rhstime(case)

`(n, t)`: the number of elements in the state vector and the minimum time in seconds for
one evaluation of the nonlinear right-hand side.
"""
function rhstime(case)
    Eω, _, transform, _ = prepare(case)
    nl = similar(Eω)
    b = @benchmarkable $transform($nl, $Eω, 0.0) seconds=BUDGET
    (length(Eω), bminimum(brun(b)).time/1e9)
end

"""
    steptime(case)

Minimum time in seconds for one `RK45.step!` of the preconditioned stepper.
"""
function steptime(case)
    Eω, linop, transform, dz = prepare(case)
    #= A fresh stepper per sample: `step!` advances z and the step size, and for a
       z-dependent linop repeated stepping would walk off the end of the fibre. =#
    b = @benchmarkable RK45.step!(s) setup=(
            s = RK45.PreconStepper($transform, $linop, $Eω, 0.0, $dz;
                                   max_dt=$dz, min_dt=$dz)
        ) evals=1 seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

"""
    proptime(case)

Minimum wall-clock time in seconds of `PROP_SAMPLES` fixed-step propagations of `case`,
including setup and output.
"""
proptime(case) = minimum(_ -> (@elapsed runcase(case, :fixed)), 1:PROP_SAMPLES)

"Format a time in seconds with a unit that keeps it readable."
fmttime(t) = t >= 1    ? @sprintf("%7.3f s ", t)   :
             t >= 1e-3 ? @sprintf("%7.3f ms", 1e3t) :
             t >= 1e-6 ? @sprintf("%7.3f µs", 1e6t) :
                         @sprintf("%7.3f ns", 1e9t)

@printf("Luna CPU benchmarks: 1 Julia thread, 1 FFTW thread, 1 BLAS thread, :estimate\n")
@printf("steps per propagation: %d\n\n", RegressionCases.NSTEPS)
@printf("%-24s %10s %12s %12s %12s\n", "case", "state", "rhs", "step", "prop")
@printf("%s\n", "-"^74)
for case in RegressionCases.CASES
    isempty(ONLY) || case.name in ONLY || continue
    n, tr = rhstime(case)
    ts = steptime(case)
    tp = proptime(case)
    @printf("%-24s %10d %12s %12s %12s\n", case.name, n,
            fmttime(tr), fmttime(ts), fmttime(tp))
    flush(stdout)
end
