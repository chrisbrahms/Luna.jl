#= Run every regression case in every mode and save the results as HDF5.

   Usage:
       julia --project=<worktree> -t 1 test/regression/run_cases.jl <outdir> [case ...]

   This script is what `test/regression/generate.jl` copies into a worktree of the base
   commit and executes there, so that the baseline is produced by that commit's Luna but
   by the *current* case definitions. It therefore has to cope with a Luna that does not
   have `Luna.set_fftw_wisdom` yet.
=#
using Luna
import FFTW
import LinearAlgebra
import Printf: @printf

length(ARGS) >= 1 || error("usage: run_cases.jl <outdir> [case ...]")
const OUTDIR = abspath(ARGS[1])
const ONLY = ARGS[2:end]

#= Deterministic FFTW/BLAS environment. Without this the baseline and the run being gated
   can get different FFT plans (and hence a different order of operations) from the shared
   wisdom file, and different reduction orders from threaded BLAS. =#
Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
if isdefined(Luna, :set_fftw_wisdom)
    Luna.set_fftw_wisdom(false)
else
    # Older commits have no switch: neutralise the wisdom cache by redefining it away.
    @eval Luna.Utils loadFFTwisdom() = nothing
    @eval Luna.Utils saveFFTwisdom() = nothing
    FFTW.forget_wisdom()
end

include(joinpath(@__DIR__, "cases.jl"))
include(joinpath(@__DIR__, "compare.jl"))
using .RegressionCases
using .RegressionCompare

isdir(OUTDIR) || mkpath(OUTDIR)

@printf("Writing regression results to %s\n", OUTDIR)
for case in RegressionCases.CASES
    isempty(ONLY) || case.name in ONLY || continue
    fpath = joinpath(OUTDIR, case.name * ".h5")
    isfile(fpath) && rm(fpath)
    for mode in RegressionCases.MODES
        t = @elapsed output = runcase(case, mode)
        saverun(OUTDIR, case.name, mode, output)
        nsteps = haskey(output["stats"], "z") ? length(output["stats"]["z"]) : -1
        @printf("  %-24s %-9s %6.1f s  %4d steps\n", case.name, mode, t, nsteps)
        flush(stdout)
    end
end
@printf("Done.\n")
