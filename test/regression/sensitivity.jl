#= One-ulp sensitivity study for the regression cases.

   Usage:
       julia --project=<worktree> -t 1 test/regression/sensitivity.jl [case ...]

   For every case and every mode, run it twice: once unperturbed and once with the initial
   frequency-domain field multiplied by `1 + eps()`, i.e. perturbed by one unit in the last
   place. The reported number is the same metric the gate uses,
   `maximum(abs, Δ)/maximum(abs, reference)`, maximised over `Eω`, `z` and every statistic.

   It measures how much a rounding-level change at the input is amplified by the
   propagation, which is the floor below which no reordering of floating-point operations
   can be expected to stay. The per-case tolerances in `tolerances.jl` are 100x these
   numbers, with a floor of 1e-12.

   For each case and mode it prints the `NREPORT` largest quantities, not just the worst
   one, because in the adaptive mode the number is usually set by `stats/dz`: the step-size
   controller's accept/reject decisions respond to a one-ulp change much more strongly than
   the field does, and the tolerance the controller forces on a case says nothing about the
   tolerance on `Eω`.
=#
using Luna
import LinearAlgebra
import Printf: @printf, @sprintf

Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Luna.set_fftw_wisdom(false)

include(joinpath(@__DIR__, "cases.jl"))
include(joinpath(@__DIR__, "compare.jl"))
using .RegressionCases
using .RegressionCompare

const ONLY = ARGS

"How many of the largest quantities to list per case and mode."
const NREPORT = 4

results = Dict{String, Dict{Symbol, Float64}}()

@printf("%-24s %-9s %12s  %s\n", "case", "mode", "sensitivity", "largest quantities")
@printf("%s\n", "-"^78)
for case in RegressionCases.CASES
    isempty(ONLY) || case.name in ONLY || continue
    results[case.name] = Dict{Symbol, Float64}()
    for mode in RegressionCases.MODES
        ref = rundict(runcase(case, mode))
        new = rundict(runcase(case, mode; perturb=eps()))
        diffs = sort(compare(ref, new); by=d -> -d.value)
        results[case.name][mode] = isempty(diffs) ? 0.0 : diffs[1].value
        top = join([@sprintf("%s %.2e", d.what, d.value)
                    for d in diffs[1:min(NREPORT, length(diffs))]], ", ")
        @printf("%-24s %-9s %12.3e  %s\n",
                case.name, mode, results[case.name][mode], top)
        flush(stdout)
    end
end

@printf("\n# Paste into tolerances.jl (100x the measured sensitivity, floor 1e-12):\n")
@printf("const TOLERANCES = Dict(\n")
for case in RegressionCases.CASES
    haskey(results, case.name) || continue
    r = results[case.name]
    tols = [max(1e-12, 100*r[mode]) for mode in RegressionCases.MODES]
    entries = join([@sprintf("%s => %.1e", repr(m), t)
                    for (m, t) in zip(RegressionCases.MODES, tols)], ", ")
    @printf("    \"%s\" => Dict(%s),\n", case.name, entries)
end
@printf(")\n")
