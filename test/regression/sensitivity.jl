#= One-ulp sensitivity study for the regression cases.

   Usage:
       julia --project=<worktree> -t 1 test/regression/sensitivity.jl [case ...]

   For every case and every mode, run it twice: once unperturbed and once with the initial
   frequency-domain field multiplied by `1 + eps()`, i.e. perturbed by one unit in the last
   place. The reported numbers are the gate's metric,
   `maximum(abs, Δ)/maximum(abs, reference)`, maximised separately over each tolerance
   class: `Eω` (the field, normalised per component and per save) and `stats` (the save grid
   and every compared statistic).

   It measures how much a rounding-level change at the input is amplified by the
   propagation, which is the floor below which no reordering of floating-point operations
   can be expected to stay. The per-case tolerances in `tolerances.jl` are 100x these
   numbers, with a floor of 1e-12.

   The comparison is the gate's, exclusions included (`RegressionCompare.skipstats`), so the
   `TOLERANCES` literal printed at the end can be pasted into `tolerances.jl` unchanged.

   For each case and mode it also prints the `NREPORT` largest quantities, so that a number
   well above 1e-12 can be attributed.
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

const CLASSES = RegressionCompare.CLASSES

# case name => mode => class => sensitivity
results = Dict{String, Dict{Symbol, Dict{Symbol, Float64}}}()

@printf("%-24s %-9s %11s %11s  %s\n",
        "case", "mode", "Eω", "stats", "largest quantities")
@printf("%s\n", "-"^104)
for case in RegressionCases.CASES
    isempty(ONLY) || case.name in ONLY || continue
    results[case.name] = Dict{Symbol, Dict{Symbol, Float64}}()
    for mode in RegressionCases.MODES
        ref = rundict(runcase(case, mode))
        new = rundict(runcase(case, mode; perturb=eps()))
        diffs = sort(compare(ref, new; skip=skipstats(mode)); by=d -> -d.value)
        results[case.name][mode] = Dict(
            c => RegressionCompare.worst(diffs, c).value for c in CLASSES)
        top = join([@sprintf("%s %.2e", RegressionCompare.describe(d), d.value)
                    for d in diffs[1:min(NREPORT, length(diffs))]], ", ")
        @printf("%-24s %-9s %11.3e %11.3e  %s\n", case.name, mode,
                results[case.name][mode][:Eω], results[case.name][mode][:stats], top)
        flush(stdout)
    end
end

"100x the measured sensitivity, floored at 1e-12."
tol(x) = max(1e-12, 100x)

@printf("\n# Paste into tolerances.jl (100x the measured sensitivity, floor 1e-12):\n")
@printf("const TOLERANCES = Dict{String, Dict{Symbol, Dict{Symbol, Float64}}}(\n")
for case in RegressionCases.CASES
    haskey(results, case.name) || continue
    @printf("    \"%s\" => Dict(\n", case.name)
    for (i, mode) in enumerate(RegressionCases.MODES)
        r = results[case.name][mode]
        last = i == length(RegressionCases.MODES)
        @printf("        %-10s => Dict(:Eω => %.1e, :stats => %.1e)%s\n",
                repr(mode), tol(r[:Eω]), tol(r[:stats]), last ? ")," : ",")
    end
end
@printf(")\n")
