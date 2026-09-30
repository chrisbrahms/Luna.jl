#= The regression gate.

   Runs every case in `test/regression/cases.jl`, in both run modes, and compares the
   result against a baseline generated from an earlier commit by
   `test/regression/generate.jl`. It always prints the observed maximum differences, even
   when it passes.

   Usage:
       julia --project=<worktree> -t 1 test/test_regression.jl

   The baseline commit is `ENV["LUNA_REGRESSION_BASE"]` if set, otherwise the merge-base of
   `HEAD` with `evanescent` (`ENV["LUNA_REGRESSION_BRANCH"]` overrides that branch). The
   baseline directory is `ENV["LUNA_REGRESSION_DIR"]` if set, otherwise
   `joinpath(Luna.Utils.cachedir(), "regression", <full sha>)`.

   Two tolerances per case and mode, one for `Eω` and one for the statistics; see
   `test/regression/tolerances.jl` and `test/regression/README.md`.

   This file is deliberately NOT part of `test/runtests.jl`: it needs a baseline, which is
   never committed.
=#
using Luna
import FFTW
import LinearAlgebra
import Printf: @printf, @sprintf
import Test: @test, @testset

#= Deterministic FFTW/BLAS environment: the same one the baseline was generated in. A
   different FFT plan or a different BLAS thread count changes the order of the
   floating-point operations and so the last bits of the result. =#
Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
# Luna.run sets its own BLAS count unless told otherwise (older commits have no switch)
isdefined(Luna, :set_blas_threads) && Luna.set_blas_threads(1)
Luna.set_fftw_wisdom(false)

include(joinpath(@__DIR__, "regression", "cases.jl"))
include(joinpath(@__DIR__, "regression", "compare.jl"))
include(joinpath(@__DIR__, "regression", "tolerances.jl"))
using .RegressionCases
using .RegressionCompare
using .RegressionTolerances

"""
    basecommit()

The commit the baseline is expected to come from: `ENV["LUNA_REGRESSION_BASE"]` if set,
otherwise the merge-base of `HEAD` with `ENV["LUNA_REGRESSION_BRANCH"]` (default
`evanescent`). Returned as a full SHA.
"""
function basecommit()
    repo = dirname(@__DIR__)
    if haskey(ENV, "LUNA_REGRESSION_BASE")
        return readchomp(`git -C $repo rev-parse $(ENV["LUNA_REGRESSION_BASE"])`)
    end
    branch = get(ENV, "LUNA_REGRESSION_BRANCH", "evanescent")
    readchomp(`git -C $repo merge-base HEAD $branch`)
end

"""
    selected() -> (cases, description)

The cases to run, as a `Vector{RegressionCases.Case}`: `ENV["LUNA_REGRESSION_ONLY"]` (a
comma-separated list of case names) if set, minus `ENV["LUNA_REGRESSION_SKIP"]`, together
with a one-line description of the selection for the header.

Both default to empty, i.e. the whole matrix. They exist because a baseline generated
before a case was added has no file for it, so the gate would report it as a failure to
load rather than as the missing baseline it is: a branch which adds a case runs the older
baselines with that case skipped, and the new case against a baseline of its own. Neither
variable is a way to make a failing case pass: the header names the lists exactly as they
were given (`only: ...`, `skipped: ...`), an unknown name is an error, and the filtering
happens before the comparison, so nothing a case does can change it.
"""
function selected()
    names(k) = haskey(ENV, k) ? String.(split(ENV[k], ',')) : String[]
    only, skip = names("LUNA_REGRESSION_ONLY"), names("LUNA_REGRESSION_SKIP")
    for n in vcat(only, skip)
        any(c -> c.name == n, RegressionCases.CASES) ||
            error("No regression case called $n (LUNA_REGRESSION_ONLY/SKIP)")
    end
    cases = isempty(only) ? RegressionCases.CASES :
            filter(c -> c.name in only, RegressionCases.CASES)
    cases = filter(c -> !(c.name in skip), cases)
    #= What the header says is what the *user asked for*, not the 21 names which survived
       it: a list of everything which ran leaves the reader to spot what is missing. =#
    parts = String[]
    isempty(only) || push!(parts, "only: " * join(only, ", "))
    isempty(skip) || push!(parts, "skipped: " * join(skip, ", "))
    (cases, isempty(parts) ? "" : "  (" * join(parts, "; ") * ")")
end

const BASE = basecommit()
const CASES, SELECTION = selected()
const BASEDIR = get(ENV, "LUNA_REGRESSION_DIR",
                    joinpath(Luna.Utils.cachedir(), "regression", BASE))

isdir(BASEDIR) || error(
    "No regression baseline at $BASEDIR.\n" *
    "Generate it with: julia --project=$(dirname(@__DIR__)) -t 1 " *
    "test/regression/generate.jl $BASE")

@printf("Regression gate\n")
@printf("  baseline commit: %s\n", BASE)
@printf("  baseline dir:    %s\n", BASEDIR)
@printf("  cases:           %d of %d%s\n", length(CASES), length(RegressionCases.CASES),
        SELECTION)
@printf("\n%-24s %-9s %11s %10s %11s %10s  %s\n",
        "case", "mode", "Eω diff", "Eω tol", "stats diff", "stats tol", "worst quantity")
@printf("%s\n", "-"^112)

worstoverall = 0.0

@testset "regression" begin
@testset "$(case.name)" for case in CASES
    baseline = try
        loadcase(BASEDIR, case.name)
    catch err
        @error "Could not load the baseline for $(case.name)" exception=err
        nothing
    end
    if isnothing(baseline)
        @test false
        continue
    end
    @testset "$mode" for mode in RegressionCases.MODES
        tols = Dict(c => tolerance(case.name, mode, c) for c in RegressionCompare.CLASSES)
        if !haskey(baseline, string(mode))
            @printf("%-24s %-9s %11s %10.1e %11s %10.1e  %s\n",
                    case.name, mode, "MISSING", tols[:Eω], "MISSING", tols[:stats],
                    "not in the baseline")
            @test false
            continue
        end
        new = rundict(runcase(case, mode))
        diffs = compare(baseline[string(mode)], new; skip=skipstats(mode))
        wE = RegressionCompare.worst(diffs, :Eω)
        wS = RegressionCompare.worst(diffs, :stats)
        w = RegressionCompare.worst(diffs)
        global worstoverall = max(worstoverall, w.value)
        @printf("%-24s %-9s %11.3e %10.1e %11.3e %10.1e  %s\n",
                case.name, mode, wE.value, tols[:Eω], wS.value, tols[:stats],
                RegressionCompare.describe(w))
        flush(stdout)
        for d in diffs
            tol = tols[d.class]
            ok = d.value <= tol
            ok || @printf("    FAIL %-44s %11.3e > %.1e\n",
                          RegressionCompare.describe(d), d.value, tol)
            @test ok
        end
        flush(stdout)
    end
end
end

@printf("%s\n", "-"^112)
@printf("largest difference over all cases and modes: %.3e\n", worstoverall)
