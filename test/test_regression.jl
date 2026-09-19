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

   This file is deliberately NOT part of `test/runtests.jl`: it needs a baseline, which is
   never committed. See `test/regression/README.md`.
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

const BASE = basecommit()
const BASEDIR = get(ENV, "LUNA_REGRESSION_DIR",
                    joinpath(Luna.Utils.cachedir(), "regression", BASE))

isdir(BASEDIR) || error(
    "No regression baseline at $BASEDIR.\n" *
    "Generate it with: julia --project=$(dirname(@__DIR__)) -t 1 " *
    "test/regression/generate.jl $BASE")

@printf("Regression gate\n")
@printf("  baseline commit: %s\n", BASE)
@printf("  baseline dir:    %s\n", BASEDIR)
@printf("\n%-24s %-9s %12s %12s  %s\n",
        "case", "mode", "max diff", "tolerance", "worst quantity")
@printf("%s\n", "-"^80)

worstoverall = 0.0

@testset "regression" begin
@testset "$(case.name)" for case in RegressionCases.CASES
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
        tol = tolerance(case.name, mode)
        if !haskey(baseline, string(mode))
            @printf("%-24s %-9s %12s %12.1e  %s\n",
                    case.name, mode, "MISSING", tol, "not in the baseline")
            @test false
            continue
        end
        new = rundict(runcase(case, mode))
        diffs = compare(baseline[string(mode)], new)
        w = RegressionCompare.worst(diffs)
        global worstoverall = max(worstoverall, isfinite(w.value) ? w.value : Inf)
        @printf("%-24s %-9s %12.3e %12.1e  %s%s\n", case.name, mode, w.value, tol, w.what,
                isempty(w.note) ? "" : " ($(w.note))")
        flush(stdout)
        for d in diffs
            ok = d.value <= tol
            ok || @printf("    FAIL %-32s %12.3e > %.1e %s\n",
                          d.what, d.value, tol, d.note)
            @test ok
        end
        flush(stdout)
    end
end
end

@printf("%s\n", "-"^80)
@printf("largest difference over all cases and modes: %.3e\n", worstoverall)
