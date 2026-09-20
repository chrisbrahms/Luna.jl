#= Generate regression baselines from an arbitrary commit.

   Usage:
       julia --project=<worktree> -t 1 test/regression/generate.jl <commit> [outdir]

   Creates (or reuses) a git worktree of `<commit>`, copies this worktree's `Manifest.toml`
   and the *current* regression scripts into it, and runs `run_cases.jl` there in a fresh
   single-threaded process. The baseline is therefore produced by `<commit>`'s Luna but by
   the case definitions on the branch being tested.

   Baselines are written to `outdir`, by default
   `joinpath(Luna.Utils.cachedir(), "regression", <full sha of commit>)`, and are never
   committed (`*.h5` is gitignored).

   The worktrees are created under `$LUNA_REGRESSION_WORKTREES` if that is set, otherwise
   next to this worktree in `../baselines`. They are left in place so that the next
   baseline generation from the same commit does not have to precompile Luna again; remove
   them with `git worktree remove <path>`.
=#
using Luna
import Printf: @printf

length(ARGS) >= 1 || error("usage: generate.jl <commit> [outdir]")

const REPO = dirname(dirname(@__DIR__))
const COMMIT = ARGS[1]

"Full SHA of `COMMIT`."
const SHA = readchomp(`git -C $REPO rev-parse $COMMIT`)

const OUTDIR = length(ARGS) >= 2 ? abspath(ARGS[2]) :
    joinpath(Luna.Utils.cachedir(), "regression", SHA)

const WTBASE = get(ENV, "LUNA_REGRESSION_WORKTREES",
                   joinpath(dirname(REPO), "baselines"))
const WT = joinpath(WTBASE, SHA[1:10])

@printf("Baseline commit: %s\n", SHA)
@printf("Worktree:        %s\n", WT)
@printf("Output:          %s\n", OUTDIR)

if isdir(WT)
    have = readchomp(`git -C $WT rev-parse HEAD`)
    have == SHA || error("$WT is at $have, not $SHA. Remove it with `git worktree remove`.")
    @printf("Reusing existing worktree.\n")
else
    isdir(WTBASE) || mkpath(WTBASE)
    run(`git -C $REPO worktree add --detach $WT $SHA`)
end

# The Manifest is gitignored, so the worktree needs this one to get the same package
# versions (and the same `dev`ed Hankel) as the branch being tested.
cp(joinpath(REPO, "Manifest.toml"), joinpath(WT, "Manifest.toml"); force=true)

# The case definitions and the runner come from *this* branch, not from the base commit.
dest = joinpath(WT, "test", "regression")
isdir(dest) || mkpath(dest)
for f in ("cases.jl", "compare.jl", "run_cases.jl")
    cp(joinpath(@__DIR__, f), joinpath(dest, f); force=true)
end

julia = first(Base.julia_cmd())
cmd = `$julia --startup-file=no --project=$WT -t 1 $(joinpath(dest, "run_cases.jl")) $OUTDIR $(ARGS[3:end])`
@printf("Running: %s\n", string(cmd))
run(cmd)
@printf("Baselines written to %s\n", OUTDIR)
