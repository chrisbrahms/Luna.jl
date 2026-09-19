# Regression gate

`test/test_regression.jl` checks that a change to Luna leaves the default CPU output
unchanged to rounding level. It compares a matrix of small propagations against a baseline
generated from an earlier commit in the same environment. The baselines are HDF5 files and
are never committed (`*.h5` is gitignored).

The gate is **not** part of `test/runtests.jl`, because it needs a baseline that does not
exist on a fresh checkout.

## Files

| File | What it is |
| --- | --- |
| `cases.jl` | The case matrix: `RegressionCases.CASES`, a `Vector{Case}`, and `runcase(case, mode)`. Self-contained, so it can be `include`d on an older commit. |
| `compare.jl` | `rundict`, `saverun`, `loadcase`, `compare`: the storage format and the difference metric. Also copied to the older commit. |
| `run_cases.jl` | Runs every case in every mode and writes the HDF5 files. Copes with a Luna that has no `set_fftw_wisdom`. |
| `generate.jl` | Makes a git worktree of a given commit, copies the Manifest and the three files above into it, and runs `run_cases.jl` there. |
| `sensitivity.jl` | The one-ulp sensitivity study the tolerances are derived from. |
| `tolerances.jl` | The per-case, per-mode tolerances. |
| `../test_regression.jl` | The gate itself. |

## Running it

Generate a baseline from the branch's base commit (once per commit; the temporary worktree
is kept so the next generation does not have to precompile Luna again):

```
julia --project=$PWD -t 1 test/regression/generate.jl evanescent
```

Then run the gate:

```
julia --project=$PWD -t 1 test/test_regression.jl
```

It prints a table of the observed maximum differences whether it passes or not.

### Where things go

- Baselines: `joinpath(Luna.Utils.cachedir(), "regression", <full sha>)`, overridden by the
  second argument of `generate.jl` and by `ENV["LUNA_REGRESSION_DIR"]` for the gate.
- Baseline worktrees: `../baselines/<sha[1:10]>` next to this worktree, overridden by
  `ENV["LUNA_REGRESSION_WORKTREES"]`. Remove them with `git worktree remove <path>`.

### Choosing the baseline commit

The gate uses, in order:

1. `ENV["LUNA_REGRESSION_BASE"]`, if set (anything `git rev-parse` accepts);
2. otherwise the merge-base of `HEAD` with `ENV["LUNA_REGRESSION_BRANCH"]` (default
   `evanescent`).

## Re-baselining onto a new base branch

`generate.jl` takes any revision `git rev-parse` accepts — a branch name, a tag, a SHA,
`HEAD`. It resolves it to a full SHA, makes (or reuses) a detached worktree of it under
`../baselines/<sha[1:10]>`, copies this worktree's `Manifest.toml` and *this branch's*
`cases.jl`, `compare.jl` and `run_cases.jl` into it, and runs `run_cases.jl` there in a
fresh `julia -t 1` process. So the baseline is always the *old* Luna running the *new* case
definitions. Verified for two different commits (the branch base and a later commit on the
branch).

When a branch's base moves — for example when `gpu/00-harness` is merged into `gpu/int-A`
and the branches after it use `gpu/int-A` as their base:

1. Make sure `cases.jl` still runs against the new base commit. It may not: `gpu/01-zmax`
   moves `zmax` out of the grids and `gpu/02-radialgrid` replaces `Hankel.QDHT` with
   `Grid.RadialGrid`, so the `Grid.RealGrid(...)`/`Grid.EnvGrid(...)` constructor calls,
   `grid.zmax`, the `Luna.run` signature and `Hankel.QDHT(R_FREE, 32, dim=3)` in `cases.jl`
   all have to be updated at that point. Once they are, `cases.jl` no longer runs against
   `evanescent`, which is expected and is why the baseline commit is an argument.
2. Generate the new baseline:

   ```
   julia --project=$PWD -t 1 test/regression/generate.jl gpu/int-A
   ```

3. Point the gate at it, either by setting `LUNA_REGRESSION_BRANCH=gpu/int-A` (the gate
   then takes the merge-base of `HEAD` with that branch) or by naming the commit outright:

   ```
   LUNA_REGRESSION_BASE=gpu/int-A julia --project=$PWD -t 1 test/test_regression.jl
   ```

4. Record the differences of the new base against the old one before throwing the old
   baseline away — a branch that legitimately changes a case (`gpu/02-radialgrid`'s radial
   deltas, for instance) has to state which case moved and by how much, and that number can
   only be measured across the re-baselining.

Baselines from different commits live in different directories, so old and new can coexist
and both gates can be run.

## The two run modes

Every case runs twice.

- `:fixed` — `min_dz = max_dz = init_dz = zmax/20`. Equal limits make `RK45.steplims!`
  return the same step every time, so the step-size controller is bypassed and every
  difference is attributable to the arithmetic rather than to a different sequence of
  steps. The 20 steps are chosen so that the step equals the `:rate` absorber's reference
  length `grid.zmax/boundary_N`, which is what `Luna.run` would otherwise cap `max_dz` at.
- `:adaptive` — the default controller, started from `init_dz = zmax/1000`. This is the
  sanity check that the controller still takes the same decisions.

## The metric

Per saved quantity: `maximum(abs, new - baseline) / maximum(abs, baseline)`, over `Eω` at
all saves, the save positions `z`, and every compared statistic. Elementwise relative
differences are meaningless in the window tapers and outside `grid.sidx`, where the field
is orders of magnitude below its peak. The reported number per case and mode is the largest
over all compared quantities.

### What is compared in which mode

`RegressionCompare.skipstats(mode)` decides. In the `:fixed` mode nothing is excluded: the
step sequence is imposed, so `stats/z` and `stats/dz` must match exactly, and checking them
also confirms that it really was imposed.

In the `:adaptive` mode `stats/z` and `stats/dz` (`RegressionCompare.STEP_STATS`) are
excluded. They record the step sequence, not the field, and the step-size controller's
accept/reject decision and its PI update respond to a one-ulp change far more strongly than
the field does — the response compounds over the steps it takes to ramp `init_dz` up to
`max_dz`. With them in the comparison the tolerance for `modeavg_env_kerr` came out at
0.73, which is no constraint on anything; with them out the adaptive tolerances are set by
`Eω` and the physical statistics, and the loosest is 5.3e-04.

`rundict`'s top-level `"z"` is a different thing — the save grid, fixed by
`Output.GridCondition` — and is compared in both modes.

Nothing is excluded from what the baseline *stores*: the HDF5 files always contain `z` and
`dz`, so this choice can be revisited without regenerating anything.

## Determinism requirements

Both the baseline and the gate run with

- `Luna.set_fftw_mode(:estimate)`,
- `Luna.set_fftw_threads(1)`,
- `LinearAlgebra.BLAS.set_num_threads(1)`,
- FFTW wisdom disabled (`Luna.set_fftw_wisdom(false)`),
- `shotnoise=false` in every case,
- a single Julia thread (`-t 1`).

Wisdom has to be off because the wisdom file is shared by every process using the same
Julia depot: `:patient` wisdom written by an unrelated run (including Luna's own
precompilation) silently changes the plan an `:estimate` run gets, and with it the order of
the floating-point operations. On an older commit, which has no `set_fftw_wisdom`,
`run_cases.jl` redefines `Luna.Utils.loadFFTwisdom`/`saveFFTwisdom` to do nothing and calls
`FFTW.forget_wisdom()`.

## Tolerances

`tolerances.jl` holds a per-case, per-mode tolerance. Each is 100x the case's one-ulp
sensitivity, with a floor of 1e-12. The sensitivity is measured by

```
julia --project=$PWD -t 1 test/regression/sensitivity.jl
```

which runs each case twice, once with the initial `Eω` multiplied by `1 + eps()`, and
prints the per-quantity breakdown and a `TOLERANCES` literal to paste into
`tolerances.jl`.

`sensitivity.jl` uses the same comparison the gate does, exclusions included, so its output
can be pasted into `tolerances.jl` unchanged.

In the `:fixed` mode the sensitivity is at rounding level for all but two cases. In the
`:adaptive` mode it is set by `Eω` and the physical statistics, mostly `peakpower`,
`peakintensity` and `energy`, and ranges from 1e-14 to 5.3e-06. The reasons for the two
`:fixed` cases above 1e-12, and for the two loosest adaptive cases, are in the header
comment of `tolerances.jl`.

Read the table `test/test_regression.jl` prints, not only its pass/fail: a case that moves
from 0 to 1e-9 in the `:fixed` mode is a real change even if it is inside the tolerance.

## Adding a case

Append a `Case` to `CASES` in `cases.jl`, add an entry to `TOLERANCES`, and regenerate the
baseline (`generate.jl` writes one HDF5 file per case, so only the new file is missing; the
gate will report it as missing until the baseline is regenerated). `cases.jl` must keep
working on the oldest commit it will be run against, so it must not use API added on the
branch under test.
