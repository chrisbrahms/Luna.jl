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
| `cases.jl` | The case matrix: `RegressionCases.CASES`, a `Vector{Case}`, and `runcase(case, mode)`. Self-contained, so it can be `include`d on an older commit; a compatibility shim at the top of it keeps it running on both sides of the Group A signature changes (see "The compatibility shim"). |
| `compare.jl` | `rundict`, `saverun`, `loadcase`, `compare`: the storage format, the difference metric and the tolerance classes. Also copied to the older commit. |
| `run_cases.jl` | Runs every case in every mode and writes the HDF5 files. Copes with a Luna that has no `set_fftw_wisdom`. |
| `generate.jl` | Makes a git worktree of a given commit, copies the Manifest and the three files above into it, and runs `run_cases.jl` there. |
| `sensitivity.jl` | The one-ulp sensitivity study the tolerances are derived from. |
| `tolerances.jl` | The per-case, per-mode, per-class tolerances. |
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

#### Which base a GPU-project branch should use

From Group B on, a branch's own base is `gpu/int-A` or a descendant of it, and the question
it has to answer is "did *this branch* change the output". That is

```
LUNA_REGRESSION_BRANCH=gpu/int-A julia --project=$PWD -t 1 test/test_regression.jl
```

which takes the merge-base of `HEAD` with `gpu/int-A`, so a branch cut from `gpu/int-A` and
a branch cut from a later integration branch each compare against their own starting point
and every case must come back exactly `0.000e+00` unless the branch states otherwise.
Generate that baseline once with `test/regression/generate.jl gpu/int-A`.

The other question — "how far has the whole project moved from where it started" — is

```
LUNA_REGRESSION_BASE=fdf8dbe3 julia --project=$PWD -t 1 test/test_regression.jl
```

which compares against `evanescent` and accumulates every branch's deltas. That comparison
is what the per-group summary records. It works because `cases.jl` still runs on
`evanescent` through the shim below; when that stops being true, this form goes away and
the per-group record becomes the chain of merge-base comparisons.

Baselines from different commits live in different directories, so both forms can be run
without regenerating anything.

### Running part of the matrix

`LUNA_REGRESSION_ONLY` and `LUNA_REGRESSION_SKIP` take comma-separated case names and
restrict the matrix (`ONLY` first, then `SKIP`). They exist for one situation: a branch
which *adds* a case has no baseline for it in any older baseline directory, and without a
filter the gate reports that as a failure to load rather than as the missing baseline it
is. Such a branch runs

```
LUNA_REGRESSION_SKIP=<new case> LUNA_REGRESSION_BASE=<old commit> julia ... test/test_regression.jl
```

for the older baselines, and generates a baseline of its own for the new case:

```
julia --project=$PWD -t 1 test/regression/generate.jl <base commit> <dir> <new case>
LUNA_REGRESSION_DIR=<dir> LUNA_REGRESSION_ONLY=<new case> LUNA_REGRESSION_BASE=<base commit> \
    julia --project=$PWD -t 1 test/test_regression.jl
```

writing it to a directory of its own so that the shared baseline directories, which the
integration branches regenerate wholesale, are left alone. The header prints the lists as
they were given (`cases: 21 of 22  (skipped: radial_field_raman)`) whenever the selection
is not the whole matrix; an unknown case name is an error, and neither variable can make a
failing case pass.

### The compatibility shim in `cases.jl`

`generate.jl` copies `cases.jl` into a worktree of the baseline commit, so the same file has
to build the same propagations on `evanescent` and on the branch under test. Group A changed
two APIs it uses: `gpu/01-zmax` took `zmax` out of the grid constructors and the grid struct
and made it a required keyword of `Luna.run` (and a positional argument of
`Boundaries.setup`), and `gpu/02-radialgrid` replaced `Hankel.QDHT` with `Grid.RadialGrid`.

The block at the top of `cases.jl` detects which of the two APIs is loaded — from the API
itself (`hasfield(Grid.RealGrid, :zmax)`, `isdefined(Grid, :RadialGrid)`), not from a commit
or a version — and the cases call `makegrid`, `runkw`, `radialgrid` and `absorber_setup`
instead of the API directly. On `evanescent` each of those reduces to exactly the call the
pre-Group-A `cases.jl` made: regenerating the `fdf8dbe3` baseline with the shimmed file
reproduces the original one bit for bit, over the 21 cases which existed then, both modes
and 502 datasets. `radial_field_raman`, added in `gpu/20-radial-device`, is the 22nd and
`rect_modal_field`, added in `gpu/26-rectmode-fix`, the 23rd; both run on `evanescent`
through the same shim, but no baseline generated before their branches contains them.

The block is marked in the file and should be deleted, with its call sites, once no baseline
in use predates `gpu/int-A`.

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

1. Make sure `cases.jl` still runs against the new base commit. Through Group A it does,
   because of the shim described above; a later branch that changes an API `cases.jl` uses
   has to extend the shim the same way, or drop the older base.
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


## The two run modes

Every case runs twice.

- `:fixed` — `min_dz = max_dz = init_dz = zmax/20`. Equal limits make `RK45.steplims!`
  return the same step every time, so the step-size controller is bypassed and every
  difference is attributable to the arithmetic rather than to a different sequence of
  steps. The 20 steps are chosen so that the step equals the `:rate` absorber's reference
  length `zmax/boundary_N`, which is what `Luna.run` would otherwise cap `max_dz` at.
- `:adaptive` — the default controller, started from `init_dz = zmax/1000`. This is the
  sanity check that the controller still takes the same decisions.

  Be aware of how little of the controller this actually exercises. The adaptive runs take
  8–58 accepted steps against the fixed mode's 20, and in all but `modeavg_field_legacy` and
  `modeavg_field_mixture` the step size reaches `max_dz` within the first few steps and stays
  there. So the controller is genuinely adapting for roughly three steps per case; after that
  the adaptive run is a fixed-step run at a different step size. It is not a substitute for
  the `:fixed` mode, and a branch that changes the error estimate or the norm should be
  checked another way as well.

## The metric

Per quantity: `maximum(abs, new - baseline) / maximum(abs, baseline)`, over `Eω`, the save
positions `z`, and every compared statistic. Elementwise relative differences are
meaningless in the window tapers and outside `grid.sidx`, where the field is orders of
magnitude below its peak.

### `Eω` is normalised per component and per save

`RegressionCompare.fielddiff` reduces the frequency axis and any transverse axes away and
forms the ratio separately for every (mode or polarisation, save) pair; the reported number
is the largest of those, and the gate prints which pair it came from.

`Eω` is `(Nω, Nz)` mode-averaged, `(Nω, Nm, Nz)` multimode or vector, `(Nω, Npol, Nk, Nz)`
in 2-D free space and `(Nω, Npol, Nkx, Nky, Nz)` in 3-D, so axis 1 is always frequency, the
last axis is always the save, axis 2 is the mode or polarisation whenever there are more
than two axes, and the rest are transverse.

One global normalisation would hide weak components. In `multimode_field_plasma` the
per-mode maxima of `|Eω|` at the last save are `[2.9e5, 1.8e1, 3.2e0, 7.2e-2]`, a spread of
4e6: a change that rewrote the weakest mode entirely would sit four million times below the
tolerance the strongest mode makes sensible. The per-component metric measures it directly —
that case's `:fixed` `Eω` sensitivity is 1.2e-13 and it is reached in mode 4.

A component/save slice whose baseline maximum is zero is skipped when the difference is zero
there too, and reported as `Inf` otherwise.

### Two tolerance classes

`RegressionCompare.classof` puts every quantity into `:Eω` (the field) or `:stats` (the save
grid `z`, the statistics and the step-count check), and `tolerances.jl` carries one tolerance
per (case, mode, class). The statistics are recorded once per accepted step, so in the
adaptive mode they inherit the step-size controller's sensitivity; sharing one tolerance with
them gave `Eω` several orders of magnitude of slack. `Eω` is the quantity every later branch
must not move, so it gets its own number.

### The number of accepted steps

Every statistic is recorded once per accepted step, so a change in the number of steps makes
every statistic fail its size check with an uninformative `Inf`. `compare` therefore checks
`length(stats["z"])` first and, if it differs, reports a single `"step count"` entry
(`changed: N steps, baseline M`) and compares no statistics. That is a hard failure in both
modes, including `:adaptive`: a different number of steps is a different propagation.

### What is compared in which mode

`RegressionCompare.skipstats(mode)` decides. In the `:fixed` mode nothing is excluded: the
step sequence is imposed, so `stats/z` and `stats/dz` must match exactly, and checking them
also confirms that it really was imposed.

In the `:adaptive` mode `stats/z` and `stats/dz` (`RegressionCompare.STEP_STATS`) are
excluded. They record the step sequence, not the field, and the step-size controller's
accept/reject decision and its PI update respond to a one-ulp change far more strongly than
the field does — the response compounds over the steps it takes to ramp `init_dz` up to
`max_dz`. With them in the comparison the tolerance for `modeavg_env_kerr` came out at
0.73, which is no constraint on anything.

Excluding them removes the *graded* part of the step sequence's sensitivity, not all of it:
the step count is still checked, in both modes, as described above.

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

The stored baselines were generated under Julia 1.13.0 with the `Manifest.toml` of the
worktree that generated them, which `generate.jl` copies into the baseline worktree so that
both sides use the same package versions. That is what makes an *exact* zero reproducible:
a different Julia, a different stdlib or a different set of JLLs -- OpenBLAS and FFTW above
all -- changes the order of the floating-point operations, and the gate then reports
rounding-level differences everywhere rather than zeros. If that happens, regenerate the
baseline in the environment the gate is being run in before concluding anything about the
branch.

## Tolerances

`tolerances.jl` holds one tolerance per case, mode and class. Each is 100x the case's
one-ulp sensitivity for that class, with a floor of 1e-12. The sensitivity is measured by

```
julia --project=$PWD -t 1 test/regression/sensitivity.jl
```

which runs each case twice, once with the initial `Eω` multiplied by `1 + eps()`, and prints
the per-quantity breakdown and a `TOLERANCES` literal to paste into `tolerances.jl`. It uses
the gate's comparison, exclusions included, so its output pastes in unchanged.

In the `:fixed` mode `Eω` is at the 1e-12 floor for every case but `multimode_field_plasma`
(6.7e-12, in the weakest mode). In the `:adaptive` mode `Eω` runs from the floor to 7.3e-05,
and the `:stats` class is looser, up to 7.3e-02, because the statistics are recorded per
step. The reasons for every number above the floor are in the header comment of
`tolerances.jl`.

The four ionising cases are the loose ones in the `:adaptive` mode, because the electron
density feeds back into the step-size controller. That feedback also fixes the energy of
`modeavg_field_plasma`: at 300 µJ its adaptive step count is not reproducible under a
one-ulp perturbation (92 steps against 98), which no tolerance can express, so the case
runs at 175 µJ, the most strongly ionising energy whose step count is stable. A case whose
step count moves under one ulp does not belong in the matrix: the step-count check is a
hard failure by design, and the only way to keep such a case green is to disable the check
for it.

Read the table `test/test_regression.jl` prints, not only its pass/fail: a case that moves
from 0 to 1e-9 in the `:fixed` mode is a real change even if it is inside the tolerance.

## Adding a case

Append a `Case` to `CASES` in `cases.jl`, add an entry to `TOLERANCES` (a `Dict` per mode
of a `Dict` per class; a missing entry falls back to the 1e-12 floor, which will usually
fail, so run `sensitivity.jl` for the new case), and regenerate the baseline (`generate.jl` writes one HDF5 file per case, so only the new file is missing; the
gate will report it as missing until the baseline is regenerated). `cases.jl` must keep
working on the oldest commit it will be run against, so it must not use API added on the
branch under test.
