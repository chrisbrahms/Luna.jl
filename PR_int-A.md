# Group A integration: the regression gate, `zmax` out of the grids, `Grid.RadialGrid`

Branch `gpu/int-A`. Base `gpu/01-zmax` (`80fe61ca`), merged with `gpu/02-radialgrid`
(`9dd6afaa`) and `gpu/00-harness` (`c0c7b76a`). All three branch from `evanescent`
(`fdf8dbe3`).

This branch is the Group A integration described in GPU_PLAN.md §5/§6. Apart from the
compatibility shim in the regression cases it contains no work of its own beyond resolving
the overlap between the three branches; `PR_01-zmax.md`, `PR_02-radialgrid.md` and
`PR_00-harness.md` are on the branch and describe what each one does.

## Merged branches

| branch | HEAD | what it does |
|---|---|---|
| `gpu/01-zmax` | `80fe61ca` | §4.13: `Grid.RealGrid`/`Grid.EnvGrid` lose the `zmax` field; `Luna.run` takes `zmax` as a required keyword and writes it to the output |
| `gpu/02-radialgrid` | `9dd6afaa` | §4.14: `Grid.RadialGrid` replaces `Hankel.QDHT` in every Luna dispatch signature; the transverse grid is written to the output |
| `gpu/00-harness` | `c0c7b76a` | the regression gate (`test/regression/`, `test/test_regression.jl`), `Luna.set_fftw_wisdom`, the `Adapt`/`GPUArraysCore` dependencies, `JLArrays` in the test target, and `benchmark/` |

Commits on this branch:

- `1021a5bf` Merge gpu/02-radialgrid into gpu/int-A (the textual merge and its conflict
  resolutions)
- `e94c2b14` Semantic sweep after the gpu/01-zmax + gpu/02-radialgrid merge
- `cb9415a8` Merge gpu/00-harness into gpu/int-A (no conflicts)
- `ff716000` Regression cases: a compatibility shim for the Group A signature changes

The merge and the sweep are separate commits so that a reviewer can see what the merge did
and what it did not do.

`Manifest.toml` is gitignored, and this worktree's copy was replaced by `gpu/00-harness`'s
so that `Adapt` and `GPUArraysCore` resolve and so that the package versions match the ones
the existing `fdf8dbe3` baseline was generated with.

## Conflicts and how they were resolved

Nine files conflicted, in thirteen hunks. All of them are radial set-up code that both
branches rewrote, and all of them have the same shape: `gpu/01-zmax` rewrote the
`Grid.RealGrid`/`Grid.EnvGrid` line to drop the leading length argument, `gpu/02-radialgrid`
rewrote the adjacent `Hankel.QDHT` line into a `Grid.RadialGrid`. Every one is resolved by
keeping both intents — the zmax-free grid constructor from 01 and the `RadialGrid` from 02.

| file | hunks | resolution |
|---|---|---|
| `examples/low_level_interface/freespace/radial.jl` | 1 | `Grid.RealGrid(800e-9, …)` + `q = Grid.RadialGrid(R, N)` |
| `…/radial_env.jl` | 1 | `Grid.EnvGrid(800e-9, …)` + `Grid.RadialGrid(R, N)` |
| `…/radial_hcf.jl` | 1 | `Grid.RealGrid(λ0, …)` + `Grid.RadialGrid(R, N)` |
| `…/radial_xypol_bbo.jl` | 1 | `Grid.RealGrid(λ0, …)` + `Grid.RadialGrid(R, N)` |
| `…/radial_xypol_env.jl` | 1 | `Grid.EnvGrid(λ0, …)` + `Grid.RadialGrid(R, N)` |
| `test/test_freespace.jl` | 3 | see below |
| `test/test_linops.jl` | 2 | 01's two grid lines + `Grid.RadialGrid(R, Nr)` / `Grid.RadialGrid(Rs, 64)` |
| `test/test_modes.jl` | 2 | 01's grid line + `q = Grid.RadialGrid(2a, 512)` |
| `test/test_boundaries.jl` | 1 | see below |

Three hunks carried content beyond the two-line pattern:

- `test/test_freespace.jl`, "evanescent channels": 01 replaced `gride.zmax/Boundaries.DEFAULT_N`
  by the literal `2e-3/Boundaries.DEFAULT_N` because the grid no longer holds the length.
  Kept, together with 02's `qe = Grid.RadialGrid(Re, 128)`.
- `test/test_freespace.jl`, the boundary-mode loop: 01's `Output.MemoryOutput(0, Ld, 3)` and
  `Luna.run(…; init_dz=0.1, boundary, zmax=Ld)` kept together with 02's
  `Grid.to_rspace(q, …)` in place of `q \ …`.
- `test/test_boundaries.jl`, "free space": 01 introduced a separate `zmax = 1e-3` local,
  which the testset's `ℓ = zmax/20` uses. Kept, together with 02's
  `q = Grid.RadialGrid(Rs, 32)`.

`src/Grid.jl`, `src/Luna.jl`, `src/Boundaries.jl`, `src/Fields.jl`, `src/LinearOps.jl`,
`src/NonlinearRHS.jl`, `src/Interface.jl` and `test/runtests.jl` merged without conflict.

**`gpu/00-harness` merged with no conflicts at all.** The only file it shares with the other
two is `src/Luna.jl`, where it adds `"fftw_wisdom" => true` to the `settings` `Dict` and a
`set_fftw_wisdom` function at the top of the module, well away from `gpu/01-zmax`'s changes
to `run`. `Project.toml`, `src/Utils.jl`, `docs/src/modules/Luna.md`, `test/test_utils.jl`,
`benchmark/` and `test/regression/` came over untouched.

Three places where both branches edited the same function and the textual merge happened to
be right were read in full:

- `Boundaries.log_setup` (01 changed its signature and its use of the length, 02 changed the
  `rdesc` line in the body). The merged function takes `zmax` as its second argument and
  dispatches `rdesc` on `sg isa Grid.RadialGrid`. It compiles and prints both the reference
  length and the transverse-collar description; `test_boundaries.jl` and `test_freespace.jl`
  exercise it.
- `Luna.run`. The order is: required-keyword check → `zmax = float(zmax)` → `max_dz` default
  → the `GridCondition` consistency check → `check_cache` → `Boundaries.setup` (which needs
  the `z0` `check_cache` may have moved) → band-limiting → `output(Grid.to_dict(grid),
  group="grid")` → `Output.hasdata(output, "zmax") || output("zmax", zmax)` → the
  `simulation_type` block, inside which 02's `output(Grid.to_dict(sg), group="spacegrid")`
  sits under the `!isnothing(sg)` guard → `RK45.solve_precon`. Both branches' writes happen
  after `check_cache`, so a resumed `HDF5Output` run finds `zmax` already present and does
  not warn, and the spacegrid group is written next to the other free-space entries.
- `Grid.to_dict`/`from_dict`: 01's generic `from_dict` (which now tolerates a legacy `zmax`
  key) and 02's `to_dict(::RadialGrid)`/`validate(::RadialGrid)` are independent.

## Semantic fixes (commit `e94c2b14`)

The textual merge leaves `gpu/02-radialgrid`'s code using the pre-01 signatures wherever 01
did not touch the same lines. The whole tree (`src/`, `test/`, `examples/`, `docs/`) was
grepped for `.zmax`, for `RealGrid(`/`EnvGrid(` calls with a leading length argument, for
`Output.MemoryOutput(0, grid.zmax, …)` and — with a paren-balancing scan, so multi-line calls
are covered — for `Luna.run(` calls with no `zmax` keyword. Three places needed changing, all
of them in code `gpu/02-radialgrid` added:

1. `test/test_freespace.jl`, the "transverse grid in the output" testset 02 added: two
   `Output.MemoryOutput(0, rgrid.zmax, 2)` and two `Luna.run(…)` calls without `zmax`. They
   now use the file's existing `L` (which is the value `rgrid.zmax` held on 02) and pass
   `zmax=L`.
2. `test/test_radialgrid.jl`, the "deprecated Hankel.QDHT entry points" testset: it built
   `Grid.RealGrid(1e-3, λ0, …)` and `Grid.EnvGrid(1e-3, λ0, …)`, which after 01 are the
   deprecated constructors and emit a warning. The leading length is dropped. The testset
   does not propagate, so nothing else in it needed the value.
3. `src/Luna.jl`: the `rcollar` entry of `run`'s docstring still called the radial transverse
   grid a `QDHT`. It now points at `Grid.RadialGrid`.

Nothing else remains. As the branch now stands, `.zmax` appears in seven places and none of
them is a grid field: a comment in `src/Luna.jl:565` explaining why `run` does
`float(zmax)`; `c.zmax` three times in `test/regression/cases.jl` and `case.zmax` twice in
`benchmark/run.jl`, which are the `zmax` field of the regression `Case` struct, not of a
grid; and one line of `test/regression/README.md`. The only calls to the deprecated grid
constructors are the ones in `test/test_grid.jl` that test the deprecation and the
docstrings in `src/Grid.jl` that document it. `test/runtests.jl` carries both branches' new
entries — `test_grid.jl` (line 23) and `test_radialgrid.jl` (line 122).

Nothing in the ~35 `examples/` lines the `gpu/01-zmax` reviewer counted on `gpu/02-radialgrid`
needed fixing after the merge: all the example files 02 touched are ones 01 touched too, so
they either conflicted (the five above) or took 01's side cleanly.

## The compatibility shim in the regression cases (commit `ff716000`)

`test/regression/generate.jl` copies `cases.jl` into a worktree of the baseline commit and
runs it there, so the baseline is the *old* Luna running the *new* case definitions. One
`cases.jl` therefore has to build the same propagations on `evanescent` and on this branch,
and Group A changed two APIs it uses:

- `gpu/01-zmax`: the grid constructors lose their leading `zmax` argument, `Grid.RealGrid`
  loses the field, `Luna.run` gains a required `zmax` keyword, and `Boundaries.setup` gains
  a `zmax` positional argument after `z0`.
- `gpu/02-radialgrid`: `Grid.RadialGrid` replaces `Hankel.QDHT`.

`gpu/00-harness` anticipated this and said `cases.jl` would have to be updated and would
then stop running against `evanescent`. It does not have to: a marked block at the top of
`cases.jl` detects which API is loaded and the cases go through it. Detection is from the
API itself, so nothing has to be edited when a further branch is merged:

```julia
const GRID_HAS_ZMAX  = hasfield(Grid.RealGrid, :zmax)   # true before gpu/01-zmax
const HAS_RADIALGRID = isdefined(Grid, :RadialGrid)     # true from gpu/02-radialgrid on
```

and four helpers built on them: `makegrid(GT, zmax, referenceλ, λ_lims, trange; kwargs...)`,
`runkw(zmax)` (which returns `(; zmax)` or an empty `NamedTuple` for splatting into
`Luna.run`), `radialgrid(R, N)` and `absorber_setup(...)` (used by `benchmark/run.jl`, which
builds an absorber outside `Luna.run`). Nine call sites in `cases.jl` (six `makegrid`, two
`radialgrid`, one `runkw`) and one in `benchmark/run.jl` use them; `compare.jl`, `run_cases.jl`, `generate.jl`, `sensitivity.jl`,
`tolerances.jl` and `test/test_regression.jl` needed no change, because they go through
`runcase`.

**The shim is numerically inert on `evanescent`, checked rather than assumed.** Regenerating
the `fdf8dbe3` baseline with the shimmed `cases.jl`, into a separate directory, reproduces
the stored baseline bit for bit: 21 files, 502 datasets, `isequal` throughout. So the Group A
gate below runs against the existing baseline and nothing was regenerated over it.

The block is marked in the file and is to be deleted, with its call sites, once no baseline
in use predates `gpu/int-A`. `test/regression/README.md` says so, and gains a section on
which base a branch should compare against.

## Regression: Group A against `evanescent`

`LUNA_REGRESSION_BASE=fdf8dbe3`, the baseline the `gpu/00-harness` implementer generated
with 21 cases, not regenerated. **460 pass, 0 fail, 88 s** (wall clock, with other agents'
jobs on the same 10-core machine, so an upper bound). Δ is the observed
`maximum(abs, Δ)/maximum(abs, baseline)` (per component and per save for `Eω`); tol is from
`tolerances.jl`. The four non-zero entries are in bold.

| case | `:fixed` Eω Δ / tol | `:fixed` stats Δ / tol | `:adaptive` Eω Δ / tol | `:adaptive` stats Δ / tol |
|---|---|---|---|---|
| `modeavg_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 5.9e-07 |
| `modeavg_field_nothg` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.3e-12 | 0 / 2.1e-04 |
| `modeavg_field_plasma` | 0 / 1.0e-12 | 0 / 2.5e-12 | 0 / 1.0e-12 | 0 / 7.3e-06 |
| `modeavg_field_raman` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.4e-09 |
| `modeavg_field_mixture` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.5e-10 |
| `modeavg_field_adk` | 0 / 1.0e-12 | 0 / 2.3e-11 | 0 / 1.0e-12 | 0 / 3.8e-05 |
| `modeavg_field_vector` | 0 / 1.0e-12 | 0 / 1.8e-09 | 0 / 1.0e-12 | 0 / 9.8e-05 |
| `modeavg_env_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 5.6e-12 | 0 / 5.3e-04 |
| `modeavg_env_raman` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.9e-12 | 0 / 4.2e-04 |
| `modeavg_env_thg` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.3e-07 |
| `gnlse_sech` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 3.6e-11 | 0 / 1.6e-09 |
| `gnlse_raman_shock` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.9e-10 | 0 / 1.2e-08 |
| `multimode_field_plasma` | 0 / 1.2e-11 | 0 / 2.9e-07 | 0 / 1.8e-07 | 0 / 7.7e-06 |
| `radial_field_kerr` | **1.901e-15** / 1.0e-12 | 0 / 1.0e-12 | **6.713e-11** / 1.2e-08 | 0 / 1.0e-12 |
| `radial_env_kerr` | **2.866e-15** / 1.0e-12 | 0 / 1.0e-12 | **4.017e-11** / 8.2e-09 | 0 / 1.0e-12 |
| `free3d_env_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 |
| `free2d_field_chi2` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 7.4e-12 | 0 / 1.0e-12 |
| `free2d_env_chi2` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 5.0e-12 | 0 / 1.0e-12 |
| `gradient_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.6e-07 | 0 / 1.6e-04 |
| `taper_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 8.7e-07 | 0 / 3.1e-05 |
| `modeavg_field_legacy` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.3e-11 | 0 / 7.8e-08 |

**Largest difference over all cases and modes: 6.713e-11**, in `radial_field_kerr`
`:adaptive`.

Every case except the two radial ones is **exactly zero** in both modes and both classes,
including the step count. That is the statement `gpu/01-zmax` and `gpu/00-harness` make about
themselves — neither changes any arithmetic — and it holds through the merge.

The two radial cases move, which is what `gpu/02-radialgrid` predicted and measured: it
replaces Hankel's `permutedims` + left GEMM by a right GEMM on a reshape in
`Boundaries.RadialCollar`, in `Fields.transform` (the input field) and in the `TransRadial`
noise setup, and folds the scalar `scaleRK` into the transform matrix. Both change the order
of the floating-point operations.

- `:fixed` (20 imposed steps, so the difference is attributable to the arithmetic alone):
  **1.901e-15** and **2.866e-15**. That is the same order as the 4.82e-15 worst case
  `PR_02-radialgrid.md` reports over its own 12 `test_freespace.jl` cases, and about one to
  three ulps. Both are three orders of magnitude inside the 1e-12 tolerance.
- `:adaptive`: **6.713e-11** and **4.017e-11**, against tolerances of 1.2e-08 and 8.2e-09.
  These are larger because the step-size controller amplifies a rounding change: the one-ulp
  sensitivity of those two cases in that mode is 1.24e-10 and 8.18e-11 (`sensitivity.jl`, in
  `PR_00-harness.md`), so the observed differences are *about half of what one ulp on the
  input field does to the same case*. The step counts are unchanged, which the gate checks
  separately and which is a hard failure if it moves.
- The `stats` class is exactly zero for both radial cases in both modes. The free-space cases
  record only `z` and `dz`; in `:fixed` those are compared and match exactly, which confirms
  the step sequence really was imposed.

This table is the Group A regression record against `evanescent`.

## Regression: `gpu/int-A` against itself

Baseline generated from this branch's `ff716000` (`test/regression/generate.jl ff716000`),
gate run with `LUNA_REGRESSION_BASE=ff716000`: **460 pass, 0 fail, every case, both modes,
both classes exactly `0.000e+00`.** Largest difference over all cases and modes:
`0.000e+00`. 88 s, again under contention.

That is the baseline Group B onwards compares against
(`LUNA_REGRESSION_BRANCH=gpu/int-A`), and it also confirms that the generator reproduces the
run environment exactly from a separate worktree and a separate process on the new
signatures.

## Tests

`julia --project=/Users/mb140/.julia/dev/Luna-gpu/gpu-int-A -t 1`, with
`Luna.set_fftw_mode(:estimate)`, `Luna.set_fftw_threads(1)` and
`LinearAlgebra.BLAS.set_num_threads(1)`. `test_processing.jl` was run with an empty `ARGS`
(the `Scans` command-line parsing reads `ARGS`). Machine: Apple M1 Pro, Julia 1.13, with
other agents' jobs running at the same time, so the times are upper bounds.

Stage 1, after the `gpu/01-zmax` + `gpu/02-radialgrid` merge:

| file | result | time |
|---|---|---|
| `test/test_grid.jl` | 88 pass | 7.3 s |
| `test/test_radialgrid.jl` | 159 pass | 9.5 s |
| `test/test_linops.jl` | 196 pass | 10.7 s |
| `test/test_modes.jl` | 587 pass | 27.2 s |
| `test/test_boundaries.jl` | 182 pass | 18.9 s |
| `test/test_fields.jl` | 179 pass | 52.7 s |
| `test/test_output.jl` | 71 pass | 11.6 s |
| `test/test_freespace.jl` | 77 pass | 308.9 s |
| `test/test_processing.jl` | 64 pass | 30.1 s |
| `test/test_interface.jl` | 301 pass | 257.3 s |

1904 assertions, no failures, errors or broken tests. Every count matches the count the
originating branch reported for that file, so the merge did not silently drop or weaken a
testset.

Rerun after the `gpu/00-harness` merge and the shim, on the `gpu/00-harness` `Manifest.toml`
(so `Adapt` and `GPUArraysCore` resolve):

| file | result | time |
|---|---|---|
| `test/test_utils.jl` | 33 pass | 3.8 s |
| `test/test_grid.jl` | 88 pass | 6.1 s |
| `test/test_radialgrid.jl` | 159 pass | 9.7 s |
| `test/test_output.jl` | 71 pass | 12.7 s |
| `test/test_interface.jl` | 301 pass | 254.3 s |
| `test/test_regression.jl` vs `fdf8dbe3` | 460 pass | 88.4 s |
| `test/test_regression.jl` vs `ff716000` (self) | 460 pass | 88.0 s |

Every time in both tables was measured with other agents' propagations running on the same
10-core machine, so all of them are upper bounds; only the pass counts and the differences
are exact.

`test_utils.jl` covers `gpu/00-harness`'s new `set_fftw_wisdom` testset and passes
unchanged at its 33 assertions.

`benchmark/run.jl` was run for the two radial cases, through the shimmed `absorber_setup`,
as a check that the benchmark still builds an absorber outside `Luna.run`:

```
case                          state          rhs         step         prop
radial_field_kerr              4128   155.541 µs     1.675 ms    64.991 ms
radial_env_kerr                4096   125.708 µs     1.535 ms    54.146 ms
```

which is within run-to-run noise of `PR_00-harness.md`'s numbers for the same two cases on
`evanescent` (156 µs / 1.68 ms / 63.3 ms and 126 µs / 1.54 ms / 56.2 ms): `Grid.RadialGrid`
did not change the per-step cost of a radial run.

All 139 tracked `.jl` files parse — `src/` (31), `test/` (43), `examples/` (62),
`benchmark/` (1), `docs/` (1), `deps/` (1) — checked with `Meta.parseall` on the file
contents, walking the result for `:error`/`:incomplete` nodes. (The stage-1 figure of 129
was `src/`, `test/` and `examples/` only, before `gpu/00-harness` added the regression and
benchmark scripts.)

Examples, run with the plotting sections cut off and (for the free-space ones) the grid sizes
reduced so each run takes seconds:

```
freespace/radial.jl                     OK (8.2 s)
freespace/radial_env.jl                 OK (3.1 s)
freespace/radial_hcf.jl                 OK (18.7 s)
freespace/radial_xypol_bbo.jl           OK (4.2 s)
freespace/radial_xypol_env.jl           OK (2.0 s)
freespace/full3D.jl                     OK (7.9 s)
basic_modeAvg.jl                        OK (5.9 s)
gradients/gradient_modeAvg.jl           OK (9.6 s)
```

## Known gaps and carried-over notes

- The regression scripts in the project scratchpad (`radial_cases.jl`,
  `reviews/scratch-02/rerun_cases.jl`) still use `Grid.RealGrid(L, …)` and `Luna.run` without
  `zmax`, so they do not run against this branch. They are not in the repository and are
  superseded by the gate; `test/regression/cases.jl` is the maintained version of what they
  did.
- `Manifest.toml` in this worktree is `gpu/00-harness`'s. It is gitignored, but it matters:
  it is the one the existing `fdf8dbe3` baseline was generated with, and the one that has
  `Adapt` and `GPUArraysCore`. A fresh worktree of this branch needs `Pkg.instantiate()` (or
  that Manifest) before `using Luna` works.
- Carried over from `gpu/00-harness`: the gate is not wired into CI; the free-space cases
  record only `z` and `dz`; the documentation build fails on six pre-existing unresolved
  `@ref`s, identically with and without this branch.
- Carried over from `gpu/02-radialgrid`: `Test.detect_ambiguities(Luna; recursive=true)`
  reports 13 new latent pairs in `Luna.setup` and `NonlinearRHS.TransRadial`, all judged
  unreachable by the reviewer.
- Carried over from `gpu/01-zmax`: the pre-existing `Scans.QueueExec` queue-file bug, and the
  two `examples/low_level_interface/` scripts (`full_modal/basic_modal_full_bothpolarisations.jl`,
  `polarisation/elliptical_env.jl`) that are broken on `evanescent` as well.
