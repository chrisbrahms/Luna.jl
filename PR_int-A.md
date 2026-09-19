# Group A integration: `zmax` out of the grids and `Grid.RadialGrid`

Branch `gpu/int-A`. Base `gpu/01-zmax` (`80fe61ca`), merged with `gpu/02-radialgrid`
(`9dd6afaa`). Both branch from `evanescent` (`fdf8dbe3`).

This branch is the Group A integration described in GPU_PLAN.md §5/§6. It contains no new
work of its own beyond resolving the overlap between the two branches; the two PR
descriptions `PR_01-zmax.md` and `PR_02-radialgrid.md` are on the branch and describe what
each one does.

## Merged branches

| branch | HEAD | what it does |
|---|---|---|
| `gpu/01-zmax` | `80fe61ca` | §4.13: `Grid.RealGrid`/`Grid.EnvGrid` lose the `zmax` field; `Luna.run` takes `zmax` as a required keyword and writes it to the output |
| `gpu/02-radialgrid` | `9dd6afaa` | §4.14: `Grid.RadialGrid` replaces `Hankel.QDHT` in every Luna dispatch signature; the transverse grid is written to the output |

Commits on this branch:

- `1021a5bf` Merge gpu/02-radialgrid into gpu/int-A (the textual merge and its conflict
  resolutions)
- `e94c2b14` Semantic sweep after the gpu/01-zmax + gpu/02-radialgrid merge

They are separate so that a reviewer can see what the merge did and what it did not do.

`gpu/00-harness` is **not** merged yet; see "Pending" below.

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

Nothing else remains: the only `.zmax` in the tree is a comment in `src/Luna.jl` explaining
why `run` does `float(zmax)`, and the only calls to the deprecated grid constructors are the
ones in `test/test_grid.jl` that test the deprecation and the docstrings in `src/Grid.jl`
that document it. `test/runtests.jl` carries both branches' new entries — `test_grid.jl`
(line 23) and `test_radialgrid.jl` (line 122).

Nothing in the ~35 `examples/` lines the `gpu/01-zmax` reviewer counted on `gpu/02-radialgrid`
needed fixing after the merge: all the example files 02 touched are ones 01 touched too, so
they either conflicted (the five above) or took 01's side cleanly.

## Tests

`julia --project=/Users/mb140/.julia/dev/Luna-gpu/gpu-int-A -t 1`, with
`Luna.set_fftw_mode(:estimate)`, `Luna.set_fftw_threads(1)` and
`LinearAlgebra.BLAS.set_num_threads(1)`. `test_processing.jl` was run with an empty `ARGS`
(the `Scans` command-line parsing reads `ARGS`). Machine: Apple M1 Pro, Julia 1.13, with
other agents' jobs running at the same time, so the times are upper bounds.

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

All 129 `.jl` files under `src/`, `test/` and `examples/` parse (`Meta.parseall` on the file
contents, walking the result for `:error`/`:incomplete` nodes).

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

## Pending

- **`gpu/00-harness` is not merged.** It carries the regression generator and
  `test/test_regression.jl`, the `settings["fftw_wisdom"]` switch, the `Adapt`/`GPUArraysCore`
  dependencies and the `benchmark/` environment. It will be merged into this branch as
  stage 2, and the regression gate re-baselined here for the changed signatures: the
  baseline generator takes the base commit as an argument, and `gpu/02-radialgrid`'s radial
  deltas against `evanescent` (worst case 4.82e-15, recorded in `PR_02-radialgrid.md`) have
  to be recorded before the re-baseline, as GPU_PLAN.md §6 requires.
- The regression scripts in the project scratchpad (`radial_cases.jl`,
  `reviews/scratch-02/rerun_cases.jl`) still use `Grid.RealGrid(L, …)` and `Luna.run` without
  `zmax`. They are not in the repository; they need updating before they are used against
  this branch.
- Carried over from `gpu/02-radialgrid`: `Test.detect_ambiguities(Luna; recursive=true)`
  reports 13 new latent pairs in `Luna.setup` and `NonlinearRHS.TransRadial`, all judged
  unreachable by the reviewer.
- Carried over from `gpu/01-zmax`: the pre-existing `Scans.QueueExec` queue-file bug, and the
  two `examples/low_level_interface/` scripts (`full_modal/basic_modal_full_bothpolarisations.jl`,
  `polarisation/elliptical_env.jl`) that are broken on `evanescent` as well.
