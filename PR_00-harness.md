# Regression gate, FFTW wisdom switch, benchmark environment

Branch `gpu/00-harness`, base `evanescent` (`fdf8dbe3`). Group A of GPU_PLAN.md.

## Motivation

The GPU project replaces most of Luna's per-step code. The agreed acceptance criterion
(GPU_PLAN.md §3) is that the default CPU path keeps producing today's output to rounding
level, gated per case against a baseline built from the base commit in the same
environment. That gate has to exist before any of the work it judges. This branch builds
it, together with the two pieces of infrastructure it depends on: a way to take the shared
FFTW wisdom cache out of the picture, and a benchmark script that gives every later branch
a before/after number.

Nothing here changes Luna's behaviour with default settings.

## What changed

### `src/Luna.jl`, `src/Utils.jl` — FFTW wisdom switch

`Luna.settings` gains `"fftw_wisdom" => true` and there is a new `Luna.set_fftw_wisdom(::Bool)`.
When the setting is `false`, `Utils.loadFFTwisdom` and `Utils.saveFFTwisdom` return
immediately, so neither `FFTW.import_wisdom` nor `FFTW.export_wisdom` is called, and
`set_fftw_wisdom(false)` additionally calls `FFTW.forget_wisdom()` once so that wisdom
already imported into the process does not reach plans made later.

The reason is in GPU_PLAN.md §2 and §4.11: the wisdom file lives in Luna's `Scratch`
directory, shared by every process using the same Julia depot — every worktree, every
concurrent agent, and Luna's own precompilation, which runs four `prop_capillary` calls at
`:patient`. `:patient` wisdom written by any of them changes the plan an `:estimate` run
gets, and with it the order of the floating-point operations in every transform. Without
this switch the regression gate compares two different FFT plans.

Defaults are unchanged.

### `docs/src/modules/Luna.md`

A "Global settings" section documenting `Luna.settings`, `set_fftw_mode`,
`set_fftw_threads` and `set_fftw_wisdom`. Two `@ref` cross-references into `Utils`, which
has no documentation page, became plain code spans so that the new `@docs` block does not
fail the documentation build. (`set_fftw_threads`'s reference to `Utils.FFTWthreads` was
already there; it was harmless only because the docstring was not in a `@docs` block.)

### `Project.toml`

`Adapt` (`"3, 4"`) and `GPUArraysCore` (`"0.1, 0.2"`) in `[deps]`, unused on this branch
and marked as such with a comment. `JLArrays` in `[extras]` and in the `test` target, for
the device-contract tests from `gpu/10-device-model` on. No Metal, no CUDA.

### `test/regression/` — the case matrix and the tooling

| File | What it is |
| --- | --- |
| `cases.jl` | `RegressionCases.CASES`, 15 cases, and `runcase(case, mode)`. |
| `compare.jl` | Storage format (HDF5) and the difference metric. |
| `run_cases.jl` | Runs everything and writes the HDF5 files. |
| `generate.jl` | Makes a worktree of a commit and runs `run_cases.jl` in it. |
| `sensitivity.jl` | The one-ulp sensitivity study. |
| `tolerances.jl` | Per-case, per-mode tolerances. |
| `README.md` | How to run all of it. |

The cases, all with `shotnoise=false` and `saveN=11`:

| Case | Geometry / physics |
| --- | --- |
| `modeavg_field_kerr` | mode-averaged field, Kerr only |
| `modeavg_field_plasma` | mode-averaged field, Kerr + PPT plasma (`PPT_options=Dict(:cache=>false)`) |
| `modeavg_field_raman` | mode-averaged field, Kerr + Raman (N₂, 0.5 bar) |
| `modeavg_field_mixture` | mode-averaged field, He/Ne mixture, low-level API |
| `modeavg_env_kerr` | mode-averaged envelope, Kerr only |
| `modeavg_env_raman` | mode-averaged envelope, Kerr + Raman |
| `gnlse_sech` | `prop_gnlse`, N = 2 sech soliton |
| `multimode_field_plasma` | 4 modes, field, Kerr + plasma |
| `radial_field_kerr` | radial free space, field |
| `radial_env_kerr` | radial free space, envelope |
| `free3d_env_kerr` | 3-D free space, envelope, `FreeGrid(1 mm, 16, 1 mm, 8)` |
| `free2d_field_chi2` | 2-D free space, field, `Chi2Field` in BBO |
| `gradient_field_kerr` | pressure gradient 2 → 0.1 bar |
| `taper_field_kerr` | core radius tapered 125 → 93.75 µm |
| `modeavg_field_legacy` | as `modeavg_field_kerr` but `boundary=:legacy` |

Every case runs in two modes:

- `:fixed` — `min_dz = max_dz = init_dz = zmax/20`. Equal limits make `RK45.steplims!`
  return the same step every time, so the controller is bypassed and a difference is
  attributable to the arithmetic. 20 steps is chosen so that the step equals the `:rate`
  absorber's reference length `grid.zmax/boundary_N`, which `Luna.run` would otherwise cap
  `max_dz` at.
- `:adaptive` — the controller, `init_dz = zmax/1000`. `zmax/1000` rather than Luna's
  absolute default `init_dz = 1e-4` so that every case actually adapts: for the 30 µm BBO
  case, `1e-4` is larger than `zmax` and would be clamped to `max_dz` on the first step.

`cases.jl` and `compare.jl` use only API that exists on the base commit, because
`generate.jl` copies them into a worktree of that commit. The baseline is therefore
produced by the base commit's Luna but by the *current* case definitions.
`run_cases.jl` copes with a Luna that has no `set_fftw_wisdom` by redefining
`Luna.Utils.loadFFTwisdom`/`saveFFTwisdom` to `nothing` and calling `FFTW.forget_wisdom()`
before any propagation.

Baselines go to `joinpath(Luna.Utils.cachedir(), "regression", <full sha>)` and are never
committed (`*.h5` is gitignored).

### `test/test_regression.jl` — the gate

Loads the baseline for `ENV["LUNA_REGRESSION_BASE"]`, or for the merge-base of `HEAD` with
`evanescent`, runs every case in both modes with the same settings, and asserts
`maximum(abs, Δ)/maximum(abs, baseline) <= tol` for `Eω`, for the save positions `z` and
for every compared statistic, per case and per mode. In the `:adaptive` mode `stats/z` and
`stats/dz` are not compared (see "What is compared in which mode"); in the `:fixed` mode
everything is. It prints the table of observed maxima whether it passes or not. It
is deliberately not in `runtests.jl`: it needs a baseline that does not exist on a fresh
checkout.

### `benchmark/`

`benchmark/Project.toml` (BenchmarkTools, Luna as a `dev` dependency) and `benchmark/run.jl`,
which for every regression case times one RHS evaluation `transform(nl, Eω, z)`, one
`RK45.step!` of the preconditioned stepper with the absorbing boundaries folded into the
operator, and one fixed-step propagation through `Luna.run`.

## Tests

| What | How | Result |
| --- | --- | --- |
| What | Result |
| --- | --- |
| `test/test_regression.jl` | **304 pass, 0 fail**. Every case, both modes, difference exactly `0.000e+00`. 66 s |
| `test/test_output.jl` | pass, 22.4 s |
| `test/test_interface.jl` | pass, 268.6 s |
| `test/test_freespace.jl` | pass, 313.9 s |
| every other `test/test_*.jl` | all 33 files exit 0; table under "Per-file test timings" |
| `test/regression/generate.jl fdf8dbe3` | baseline written, 15 cases × 2 modes |
| `test/regression/generate.jl HEAD` | baseline written from a second, unrelated commit — see "Re-baselining" |
| `LUNA_REGRESSION_BASE=211ab5ed test/test_regression.jl` | **304 pass, 0 fail**, all `0.000e+00` against that second baseline |
| `test/regression/sensitivity.jl` | table below |
| `benchmark/run.jl` | table below |
| `include("docs/make.jl")` | fails, identically to the base commit — see "Known gaps" 5 |

Commands:

```
julia --project=. -t 1 test/regression/generate.jl fdf8dbe3   # once per baseline commit
julia --project=. -t 1 test/test_regression.jl
julia --project=. -t 1 test/regression/sensitivity.jl
julia --project=benchmark -t 1 benchmark/run.jl

julia --project=. -t 1 -e 'using Luna, LinearAlgebra
  Luna.set_fftw_mode(:estimate); Luna.set_fftw_threads(1)
  LinearAlgebra.BLAS.set_num_threads(1)
  include("test/test_output.jl")'
```

`test/test_regression.jl` is not in `runtests.jl`, so `Pkg.test()` is unaffected.

### Re-baselining onto a new base

`generate.jl` takes any revision `git rev-parse` accepts. It resolves it to a full SHA,
makes (or reuses) a detached worktree of it under `../baselines/<sha[1:10]>`, copies this
worktree's `Manifest.toml` and *this branch's* `cases.jl`, `compare.jl` and `run_cases.jl`
into it, and runs `run_cases.jl` there in a fresh `julia -t 1` process: the old Luna, the
new case definitions. Exercised here for two different commits — the branch base
`fdf8dbe3` and the branch head `211ab5ed` — which land in separate directories and can both
be gated against.

To move the gate onto `gpu/int-A` once this branch is merged into it:

```
julia --project=. -t 1 test/regression/generate.jl gpu/int-A
LUNA_REGRESSION_BASE=gpu/int-A julia --project=. -t 1 test/test_regression.jl
```

`LUNA_REGRESSION_BRANCH=gpu/int-A` does the same thing via the merge-base. The step that
needs hands is **`cases.jl` itself**: `gpu/01-zmax` moves `zmax` out of the grids and
`gpu/02-radialgrid` replaces `Hankel.QDHT`, so the `Grid.RealGrid`/`Grid.EnvGrid`
constructor calls, the reads of `grid.zmax`, the `Luna.run` signature and
`Hankel.QDHT(R_FREE, 32, dim=3)` in `cases.jl` have to be updated as part of the `int-A`
merge. After that `cases.jl` no longer runs against `evanescent`, which is expected — that
is why the baseline commit is an argument rather than a constant. Record the deltas of the
new base against the old baseline before discarding it; `gpu/02-radialgrid`'s radial cases
are the ones expected to move. The full procedure is in `test/regression/README.md`.

## Measurements

Machine: Apple M1 Pro, 10 cores, macOS, Julia 1.13.0. Everything below with `-t 1`,
`Luna.set_fftw_mode(:estimate)`, `Luna.set_fftw_threads(1)`,
`LinearAlgebra.BLAS.set_num_threads(1)`, `Luna.set_fftw_wisdom(false)`. Other agents were
running propagations on the same machine throughout, so the wall-clock numbers are upper
bounds and repeatable only to some tens of percent.

### Regression gate against `fdf8dbe3`

Every case, both modes: `0.000e+00`, in every run. Largest difference over all cases and
modes: `0.000e+00`. Per-case deltas, `:fixed`/`:adaptive`, against the tolerances in
`tolerances.jl`:

| Case | `:fixed` Δ | tol | `:adaptive` Δ | tol |
| --- | --- | --- | --- | --- |
| `modeavg_field_kerr` | 0 | 1.0e-12 | 0 | 5.9e-07 |
| `modeavg_field_plasma` | 0 | 2.5e-12 | 0 | 7.3e-06 |
| `modeavg_field_raman` | 0 | 1.0e-12 | 0 | 2.4e-09 |
| `modeavg_field_mixture` | 0 | 1.0e-12 | 0 | 1.5e-10 |
| `modeavg_env_kerr` | 0 | 1.0e-12 | 0 | 5.3e-04 |
| `modeavg_env_raman` | 0 | 1.0e-12 | 0 | 4.2e-04 |
| `gnlse_sech` | 0 | 1.0e-12 | 0 | 1.6e-09 |
| `multimode_field_plasma` | 0 | 2.9e-07 | 0 | 7.7e-06 |
| `radial_field_kerr` | 0 | 1.0e-12 | 0 | 1.2e-08 |
| `radial_env_kerr` | 0 | 1.0e-12 | 0 | 8.2e-09 |
| `free3d_env_kerr` | 0 | 1.0e-12 | 0 | 1.0e-12 |
| `free2d_field_chi2` | 0 | 1.0e-12 | 0 | 7.4e-12 |
| `gradient_field_kerr` | 0 | 1.0e-12 | 0 | 1.6e-04 |
| `taper_field_kerr` | 0 | 1.0e-12 | 0 | 3.1e-05 |
| `modeavg_field_legacy` | 0 | 1.0e-12 | 0 | 7.8e-08 |

Exact zero, not "below tolerance", is a stronger statement than the brief asks for and it
also confirms that the baseline generator reproduces the run environment exactly: the
baseline came from a separate Julia process, in a separate worktree, from a separate
checkout. The same holds against a baseline generated from a second commit (`211ab5ed`, the
branch head): 326/326, all `0.000e+00`.

### One-ulp sensitivity

Initial `Eω` multiplied by `1 + eps()`; the number is the gate's metric maximised over
`Eω`, the save positions `z` and every compared statistic — which in the `:adaptive` mode
excludes `stats/z` and `stats/dz`, see "What is compared in which mode" below.

| Case | `:fixed` | driven by | `:adaptive` | driven by |
| --- | --- | --- | --- | --- |
| `modeavg_field_kerr` | 1.59e-15 | `peakpower` | 5.89e-09 | `peakintensity` |
| `modeavg_field_plasma` | 2.50e-14 | `peak_ionisation_rate` | 7.29e-08 | `peak_ionisation_rate` |
| `modeavg_field_raman` | 1.59e-15 | `peakpower` | 2.38e-11 | `peakintensity` |
| `modeavg_field_mixture` | 1.32e-15 | `Eω` | 1.48e-12 | `peakpower` |
| `modeavg_env_kerr` | 1.09e-15 | `peakintensity` | 5.34e-06 | `peakpower` |
| `modeavg_env_raman` | 1.87e-15 | `peakintensity` | 4.24e-06 | `peakintensity` |
| `gnlse_sech` | 2.11e-15 | `peakintensity` | 1.64e-11 | `peakintensity` |
| `multimode_field_plasma` | 2.95e-09 | `transverse_integral_error_rel` | 7.74e-08 | `peak_ionisation_rate` |
| `radial_field_kerr` | 2.05e-15 | `Eω` | 1.24e-10 | `Eω` |
| `radial_env_kerr` | 2.25e-15 | `Eω` | 8.18e-11 | `Eω` |
| `free3d_env_kerr` | 9.30e-16 | `Eω` | 1.25e-15 | `Eω` |
| `free2d_field_chi2` | 8.48e-16 | `Eω` | 7.37e-14 | `Eω` |
| `gradient_field_kerr` | 1.59e-15 | `peakpower` | 1.61e-06 | `zdw` |
| `taper_field_kerr` | 1.27e-15 | `Eω` | 3.06e-07 | `zdw` |
| `modeavg_field_legacy` | 1.87e-15 | `peakintensity` | 7.77e-10 | `peakpower` |

Tolerances in `tolerances.jl` are 100× these, floored at 1e-12.

### What is compared in which mode

`RegressionCompare.skipstats(mode)` decides. In the `:fixed` mode nothing is excluded: the
step sequence is imposed, so `stats/z` and `stats/dz` must match exactly, and checking them
also confirms that it really was imposed.

In the `:adaptive` mode `stats/z` and `stats/dz` (`RegressionCompare.STEP_STATS`) are
excluded. They record the step sequence, not the field. The controller's accept/reject
decision and its PI update respond to a one-ulp change far more strongly than the field
does, and the response compounds over the steps it takes to ramp `init_dz` up to `max_dz`.
With them in, the `:adaptive` sensitivity was set by `dz` in eight of the fifteen cases and
by `z` in a ninth, and the resulting tolerances were 7.3e-01 for `modeavg_env_kerr` and
5.9e-01 for `modeavg_env_raman` — no constraint on anything, applied to `Eω` as well.
Excluding them drops the loosest adaptive tolerance by a factor of about 1400, to 5.3e-04,
and leaves the numbers set by `Eω` and the physical statistics:

| | before | after |
| --- | --- | --- |
| loosest adaptive tolerance | 7.3e-01 (`modeavg_env_kerr`) | 5.3e-04 (`modeavg_env_kerr`) |
| next loosest | 5.9e-01 (`modeavg_env_raman`) | 4.2e-04 (`modeavg_env_raman`) |
| median adaptive tolerance | 1.6e-06 | 6.0e-08 |
| assertions in the gate | 326 | 304 |

`rundict`'s top-level `"z"` is a different thing — the save grid, fixed by
`Output.GridCondition` — and is compared in both modes. Nothing is excluded from what the
baseline *stores*, so the choice can be revisited without regenerating anything.

The two loosest adaptive cases, `gradient_field_kerr` (1.6e-04) and `taper_field_kerr`
(3.1e-05), are now set by `stats/zdw` and, for the gradient, `stats/density` and
`stats/pressure`. Those are properties of the medium at whatever z the stepper landed on,
so they inherit part of the step sequence's sensitivity indirectly, but they are not
step-sequence records and would catch a real change in the density or taper function, so
they stay in. `Eω` in those two cases is 2.6e-09 and 8.6e-09.

Cases above 1e-12 in the `:fixed` mode, as the brief asks:

- **`modeavg_field_plasma`, 2.50e-14, `stats/peak_ionisation_rate`.** The PPT rate is
  exponential in the field amplitude, so it amplifies a one-ulp change in the field by
  about an order of magnitude. `Eω` itself is at 1.4e-15.
- **`multimode_field_plasma`, 2.95e-09, `stats/transverse_integral_error_rel` and `_abs`.**
  These are `HCubature`'s own error estimates for the modal overlap integral in
  `TransModal`. The adaptive quadrature subdivides differently when the integrand changes
  in its last bits, so the error *estimate* moves by far more than the integral does. The
  next largest quantity in that case is `stats/mode_reconstruction_error` at 2.7e-13 and
  everything else is at rounding level.

In the `:adaptive` mode the number is set by `stats/dz` in eight of the fifteen cases and by
`stats/z` (the running sum of `dz`) in a ninth; the four free-space cases record no
statistics, so theirs is `Eω` itself, and it is at 1e-10 or below. The step-size
controller's accept/reject decision and its PI update respond to a one-ulp change far more
strongly than the field does, and the response compounds over the steps it takes to ramp
`init_dz` up to `max_dz`. `Eω` is three to six orders of magnitude tighter than the
whole-case number: for `modeavg_field_kerr` the case number is 8.2e-06 while
`peakintensity`, the largest non-step quantity, is 5.9e-09.

**This is the weak point of the gate as specified** — see "Open questions" below.

### Benchmarks (`benchmark/run.jl`)

`state` is the number of elements in `Eω`. `rhs` and `step` are BenchmarkTools minima,
`prop` the minimum of three fixed-step (20-step) propagations including setup and output.

| Case | state | rhs | step | prop |
| --- | ---: | ---: | ---: | ---: |
| `modeavg_field_kerr` | 2049 | 41.2 µs | 495 µs | 36.7 ms |
| `modeavg_field_plasma` | 2049 | 112 µs | 927 µs | 5.60 s |
| `modeavg_field_raman` | 2049 | 251 µs | 1.76 ms | 72.9 ms |
| `modeavg_field_mixture` | 2049 | 42.5 µs | 580 µs | 39.5 ms |
| `modeavg_env_kerr` | 2048 | 21.1 µs | 366 µs | 31.1 ms |
| `modeavg_env_raman` | 2048 | 91.8 µs | 815 µs | 45.5 ms |
| `gnlse_sech` | 2048 | 20.7 µs | 516 µs | 31.4 ms |
| `multimode_field_plasma` | 8196 | 4.86 ms | 30.5 ms | 6.47 s |
| `radial_field_kerr` | 4128 | 156 µs | 1.68 ms | 63.3 ms |
| `radial_env_kerr` | 4096 | 126 µs | 1.54 ms | 56.2 ms |
| `free3d_env_kerr` | 16384 | 319 µs | 5.19 ms | 133 ms |
| `free2d_field_chi2` | 16448 | 1.87 ms | 14.4 ms | 508 ms |
| `gradient_field_kerr` | 2049 | 80.8 µs | 1.16 ms | 235 ms |
| `taper_field_kerr` | 2049 | 309 µs | 3.65 ms | 101 ms |
| `modeavg_field_legacy` | 2049 | 41.1 µs | 496 µs | 35.2 ms |

Notes on the numbers:

- `step` is close to 6 × `rhs` plus the propagator, as expected for a 6-stage FSAL
  evaluation.
- The two plasma cases' `prop` is dominated by setup, not by stepping: they pass
  `PPT_options=Dict(:cache=>false)` (so that the case does not write to the PPT rate cache
  shared with every other worktree), which makes `makePPTcache` recompute the rate table on
  every `prepare()`, about 5 s. 20 steps of `modeavg_field_plasma` is 19 ms.
- `gradient_field_kerr` and `taper_field_kerr` are slower per step than the constant case
  because their operator is z-dependent: it is rebuilt at every stage.
- `multimode_field_plasma`'s RHS is two orders of magnitude above the mode-averaged cases:
  it is the adaptive `HCubature` modal integral at every evaluation, and it is the main
  target of §4.4.

### Per-file test timings

Fresh `julia --project=. -t 1` process per file, `:estimate`, one FFTW and one BLAS
thread, `using Luna` (≈5 s with a warm precompile cache) included in every number. All 33
files exit 0. Measured while two other agents were running propagations and a full
`Pkg.test()` on the same 10-core machine, so these are upper bounds: the short files are
dominated by load time and are reliable, the long ones are probably inflated, perhaps by a
factor of two.

Total 2830 s ≈ 47 min, of which `test_scans.jl`, `test_polarisation_field.jl`,
`test_freespace.jl`, `test_interface.jl` and `test_multimode.jl` are 1839 s (65%).

| File | s | File | s | File | s |
| --- | ---: | --- | ---: | --- | ---: |
| `test_scans.jl` | 714.9 | `test_boundaries.jl` | 33.0 | `test_maths.jl` | 15.5 |
| `test_polarisation_field.jl` | 695.1 | `test_gnlse.jl` | 29.0 | `test_polarisation_env.jl` | 13.5 |
| `test_freespace.jl` | 313.9 | `test_tapers.jl` | 24.1 | `test_mixtures.jl` | 13.1 |
| `test_interface.jl` | 268.6 | `test_noise.jl` | 22.7 | `test_physdata.jl` | 12.3 |
| `test_multimode.jl` | 246.1 | `test_output.jl` | 22.4 | `test_rk45.jl` | 11.9 |
| `test_fields.jl` | 66.5 | `test_linops.jl` | 18.2 | `test_capillary.jl` | 11.8 |
| `test_vectorplasma.jl` | 43.9 | `test_gradient.jl` | 17.5 | `test_stats.jl` | 11.4 |
| `test_modes.jl` | 37.0 | `test_kerr.jl` | 10.3 | `test_chi2.jl` | 9.9 |
| `test_processing.jl` | 35.4 | `test_linearprop.jl` | 9.6 | `test_antiresonant.jl` | 9.3 |
| `test_ionisation.jl` | 33.8 | `test_utils.jl` | 9.1 | `test_rect_modes.jl` | 8.6 |
| | | `test_raman.jl` | 8.2 | `test_tools.jl` | 7.5 |
| | | `test_polarisation.jl` | 5.7 | | |

`test/test_regression.jl` itself is 69 s, plus the cost of generating the baseline: about
2 min of propagations, and a full precompilation of Luna in the baseline worktree the first
time a given commit is used.

For the branch briefs: a branch touching `NonlinearRHS`, `LinearOps` or `Luna.run` wants
`test_freespace.jl`, `test_multimode.jl`, `test_interface.jl` and the regression gate,
which is about 15 min. `test_scans.jl` and `test_polarisation_field.jl` are worth running
only on branches that touch `Scans.jl` or polarisation, and on the integration branches.

### FFTW threading on Julia 1.13

GPU_PLAN.md §2 records a report that FFTW.jl's task-based threading callback
(`FFTW/src/providers.jl:56-78`) segfaults on Julia ≥ 1.12, and that the fork deregisters
it. **This did not reproduce here.** FFTW.jl v1.10.0, `fftw_provider == "fftw"`, Julia
1.13.0, `julia -t 4`, `Luna.set_fftw_threads(4)` (so `Utils.FFTWthreads() == 4` and the
`spawnloop` callback is registered, since it is registered whenever `nthreads() > 1`):

- a mode-averaged capillary propagation (8192-point grid) ran to completion, 23 steps, no
  crash;
- a 3-D free-space envelope propagation on `FreeGrid(2 mm, 64, 2 mm, 64)` with an
  8192-point time grid — the largest multi-axis plans Luna makes — ran three times at
  `:estimate` (20.1, 16.7, 18.1 s) and three times at `:patient` (28.8, 28.6, 24.7 s), six
  propagations, no crash.

No workaround adopted, as instructed. If it does appear later it will most likely be
platform- or FFTW-build-dependent, so it is worth re-checking on the Linux CI runners.

## Known gaps and deviations from GPU_PLAN.md

1. **The free-space cases record no statistics.** The brief asks for "default statistics
   on" for every case. There is no `Stats.default` method for free-space geometries —
   `Stats.default` dispatches on `Modes.AbstractMode` and `Modes.ModeCollection` only — and
   the individual stats functions do not all work on a 3-D `Eω`: `Stats.ω0`'s `squeeze`
   has methods for 1- and 2-dimensional arrays only, and `Stats.energy` indexes `Eω[:, i]`.
   The four free-space cases (`radial_field_kerr`, `radial_env_kerr`, `free3d_env_kerr`,
   `free2d_field_chi2`) therefore compare `Eω` and `z` only, which is what the free-space
   examples do (their `Stats.collect_stats` lines are commented out). Adding free-space
   statistics is a change to `Stats.jl` and out of scope here; it would strengthen the gate
   and is worth a separate issue.

2. **The `:adaptive` mode does not compare `stats/z` and `stats/dz`.** This is a deviation
   from the brief, made deliberately and on instruction. With them in, the rule "one
   tolerance per case per mode, 100× the measured sensitivity" produced adaptive tolerances
   of 0.73 and 0.59 for the two envelope Kerr cases, which constrained nothing. They are
   still compared, and must be exact, in the `:fixed` mode. See "What is compared in which
   mode". The residual looseness is `stats/zdw`/`density`/`pressure` in the gradient and
   taper cases, which are medium properties sampled at the step positions; they are kept
   because they would catch a real change in the density or taper function.

3. **`Luna.prop_capillary_args` does not exist** under that name: the `*_args` functions
   live in `Luna.Interface` and are not re-exported (only `prop_capillary` and `prop_gnlse`
   are, `Luna.jl:108-109`). `cases.jl` calls `Interface.prop_capillary_args` and
   `Interface.prop_gnlse_args`. Not a bug, just a correction to the brief.

4. **`Interface.boundary_kwargs` is duplicated in `cases.jl`.** The absorbing-boundary
   options have to go both to `prop_capillary_args` (which records them in `saveargs`) and
   to `Luna.run` (which uses them); `prop_capillary` does that with a private helper.
   `cases.jl` repeats the four-key list rather than calling the helper, so that the file
   keeps working on older commits where the helper may differ. If the key list grows, both
   copies need updating.

5. **The documentation build fails, identically, with and without this branch.**
   `include("docs/make.jl")` terminates at `[:cross_references]` with six unresolved
   `@ref`s — three `LinearOps.βz`, `Luna.PhysData.crystal_internal_angle`, `norm_free`,
   `LinearOps.make_const_linop`, in `docs/src/modules/Boundaries.md` and
   `NonlinearRHS.md`. The same six, and only those six, appear on the base commit
   `fdf8dbe3` in its own worktree, so they are pre-existing and not caused by this branch;
   the new `@docs` block in `Luna.md` produces none. Not fixed here (COMMON.md: report
   bugs, do not fix them outside the brief) — worth an issue, because it means the docs
   have not built on `evanescent` for a while.

6. **`benchmark/Project.toml` has no `[sources]` entry.** `Pkg.develop(path=".")` adds one
   with an absolute path, which is machine-specific; the file says so and the committed
   version does not contain it. Consequence: a bare `Pkg.instantiate()` in `benchmark/`
   without the `develop` step first would resolve Luna from the registry. The documented
   setup line does the `develop` first.

7. **`benchmark/` has no `SUITE`/PkgBenchmark layout.** `CLAUDE.md` in the main working
   copy describes a `benchmark/benchmarks.jl` with a `SUITE` indexed by integrator and a
   `workprecision.jl`; neither exists on `evanescent` (GPU_PLAN.md §1 notes this), and the
   brief asks only for `run.jl`. `CLAUDE.md` is corrected in `gpu/32-docs`.

## Open questions

- **Per-quantity tolerances.** Still one tolerance per (case, mode), as the brief
  specifies. Excluding the step diagnostics from the `:adaptive` comparison has made that
  good enough — the loosest tolerance is 5.3e-04 and the median 6.0e-08 — so a `Dict`
  keyed by quantity as well is no longer needed. It would still be the way to tighten
  `gradient_field_kerr` and `taper_field_kerr`, whose numbers are set by `stats/zdw` rather
  than by the field.
- **Re-baselining.** `gpu/int-A` has to re-baseline for `gpu/01-zmax`'s changed signatures.
  `generate.jl` takes the commit as an argument and `cases.jl` is copied from the branch
  under test, so the mechanism is there, but `cases.jl` itself will need its `Luna.run` and
  grid constructor calls updated at that point, and it then stops being runnable against
  `evanescent`. The plan already anticipates this (§6, "regression gate re-baselined here
  for the changed signatures").
- **CI.** The gate is not wired into CI, because a baseline has to be generated first and
  takes a worktree plus a full precompilation. Worth deciding before `gpu/int-A` whether
  the integration branches run it by hand (as now) or whether CI generates a baseline from
  the merge-base.
