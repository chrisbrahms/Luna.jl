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
| `cases.jl` | `RegressionCases.CASES`, 20 cases, and `runcase(case, mode)`. |
| `compare.jl` | Storage format (HDF5), the difference metric and the tolerance classes. |
| `run_cases.jl` | Runs everything and writes the HDF5 files. |
| `generate.jl` | Makes a worktree of a commit and runs `run_cases.jl` in it. |
| `sensitivity.jl` | The one-ulp sensitivity study. |
| `tolerances.jl` | Per-case, per-mode, per-class tolerances. |
| `README.md` | How to run all of it. |

The cases, all with `shotnoise=false` and `saveN=11`:

| Case | Geometry / physics |
| --- | --- |
| `modeavg_field_kerr` | mode-averaged field, Kerr only |
| `modeavg_field_plasma` | mode-averaged field, Kerr + PPT plasma (`PPT_options=Dict(:cache=>false)`) |
| `modeavg_field_raman` | mode-averaged field, Kerr + Raman (N₂, 0.5 bar) |
| `modeavg_field_mixture` | mode-averaged field, He/Ne mixture, low-level API |
| `modeavg_field_adk` | mode-averaged field, Kerr + **ADK** plasma (`plasma=:ADK`) |
| `modeavg_field_vector` | **elliptically polarised** (ε = 0.5) field, vector Kerr + plasma, two polarisation components |
| `modeavg_env_kerr` | mode-averaged envelope, Kerr only |
| `modeavg_env_raman` | mode-averaged envelope, Kerr + Raman |
| `modeavg_env_thg` | mode-averaged envelope, **`Kerr_env_thg`** (`thg=true`) |
| `gnlse_sech` | `prop_gnlse`, N = 2 sech soliton, no Raman, no shock |
| `gnlse_raman_shock` | `prop_gnlse`, **`raman=true, shock=true`** |
| `multimode_field_plasma` | 4 modes, field, Kerr + plasma |
| `radial_field_kerr` | radial free space, field |
| `radial_env_kerr` | radial free space, envelope |
| `free3d_env_kerr` | 3-D free space, envelope, `FreeGrid(1 mm, 16, 1 mm, 8)` |
| `free2d_field_chi2` | 2-D free space, field, `Chi2Field` in BBO |
| `free2d_env_chi2` | 2-D free space, envelope, **`Chi2Env`** in BBO |
| `gradient_field_kerr` | pressure gradient 2 → 0.1 bar |
| `taper_field_kerr` | core radius tapered 125 → 93.75 µm |
| `modeavg_field_legacy` | as `modeavg_field_kerr` but `boundary=:legacy` |

The five in bold were added after review round 1 (finding 8) so that the gate covers the
responses Groups C and D rewrite: `Kerr_env_thg`, `Chi2Env`, ADK ionisation, the vector
response path and `prop_gnlse`'s Raman and self-steepening branches.

All four free-space cases now record `z` and `dz` through `Stats.collect_stats(grid, Eω)`,
so their step sequence is checked too. There is still no `Stats.default` for free-space
geometries (known gap 1).

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
for every compared statistic, per case and per mode. There are two tolerances per case and
mode, one for `Eω` and one for the statistics (see "Two tolerance classes"). In the
`:adaptive` mode `stats/z` and `stats/dz` are not compared, but the number of accepted steps
is checked in both modes. It prints a table of observed maxima whether it passes or not,
with the (component, save) index the `Eω` maximum came from. It is deliberately not in
`runtests.jl`: it needs a baseline that does not exist on a fresh checkout.

### `benchmark/`

`benchmark/Project.toml` (BenchmarkTools, and Luna through a relative `[sources]` entry, so
a bare `Pkg.instantiate()` in `benchmark/` is enough) and `benchmark/run.jl`,
which for every regression case times one RHS evaluation `transform(nl, Eω, z)`, one
`RK45.step!` of the preconditioned stepper with the absorbing boundaries folded into the
operator, and one fixed-step propagation through `Luna.run`.

## Tests

| What | Result |
| --- | --- |
| `test/test_regression.jl` | **436 pass, 0 fail**. Every case, both modes, both classes, difference exactly `0.000e+00`. 90 s |
| `test/test_utils.jl` | **33 pass, 0 fail**, 3.7 s (includes the new `set_fftw_wisdom` testset) |
| `test/test_output.jl` | pass, 22.4 s |
| `test/test_interface.jl` | pass, 268.6 s |
| `test/test_freespace.jl` | pass, 313.9 s |
| every other `test/test_*.jl` | all 33 files exit 0; table under "Per-file test timings" |
| `test/regression/generate.jl fdf8dbe3` | baseline written, 20 cases × 2 modes |
| `test/regression/generate.jl <other commit>` | exercised for two further commits — see "Re-baselining" |
| `test/regression/sensitivity.jl` | table below |
| `benchmark/run.jl` | table below |
| `include("docs/make.jl")` | fails, identically to the base commit — see "Known gaps" 5 |

The per-file timings below predate review round 1; only `test_utils.jl` changed since (it
gained a testset that adds 0.1 s), and `test_regression.jl` went from 66 s to 90 s with the
five new cases.

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
new case definitions. Exercised for three different commits: the branch base `fdf8dbe3`, an
intermediate commit `211ab5ed`, and (by the reviewer, independently, into a redirected
output directory) `da01cb1a`. Each lands in its own directory and can be gated against with
`LUNA_REGRESSION_BASE`.

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

### The metric

Per quantity, `maximum(abs, Δ)/maximum(abs, baseline)`. Elementwise relative differences are
meaningless in the window tapers and outside `grid.sidx`.

**`Eω` is normalised per component and per save.** `RegressionCompare.fielddiff` reduces the
frequency axis and any transverse axes away and forms the ratio separately for each (mode or
polarisation, save) pair; the reported value is the largest, and the gate prints which pair
it came from. `Eω` is `(Nω, Nz)`, `(Nω, Nm, Nz)`, `(Nω, Npol, Nk, Nz)` or
`(Nω, Npol, Nkx, Nky, Nz)`, so axis 1 is frequency, the last axis is the save, axis 2 is the
component when there are more than two axes, and the rest are transverse.

One global normalisation hid weak components. In `multimode_field_plasma` the per-mode
maxima of `|Eω|` at the last save are `[2.9e5, 1.8e1, 3.2e0, 7.2e-2]`, a spread of 4e6: a
change that rewrote mode 4 entirely would have sat four million times below a tolerance set
by mode 1. That case is the only multimode one in the matrix and `gpu/22-modal` rewrites the
transform it exercises.

**Two tolerance classes.** `RegressionCompare.classof` puts each quantity into `:Eω` or
`:stats` (the save grid `z`, the statistics, the step-count check), and `tolerances.jl`
carries one tolerance per (case, mode, class). The statistics are recorded once per accepted
step and so inherit the step-size controller's sensitivity in the adaptive mode; sharing one
tolerance with them gave `Eω` up to four orders of magnitude of slack.

**The step count is checked explicitly.** Every statistic is recorded per accepted step, so a
change in the number of steps would fail every statistic's size check with an uninformative
`Inf`. `compare` checks `length(stats["z"])` first and, if it differs, reports a single
`"step count"` entry (`changed: N steps, baseline M`) and compares no statistics. Hard
failure in both modes.

**`stats/z` and `stats/dz` are excluded in `:adaptive` only.** They record the step sequence,
not the field, and the controller responds to a one-ulp change far more strongly than the
field does; with them in, the tolerance for `modeavg_env_kerr` came out at 0.73. In `:fixed`
they are compared and must be exact, which is also the check that the step sequence really
was imposed. Nothing is excluded from what the baseline *stores*, so the choice can be
revisited without regenerating.

### Regression gate against `fdf8dbe3`

**436 pass, 0 fail.** Every case, every mode, every class: `0.000e+00`. Largest difference
over all cases and modes: `0.000e+00`. The tolerances are from `tolerances.jl`; Δ is the
observed difference.

| Case | `:fixed` Eω Δ / tol | `:fixed` stats Δ / tol | `:adaptive` Eω Δ / tol | `:adaptive` stats Δ / tol |
| --- | --- | --- | --- | --- |
| `modeavg_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 5.9e-07 |
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
| `radial_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.2e-08 | 0 / 1.0e-12 |
| `radial_env_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 8.2e-09 | 0 / 1.0e-12 |
| `free3d_env_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 |
| `free2d_field_chi2` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 7.4e-12 | 0 / 1.0e-12 |
| `free2d_env_chi2` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 5.0e-12 | 0 / 1.0e-12 |
| `gradient_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.6e-07 | 0 / 1.6e-04 |
| `taper_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 8.7e-07 | 0 / 3.1e-05 |
| `modeavg_field_legacy` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.3e-11 | 0 / 7.8e-08 |

Exact zero, not "below tolerance", is a stronger statement than the brief asks for, and it
also confirms that the baseline generator reproduces the run environment exactly: the
baseline came from a separate Julia process, in a separate worktree, from a separate
checkout. The same held for the two other baseline commits it has been run against
(`211ab5ed` and, by the reviewer, `da01cb1a`), before the five new cases were added.

### One-ulp sensitivity

Initial `Eω` multiplied by `1 + eps()`; the numbers are the gate's metric maximised over
each tolerance class, using the gate's own comparison (`:adaptive` excludes `stats/z` and
`stats/dz`). Tolerances are 100× these, floored at 1e-12.

| Case | `:fixed` Eω | `:fixed` stats | driven by | `:adaptive` Eω | `:adaptive` stats | driven by |
| --- | --- | --- | --- | --- | --- | --- |
| `modeavg_field_kerr` | 8.60e-16 | 1.59e-15 | `peakpower` | 1.45e-15 | 5.89e-09 | `peakintensity` |
| `modeavg_field_plasma` | 1.40e-15 | 2.50e-14 | `peak_ionisation_rate` | 2.20e-15 | 7.29e-08 | `peak_ionisation_rate` |
| `modeavg_field_raman` | 1.41e-15 | 1.59e-15 | `peakpower` | 2.26e-15 | 2.38e-11 | `peakintensity` |
| `modeavg_field_mixture` | 1.32e-15 | 9.52e-16 | `Eω` | 2.86e-15 | 1.48e-12 | `peakpower` |
| `modeavg_field_adk` | 1.40e-15 | 2.27e-13 | `peak_ionisation_rate` | 2.20e-15 | 3.81e-07 | `peak_ionisation_rate` |
| `modeavg_field_vector` | 1.24e-15 | 1.84e-11 | `transverse_integral_error_rel` | 1.79e-15 | 9.78e-07 | `peak_ionisation_rate` |
| `modeavg_env_kerr` | 8.84e-16 | 1.09e-15 | `peakintensity` | 5.61e-14 | 5.34e-06 | `peakpower` |
| `modeavg_env_raman` | 1.15e-15 | 1.87e-15 | `peakintensity` | 2.92e-14 | 4.24e-06 | `peakintensity` |
| `modeavg_env_thg` | 1.58e-15 | 1.79e-15 | `peakpower` | 1.92e-15 | 2.26e-09 | `peakpower` |
| `gnlse_sech` | 8.14e-16 | 2.11e-15 | `peakintensity` | 3.61e-13 | 1.64e-11 | `peakintensity` |
| `gnlse_raman_shock` | 9.65e-16 | 1.87e-15 | `peakintensity` | 1.90e-12 | 1.25e-10 | `peakintensity` |
| `multimode_field_plasma` | 1.16e-13 | 2.95e-09 | `transverse_integral_error_rel` | 1.80e-09 | 7.74e-08 | `peak_ionisation_rate` |
| `radial_field_kerr` | 2.05e-15 | 0 | `Eω` | 1.24e-10 | 0 | `Eω` |
| `radial_env_kerr` | 2.25e-15 | 0 | `Eω` | 8.18e-11 | 0 | `Eω` |
| `free3d_env_kerr` | 9.30e-16 | 0 | `Eω` | 1.25e-15 | 0 | `Eω` |
| `free2d_field_chi2` | 1.24e-15 | 0 | `Eω` | 7.38e-14 | 0 | `Eω` |
| `free2d_env_chi2` | 1.13e-15 | 0 | `Eω` | 5.04e-14 | 0 | `Eω` |
| `gradient_field_kerr` | 1.24e-15 | 1.59e-15 | `peakpower` | 2.58e-09 | 1.61e-06 | `zdw` |
| `taper_field_kerr` | 1.28e-15 | 1.10e-15 | `Eω` | 8.72e-09 | 3.06e-07 | `zdw` |
| `modeavg_field_legacy` | 1.02e-15 | 1.87e-15 | `peakintensity` | 1.29e-13 | 7.77e-10 | `peakpower` |

"driven by" is the quantity that sets the `stats` number; the free-space cases record only
`z` and `dz`, which are exact in `:fixed` and excluded in `:adaptive`, so their `stats`
sensitivity is zero and the case is carried entirely by `Eω`.

Every `:fixed` number above the 1e-12 floor, as the brief asks:

- **`multimode_field_plasma`, Eω 1.16e-13, in mode 4 of 4 at save 8.** With the old global
  normalisation this was invisible; the per-component metric measures it. The four modes'
  maxima span 4e6, so mode 4 is where any modal-transform change will show first.
- **`modeavg_field_plasma` 2.50e-14 and `modeavg_field_adk` 2.27e-13, `peak_ionisation_rate`.**
  Both ionisation rates are exponential in the field amplitude. ADK amplifies about ten
  times more than PPT, which goes through a spline.
- **`modeavg_field_vector` 1.84e-11 and `multimode_field_plasma` 2.95e-09,
  `transverse_integral_error_rel`/`_abs`.** `HCubature`'s own error estimates for the modal
  overlap integral: the adaptive quadrature subdivides differently when the integrand moves
  in its last bits, so the error *estimate* moves far more than the integral. In
  `multimode_field_plasma` the next largest statistic is `mode_reconstruction_error` at
  2.74e-13.

In the `:adaptive` mode the `stats` class is loose by design (up to 5.34e-06 raw, 5.3e-04 as
a tolerance) because every statistic is recorded per accepted step. The `Eω` class is two to
four orders of magnitude tighter, which is the point of splitting the classes; the loosest
are the two z-dependent cases, `taper_field_kerr` 8.72e-09 and `gradient_field_kerr`
2.58e-09, where the operator is rebuilt at every stage so the step positions feed back into
the field itself.

### What review round 1 changed here

Before the changes below, `Eω` and the statistics shared one tolerance per (case, mode), and
`Eω` was normalised by the maximum over the whole array including the save axis.

| | before | after |
| --- | --- | --- |
| loosest `Eω` tolerance, `:adaptive` | 5.3e-04 (shared with `stats`) | 8.7e-07 |
| `Eω` tolerance, `gradient_field_kerr` `:adaptive` | 1.6e-04 | 2.6e-07 |
| `Eω` sensitivity, `multimode_field_plasma` `:fixed` | not separately measured; whole-array metric dominated by mode 1 | 1.16e-13, located in mode 4 |
| assertions in the gate | 304 | 436 |
| cases | 15 | 20 |

### Benchmarks (`benchmark/run.jl`)

`state` is the number of elements in `Eω`. `rhs` and `step` are BenchmarkTools minima,
`prop` the minimum of three fixed-step (20-step) propagations including setup and output.

| Case | state | rhs | step | prop |
| --- | ---: | ---: | ---: | ---: |
| `modeavg_field_kerr` | 2049 | 41.2 µs | 495 µs | 36.7 ms |
| `modeavg_field_plasma` | 2049 | 112 µs | 927 µs | 5.60 s |
| `modeavg_field_raman` | 2049 | 251 µs | 1.76 ms | 72.9 ms |
| `modeavg_field_mixture` | 2049 | 42.5 µs | 580 µs | 39.5 ms |
| `modeavg_field_adk` | 2049 | 133 µs | 1.04 ms | 55.0 ms |
| `modeavg_field_vector` | 4098 | 5.03 ms | 31.2 ms | 6.09 s |
| `modeavg_env_kerr` | 2048 | 21.1 µs | 366 µs | 31.1 ms |
| `modeavg_env_raman` | 2048 | 91.8 µs | 815 µs | 45.5 ms |
| `modeavg_env_thg` | 2048 | 40.3 µs | 467 µs | 35.6 ms |
| `gnlse_sech` | 2048 | 20.7 µs | 516 µs | 31.4 ms |
| `gnlse_raman_shock` | 2048 | 77.4 µs | 864 µs | 42.2 ms |
| `multimode_field_plasma` | 8196 | 4.86 ms | 30.5 ms | 6.47 s |
| `radial_field_kerr` | 4128 | 156 µs | 1.68 ms | 63.3 ms |
| `radial_env_kerr` | 4096 | 126 µs | 1.54 ms | 56.2 ms |
| `free3d_env_kerr` | 16384 | 319 µs | 5.19 ms | 133 ms |
| `free2d_field_chi2` | 16448 | 1.87 ms | 14.4 ms | 508 ms |
| `free2d_env_chi2` | 32768 | 3.28 ms | 25.6 ms | 828 ms |
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
- `modeavg_field_vector` costs as much as the four-mode case despite having two components:
  elliptical polarisation turns a mode-averaged run into a two-mode `TransModal`, so it pays
  the adaptive modal integral. Its `prop` is setup-dominated for the same PPT reason as the
  other plasma cases.
- `free2d_env_chi2` has twice the state of `free2d_field_chi2` because the envelope grid is
  built with `thg=true`, and costs about 1.8x as much per step.
- The ten cases whose numbers are unchanged from the pre-review run were re-measured for
  two of them (`modeavg_env_kerr` 21.2 µs / 358 µs / 32.5 ms, `free2d_field_chi2` 1.85 ms /
  14.1 ms / 496 ms) to confirm the boundary-keyword fix changed nothing.
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

## Changes made after review round 1

Review: `reviews/gpu-00-harness-1.md`, verdict "approve with minor fixes". All eleven
actionable findings addressed; findings 12 and 13 were report-only and are reflected in the
FFTW and documentation sections above.

| # | Finding | What changed |
| --- | --- | --- |
| 1 | `Eω` normalised by the global maximum, so weak modes were unchecked | `RegressionCompare.fielddiff`: normalise per component and per save, report the index. `multimode_field_plasma`'s `:fixed` `Eω` sensitivity is now measurable at 1.16e-13, in mode 4 |
| 2 | one tolerance per (case, mode) left `Eω` orders of magnitude of slack | two tolerance classes, `:Eω` and `:stats`, keyed per case and mode; `sensitivity.jl` measures and prints both |
| 3 | a change in the step count would fail every statistic with an uninformative `Inf` | `compare` checks `length(stats["z"])` first and reports one named `"step count"` failure; hard failure in both modes. README now states how few controller steps the adaptive mode really exercises |
| 4 | `loadFFTwisdom`'s early return also skipped `FFTW.set_num_threads` | moved above the guard, docstring updated |
| 5 | `Scans` workers do not inherit the setting | said so in the `set_fftw_wisdom` docstring, for all three settings |
| 6 | no test for `set_fftw_wisdom` | new testset in `test/test_utils.jl`: default, both functions return `nothing`, cache mtime unchanged, no pidlock, restore |
| 7 | free-space cases recorded no statistics, so their step sequence was unchecked | all five now use `Stats.collect_stats(grid, Eω)`, recording `z` and `dz` |
| 8 | response coverage gaps before Groups C and D | five new cases: `gnlse_raman_shock`, `modeavg_env_thg`, `free2d_env_chi2`, `modeavg_field_adk`, `modeavg_field_vector` |
| 9 | `benchmark/run.jl` dropped all boundary keywords but `:boundary` | all four mapped onto `Boundaries.setup`; unused `linop`/`dz` bindings removed |
| 10 | `Project.toml` and PR staleness | `JLArrays = "0.1, 0.2, 0.3"` (current release 0.3.3); stale "branch head", the 326-vs-304 paragraph and the duplicated table header fixed |
| 11 | nits | `-` printed when the maximum is zero; soliton comment corrected (N = 2 is second order, `GNLSE_LENGTH` is 0.2 soliton periods); `runcase` now suppresses `Info` and below rather than everything, so warnings reach the gate output; `benchmark/Project.toml` carries a relative `[sources]` entry |

On finding 11's `jldoctest` note: `set_fftw_mode`'s docstring contains a doctest that
Documenter will execute now that the function is in a `@docs` block, mutating
`settings["fftw_flag"]` during the docs build. It should pass (it asserts `0x00000020`,
`FFTW.PATIENT`, which is the default). It could not be confirmed, because the docs build
fails earlier on the pre-existing cross-reference errors in known gap 5. No change made.

Everything in this section is in the commits after `1effca40`. The gate is still exactly
zero, on a baseline regenerated from `fdf8dbe3` with the new `cases.jl` and `compare.jl`.

## Known gaps and deviations from GPU_PLAN.md

1. **The free-space cases record only `z` and `dz`, not physical statistics.** The brief
   asks for "default statistics on" for every case. There is no `Stats.default` method for
   free-space geometries — it dispatches on `Modes.AbstractMode` and `Modes.ModeCollection`
   only — and the individual functions do not work on a 3- or 4-dimensional `Eω`:
   `Stats.ω0`'s `squeeze` has 1- and 2-dimensional methods only, and `Stats.energy` and
   `Stats.peakpower` index `Eω[:, i]`. Since review round 1 the five free-space cases do use
   `Stats.collect_stats(grid, Eω)`, which appends `Stats.zdz!` and works on all three
   geometries, so they now record `z` and `dz` — enough for the fixed-step check that the
   step sequence was imposed and for the step-count check in the adaptive mode. Energy,
   peak power and beam size are still not checked there. Adding them is a change to
   `Stats.jl`, out of scope here, and worth a separate issue.

2. **The `:adaptive` mode does not compare `stats/z` and `stats/dz`, and the tolerances are
   keyed by class rather than by quantity.** Both are deviations from the brief, made
   deliberately and on instruction. With the step diagnostics in, the rule "one tolerance
   per case per mode, 100× the measured sensitivity" produced adaptive tolerances of 0.73
   and 0.59 for the two envelope Kerr cases, which constrained nothing. They are still
   compared, and must be exact, in `:fixed`, and the step *count* is a hard failure in both
   modes. The residual looseness in the `:stats` class is `stats/zdw`/`density`/`pressure`
   in the gradient and taper cases, which are medium properties sampled at the step
   positions; they are kept because they would catch a real change in the density or taper
   function, and they no longer contaminate the `Eω` tolerance.

   Still outstanding: the `:adaptive` `Eω` tolerances for `taper_field_kerr` (8.7e-07) and
   `gradient_field_kerr` (2.6e-07) are two to three orders of magnitude looser than the
   rest. That is a real property of those cases — the operator is rebuilt at every stage, so
   the step positions feed into the field — not an artefact of the metric, and the `:fixed`
   mode holds both to 1e-12.

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

6. **`benchmark/Project.toml` now carries a relative `[sources]` entry** (`Luna = {path = ".."}`),
   so a bare `Pkg.instantiate()` in `benchmark/` is enough and nothing writes an absolute
   path into a tracked file. That form needs Julia >= 1.11, which is why the benchmark
   environment's `julia` compat is 1.11 while Luna's stays at 1.9; the benchmark environment
   is a developer tool and is not part of what `Pkg.add("Luna")` resolves. Verified by
   deleting `benchmark/Manifest.toml` and instantiating from scratch. The older documented
   setup line does the `develop` first.

7. **`benchmark/` has no `SUITE`/PkgBenchmark layout.** `CLAUDE.md` in the main working
   copy describes a `benchmark/benchmarks.jl` with a `SUITE` indexed by integrator and a
   `workprecision.jl`; neither exists on `evanescent` (GPU_PLAN.md §1 notes this), and the
   brief asks only for `run.jl`. `CLAUDE.md` is corrected in `gpu/32-docs`.

## Open questions

- **Per-quantity tolerances.** Two classes now, not one number and not one per quantity.
  That is enough for `Eω`, which is what matters; the remaining candidate for a third class
  would be the `HCubature` error estimates (`transverse_integral_error_*`), which set the
  `:fixed` `stats` tolerance for the two modal cases at 1.8e-09 and 2.9e-07 while everything
  else in those cases is at rounding level. Worth doing if `gpu/22-modal` turns out to move
  them.
- **Response coverage.** Finding 8 is addressed for `Kerr_env_thg`, `Chi2Env`, ADK, the
  vector path and GNLSE Raman/shock. `Kerr_field_nothg` is still uncovered: it is what
  `prop_capillary` selects for a `RealGrid` with `thg=false`, and adding a case for it is
  one line if Group D wants it.
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
