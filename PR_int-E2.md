# `gpu/int-E2`: Group E2 integration

Base: `gpu/27-linop-integral` b5d69c58, itself cut from `gpu/int-E` 5ec21f77. This branch
merges the other two Group E2 branches into it, makes the three consistent, and records the
project's final regression position against `evanescent`.

Group E2 is the last implementation group in `GPU_PLAN.md` before Group F (documentation
consolidation and the two conditional branches). With it the plan's §6 sequence is complete.

## Commits

| commit | what |
| --- | --- |
| `1ffbc8fd` | merge of `gpu/26-rectmode-fix` |
| `b7f8fdf8` | merge of `gpu/25-modal-error-stat` |
| `00aee73e` | the semantic sweep (`Stats.default` argument order, the `test_device.jl` guard, `docs/src/gpu.md`) |
| `a7d8425d` | re-measured tolerances for `gradient_field_kerr` and `taper_field_kerr`, and the recorded decision about the free-space cases' statistics |

## Merged branches

| branch | HEAD | merge commit | what it brings |
| --- | --- | --- | --- |
| `gpu/26-rectmode-fix` | b40c5b0e | 1ffbc8fd | `x1 >= ul[2]` → `x2 >= ul[2]` in `NonlinearRHS._points!`, the `rect_modal_field` regression case and its tolerances |
| `gpu/25-modal-error-stat` | 0495f948 | b7f8fdf8 | the Kronrod transverse-integral-error statistic, `Stats.default` for the radial and free-space transforms, explicit transverse reductions |

Both were cut from `gpu/int-E` 5ec21f77 and merged with `--no-ff` in the order the brief
gives. `gpu/27-linop-integral` is the base rather than a merge, so its commits are already
in the history.

## Conflicts

**None.** Both merges completed with the `ort` strategy and no conflicted paths; the only
file either merge had to combine hunks in was `src/NonlinearRHS.jl` (26's `_points!` fix
against 25's two docstrings, far apart in the file). The places the plan expected trouble:

- `src/Interface.jl`: 27's `tabvalues`/`aefftol` block ends where 25's `Stats.default(...;
  gas=gas, stats_kwargs...)` call begins — adjacent, not overlapping. The merged
  `prop_capillary_args` keeps both, which I read back to confirm (`src/Interface.jl:605-621`).
- `test/regression/cases.jl`: 26 appends a case and its setup function at the end, 25 edits
  a comment block about the free-space statistics in the middle.
- `test/regression/tolerances.jl`: only 26 touches it.
- `test/test_device.jl`, `test/test_metal.jl`, `docs/src/gpu.md`: only 25 touches them.

Both reviewers had already run `git merge-tree --write-tree` for these pairs and reported a
clean tree; that is what happened.

## The semantic sweep

Commit `00aee73e`. Three things the merge left inconsistent rather than conflicting.

### One argument order for `Stats.default`

`gpu/25` added

```julia
Stats.default(grid, Eω, transform::Union{TransRadial, TransFree, TransFree2D}, linop; ...)
```

while the two modal methods are

```julia
Stats.default(grid, Eω, mode::Modes.AbstractMode,   linop, transform; ...)
Stats.default(grid, Eω, modes::Modes.ModeCollection, linop, transform; ...)
```

so the same two objects appeared in opposite orders depending on the geometry, and a caller
who wrote the modal order for a free-space run got a `MethodError`. The free-space method
now takes `(grid, Eω, linop, transform)`: the modal order with the mode left out, which is
also the order `Luna.run(Eω, grid, linop, transform, FT, output)` takes them in. Dispatch is
on the fourth positional argument, so there is no ambiguity with the modal methods.

Changed with it: the twelve free-space examples under
`examples/low_level_interface/freespace/`, `test/test_stats.jl`, `test/test_device.jl`,
`test/test_metal.jl`, and the two `!!! note "Statistics"` blocks `gpu/25` added to
`Luna.setup`'s docstrings. This is a source-incompatible change to a method one branch old
and never released; nothing outside the worktree can be calling it.

### The `test_device.jl` guard gap

`test/test_device.jl` guards its JLArray tests with

```julia
have_jlarrays = try; @eval import JLArrays; true; catch; false; end
if !have_jlarrays
    @warn "..."
else
    ...
end # have_jlarrays
```

Five JLArray testsets added by `gpu/21-free-device` sit *after* that block's `end`, and
`gpu/25` added three more next to them, so in an environment without JLArrays — which is
what `julia --project=<worktree>` gives, since `JLArrays` is a test-only dependency — the
file ended with eight errors instead of skipping them. Pre-existing, noted by the `gpu/25`
reviewer. The eight now sit in a second `if have_jlarrays ... end` block; the two host-only
`Float32` testsets and the two buffer-count testsets between the blocks still run.

### The documentation

`docs/src/gpu.md`'s "What runs where" listed "two pieces of per-step work [that] are still
host code inside otherwise device-capable runs": the crystal-optics normalisation and
`fwhm_r`. After 25 and 27 there are three, and the statistics one covers more than `fwhm_r`.
It now names the host-only statistics (`fwhm_r` and the rest of `Stats.beam_profile`, the
modal reconstruction error, the transverse quadrature error), says that one such member
makes its whole set host-only, and names `linop_integral=:quadrature`'s host evaluation of
the operator — so the section agrees with the statistics text and with "Tapers and pressure
gradients" instead of predating both.

The paragraphs about the absorbing boundaries and the default statistics were at the end of
"### Memory" and are not about memory; they now have their own "### Boundaries and
statistics" heading, which is also what "What runs where" now points at.

`docs/src/developer/device_model.md` needed nothing: 25's "Statistics"/"Radial and
free-space states" sections and 27's "The integrated linear operator" section were written
against each other's branches and already agree.

### Recorded, not changed

`gpu/25`'s review suggested adding the default statistics to the two radial gate cases.
Not done, and the reason is in `test/regression/cases.jl` next to `freestats`: `cases.jl`
has to run unchanged on every baseline commit, and none of the three in use (`b641025b`,
`782f55d1`, `fdf8dbe3`) has a `Stats.default` for a radial transform. Recording those
statistics would make the gate impossible to run against anything before `gpu/25`, which is
what the project's regression record is made of. It becomes available once the oldest
baseline in use is `gpu/int-E2` or later.

The detached baseline worktree `Luna-gpu/baselines/gpu26-5ec21f77`, which `gpu/26` left
behind, has been removed (`git worktree remove`). Its one baseline file, the pre-fix
`rect_modal_field`, is kept in the scratchpad and is what the int-E comparison below uses.

## The regression gate

M1 Pro, Julia 1.13.0, `--project=<worktree>`, `-t 1`, `Luna.set_fftw_mode(:estimate)`,
one FFTW thread, one BLAS thread, `Luna.set_fftw_wisdom(false)`.

Baselines were regenerated from the merge commit `00aee73e` over the whole 23-case matrix,
including `rect_modal_field`. The three older baseline directories
(`b641025b`, `782f55d1`, `fdf8dbe3`) predate that case, so they are run with
`LUNA_REGRESSION_SKIP=rect_modal_field`; a `rect_modal_field` baseline was generated
separately from each of the three (the `gpu/int-E` one is `gpu/26`'s, kept from its scratch
worktree) and gated with `LUNA_REGRESSION_ONLY`/`LUNA_REGRESSION_DIR`.

### 1. Against itself (`00aee73e`)

**506 pass, 0 fail, 0 error**, 23 cases in both modes: every case, both classes,
**exactly `0.000e+00`**. The run is reproducible in this environment.

### 2. Against `gpu/int-E` `b641025b`: what Group E2 changed

**443 pass, 23 fail** over the 22 older cases, plus the `rect_modal_field` comparison run
on its own. Every case other than the three below is exactly `0.000e+00` in both classes
and both modes, and no case changed its accepted step count except `rect_modal_field`.

| case | mode | Eω diff | Eω tol | stats diff | stats tol | worst quantity |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `gradient_field_kerr` | fixed | 7.968e-03 | 1.0e-12 | 6.733e-06 | 1.0e-12 | `Eω [save 11/11]` |
| `gradient_field_kerr` | adaptive | 7.275e-03 | 2.5e-11 | 4.899e-04 | 1.6e-04 | `Eω [save 11/11]` |
| `taper_field_kerr` | fixed | 9.228e-02 | 1.0e-12 | 5.085e-04 | 1.0e-12 | `Eω [save 11/11]` |
| `taper_field_kerr` | adaptive | 9.019e-02 | 1.0e-12 | 1.175e-03 | 3.5e-05 | `Eω [save 11/11]` |
| `rect_modal_field` | fixed | 2.095e-01 | 1.0e-12 | 9.154e-01 | 2.6e-12 | `stats/transverse_integral_error_rel` |
| `rect_modal_field` | adaptive | 2.095e-01 | 2.1e-12 | `Inf` | 7.3e-09 | step count (73 steps, baseline 74) |

These six rows are the two intended physics corrections and nothing else. The gradient and
taper numbers are `gpu/27`'s propagator change; the `rect_modal_field` numbers reproduce
`PR_26-rectmode-fix.md`'s pre-fix/post-fix table to every digit, which also shows that
`gpu/25` and `gpu/27` leave that case alone.

### 3. Against `gpu/int-A` `782f55d1`: Groups B to E2

**443 pass, 23 fail.** The six rows above, unchanged to every digit, plus the rows Groups B
to E recorded (`PR_int-E.md`, "Against `gpu/int-A`"): the four ionising cases and
`multimode_field_plasma` at the `peak_ionisation_rate`/`transverse_integral_error_rel`
level (1.1e-15 to 2.4e-04, all inside tolerance), and the two `free2d_*_chi2` cases at
7.7e-16 to 4.8e-14. The three radial cases are exactly `0.000e+00` here.

That the gradient, taper and `rect_modal_field` numbers are identical against `b641025b`
and `782f55d1` is the statement that Groups B, C, D and E left those three cases untouched.

### 4. Against `evanescent` `fdf8dbe3`: the whole project

**443 pass, 23 fail**, and this is the project's final regression position. The same rows
again, plus the radial rows `gpu/02-radialgrid` recorded when `Grid.RadialGrid` replaced
`Hankel.QDHT`:

| case | mode | Eω diff | Eω tol | stats diff | stats tol | worst quantity |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `radial_field_kerr` | fixed | 1.901e-15 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 1/1, save 10/11]` |
| `radial_field_kerr` | adaptive | 6.713e-11 | 1.2e-08 | 0.000e+00 | 1.0e-12 | `Eω [component 1/1, save 2/11]` |
| `radial_env_kerr` | fixed | 2.866e-15 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 1/1, save 11/11]` |
| `radial_env_kerr` | adaptive | 4.017e-11 | 8.2e-09 | 0.000e+00 | 1.0e-12 | `Eω [component 1/1, save 2/11]` |
| `radial_field_raman` | fixed | 2.134e-15 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 1/1, save 11/11]` |
| `radial_field_raman` | adaptive | 1.951e-11 | 1.1e-09 | 0.000e+00 | 1.0e-12 | `Eω [component 1/1, save 2/11]` |
| `modeavg_field_plasma` | fixed | 1.131e-15 | 1.0e-12 | 4.769e-15 | 1.0e-12 | `stats/peak_ionisation_rate` |
| `modeavg_field_plasma` | adaptive | 5.616e-09 | 9.2e-07 | 7.103e-06 | 1.2e-03 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | fixed | 1.245e-15 | 1.0e-12 | 5.202e-15 | 1.0e-12 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | adaptive | 2.283e-13 | 4.8e-10 | 8.550e-08 | 1.8e-04 | `stats/peak_ionisation_rate` |
| `modeavg_field_vector` | fixed | 1.133e-15 | 1.0e-12 | 6.029e-13 | 8.0e-11 | `stats/transverse_integral_error_rel` |
| `modeavg_field_vector` | adaptive | 2.271e-10 | 6.3e-09 | 3.809e-05 | 1.1e-03 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | fixed | 5.789e-14 | 6.7e-12 | 4.768e-13 | 5.8e-11 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | adaptive | 2.383e-07 | 7.3e-05 | 2.319e-04 | 7.3e-02 | `stats/transverse_integral_error_rel` |
| `free2d_field_chi2` | fixed | 7.683e-16 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 2/2, save 8/11]` |
| `free2d_field_chi2` | adaptive | 4.841e-14 | 7.4e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 2/2, save 2/11]` |
| `free2d_env_chi2` | fixed | 8.348e-16 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 1/2, save 11/11]` |
| `free2d_env_chi2` | adaptive | 3.151e-14 | 5.0e-12 | 0.000e+00 | 1.0e-12 | `Eω [component 2/2, save 2/11]` |
| `gradient_field_kerr` | fixed | 7.968e-03 | 1.0e-12 | 6.733e-06 | 1.0e-12 | `Eω [save 11/11]` |
| `gradient_field_kerr` | adaptive | 7.275e-03 | 2.5e-11 | 4.899e-04 | 1.6e-04 | `Eω [save 11/11]` |
| `taper_field_kerr` | fixed | 9.228e-02 | 1.0e-12 | 5.085e-04 | 1.0e-12 | `Eω [save 11/11]` |
| `taper_field_kerr` | adaptive | 9.019e-02 | 1.0e-12 | 1.175e-03 | 3.5e-05 | `Eω [save 11/11]` |
| `rect_modal_field` | fixed | 2.095e-01 | 1.0e-12 | 9.154e-01 | 2.6e-12 | `stats/transverse_integral_error_rel` |
| `rect_modal_field` | adaptive | 2.095e-01 | 2.1e-12 | `Inf` | 7.3e-09 | step count (73 steps, baseline 74) |

Every other case is exactly `0.000e+00` against `evanescent` in both classes and both modes.
Apart from the six rows of the two intended corrections, the largest difference anywhere in
the matrix is **2.383e-07** (`multimode_field_plasma`, `:adaptive` `Eω`, against a 7.3e-05
tolerance) and the largest in the `:fixed` mode, where the step sequence is imposed, is
**5.789e-14**. No case except `rect_modal_field` changed its accepted step count over the
whole project.

### Tolerances for the two moved cases

Commit `a7d8425d`. A fresh one-ulp sensitivity run on this branch
(`test/regression/sensitivity.jl gradient_field_kerr taper_field_kerr`):

| case | mode | `Eω` sensitivity | `stats` sensitivity | largest quantity |
| --- | --- | ---: | ---: | --- |
| `gradient_field_kerr` | fixed | 1.050e-15 | 1.246e-15 | `stats/peakintensity` |
| `gradient_field_kerr` | adaptive | 2.462e-13 | 1.570e-06 | `stats/zdw` |
| `taper_field_kerr` | fixed | 1.605e-15 | 1.612e-15 | `stats/peakintensity` |
| `taper_field_kerr` | adaptive | 7.566e-15 | 3.467e-07 | `stats/peakintensity` |

so the tolerances (100x, floor 1e-12) become

| case | mode | `Eω` before | `Eω` now | `stats` before | `stats` now |
| --- | --- | ---: | ---: | ---: | ---: |
| `gradient_field_kerr` | fixed | 1.0e-12 | 1.0e-12 | 1.0e-12 | 1.0e-12 |
| `gradient_field_kerr` | adaptive | 2.6e-07 | **2.5e-11** | 1.6e-04 | 1.6e-04 |
| `taper_field_kerr` | fixed | 1.0e-12 | 1.0e-12 | 1.0e-12 | 1.0e-12 |
| `taper_field_kerr` | adaptive | 8.7e-07 | **1.0e-12** | 3.1e-05 | **3.5e-05** |

The `:adaptive` `Eω` tolerances fall by four and five orders of magnitude, which is a
statement about the new propagator rather than about the cases: the one-point rule rebuilt
the operator at every stage, so where the step-size controller landed fed back into the
field, and a one-ulp perturbation of the input moved it by 1e-9. With Φ read off a table
built once at setup, that feedback is gone and the adaptive runs are as insensitive as the
fixed-step ones. The statistics tolerances barely move: they are still recorded once per
accepted step. These two cases are no longer the loosest in the matrix.

## Tests

M1 Pro, Julia 1.13.0, `-t 1` unless stated, `Luna.set_fftw_mode(:estimate)`, one FFTW
thread, one BLAS thread, `Luna.set_fftw_wisdom(false)`. Several other Julia processes were
on the machine throughout, so the wall times are upper bounds.

### CPU test files

`--project=<worktree>`, each file `include`d into one session.

| file | pass | fail | error |
| --- | ---: | ---: | ---: |
| `test_interface.jl` | 369 | 0 | 0 |
| `test_linops.jl` | 354 | 0 | 0 |
| `test_modes.jl` | 724 | 0 | 0 |
| `test_multimode.jl` | 15 | 0 | 0 |
| `test_rect_modes.jl` | 89 | 0 | 0 |
| `test_freespace.jl` | 77 | 0 | 0 |
| `test_boundaries.jl` | 182 | 0 | 0 |
| `test_gradient.jl` | 7 | 0 | 0 |
| `test_tapers.jl` | 2 | 0 | 0 |
| `test_output.jl` | 117 | 0 | 0 |
| `test_rk45.jl` | 51 | 0 | 0 |
| `test_stats.jl` | 167 | 0 | 0 |
| **total** | **2154** | **0** | **0** |

`test_rk45.jl` writes its assertions at the top level rather than inside a `@testset`, so
it prints no summary of its own; the 51 is from a run with the file wrapped in one.

### Device tests

| what | result |
| --- | --- |
| `test_device.jl`, JLArrays, `-t 1` | **1432 pass, 0 fail, 0 error**, 76 testsets |
| `test_device.jl`, JLArrays, `-t 4` | **1432 pass, 0 fail, 0 error**, 76 testsets — identical |
| `test_metal.jl`, Metal on this machine's GPU | **790 pass, 0 fail, 0 error**, 31 testsets |

The Metal run was the only Metal process on the machine while it ran.

`test_device.jl` under `--project=<worktree>`, i.e. with no JLArrays, now skips the device
tests with one warning instead of erroring; that is the guard fix above.

### Documentation

`docs/make.jl` in a stacked environment with Documenter: **four** unresolved cross
references, all pre-existing and none from this branch —
`loadFFTwisdom` and `saveFFTwisdom` in `modules/Luna.md`, `AbstractOutput` in
`modules/Output.md`, `Luna.PhysData.crystal_internal_angle` in `modules/NonlinearRHS.md`.
`gpu/int-E` reported six; the two `LinearOps.βz` references in `modules/Boundaries.md`
were fixed on `gpu/27`. The build still terminates on them before rendering, as it did on
`gpu/int-E`; nothing here made it worse.

### Parsing and examples

Every `.jl` file under `src/`, `test/` and `examples/` parses: **139 files, all OK.**

One example per geometry through `Luna.run` (the propagation part of each file, cut before
the first plotting call, which needs a display), plus the gradient, taper and rectangular
examples:

| example | geometry | result |
| --- | --- | --- |
| `low_level_interface/basic_modeAvg.jl` | mode-averaged, field | OK, 8.2 s |
| `low_level_interface/basic_modeAvg_env.jl` | mode-averaged, envelope | OK, 2.3 s |
| `low_level_interface/basic_modal.jl` | multimode, adaptive cubature | OK, 158.7 s |
| `low_level_interface/full_modal/basic_modal_full.jl` | multimode, `full=true` | OK, 436.5 s |
| `low_level_interface/freespace/radial.jl` | radial | OK, 33.4 s |
| `low_level_interface/freespace/free2D.jl` | 2-D Cartesian | OK, 5.0 s |
| `low_level_interface/freespace/full3D_env.jl` | 3-D Cartesian | OK, 89.1 s |
| `low_level_interface/freespace/free2D_bbo.jl` | 2-D Cartesian, χ⁽²⁾ | OK, 16.4 s |
| `low_level_interface/gnlse/simplescg_modeAvg_env.jl` | GNLSE | OK, 23.3 s |
| `low_level_interface/gradients/gradient_modeAvg.jl` | pressure gradient | OK, 8.0 s |
| `low_level_interface/gradients/gradient_modal.jl` | pressure gradient, multimode | OK, 29.7 s |
| `low_level_interface/tapers/taper_modeAvg.jl` | taper | OK, 10.2 s |
| `low_level_interface/rectangular/rectangular_modal.jl` | rectangular, `a > b` | OK, 167.3 s |
| `simple_interface/gnlse_sol.jl` | `prop_gnlse` | OK, 2.6 s |

The five free-space examples exercise the new `Stats.default(grid, Eω, linop, transform)`
argument order; the two gradient examples and the taper example exercise
`linop_integral=:auto`; `rectangular_modal.jl` is the 50 x 10 µm guide, so it is the example
whose answer the `gpu/26` fix changes.

The example runner does not call `Luna.set_fftw_wisdom(false)` — it is the script `gpu/int-E`
used, unchanged — so these runs read and wrote the shared FFTW wisdom cache. It does not
touch the regression record: the gate, the baseline generation and the sensitivity run all
disable wisdom themselves.


## State of the project

Group E2 is the last implementation group. What the branch leaves behind:

### What runs where

Device-capable (Metal tested on this machine, JLArrays in `test_device.jl`, CUDA untested
for want of a host): the mode-averaged transform, `TransRadial`, `TransFree2D`,
`TransFree`, and `TransModalFixed` (`modal_integral=:fixed`); every nonlinear response
Luna ships (Kerr with and without THG, the two χ⁽²⁾ responses, plasma, Raman); the
absorbing boundaries and the transverse collars; the stepper; the default statistics,
behind the measured `Stats.STATS_DEVICE_MINLEN` threshold; and a z-dependent linear
operator, through the tabulated integral.

Host only: the adaptive multimode transform `TransModal` (the CPU default), `prop_gnlse`,
the crystal-optics normalisation's refractive indices, `Stats.beam_profile` and the two
transverse-integral error statistics, and `linop_integral=:quadrature`'s evaluation of the
operator.

### The two intended physics corrections

Both are exceptions to `GPU_PLAN.md` §9 decision 1 ("identical to machine precision"), both
were taken with the user (decisions 6 and 7), and both are recorded in the regression table
above rather than hidden inside a tolerance.

1. **`gpu/27-linop-integral`**: the one-point propagator `exp(L(t2)·(t2 − t1))` for a
   z-dependent linear operator is replaced by `exp(Φ(t2) − Φ(t1))`. Every taper and every
   pressure gradient moves by the amount the one-point rule was wrong — first order in the
   step, and invisible to the step-size controller because it cancels out of the embedded
   error estimate. Uniform fibres are bit-for-bit unchanged.
2. **`gpu/26-rectmode-fix`**: the Cartesian in-domain test of the adaptive modal integral
   tested `x1` against the upper limit of the second coordinate. Every `RectMode`
   propagation with `a > b` and `full=true` moves. Square guides and polar domains are
   bit-for-bit unchanged.

`examples/low_level_interface/rectangular/rectangular_modal.jl` is a 50 × 10 µm guide, so
it is one of the runs that now gives a different, correct answer.

### Open for Group F

- **`gpu/32-docs`**: `docs/src/gpu.md` and `docs/src/developer/device_model.md` still carry
  the "Work in progress" admonition and per-branch measurement tables. They are internally
  consistent as of this branch but want consolidating, together with `CLAUDE.md` and the
  README.
- **`gpu/30-neff-grid`** (conditional): §4.5 layer 3, the tabulated `neff` grid, only if a
  benchmark justifies it. `gpu/23`/`gpu/27`'s tabulation of Φ, β and `Aeff` already removes
  the per-stage host evaluation from a tapered or graded device run, which was the reason
  layer 3 was proposed; the decision is now whether anything is left to gain.
- **`gpu/31-ka-kernels`** (optional): KernelAbstractions pointwise kernels on the CPU, only
  if profiling shows the broadcasts dominate.
- CUDA has never been run. `LUNA_TEST_CUDA=1` on a CUDA host is the remaining item in
  §10.3 before any of this is opened as an upstream PR.

### Known gaps carried forward

Unchanged from `PR_int-E.md` except where a Group E2 branch closed one:

- Closed by `gpu/25`: `Stats.default` no longer raises a `MethodError` for a
  `TransModalFixed`; the default statistics now work on a radial or free-space state (the
  `Stats.squeeze` rank limitation is gone).
- Closed by `gpu/26`: the `x1 >= ul[2]` typo.
- Closed by `gpu/27`: the one-point propagator; `tabulate_linop` is no longer a change of
  discretisation a user has to opt into, because the integral is now the definition.
- Still open: `Stats.electrondensity`'s mode-averaged host branch scales the shared `Et` in
  place, so a user statistic after it in a set sees a scaled field (pre-existing);
  `Stats.peakintensity` of more than one mode, `fwhm_r` and the two transverse-integral
  error statistics have no device form; the crystal-optics normalisation stages through the
  host; a tapered mode collection with `modal_integral=:fixed` re-evaluates the mode fields
  on the host whenever `z` moves; a free-space energy functional has no device form.
