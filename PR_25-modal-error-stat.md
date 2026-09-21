# The Kronrod error statistic, and `Stats.default` for radial and free-space states

Branch `gpu/25-modal-error-stat`, base `gpu/int-E` (`5ec21f77`). GPU_PLAN.md §4.8, §6
Group E2, §9 decision 8, §11.

Two gaps `gpu/int-E` recorded, both in `Stats.jl`:

- `Stats.default` raised a `MethodError` on a `NonlinearRHS.TransModalFixed` unless
  `mode_error=false` was passed, so a `modal_integral=:fixed` run recorded no
  transverse-integral diagnostic at all. `NonlinearRHS.integral_error!` existed and was
  tested; nothing called it per step.
- `Stats.default` had no method for a radial or free-space transform, and the individual
  statistics could not have run on such a state anyway: `Stats.squeeze` had a method for a
  vector and one for a matrix. That is why the free-space examples recorded no statistics
  at all and carried a commented-out `Stats.collect_stats` line, and why the regression
  gate's free-space cases record only `z` and `dz`.

**The default CPU path is unchanged: the regression gate is exactly `0.000e+00` on all 22
cases, in both modes and both classes, against `b641025b` (the `gpu/int-E` merge commit).**

## What changed

### `src/Stats.jl`

**The Kronrod statistic.** `Stats.transverse_integral_error(t::TransModalFixed)` records,
per accepted step, the fixed quadrature rule's own embedded error estimate:

| dataset | what it is |
| --- | --- |
| `transverse_points` | the rule's node count, fixed for the propagation |
| `transverse_integral_error_abs` | the root-mean-square over `(ω, mode)` of `NonlinearRHS.integral_error!`, i.e. the embedded coarse rule minus the rule the propagation uses |
| `transverse_integral_error_rel` | the same over the root-mean-square of the polarisation itself |

Those are the three datasets `Stats.mode_reconstruction_error` records besides
`mode_reconstruction_error` itself, which the fixed rule has no way to compute: it needs
the adaptive transform's single-point machinery. Both error datasets are `NaN` where the
rule has no embedded coarse rule (`NonlinearRHS.has_error_estimate`, i.e. without
`modal_kronrod=true` on a radial or Cartesian rule); the node count is recorded either way.
The transform is evaluated once per call on the accepted state, as
`mode_reconstruction_error` does, so the estimate belongs to the state which was saved
rather than to whichever stage the stepper last called.

`Stats.default`'s `mode_error=true` now picks the diagnostic the transform has
(`_mode_error_stat`, one method per transform type), so it no longer errors on a
`TransModalFixed`.

The statistic is host-only, like `fwhm_r` and `peakintensity` of more than one mode, which
already make every multimode set host-only. The *transform* it evaluates may be on a
device, so it stages the host state into the transform's units and array type and scales
the absolute error back by `E_ref`; that path is tested on JLArray.

**The transverse reductions.** `Stats.squeeze`, with its `Array{T,1}` and `Array{T,2}`
methods, is replaced by two functions which do not care about the rank:

- `_dropfreq(x)` drops the frequency axis after a reduction along it — a scalar for a
  mode-averaged state, one value per column otherwise;
- `_bycolumn(f, Eω)` applies an energy functional to each slice of the second axis, which
  is the mode axis of a modal state and the polarisation axis of a free-space one, with
  everything past it left to the functional.

For a state of rank 1 or 2 both are the expressions they replace, element for element,
which is why the gate is exactly zero.

**Radial and free-space statistics.** A new
`Stats.default(grid, Eω, transform::Union{TransRadial, TransFree, TransFree2D}, linop;
windows, gas, userfuns, collar, beam_profile, stats_device)`. `linop` is accepted for
symmetry with the modal methods and is not used: a free-space operator has no mode and so
no zero-dispersion wavelength to record.

A radial or free-space state is `(nω, npol, nk...)` and lives in transverse *reciprocal*
space, which decides the shape of everything below.

- **`energy`** — the total energy, from `Fields.energyfuncs(grid, spacegrid)[2]`, one
  number per polarisation component. `Stats.energy`, `energy_window` and `energy_λ` take
  an optional transverse grid; the device form of all three is now a `Stats.EnergyWeights`
  — a prefactor and one weight vector per *weighted* axis — reduced over every axis but
  the polarisation one. A radial grid weights the frequency axis (`SimpsonEven` or a plain
  sum, as before) and the reciprocal radial axis (`Grid.integrate_k`'s weights); the
  Cartesian grids weight nothing, because their functional is a plain sum over every axis
  including the frequency one. The weights are still *checked* against the functional the
  caller passed, now on a probe state of the run's own shape, so a functional which is not
  the grid's own gives a host-only statistic rather than a silently different quantity.
- **`peakintensity`, `ω0`, `fwhm_t_min`, `fwhm_t_max` on the propagation axis.**
  `Stats.onaxis(f, spacegrid)` wraps any statistic so that it is evaluated there, because
  a transverse average of a duration is not a duration. The field on axis is a fixed
  linear combination of the transverse samples — `Grid.onaxis`'s weights on a radial grid,
  and on a Cartesian grid the inverse DFT evaluated at `x = 0`, whose weights are exactly
  `±1/N` on a centred axis (`cispi` of an integer is exactly `±1`) — so the projection is
  one weighted reduction over the transverse axes, and the inverse time transform which
  follows is of the `(nω, npol)` result rather than of the state. `ω0`, `fwhm_t`,
  `peakintensity` and `electrondensity` then apply unchanged, device branches included.
  Two new constructors take a field which is already in V/m, which is what a projected
  free-space state is: `peakintensity(grid)` (`c ε₀/2 max_t Σ_pol |E|²`) and
  `electrondensity(grid, ionrate, dfun)` (two components drive the rate through their
  quadrature sum, as in the multimode method).
- **`fwhm_r` (or `fwhm_x`/`fwhm_y`) and `collar_energy_fraction`** —
  `Stats.beam_profile(grid, spacegrid; collar)`, from the transverse fluence profile
  `Σ_{ω,pol} |E(ω, pol, r)|²` in **real** space. On a radial grid the profile is mirrored
  about the axis (`Grid.rsymmetric`), with the on-axis sample taken from the state; on a
  Cartesian grid the width is measured on the cuts through the axis. The collar fraction is
  the part of the transverse integral of the profile which lies where the transverse
  absorber is active (`Boundaries.rprofile`); `collar` defaults to
  `Boundaries.DEFAULT_RCOLLAR`, which is `Luna.run`'s default `rcollar`, and `nothing`
  skips the dataset.
- `density`, `pressure`, `electrondensity`/`peak_ionisation_rate` on axis with plasma,
  `z` and `dz`, as for the modal sets.

`Stats.collect_stats` now plans the analytic transform only when a statistic reads it
(`needs_time`). The free-space set does not — its on-axis statistics transform the
projected field instead — and on a `RealGrid` that buffer is four times the state, which
on a 3-D grid is the whole memory budget.

### `src/Interface.jl`

`_statskwargs`, which forced `mode_error=false` for a `TransModalFixed`, is removed;
`prop_capillary_args` passes `stats_kwargs` straight through. A `modal_integral=:fixed`
run now records `transverse_points` and the two error datasets.

### `src/Plotting.jl`

`Plotting.stats` plots `fwhm_x`, `fwhm_y` and `collar_energy_fraction` when they are
present, so a free-space output plots as a modal one does.

### `examples/low_level_interface/freespace/`

All twelve examples build the statistics with `Stats.default(grid, Eω, transform, linop)`
and pass them to `Output.MemoryOutput`, in place of the commented-out
`Stats.collect_stats(grid, Eω, Stats.ω0(grid))` line which would have thrown. Each was run
with its grids and lengths shrunk so a run takes seconds; all twelve complete and record
11 to 13 statistics -- 13 for the 3-D χ⁽²⁾ case, which has two polarisation components and
two transverse widths.

### `docs/`

`docs/src/developer/device_model.md` gains a "Radial and free-space states" subsection
under "Statistics" (the energy weights, the on-axis projection, why the beam profile has
no device form) and its tables are corrected. `docs/src/gpu.md` says that the free-space
sets are host-only because of the beam profile and what `beam_profile=false` buys, and its
`modal_integral=:fixed` paragraph now describes the statistic instead of promising it.
`src/Luna.jl`'s and `src/NonlinearRHS.jl`'s notes about `mode_error=false` are replaced.

## Tests

### `test/test_stats.jl` (1 testset before, 10 after)

- **the transverse integral error statistic** (23): the statistic against
  `NonlinearRHS.integral_error!` of the same state, with `kronrod=true` (finite, and
  `< 1e-3` relative) and without (`NaN`, with the node count still recorded); that
  `Stats.default` puts it in the set for a `TransModalFixed` and
  `mode_reconstruction_error` for a `TransModal`.
- **`<geometry>` statistics, `<grid type>`** (6 testsets, 22–23 each): every statistic of
  the radial, 2-D and 3-D free-space sets against a hand computation on the same state —
  the energy functional applied directly, the on-axis field from `FFTW.ifft`/`Grid.onaxis`
  and its analytic transform, the fluence profile and its FWHM, the collar mask. Also that
  the set is host-only exactly because of `BeamProfile`, and that dropping it changes
  nothing about the rest.
- **on-axis electron density** (3): against `Ionisation.IonRateADK` and
  `Maths.cumtrapz!` applied to the on-axis field by hand.
- **a radial propagation with the default statistics** (6): end to end, energy conserved
  to 1e-3 by a Kerr-only run, the beam diffracting, nothing `NaN`.

### `test/test_device.jl` (JLArray)

- **radial and free-space statistics on JLArray** (190): the sets on JLArray against the
  host, key by key and **shape by shape**, for eight cases — radial field-resolved,
  envelope, two polarisation components, plasma and an energy window; 2-D free space,
  field-resolved and envelope; 3-D free space. Plus that `BeamProfile` is what makes a set
  host-only.
- **the free-space device statistics undo the unit scaling** (36): the same set built for
  `E_ref = 1024` and given the state divided by 1024 gives the physical answers back, for
  all three geometries.
- **the transverse integral error on a device transform** (9): a host-built statistic with
  a `TransModalFixed` on JLArray behind it, against the same thing on the host.

The shape assertion earned its place: a reduction over several axes leaves them as
singletons, so the device energy came back as a `1×npol×1` array where the host branch
gives a vector. `Stats._weightedenergy` now `vec`s it.

### `test/test_metal.jl` (hardware)

- **radial and free-space statistics on Metal** (new): a radial field-resolved, a radial
  envelope and a 3-D free-space propagation with `Stats.default` and fixed steps, on Metal
  against **the same propagation on the CPU in `Float32`** — the same precision, so what is
  compared is the device arithmetic and the device branches against the host branches, not
  `Float32` against `Float64`. Every recorded statistic, key by key, at 1e-3. Then the
  whole set as a user gets it (`beam_profile=true`), which is host-only, against the same
  host reference.

  Two things about that comparison, both measured on the CPU (`Float32` against `Float64`
  on the same code path, 2 mm propagation, fixed steps):

  - the whole set differs between the two precisions by 1e-7 to 1e-6 relative, so the 1e-3
    the Metal comparison uses is comparing the device arithmetic and not the precision;
  - `collar_energy_fraction` is the exception and is compared **absolutely**. It is the
    ratio of the fluence in the collar to the total, and with the beam this far inside the
    aperture the numerator is at the rounding level: the same run records 1.6e-18 of its
    energy in the collar in `Float64` and 1.1e-15 to 1.3e-14 in `Float32`. Both mean
    "nothing has reached the collar"; neither is a number to compare relatively. The
    `beam_profile` docstring says so. Where the collar fraction is real — the 3-D case,
    which loses 2.8 % of the beam to the transverse absorber — the two precisions agree to
    2 % (4.60e-12 against 4.53e-12).

## The regression record

M1 Pro, Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread, wisdom off,
through the timing wrapper.

| baseline | result | largest difference |
| --- | --- | --- |
| `b641025b` (`gpu/int-E`, the merge commit) | **466 pass, 0 fail** | **`0.000e+00` on every row**, both classes, both modes |
| `fdf8dbe3` (`evanescent`) | **466 pass, 0 fail** | 2.319e-04 (`multimode_field_plasma` adaptive statistics) |

Against `evanescent` every non-zero row is one `gpu/int-E` already recorded, to every
digit: the plasma rows from `gpu/13`'s scan against the old serial `cumtrapz!`, the χ⁽²⁾
and radial rows from `gpu/15` and `gpu/02`, and the two modal rows from `gpu/22`'s
accumulation order. No case changed its number of accepted steps against either baseline.

The gate's own cases are unchanged. The free-space cases still use `RegressionCases.freestats`
(`z` and `dz` only) rather than the new `Stats.default`: adding the new statistics to them
would need the per-case tolerances re-measured by the one-ulp sensitivity study, which the
brief makes conditional on exactly that and which is a baseline-regeneration job rather than
a statistics one. The comment in `cases.jl` which said there was no such method is updated
to say why they are still not used there.

## Tests run

All on M1 Pro / Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread,
wisdom off, through the timing wrapper. Environments as each test file's header describes:
a stacked environment with `JLArrays` for `test_device.jl`, one with `Metal` for
`test_metal.jl`.

| what | result |
| --- | --- |
| `test/test_regression.jl` vs `b641025b` | **466 pass, 0 fail**, `0.000e+00` everywhere |
| `test/test_regression.jl` vs `fdf8dbe3` | **466 pass, 0 fail**, the `gpu/int-E` record to every digit |
| `test/test_device.jl` (JLArrays) | **1383 pass, 0 fail, 71 testsets** (1148/68 on the base) |
| `test/test_metal.jl` (Metal, hardware) | METALRES |
| `test/test_stats.jl` | **167 pass, 0 fail, 10 testsets** (1/1 on the base) |
| `test/test_multimode.jl` | 5 pass, 0 fail, 4 testsets |
| `test/test_freespace.jl` | 77 pass, 0 fail, 43 testsets |
| `test/test_interface.jl` | **292 pass, 0 fail, 12 testsets** (286/12 on the base) |
| `test/test_output.jl` | 117 pass, 0 fail, 8 testsets |
| `test/test_modes.jl` | 724 pass, 0 fail, 8 testsets |
| `docs/make.jl` | the six unresolved cross-references `gpu/int-E` recorded, and no new ones |
| the twelve free-space examples | **12 of 12 run**, each recording 11-13 statistics |

## Known gaps and open questions

- **A `modal_integral=:fixed` run now pays for a statistic it did not before.**
  `Stats.transverse_integral_error` evaluates the transform once per accepted step, which
  is one of the six right-hand-side evaluations a step makes, plus one more matrix product
  for the embedded rule. That is exactly what `Stats.mode_reconstruction_error` costs an
  adaptive run, and it is the price of the diagnostic; `stats_kwargs=Dict(:mode_error =>
  false)` turns it off, and `stats_period` amortises it.
- **The free-space default set is host-only**, because `Stats.beam_profile` has no device
  form: it applies the inverse transverse transform to the whole state, which is an `N×N`
  matrix product on a radial grid (and would need a complex copy of the Hankel matrix on
  the device) and an FFT over the transverse axes on a Cartesian one. `beam_profile=false`
  leaves a set every member of which is device-capable, which is how the device branches
  are reached in the tests. Given `gpu/int-E`'s re-measurement — the device path does not
  pay at any size Luna produces, and is two orders of magnitude worse on a radial state —
  giving the beam profile a device form would not change which path `:auto` takes today.
- **The on-axis projection is repeated per wrapper.** A set with three on-axis statistics
  reduces the state three times. Sharing it would mean caching the result against `z`,
  which a repeated `z` would silently get wrong, and the reductions are a small part of the
  step which produced the state. A two-phase statistics protocol (the one `gpu/24` already
  wants for the batched transfers) would give the projection a natural home.
- **`Stats.beam_profile` needs one more host buffer the size of the state**, for the field
  in transverse real space. On a large 3-D grid that matters; `beam_profile=false` or
  `Output.nostats` is the lever, as it already is for the per-step cost.
- **`fwhm_x`/`fwhm_y` rather than `fwhm_r` on a Cartesian transverse grid.** The brief says
  `fwhm_r`; a Cartesian grid has no radius, and a 3-D beam has two widths, so the keys are
  per axis and `Plotting.stats` plots them.
- **The collar width is a keyword, not read from the boundaries.** `Stats.default` is built
  before `Luna.run` decides anything about the absorbers, so `collar` defaults to
  `Boundaries.DEFAULT_RCOLLAR` — the same default `Luna.run` uses — and has to be set by
  hand if `rcollar` is. On a Cartesian grid it is ignored, because `Boundaries.rprofile`
  ignores it there too and uses the grid's own window.
- **The energy statistic's probe check is a runtime check, not a type-level one**, as on
  `gpu/24`. It now runs on a probe state of the run's own shape, so it also catches a
  functional built for a different transverse grid.
- **`Stats.peakpower` has no free-space form.** The instantaneous power of a free-space
  field is the transverse integral of the intensity, which needs the time-domain field of
  the whole state — the transform the free-space set exists to avoid. It is not in the set;
  `energy` and the on-axis peak intensity are.
- **`electrondensity` is field-resolved only**, as it already was: the constructors are
  typed on `Grid.RealGrid`.

## Deviations from GPU_PLAN.md and the brief

1. **`beam_profile` and the collar fraction are one statistic**, not two. Both are read off
   the same real-space transverse profile, which is the one expensive thing in the set, and
   the statistics protocol has no way for two statistics to share a per-call intermediate.
2. **`Stats.default`'s free-space method takes `collar` and `beam_profile` keywords** which
   the brief does not name. The first is how the collar fraction is reachable at all (see
   above); the second is what makes the device branches reachable.
3. **The gate's free-space cases are not changed** to use the new `Stats.default`. The
   brief makes that conditional on the tolerances staying measured; they would have to be
   re-measured by the sensitivity study, which is a baseline job.
4. **`Stats.collect_stats` no longer plans the analytic transform when nothing reads it.**
   `gpu/24`'s review recorded that it was planned and allocated either way "which is what
   Luna has always done"; that is now false, deliberately, because the free-space set is the
   first one which reads no time-domain field and the buffer is four times the state.
5. **`src/Plotting.jl` is a scope extension** — three lines, so that the new keys plot.
