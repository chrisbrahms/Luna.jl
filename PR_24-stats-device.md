# The default statistics on the device

Branch `gpu/24-stats-device`, base `gpu/int-D` (`90826dc4`). GPU_PLAN.md §4.2, §4.8
(second paragraph), §6 Group E, §11.

`gpu/11-boundaries-output` made a device run copy the whole state to the host on every
accepted step to compute the default statistics, and warn once that it was doing so. This
branch removes that copy: the statistics Luna ships are evaluated on the state where the
stepper holds it. A statistics function a user writes still gets the host copy, and the
warning now names it.

**The default CPU path is unchanged: the regression gate is exactly `0.000e+00` for all
21 cases, in both modes and both classes, against `fa556e6f` (source-identical to the
base) and against `fdf8dbe3` (`evanescent`).**

## Motivation

`Stats.jl` was closures over host `Vector`s and host reductions. On a device that meant
`Luna.ScaledOutput` copied the state down and unscaled it on every step whose statistics
fired — for the mode-averaged Kerr propagation `gpu/11` measured, the statistics were the
single largest per-step cost after the right-hand side, and `stats_period` existed mainly
to amortise the copy. The statistics themselves are reductions: they belong on the device.

## What changed

### `src/Stats.jl`

The statistics in the two `Stats.default` sets are callable structs instead of closures,
so that they can carry three traits:

| trait | what it says | fallback |
| --- | --- | --- |
| `Stats.device_capable(f)` | `f` can be evaluated on the state as the stepper holds it | `false` |
| `Stats.needs_time(f)` | `f` reads the time-domain field `Et` | `true` |
| `Stats.statlabel(f)` | what to call `f` in a diagnostic | the type's name |

The fallbacks are what a user-written closure gets, so a user statistic behaves exactly as
before: it is handed a host copy in physical units, and it makes `ScaledOutput` produce
one.

Each struct has two branches:

- the **host branch** is the code Luna has always run, transcribed unchanged;
- the **device branch** is reductions and broadcasts over the state.

Which one runs is decided by `Utils.isdevice`, not by precision. A `Float32` run on the
host therefore also takes the host branch, and nothing about the CPU output moves — which
is why the gate is exactly zero rather than merely within tolerance.

The device branches:

| statistic | device form |
| --- | --- |
| `ω0` | the two sums `Maths.moment` takes as two fused reductions, instead of materialising `abs2.(Eω)`. Scale-invariant |
| `energy`, `energy_λ`, `energy_window` | the frequency integral as one weighted reduction; a window folds into the weights squared |
| `peakpower` | one `maximum` (plus one over the polarisation sum for several modes) |
| `peakintensity` (effective area) | one `maximum` |
| `peakintensity` (on axis) | one `maximum`; for a single mode the transverse field is a scalar factor which comes out of it |
| `fwhm_t` | `abs2` into a device buffer and one copy of *that* to the host; the root-finding is then the same on the same samples |
| `electrondensity` (mode-averaged) | one broadcast of `Ionisation.ratekernel` off the state, then two reductions |
| `density`, `pressure`, `core_radius`, `zdw`, `zdz!` | unchanged: they read only `z` |
| `fwhm_r`, `mode_reconstruction_error`, `peakintensity` of more than one mode | **host-only**, `device_capable == false` |

Three of those need a word.

**The energy functional.** `Fields.energyfuncs(grid)[2]` integrates the spectral power
density with `NumericalIntegration`'s `SimpsonEven` on a `RealGrid` and with a plain `sum`
on an `EnvGrid`. `SimpsonEven` is the alternative extended Simpson rule: every sample with
weight 1 except the first four and last four, which carry 17/48, 59/48, 43/48 and 49/48.
As a weight vector that is one fused reduction on any array type. But `Stats.energy` takes
the functional as an *argument*, so the weights could be the wrong quantity for the
functional the caller passed (`Fields.energyfuncs(grid, spacegrid)[2]`, or their own). They
are therefore checked against it at construction, on a deterministic probe spectrum; if
they disagree the statistic is host-only. Nothing is silently replaced.

**The electron density.** Only `frac[end]` of the cumulative ionisation integral is ever
read, and that end point is the trapezoid rule over the whole window. The device branch is
therefore one weighted reduction, not a prefix scan followed by a scalar read of the last
element — which a device array refuses anyway. The intensity conversion
(`1/sqrt(ε₀c·aeff/2)`) and the unit scaling are folded into the field reference
`Ionisation.ratekernel` multiplies each sample by, so the field itself is never rescaled.
The rate is converted with `Ionisation.device_rate` at `prepare` time; a rate with no
device kernel (the direct `IonRatePPT`, a user's own) leaves the statistic host-only, as
does `oversampling != 1`.

**The unit scaling.** On a device the state is scaled (`e = E/E_ref`), because nothing has
unscaled it — `ScaledOutput` unscales only into a host buffer, which is the copy this
branch exists to avoid. The device branches apply `E_ref` themselves, always to the
*scalar result* of a reduction and in `Float64`, never to the field: in a `Float32` run
the field is deliberately of order one, and putting SI units back on it before reducing
would throw away the reason the scaling exists (the spectral field of a 300 µJ pulse is
around 1e13 V/m·s; its square summed over the grid is within a factor of 1e7 of `Float32`'s
range). `Stats.StatsContext` carries `E_ref`; `Stats.collect_stats` takes it as a keyword
(default `1.0`) and `Stats.default` reads it from the transform through
`Luna.runscaling`. `Stats.default`'s signature is unchanged.

Two more mechanical changes in the same file:

- `Stats.prepare(f, ctx)` is the second construction phase: `collect_stats` builds a
  `StatsContext` from the state prototype and the scaling and calls it on every statistic,
  which is where device buffers, mirrors and the converted ionisation rate are allocated.
  The fallback returns `f` unchanged, so a closure passes straight through.
- `Stats.plan_analytic` builds its buffers with `similar` and plans through
  `Utils.plan_ft`/`Utils.plan_ift`, so the analytic signal is one inverse FFT on whatever
  array type the state lives on, shared by every statistic which needs it. `inv` of a
  forward complex plan is the same `ScaledPlan` FFTW's `plan_ifft` returns, and the
  rFFT-to-FFT copy (a scalar loop, `copyto_fft!`) is now one broadcast writing the same
  `2*Eω[i]`, so the host arithmetic does not move. `collect_stats` skips the transform
  entirely when no statistic reads `Et`.

### `src/Device.jl`

- `Luna.stats_device_capable(x)` and `Luna.stats_host_list(x)`: the hooks `ScaledOutput`
  reaches the traits through. They answer for an output handler (`MemoryOutput`,
  `HDF5Output`), unwrap `Output.PeriodicStats`, are trivially `true`/empty for
  `Output.nostats`, and fall back to `false` and the type's name. `Stats` adds the methods
  for its own collector. `Device.jl` comes before `Stats.jl` in the include list, so the
  generic functions live here and the specific methods there.
- `ScaledOutput` gains a `devstats` field, decided once at construction as
  `isdevice(y) && stats_device_capable(o)`: the statistics function is fixed for the
  propagation and so is where the state lives. When it is true nothing copies the state to
  the host on any step. `needs_host_y` still gates the copy per step for a set which does
  need it, so `stats_period` keeps skipping it on a step `PeriodicStats` will not fire on.
- `devstats` is only ever true for a state **on a device**. A scaled host run (`Float32` on
  the CPU) still copies and unscales, because the statistics would otherwise take their
  host branches on a scaled field.
- The one-time warning names the statistics which have no device form, and no longer points
  forward to this branch.

`Luna.run` is unchanged: it already constructs `ScaledOutput`, and the decision is inside
the constructor. (`gpu/23-tabulated-linop` edits `Luna.run` at the operator-building lines;
this branch touches none of them.)

### `src/Interface.jl`

`prop_capillary_args` passes the state itself to `Stats.default` again. `gpu/11` had to
pass a host-shaped *template* because the `EnvGrid` `plan_analytic` planned an FFTW
transform directly on a copy of it, which fails on a device array; `plan_analytic` now
plans for the array type it is given. The `stats_period` docstring says what raising it
buys now.

### `docs/`

`docs/src/developer/device_model.md` gains a "Statistics" section (the traits, the two
branches, the device forms, the scaling and the `ScaledOutput` decision).
`docs/src/gpu.md`'s "What runs where" and `stats_period` sections say that the default
statistics run on the device and that a statistics function you write does not.
`docs/src/modules/Luna.md` lists the two new hooks.

## Tests

### `test/test_device.jl` (JLArray)

Four new testsets, plus the existing statistics assertions tightened:

- **boundaries and default statistics on JLArray** (was 5 assertions, now 32): every
  recorded statistic of the propagation, not only the energy, compared key by key at
  1e-10.
- **every default statistic on JLArray** (74): the default set evaluated directly on one
  state, host against JLArray, for five cases — field-resolved Kerr, envelope Kerr, the
  on-axis peak intensity, an energy window, and argon at 150 µJ with a tabulated PPT rate
  so that `electrondensity` is a real comparison (the case ionises ≳1e-4 of the gas, which
  is asserted). Comparing one state rather than two propagations means a failure names the
  quantity, not the run.
- **the device statistics undo the unit scaling** (14): on JLArray in `Float64` `E_ref` is
  1, so it is set by hand — the same statistics built for `E_ref = 1024` and given the
  state divided by 1024 must return the physical answers. That is exact arithmetic on the
  scaling, independent of precision, which a `Float32` comparison cannot separate out.
- **device_capable and the host statistics list** (10): the traits, their propagation
  through `MemoryOutput` and `PeriodicStats`, a user closure making the set host-only and
  appearing in the list by name, and the multimode default set reporting exactly
  `FWHMr`, `ModeReconstructionError` and `PeakIntensityModes`.
- **the default statistics skip the device-to-host copy** (10): the exit criterion.
  `ScaledOutput.ybuf` is the only place a device-to-host copy of the state lands and
  `_tohost_unscale!` is the only thing which writes it, so a sentinel written into it and
  found unchanged after three calls is the copy count — the same instrument gpu/11's
  `stats_period` test uses. Adding a user statistic brings the copy and the warning back.

### `test/test_metal.jl` (hardware)

- **default statistics on Metal** (new): `prop_capillary` with the default statistics on
  Metal against `prop_capillary` on the CPU in `Float32` — the same precision, so what is
  compared is the device arithmetic and the device branches against the host branches, not
  `Float32` against `Float64`. Every recorded statistic, key by key, at 1e-3. Then the copy
  count, as in the JLArray test, on a real `MtlArray` state with a real `E_ref`; and the
  user-statistic fallback with its warning.
- `metalcase`/`metalgradientcase` pass the state itself to `Stats.default` instead of a
  host template, so every existing Metal statistics assertion now exercises the device
  branches.

### CPU

`test/test_stats.jl`, `test/test_output.jl`, `test/test_polarisation_env.jl`,
`test/test_mixtures.jl`, `test/test_vectorplasma.jl`, `test/test_modes.jl` and
`test/test_interface.jl` are unchanged and pass.

## The regression record

M1 Pro, Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread, wisdom off,
`shotnoise=false`, through the timing wrapper.

| baseline | result | largest difference |
| --- | --- | --- |
| `fa556e6f` (source-identical to the base `90826dc4`) | **460 pass, 0 fail** | **`0.000e+00` on every row**, both classes, both modes |
| `fdf8dbe3` (`evanescent`) | **460 pass, 0 fail** | 7.240e-05 (`multimode_field_plasma` adaptive statistics) |

The `fa556e6f` gate is exactly zero on all 21 cases: the statistics class as well as the
`Eω` class, in both the fixed and the adaptive mode. Nothing about the default CPU path
moved, which is what the two-branch design is for -- the host branches are the same
expressions as before and are selected by `Utils.isdevice`, so a `Float32` CPU run is
unchanged too.

Against `evanescent` every non-zero row is one `gpu/int-D` already recorded, to every
digit (the plasma rows from `gpu/13`'s scan against the old serial `cumtrapz!`, the χ⁽²⁾
rows from `gpu/15`'s FMA contraction, the radial ulps from `gpu/02`). No case changed its
number of accepted steps against either baseline.

## Tests run

All on M1 Pro / Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread,
wisdom off, through the timing wrapper. Environments as each test file's header
describes: a stacked environment with `JLArrays` for `test_device.jl`, one with `Metal`
for `test_metal.jl`.

| what | result |
| --- | --- |
| `test/test_regression.jl` vs `fa556e6f` | **460 pass, 0 fail**, `0.000e+00` everywhere |
| `test/test_regression.jl` vs `fdf8dbe3` | **460 pass, 0 fail**, the `gpu/int-D` record |
| `test/test_device.jl` (JLArrays) | **741 pass, 0 fail**, 35 testsets (611/31 on the base) |
| `test/test_metal.jl` (Metal 1.11, hardware) | **415 pass, 0 fail**, 17 testsets (386/16 on the base) |
| `test/test_stats.jl` | 1 pass, 0 fail |
| `test/test_output.jl` | 117 pass, 0 fail, 8 testsets |
| `test/test_mixtures.jl` | 2049 pass, 0 fail, 2 testsets |
| `test/test_vectorplasma.jl` | 2 pass, 0 fail |
| `test/test_polarisation_env.jl` | 4 pass, 0 fail |
| `test/test_modes.jl` | 587 pass, 0 fail |
| `test/test_interface.jl` | 340 pass, 0 fail, 11 testsets |

## Benchmarks

`benchmark/stats.jl` (new) times one call of the statistics function, which is what
`Luna.run` does through the output handler on every accepted step. Two columns:
**host+copy**, the state copied to the host and unscaled and then the host branches --
what every device run did before this branch, and what a set containing a user statistic
still does; and **device**, the statistics on the state where it is. The **step** column
is one RK45 step of the same propagation, for scale. M1 Pro, `-t 1`, one FFTW thread, one
BLAS thread, `:estimate`, no wisdom. Argon at 1 bar, mode-averaged, one column.

responses: Kerr

| device | trange | state | host+copy | device | step |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPU Float64 | 400 fs | 1025 | 134.0 µs | 133.7 µs | 218.8 µs |
| CPU Float32 | 400 fs | 1025 | 135.5 µs | 133.5 µs | 182.3 µs |
| Metal Float32 | 400 fs | 1025 | 301.5 µs | **2.595 ms** | 2.00 ms |
| CPU Float64 | 1600 fs | 4097 | 280.6 µs | 279.6 µs | 932.7 µs |
| CPU Float32 | 1600 fs | 4097 | 280.6 µs | 270.2 µs | 760.7 µs |
| Metal Float32 | 1600 fs | 4097 | 526.5 µs | **3.180 ms** | 2.54 ms |
| CPU Float64 | 6400 fs | 16385 | 1.034 ms | 1.017 ms | 5.12 ms |
| CPU Float32 | 6400 fs | 16385 | 1.026 ms | 874.8 µs | 3.43 ms |
| Metal Float32 | 6400 fs | 16385 | 1.281 ms | **3.646 ms** | 3.27 ms |

responses: Kerr + plasma (tabulated rate), i.e. with `Stats.electrondensity`

| device | trange | state | host+copy | device | step |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPU Float64 | 400 fs | 1025 | 156.8 µs | 157.0 µs | 521.2 µs |
| CPU Float32 | 400 fs | 1025 | 156.8 µs | 152.9 µs | 458.1 µs |
| Metal Float32 | 400 fs | 1025 | 330.1 µs | **3.529 ms** | 5.20 ms |
| CPU Float64 | 6400 fs | 16385 | 1.324 ms | 1.316 ms | 9.73 ms |
| CPU Float32 | 6400 fs | 16385 | 1.312 ms | 1.160 ms | 7.67 ms |
| Metal Float32 | 6400 fs | 16385 | 1.590 ms | **4.979 ms** | 5.75 ms |

**On the CPU the two are the same, as they must be** -- the device column takes the host
branch there, so the small differences are the copy the host column makes and noise.

**On Metal the device statistics are 2 to 9 times slower than the host copy.** That is
the opposite of what this branch set out to do, and it is not a defect in any one
statistic. Attribution (same machine, mode-averaged Kerr, one column):

| | analytic transform | one reduction ending in a scalar read | `abs2` broadcast + full copy of it to the host | whole default set |
| --- | ---: | ---: | ---: | ---: |
| CPU Float64, 400 fs | 5.5 µs | 0.50 µs | 0.79 µs | 137.7 µs |
| Metal, 400 fs | 207.6 µs | 389.0 µs | 379.6 µs | 2.789 ms |
| CPU Float64, 6400 fs | 284.5 µs | 8.1 µs | 13.3 µs | 1.018 ms |
| Metal, 6400 fs | 247.4 µs | 431.6 µs | 408.7 µs | 3.647 ms |

Every number in the Metal row is **independent of the grid size**: what is being measured
is the round trip, not the arithmetic. One reduction that ends in a device-to-host scalar
read costs ~400 µs on Metal, and the default set does six of them (`ω0` twice, `energy`,
`peakpower`, `peakintensity`, and `fwhm_t`'s copy), plus one MPSGraph inverse FFT at
~210-250 µs: 6 x 400 µs + 250 µs is the 2.8-3.6 ms measured. The host path instead makes
*one* transfer, of 8-128 kB, and then does the whole set in FFTW and host reductions.

So for the one geometry Luna can currently put on a GPU -- mode-averaged, a single column
-- this branch meets its exit criterion (no per-step copy of the state) but costs more
wall time than the copy it removes. What to do about it is in "Known gaps" below.

## Known gaps and open questions

- **The device statistics are slower than the host copy on Metal, for a single-column
  mode-averaged run** (the benchmark above). The cost is six device-to-host round trips
  per call at ~400 µs each, not the reductions. Two things would change it, neither in
  this branch's scope:
  1. **Batch the transfers.** Every scalar reduction could write into one small device
     buffer which is transferred once per call, which by the numbers above would take the
     set from ~2.8 ms to ~0.6 ms on Metal. That needs a two-phase statistics protocol
     (reduce on the device, finish on the host) rather than the one-phase
     `(d, Eω, Et, z, dz)` contract, so it is a design change, not a tweak. `peakpower`
     and `peakintensity` also recompute the same `maximum(abs2, Et)` that `fwhm_t`
     already brings to the host, which the same mechanism would remove.
  2. **The multi-column transforms** (`gpu/20`, `gpu/21`, `gpu/22`). There the state is
     100-1000x larger, so the copy the host path makes stops being free while the round
     trips stay at ~400 µs. This branch's arithmetic is written for those shapes already
     (every reduction has a `dims=1` form).
  Until then, `stats_period` is the lever on Metal, and it now skips the whole cost
  rather than only the copy.
- **`Stats.electrondensity`'s mode-averaged host branch mutates the shared `Et`.**
  `Maths.oversample(t, Et; factor=1)` returns its argument, so `@. Eto /= sqrt(...)`
  divides the analytic field the whole statistics set shares. Any statistic *after* it in
  the set sees a scaled `Et` -- which in the default sets is nothing (`energy_λ` reads
  `Eω`), but a `userfuns` entry would. This is pre-existing (`evanescent` does the same)
  and is not fixed here, per the project's rule about unrelated bugs. The device branch
  does not mutate `Et`, because it folds the conversion into the rate's field reference
  instead; the difference is unreachable, since a user statistic forces the host copy and
  therefore the host branch for everything in the set.
- **`Stats.peakpower(grid, Eω, window; label)`** (the windowed peak power) is still a
  closure with its own transform and buffers, so it is host-only and makes the set it is
  in host-only. It is not in either default set; `examples/low_level_interface/basic_modeAvg.jl`
  passes it through `userfuns`. Giving it a device form is the same work as the others
  and was left out to keep the diff to the default sets.
- **`peakintensity` of more than one mode has no device form.** The projection is a GEMM
  over the mode axis, not a scalar factor. It only appears in the multimode default set,
  whose transform is host-only anyway and which also contains `fwhm_r` and
  `mode_reconstruction_error`.
- **`ω0` and `energy` reduce in the state's precision on a device.** In `Float32` the
  relative error of a tree reduction over n terms grows like sqrt(n)*eps; at n = 16385
  that is ~5e-6, well inside the 1e-3 the Metal comparison uses but worth knowing if
  someone reads `ω0` to six digits from a `Float32` run. The prefactors and `E_ref` are
  applied in `Float64` afterwards, so only the sum itself is single precision.
- **`Stats._specof` uses `Base.typename(typeof(x)).wrapper`** to get `MtlArray` from
  `MtlArray{ComplexF32, 1, Metal.PrivateStorage}`. That is the array type a `DeviceSpec`
  holds, and the statistics are the only place in Luna which has the array but not the
  spec. It runs once per propagation, at construction. A `Luna.specof(x)` in `Device.jl`
  would be the tidier home for it; it was kept private here to keep the diff off a file
  the other Group E branches are also editing.
- **`Stats.energy`'s probe.** The device weights are validated against the energy
  functional the caller passed on one deterministic probe spectrum. Any energy functional
  is a positive linear functional of `abs2.(Eω)`, and two such functionals that agree on a
  spectrum with no zeros are the same to the precision of the check, so this is sound
  rather than merely likely -- but it is a runtime check, not a type-level one.

## Deviations from GPU_PLAN.md and the brief

1. **The electron density does not use the batched scan.** The brief says "electron
   density (via the batched scan from gpu/13)". `Maths.cumtrapz_scan!` produces the whole
   cumulative integral, of which `Stats.electrondensity` reads only the last element --
   and reading the last element of a device array is scalar indexing, which
   `allowscalar(false)` refuses. The last element of the cumulative trapezoid is the
   trapezoid rule over the whole window, so the device branch is one weighted reduction
   instead: the same quantity, one pass instead of two, and no scalar index. The host
   branch keeps `Maths.cumtrapz!`, which is what keeps the gate exactly zero.
2. **The per-step copy is skipped only when the state is on a device.** The brief says
   "the per-step host copy happens only when a non-device-capable statistic is present".
   `ScaledOutput` also requires the state to be on a device: a scaled *host* run
   (`Float32` on the CPU) still copies and unscales, because the statistics would
   otherwise take their host branches on a scaled field and be wrong by `E_ref^2`. The
   copy is host-to-host there and costs almost nothing.
3. **`collect_stats` gained an `Eref` keyword** (defaulting to `1.0`, so the signature is
   compatible) and prepares its statistics twice when the set turns out not to be
   device-capable. `Stats.default` and `Interface` are unchanged in signature, as the
   brief requires.
4. **`Interface.prop_capillary_args` passes the state, not a host template**, undoing the
   workaround `gpu/11-boundaries-output` added for the `EnvGrid` `plan_analytic`. That is
   a change to a line the brief did not name, but it is the line that decides what the
   statistics are built for and there is no way to build device statistics without it.
5. **`benchmark/stats.jl` is a new file** rather than a section of `benchmark/device.jl`,
   so that the statistics can be timed without re-running the response benchmarks.
