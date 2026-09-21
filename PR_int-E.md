# `gpu/int-E`: Group E integration

Base `gpu/21-free-device` (`d7aeb45d`, which contains `gpu/20-radial-device` at
`73e8dba6`). Merges `gpu/23-tabulated-linop` (`4f389f80`), `gpu/24-stats-device`
(`f7ab7fe0`) and `gpu/22-modal` (`7187c795`), largest last, each as its own merge commit,
followed by one commit of integration work.

After this branch every geometry Luna has runs on a device except the *adaptive*
multimode transverse integral and `prop_gnlse`:

| geometry | transform | device | reached through |
| --- | --- | --- | --- |
| mode-averaged | `TransModeAvg` | yes | `prop_capillary`, low level |
| radial free space | `TransRadial` | yes | low level only |
| 2-D Cartesian free space | `TransFree2D` | yes | low level only |
| 3-D free space | `TransFree` | yes | low level only |
| multimode, fixed quadrature | `TransModalFixed` | yes | `prop_capillary(…; modal_integral=:fixed)`, low level |
| multimode, adaptive cubature | `TransModal` | **no** (host `Float64` only) | the default |
| GNLSE | `TransModeAvg` + `NormModeAvgGNLSE` | **no** | `prop_gnlse` |

## Commits

```
844859e7  Merge gpu/23-tabulated-linop into gpu/int-E
0e09c4bc  Merge gpu/24-stats-device into gpu/int-E
b641025b  Merge gpu/22-modal into gpu/int-E
196c4823  int-E: the semantic sweep, the TransRadial aliasing and the Boundaries element types
```

## Conflicts and how they were resolved

### `gpu/23-tabulated-linop` (one conflict)

- **`test/test_device.jl`**, the seam after the `zcase` helper. A pure add/add: `gpu/23`
  adds the `gradientcase`/`tapercase` one-liners which wrap `zcase`, `gpu/20` adds
  `radialdiff`, `radialcase` and the "Hankel step as one GEMM" testset at the same point.
  Both kept, `gpu/23`'s first so they sit next to the function they wrap.

Everything else auto-merged, as the `gpu/23` review's trial merges predicted, including
`Luna.run`: `gpu/23`'s tabulation block sits after `gpu/20`/`gpu/21`'s
`upload_like(Eω, linop)` and before the metadata writes, which is the order `gpu/23`
tested.

### `gpu/24-stats-device` (two conflicts)

- **`src/Interface.jl`**, the `setup`/`Stats.default` block of `prop_capillary_args`.
  `gpu/23` added the `aefftol`/`aeffspan` keywords to the `setup` call; `gpu/24` replaced
  the `Stats.default` line below it (a host-shaped template → the state itself) and its
  comment. Kept `gpu/23`'s keywords and `gpu/24`'s call and comment.
- **`test/test_metal.jl`**, `metalgradientcase`'s keyword list: `gpu/23` added
  `tabulate_linop`, `linop_tol`, `nsteps` and `rtol`; `gpu/24` added `stats_device`. All
  five kept; the two branches' changes to the body auto-merged.

### `gpu/22-modal` (six conflicts)

- **`src/Luna.jl`, `runscaling`.** `gpu/20`/`gpu/21` added `TransRadial`, `TransFree` and
  `TransFree2D`; `gpu/22` added `TransModalFixed`. All five methods kept. The docstring
  now says what is true after the merge: every transform carries a scaling except
  `TransModal`, which cannot, because its cubature driver is host scalar code returning
  `Vector{Float64}`.
- **`src/Luna.jl`, the `run` comment** on which transforms can produce a device `Eω`:
  rewritten as "every transform except `TransModal`".
- **`src/Interface.jl`, `prop_capillary_args`.** The `setup` call now carries `gpu/23`'s
  `aefftol`/`aeffspan` *and* `gpu/22`'s `modal_integral`/`nr`/`nθ`/`kronrod`, and the
  `Stats.default` call is `gpu/24`'s line with `gpu/22`'s
  `_statskwargs(transform, stats_kwargs)`, which turns off the mode-reconstruction-error
  statistic for `modal_integral=:fixed`.
- **`src/Interface.jl`, the four `setup` methods.** `prop_capillary_args` passes all six
  keywords to whichever method matches, so all four methods now accept all six; the
  mode-averaged pair ignores the modal ones and the multimode pair ignores
  `aefftol`/`aeffspan`. The multimode methods lose their `_cpu_only!` call, as on
  `gpu/22`: `Luna.setup_modal` refuses `:adaptive` on a device itself and names
  `modal_integral=:fixed` as the fix.
- **`src/Interface.jl`, `_cpu_only!`.** `gpu/20`, `gpu/21` and `gpu/22` each rewrote the
  comment block. It is now `prop_gnlse`'s alone; the comment says so and says that the
  free-space device paths are reachable through the low-level interface only. Its error
  message no longer points at mode-averaged propagation as the only device-capable
  alternative.
- **`docs/src/gpu.md`** (admonition, "What runs where", "Performance") and
  **`docs/src/developer/device_model.md`** (the benchmark list, where `gpu/22` appends its
  "The modal transforms" section): see "Documentation" below. Neither side was taken as
  it stood.
- **`test/test_device.jl`**: the import list (`Boundaries` from `gpu/20`, `RectModes` and
  `Random.MersenneTwister` from `gpu/22`) and two add/add helper/testset seams, all kept.
- **`test/test_metal.jl`**: one add/add seam between `gpu/24`'s device-statistics testset
  and `gpu/22`'s modal case, both kept.

`using Luna` compiles after each of the three merges (checked between them), and
`test/test_boundaries.jl`'s host half runs.

## The semantic sweep

- **One `runscaling`** covering `TransModeAvg`, `TransRadial`, `TransFree`, `TransFree2D`
  and `TransModalFixed`. `TransModal` falls through to `UNIT_SCALING`, which is correct
  rather than a gap: it cannot produce a scaled state.
- **`Luna.run`** has `gpu/23`'s tabulation and `gpu/24`'s `ScaledOutput`/`devstats` both
  present and in the order each branch tested: the output is wrapped first (before
  `check_cache` and `Boundaries.setup`, which then see host, physical-unit data), the
  operator is uploaded after `Boundaries.setup`, and the tables are built last, over
  `[z0, zmax + max(max_dz, init_dz)]` with the absorber's `max_dz`.
- **The `Interface` device sentinel** resolves to `Luna.device_request()` for a
  mode-averaged call and for a multimode call with `modal_integral=:fixed`
  (`hasdevicepath`), and to the host otherwise. Radial and free-space device paths are
  reachable through the low-level interface only, because `prop_capillary`/`prop_gnlse`
  never build one.
- **`_cpu_only!`** is `prop_gnlse`'s only. Multimode `:adaptive` is refused by
  `Luna.setup_modal`, which knows which transverse integral was asked for.
- **The statistics for a multimode `:fixed` run** still rebuild for the host, as
  `gpu/24`'s review asked to be checked: `Stats.default`'s multimode set contains
  `FWHMr` and `PeakIntensityModes`, neither of which has a device form, so
  `Stats.collect_stats` builds the whole set for a host copy and `ScaledOutput.devstats`
  is false. `_statskwargs` removes `mode_error` for `TransModalFixed`, so the set is
  built at all.
- **Documentation.** `docs/src/gpu.md`'s admonition and "What runs where" are rewritten
  for the post-E state rather than taking either side of the conflict, including the two
  pieces of per-step host work that remain inside otherwise device-capable runs (the
  crystal-optics normalisation stages its refractive indices through the host, and
  `fwhm_r` has no device form). "Performance" keeps all four benchmark scripts and
  `gpu/22`'s multimode table. `device_model.md` keeps all four benchmark entries and
  `gpu/22`'s new section.

## `TransRadial`: `Pωo === Eωo`

`gpu/21` made the argument for `TransFree`/`TransFree2D` and left `TransRadial` out of
scope. It transfers unchanged and is now done: `to_time!` writes the field into `Eωo` and
transforms out of it into `Eto_k`; the two Hankel steps and the response protocol work on
`Eto_k`, `Eto_r`, `Pto_r` and `Pto_k`; `to_freq!` then writes the polarisation into the
same buffer. The two never appear as the input and the output of one FFT call, which a
device plan would reject.

One `(nωo, npol, nr)` complex array per radial run. `TransRadial` now holds **five**
field-sized buffers, not six, and seven with the modified shot-noise model instead of
eight. `test/test_device.jl` asserts the count by `objectid`, as it does for the two
Cartesian transforms, for a field and an envelope grid.

## `Boundaries.clampdecay`, `addloss`, `addloss_k`

`Luna.run` calls `Boundaries.setup` *before* it uploads the operator, so in a propagation
these three always see a host `Float64` operator and nothing about the default path
changes. A caller who uploads the operator first, or builds it on a device — which
`benchmark/free.jl` does — used to get a `ComplexF32` device operator promoted back to
`ComplexF64` on the host, which Metal's kernel compiler rejects outright.

`Boundaries._likeop` now converts a scalar `ratemax` with `Luna.scalar` and an array `α`
with `Luna.upload_like`, so the result keeps the operator's array type and element type.
The closure forms cannot know the operator's type until they are called, so they build the
converted array on the first call and keep it. `α./2` is formed once at construction
rather than per call — the same arithmetic, one allocation per propagation instead of one
per stage.

On the host `Float64` path `upload_like` returns its argument and `scalar` is a no-op, so
the path the gate exercises is bit-identical; the gate against the merge commit is exactly
zero on all 22 cases.

New testset in `test/test_boundaries.jl`, "the operator absorbers on an uploaded
operator" (19 assertions, JLArrays): the operator uploaded to a `JLArray{ComplexF32, 3}`
first, all three functions in both their array and their closure form, each checked for
array type, element type, agreement with the host result to 1e-6, and — for the closures —
that the second call reuses the conversion rather than rebuilding it.

## The statistics threshold, re-measured

`gpu/24` left `Stats.STATS_DEVICE_MINLEN` at `2^22` elements with the note that no state
Luna produced could reach it, and asked `gpu/int-E` to re-measure once a device-capable
multi-column transform existed. Measured on real radial (`gpu/20`) and 3-D Cartesian
(`gpu/21`) Metal states. M1 Pro, Metal 1.11, Julia 1.13.0, `-t 1`, one FFTW thread, one
BLAS thread, `:estimate`, no wisdom, `allowscalar(false)`. One call of the statistics
function, which is what `Luna.run` does through the output handler on every accepted step;
"one step" is one RK45 step of the same propagation, for scale.

The set is `ω0`, `energy`, `peakpower`, `fwhm_t` and `density` — the device-capable part
of the default mode-averaged set. Two notes on how it had to be built:

- `energy` uses the grid's own spectral functional (`Fields.energyfuncs(grid)[2]`), not the
  free-space one. A free-space energy functional integrates over the transverse grid as
  well, so it has no weight-vector form, `Stats._devenergy` refuses it, and the statistic —
  and with it the whole set — would be host-only. The arithmetic and the number of round
  trips are the same either way, which is what is being timed.
- The state is handed to the statistics as a native `(nω, ncols)` device array holding the
  same numbers. `Stats.squeeze` has methods for a vector and a matrix only, so the default
  statistics **do not run on a 3-D `(nω, npol, nr)` state at all** — on the host either,
  and on `evanescent` too (see "Known gaps"). The flattened copy has the same element
  count and the same columns, so every reduction, round trip and transfer is the one a
  radial or 3-D run would pay.

| state | elements | host + copy | device | one step | `:auto` chose |
| --- | ---: | ---: | ---: | ---: | --- |
| radial, `N = 256` | 33024 | 2.03 ms | 333 ms | 3.49 ms | host |
| radial, `N = 1024` | 132096 | 7.65 ms | 1.44 s | 10.5 ms | host |
| 3-D, 64 x 64 | 524288 | 662 ms | 615 ms | 9.28 ms | host |
| 3-D, 128 x 128 | 2097152 | 5.23 s | 4.88 s | 26.1 ms | host |

**Decision: `STATS_DEVICE_MINLEN` stays at `2^22`.** The device path does not pay at any
size measured. It is 160x worse on a radial state, and on the two 3-D states — the only
ones within a factor of two of the threshold — it is 7 % and 4 % better, which is not
worth choosing a path that is two orders of magnitude worse one geometry over.

Per-statistic attribution on the same real radial `N = 256` state (host branch against
device branch, one call each):

| statistic | host | device |
| --- | ---: | ---: |
| `ω0` | 84.6 µs | 791 µs |
| `energy` | 46.6 µs | 368 µs |
| `peakpower` | 284 µs | 961 µs |
| `fwhm_t` | 1.38 ms | 303 ms |
| whole set | 1.67 ms | 347 ms |

Two things come out of it.

1. **`fwhm_t` is the whole radial gap** (303 ms of 333). Its device branch reduces the
   *scaled* field: `f.Pd .= abs2.(Et)` runs on the state as the stepper holds it, so the
   columns far off axis, where the physical field is 20 orders below the peak, underflow
   to exactly zero in `Float32`, and the host root-finding which follows is far slower on a
   flat-zero column than on the small but non-zero numbers the host path hands it. The
   *values* are unaffected (an FWHM is scale-invariant, and a zero column has no FWHM
   either way); the cost is not. Recorded as a gap, not fixed here.
2. **Neither path is usable on a large free-space grid.** Both are dominated by the same
   per-column host root-finding in `fwhm_t`: 0.6 s at 4096 columns and 5 s at 16384,
   against a 9 ms and a 26 ms step. `stats_period`, or `Output.nostats`, is the answer
   there, and `docs/src/gpu.md` now says so.

The measurement script is `scratchpad/statsthresh.jl` and the attribution
`scratchpad/statsattrib2.jl`; neither is committed, because `benchmark/stats.jl` already
covers the mode-averaged case and these two build states no benchmark environment can
reach without a GPU.

## The regression gate

M1 Pro, Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread, wisdom off,
`shotnoise=false`. 22 cases in two modes each. `radial_field_raman` (added by `gpu/20`)
had no baseline before this branch; it was generated for `fdf8dbe3`, `782f55d1` and
`fa556e6f` from their own worktrees with the current `cases.jl`, so all four baselines
below hold the full matrix and the gate runs 22 of 22 with no filter.

| baseline | what it answers | result | largest difference |
| --- | --- | --- | --- |
| `b641025b` (this branch's last merge commit) | did the integration work after the merges change anything | **466 pass, 0 fail** | **0.000e+00**, every case, both modes, both classes |
| `fa556e6f` (source-identical to `gpu/int-D`) | what Group E did | **466 pass, 0 fail** | 3.043e-04 (`multimode_field_plasma` adaptive statistics) |
| `782f55d1` (`gpu/int-A`) | Groups B-E | *see below* | |
| `fdf8dbe3` (`evanescent`) | the whole project so far | *see below* | |

### Against `fa556e6f`: the Group E table

Every row is exactly `0.000e+00` except the two which go through a modal transform:

| case | mode | Eω diff | Eω tol | stats diff | stats tol | worst quantity |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `modeavg_field_vector` | fixed | 1.452e-15 | 1.0e-12 | 4.676e-13 | 8.0e-11 | `stats/transverse_integral_error_rel` |
| `modeavg_field_vector` | adaptive | 3.086e-11 | 6.3e-09 | 5.159e-06 | 1.1e-03 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | fixed | 6.913e-14 | 6.7e-12 | 5.417e-13 | 5.8e-11 | `stats/transverse_integral_error_abs` |
| `multimode_field_plasma` | adaptive | 3.146e-07 | 7.3e-05 | 3.043e-04 | 7.3e-02 | `stats/transverse_integral_error_rel` |

These are `gpu/22`'s rows, reproduced here to every digit: the batched column evaluator
changes the order of accumulation in the transverse integral, and the worst quantity in
each is `HCubature`'s own error *estimate*, which moves by far more than the integral does
because the adaptive rule subdivides differently when the integrand changes in its last
bits. The accepted step counts are unchanged. Nothing from `gpu/23` (tabulation is
opt-in), `gpu/24` (two branches, host branch transcribed) or the integration work
(`gpu/20`'s and `gpu/21`'s device paths are not on the default CPU path) moves any row,
which is what the `b641025b` row says independently.

### Against `782f55d1` (`gpu/int-A`): Groups B to E

**466 pass, 0 fail**, largest difference 2.319e-04. The non-zero rows:

| case | mode | Eω diff | Eω tol | stats diff | stats tol | worst quantity |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `modeavg_field_plasma` | fixed | 1.131e-15 | 1.0e-12 | 4.769e-15 | 1.0e-12 | `stats/peak_ionisation_rate` |
| `modeavg_field_plasma` | adaptive | 5.616e-09 | 9.2e-07 | 7.103e-06 | 1.2e-03 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | fixed | 1.245e-15 | 1.0e-12 | 5.202e-15 | 1.0e-12 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | adaptive | 2.283e-13 | 4.8e-10 | 8.550e-08 | 1.8e-04 | `stats/peak_ionisation_rate` |
| `modeavg_field_vector` | fixed | 1.133e-15 | 1.0e-12 | 6.029e-13 | 8.0e-11 | `stats/transverse_integral_error_rel` |
| `modeavg_field_vector` | adaptive | 2.271e-10 | 6.3e-09 | 3.809e-05 | 1.1e-03 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | fixed | 5.789e-14 | 6.7e-12 | 4.768e-13 | 5.8e-11 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | adaptive | 2.383e-07 | 7.3e-05 | 2.319e-04 | 7.3e-02 | `stats/transverse_integral_error_rel` |
| `free2d_field_chi2` | fixed | 7.683e-16 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_field_chi2` | adaptive | 4.841e-14 | 7.4e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_env_chi2` | fixed | 8.348e-16 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_env_chi2` | adaptive | 3.151e-14 | 5.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |

The plasma rows are `gpu/13`'s scan against the old serial `cumtrapz!`, the χ⁽²⁾ rows are
`gpu/15`'s FMA contraction, and the two modal rows are `gpu/22`'s accumulation order added
to what `gpu/int-D` already recorded. Every row is inside its tolerance and no case changed
its accepted step count.

## Known gaps and open questions

Carried forward from the branch PRs unless marked new.

- **New: the default statistics do not run on a 3-D `(nω, npol, nr)` state.**
  `Stats.squeeze` has methods for `Array{T,1}` and `Array{T,2}` only, so `Stats.ω0` (and
  anything else which squeezes) raises a `MethodError` on a radial or free-space state,
  on the host as well as on a device. This is pre-existing — `evanescent` has the same two
  methods — and it is why the free-space examples pass `Output.nostats` and comment out the
  `Stats.collect_stats` line. Nothing on this branch made it worse; it is recorded because
  the statistics threshold measurement is the first thing in the project to try.
- **New: `Stats.fwhm_t`'s device branch reduces the scaled field**, so on a multi-column
  state the columns which underflow to zero in `Float32` make the host root-finding which
  follows far slower than on the host path (the table above). The values are unaffected.
- **`Stats.electrondensity`'s mode-averaged host branch scales the shared `Et` in place**
  (`Maths.oversample(t, Et; factor=1)` returns its argument), so a user statistic *after*
  it in a set sees a scaled field. Pre-existing on `evanescent`, not fixed on `gpu/24`,
  not fixed here.
- **`x1 >= ul[2]` in the Cartesian branch of the adaptive modal point rule** looks like a
  typo for `x2 >= ul[2]` and drops a strip of the guide for an `a > b` `RectModes.RectMode`
  with `full=true`. Pre-existing; `gpu/22`'s review worked out what it costs
  (`PR_22-modal.md`, known gap 6). Fixing it changes the answer of every such run, so it
  wants its own commit with a regression case. User decision.
- **`Stats.default` raises a `MethodError` for a `TransModalFixed` at the low level**
  unless `mode_error=false` is passed; `prop_capillary` is protected by
  `Interface._statskwargs`. `gpu/25-modal-error-stat` fixes it together with the Kronrod
  error statistic (`NonlinearRHS.integral_error!` exists and is tested; nothing calls it
  per step).
- **`Stats.peakintensity` of more than one mode, `fwhm_r` and the mode reconstruction
  error have no device form**, so a multimode statistics set is always computed on a host
  copy, whichever transverse integral produced it.
- **The crystal-optics normalisation evaluates its refractive indices on the host** and
  stages the result; the propagation itself stays on the device.
- **A tapered mode collection with `modal_integral=:fixed` re-evaluates the mode fields on
  the host** whenever `z` moves. The adaptive rule does the same at every point of every
  round, so it is not a regression; `Modes.scale_invariant` from the fork, or `gpu/23`'s
  tabulation, would fix it.
- **`tabulate_linop` is a change of discretisation, not of the equation.** A result
  computed with it is not comparable element by element with one computed without it
  (`PR_23-tabulated-linop.md`, "What it does to the answer"). Whether it should become the
  default is the user's decision at the Group E stop.
- **A free-space energy functional has no device form.** `Stats._devenergy` refuses
  `Fields.energyfuncs(grid, spacegrid)[2]` because it integrates over the transverse grid
  as well, so `Stats.energy` built from it is host-only and makes its whole set host-only.
  Given the threshold measurement above this costs nothing today.

## Base for Group F and E2

`gpu/int-E` is the base for `gpu/25-modal-error-stat` (Group E2) and for Group F
(`gpu/30-neff-grid`, `gpu/31-ka-kernels`, `gpu/32-docs`). Two things for them:

- The regression baselines for all four commits used here are in
  `Luna.Utils.cachedir()/regression/<sha>` and now include `radial_field_raman`, so a
  Group F branch can run the gate against `fa556e6f`, `782f55d1` or `fdf8dbe3` with no
  filter. A branch cut from `gpu/int-E` should generate its own baseline from
  `test/regression/generate.jl gpu/int-E` and require exactly `0.000e+00`.
- `docs/src/gpu.md` and `docs/src/developer/device_model.md` are consistent as of this
  branch but still carry the "Work in progress" admonition and per-branch measurement
  tables; `gpu/32-docs` consolidates them.

### Against `fdf8dbe3` (`evanescent`): the whole project so far

**466 pass, 0 fail**, largest difference 2.319e-04 — the same number and, row for row, the
same values as against `gpu/int-A`, plus the radial rows `gpu/02-radialgrid` recorded:

| case | mode | Eω diff | Eω tol | stats diff | stats tol | worst quantity |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `modeavg_field_plasma` | fixed | 1.131e-15 | 1.0e-12 | 4.769e-15 | 1.0e-12 | `stats/peak_ionisation_rate` |
| `modeavg_field_plasma` | adaptive | 5.616e-09 | 9.2e-07 | 7.103e-06 | 1.2e-03 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | fixed | 1.245e-15 | 1.0e-12 | 5.202e-15 | 1.0e-12 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | adaptive | 2.283e-13 | 4.8e-10 | 8.550e-08 | 1.8e-04 | `stats/peak_ionisation_rate` |
| `modeavg_field_vector` | fixed | 1.133e-15 | 1.0e-12 | 6.029e-13 | 8.0e-11 | `stats/transverse_integral_error_rel` |
| `modeavg_field_vector` | adaptive | 2.271e-10 | 6.3e-09 | 3.809e-05 | 1.1e-03 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | fixed | 5.789e-14 | 6.7e-12 | 4.768e-13 | 5.8e-11 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | adaptive | 2.383e-07 | 7.3e-05 | 2.319e-04 | 7.3e-02 | `stats/transverse_integral_error_rel` |
| `radial_field_kerr` | fixed | 1.901e-15 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `radial_field_kerr` | adaptive | 6.713e-11 | 1.2e-08 | 0.000e+00 | 1.0e-12 | `Eω` |
| `radial_env_kerr` | fixed | 2.866e-15 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `radial_env_kerr` | adaptive | 4.017e-11 | 8.2e-09 | 0.000e+00 | 1.0e-12 | `Eω` |
| `radial_field_raman` | fixed | 2.134e-15 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `radial_field_raman` | adaptive | 1.951e-11 | 1.1e-09 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_field_chi2` | fixed | 7.683e-16 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_field_chi2` | adaptive | 4.841e-14 | 7.4e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_env_chi2` | fixed | 8.348e-16 | 1.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |
| `free2d_env_chi2` | adaptive | 3.151e-14 | 5.0e-12 | 0.000e+00 | 1.0e-12 | `Eω` |

`radial_field_raman` is new to the matrix (`gpu/20`) and behaves like the other two radial
cases: it moves against `evanescent` by the `Grid.RadialGrid` ulps `gpu/02` recorded and by
nothing against `gpu/int-A` or `gpu/int-D`. Every row is inside its tolerance; no case
changed its accepted step count against any of the four baselines.

## Tests

M1 Pro, Julia 1.13.0, `-t 1` unless stated, `Luna.set_fftw_mode(:estimate)`,
`set_fftw_threads(1)`, `BLAS.set_num_threads(1)`, `set_fftw_wisdom(false)`, through the
timing wrapper. `test_device.jl` and `test_boundaries.jl` need `JLArrays`, `test_metal.jl`
needs `Metal`, both from a stacked environment as each file's header describes.

| file | result |
| --- | --- |
| `test/test_regression.jl` vs `b641025b` | **466 pass, 0 fail**, `0.000e+00` on every row |
| `test/test_regression.jl` vs `fa556e6f` | **466 pass, 0 fail** |
| `test/test_regression.jl` vs `782f55d1` | **466 pass, 0 fail** |
| `test/test_regression.jl` vs `fdf8dbe3` | **466 pass, 0 fail** |
| `test/test_device.jl` (JLArrays), `-t 1` | **1148 pass, 0 fail, 68 testsets** |
| `test/test_device.jl` (JLArrays), `-t 4` | **1148 pass, 0 fail, 68 testsets** |
| `test/test_metal.jl` (Metal 1.11, hardware) | **721 pass, 0 fail, 30 testsets** |
| `test/test_interface.jl` | 360 pass, 0 fail, 13 testsets |
| `test/test_linops.jl` | 270 pass, 0 fail, 5 testsets |
| `test/test_modes.jl` | 724 pass, 0 fail, 8 testsets |
| `test/test_multimode.jl` | 6 pass, 0 fail, 5 testsets |
| `test/test_freespace.jl` | 77 pass, 0 fail, 43 testsets |
| `test/test_boundaries.jl` (JLArrays) | 203 pass, 0 fail, 12 testsets |
| `test/test_stats.jl` | 1 pass, 0 fail |
| `test/test_output.jl` | 117 pass, 0 fail, 8 testsets |

`test_device.jl` gives the same 1148 at one and at four Julia threads, which is the point
of running it twice: nothing in the device paths depends on the thread count.

### Tests added on this branch

- `test/test_device.jl`, **"TransRadial holds five field-sized buffers"** (12): `Pωo ===
  Eωo` by identity, the field-sized buffer count for a field and an envelope grid, that
  the four time-domain buffers are four distinct arrays, and seven buffers with the
  modified shot-noise model.
- `test/test_boundaries.jl`, **"the operator absorbers on an uploaded operator"** (19):
  described under "`Boundaries.clampdecay`, `addloss`, `addloss_k`" above.

### Parsing and examples

All 139 `.jl` files under `src/`, `test/` and `examples/` parse.

One example per geometry run through `Luna.run` (the propagation part of each, up to the
first plotting call):

| example | geometry | time |
| --- | --- | ---: |
| `low_level_interface/basic_modeAvg.jl` | mode-averaged, field | 7.8 s |
| `low_level_interface/basic_modeAvg_env.jl` | mode-averaged, envelope | 2.0 s |
| `low_level_interface/basic_modal.jl` | multimode, radial integral | 159.3 s |
| `low_level_interface/full_modal/basic_modal_full.jl` | multimode, full 2-D integral | 456.4 s |
| `low_level_interface/rectangular/rectangular_modal.jl` | multimode, Cartesian domain | 170.4 s |
| `low_level_interface/freespace/radial.jl` | radial free space | 32.3 s |
| `low_level_interface/freespace/free2D.jl` | 2-D Cartesian free space | 9.8 s |
| `low_level_interface/freespace/free2D_bbo.jl` | 2-D Cartesian, χ⁽²⁾ | 16.4 s |
| `low_level_interface/freespace/full3D_env.jl` | 3-D free space | 82.5 s |
| `low_level_interface/gradients/gradient_modeAvg.jl` | pressure gradient | 10.2 s |
| `low_level_interface/tapers/taper_modeAvg.jl` | taper | 20.1 s |
| `low_level_interface/gnlse/simplescg_modeAvg_env.jl` | GNLSE, low level | 23.2 s |
| `simple_interface/gnlse_sol.jl` | GNLSE, simple interface | 2.4 s |

(`free2D.jl` calls `pygui(true)` on line 5, so the runner's "stop at the first plotting
call" rule cut the whole file; it was re-run with `pygui` treated like the `PyPlot` import
and takes 9.8 s.)

### Documentation

`include("docs/make.jl")` still ends with `makedocs` refusing to render on unresolved
cross-references, as it does on `gpu/int-D`. **Six** unresolved references here against
`gpu/int-D`'s **nine**, and none of the six is new:

| reference | on `gpu/int-D` | here |
| --- | ---: | ---: |
| `LinearOps.βz` | 3 | 2 |
| `AbstractOutput` | 1 | 1 |
| `loadFFTwisdom`, `saveFFTwisdom` | 2 | 2 |
| `Luna.PhysData.crystal_internal_angle` | 1 | 1 |
| `LinearOps.make_const_linop`, `norm_free` | 2 | 0 (fixed by `gpu/23`) |

Two references *were* new on the first build of this branch and are fixed in `90fed750`:
`Boundaries.addloss` (the element-type commit had put two helper definitions between its
docstring and its first method) and `Luna.prop_capillary` in `LinearOps.jl` (which does not
resolve from inside `module LinearOps`).
