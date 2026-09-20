# The output and precision boundary; absorbing boundaries on a device

Branch `gpu/11-boundaries-output`, base `gpu/10-device-model` (`c1ffc410`). GPU_PLAN.md
§3, §4.1, §4.2, §4.7, §4.8, §4.11, §6 Group B.

This is the second branch of Group B. It makes the absorbing boundaries
(`Boundaries.RateAbsorber`/`LegacyAbsorber`, the transverse collars) and the output
device/precision boundary (`Luna.ScaledOutput`) work for the mode-averaged transform on
any array type, and threads `device`/`precision`/`stats_period` through the simple
interface. After this branch, `using Metal; prop_capillary(...)` with Kerr only
(constant and gradient pressure, `boundary=:rate`, default statistics) runs on the GPU.

**The default CPU Float64 path is unchanged: the regression gate is exactly `0.000e+00`
for all 21 cases, in both modes and both classes, against both `gpu/10-device-model` and
`gpu/int-A`** -- not merely within tolerance. See "Why the default CPU path did not
move".

## Motivation

`gpu/10-device-model` made mode-averaged Kerr propagation itself run on any array type,
but stopped at the edges: the absorbing boundaries were scalar loops over host arrays,
`Luna.run` refused anything but `boundary=:none` on a device, `Output.jl` allocated
`ComplexF64` unconditionally and never unscaled a `Float32` run, and `Interface.jl`
hardcoded `device=Luna.HostSpec()` at every call site so the simple interface never saw
the device model at all. This branch closes those gaps for the mode-averaged transform.

## What changed

### `src/Output.jl`

Stays device-unaware -- it still only knows how to save an array and a dictionary.

- `Output.willsave(o, y, t, dt)`: whether calling `o` right now would save at least one
  data point, without saving anything. `MemoryOutput`/`HDF5Output` methods for a
  `GridCondition` save condition (stateless, safe to call speculatively); the fallback is
  the conservative `true`.
- `Output.PeriodicStats(statsfun, period)`: evaluate a statistics function only every
  `period`-th accepted step (first call always fires), returning `nothing` in between.
  `MemoryOutput`/`HDF5Output`'s call operators now skip a `nothing` statistics result
  instead of erroring on it; `HDF5Output`'s cachehash/`append_stats!` tolerate a step
  with nothing collected since the last save.
- `MemoryOutput`/`HDF5Output` allocate their solution array with `eltype(y)` instead of a
  hardcoded `ComplexF64`, so a `Float32` run saves `Float32`.
- `HDF5Output`'s `cachehash` includes `eltype(y)`, so resuming a run in a different
  precision than the cache was written in is refused rather than misinterpreted.

### `src/Device.jl` (included after `Output.jl` now, not right after `Utils.jl`)

`Luna.ScaledOutput(o, y, Eref)`: the output/precision boundary. `Luna.run` wraps the
output handler in one whenever the state is on a device or the run is scaled
(`E_ref != 1`, i.e. every `Float32` run, host or device); on the default CPU `Float64`
path the handler is passed through unchanged. Two reusable host buffers in the state's
element type: `ybuf` for the unscaled per-step `y` (statistics, an `HDF5Output`'s resume
cache -- filled only when needed, the cache case gated by `willsave`), `ibuf` for a saved
field (filled lazily inside the closure `Luna.run`'s `stepfun` already passes as `yfun`,
so a non-saving step costs nothing). Warns once per propagation when a device run forces
a host copy for statistics. The stepper's own arrays are never modified.

`Device.jl` moved after `Output.jl` in `Luna.jl`'s include list so that `ScaledOutput`
can dispatch on `Output.MemoryOutput`/`HDF5Output` and extend `Output.willsave`/
`check_cache`/`hasdata`. Nothing between `Utils.jl` and `Output.jl` needed `Device.jl`.

### `src/Luna.jl`

- The restriction to `boundary=:none` on a device is removed.
- `Et` (the boundary/collar scratch buffer) is built with `similar(Eω, ...)`, sized and
  typed from the grid and the state, instead of an actual inverse FFT via `\`, which the
  generic device planners do not support. Its contents never mattered -- every absorber
  overwrites it before reading it back.
- `check_cache`'s resumed field is rescaled (divided by the new run's `E_ref`) as well as
  reuploaded, since the cache is always written unscaled and in physical units.
- A new `runscaling(transform)` returns the `UnitScaling` a transform carries
  (`NonlinearRHS.TransModeAvg` does; every other transform defaults to the identity),
  used to decide whether `ScaledOutput` is needed.

### `src/Boundaries.jl`

`RateAbsorber`, `LegacyAbsorber` and the transverse collars (`RadialCollar`,
`CartesianCollar`) rewritten as broadcasts and reductions over mirrored arrays, on every
backend:

- Both absorbers hold an explicit inverse plan (`IFT = Utils.plan_ift(FT)`) and use
  `mul!` instead of `ldiv!`, matching `NonlinearRHS.to_time!`'s discipline and required
  for the generic device planners (which do not implement `\`).
- The temporal collar's per-step factor and its application to `Et` are plain broadcasts
  over the whole time axis rather than a loop restricted to the indices where `αt > 0`:
  the absorber profile is exactly 1 in the flat interior, so multiplying there is an
  exact no-op -- this is why the gate is `0.000e+00` for `:rate` and `:legacy`, not
  merely within tolerance.
- The removed-energy bookkeeping (the one-time "boundary has removed X%" warning only --
  not part of `Eω` or any recorded statistic) is one fused reduction (`_absorbed`,
  shared by all three structs) instead of a scalar accumulation inside the mutating loop.
- `RadialCollar`'s `αr`/`weight` are stored in the state's real precision (previously
  always `Float64`); `CartesianCollar`'s `αxy` is mirrored with `Luna.upload_like`. Both
  stay host-only: `TransRadial`/`TransFree*` have no device path yet (`gpu/21`).
- `LegacyAbsorber`'s `ωwin`/`twin` multiplies use mirrors of the grid's windows
  (`upload_like`) instead of the grid's own `Float64` vectors directly.

`ridcs`/`tidcs` index vectors are removed rather than converted to masks (a deviation
from the plan's wording; see "Deviations" below).

### `src/Interface.jl`

`prop_capillary`/`prop_gnlse` gain `device` (default `nothing`, a sentinel meaning "not
specified" -- see review round 1, finding 2, below for why it is not simply
`Luna.device_request()`), `precision` (default `nothing`) and `stats_period` (default
`1`) keywords, replacing the `device=Luna.HostSpec()` `gpu/10` hardcoded at every call
site. An untouched call reproduces today's behaviour: `prop_capillary_args` resolves the
sentinel to `Luna.device_request()` (following a loaded GPU package) only when the
propagation is mode-averaged *and* every response it was built with is device-capable
(`Nonlinear.device_capable`); otherwise, and always for `prop_gnlse`, it resolves to the
CPU regardless of `Luna.settings["device"]`, exactly as before this keyword existed.

Multimode and radial propagation (`modes` a collection) and `prop_gnlse` are not device-
or reduced-precision-capable yet: their `setup`/`prop_gnlse_args` methods accept and
validate `device`/`precision` with a new `_cpu_only!` helper, erroring with a message
naming the actual limitation rather than an unrelated `MethodError` or a silent fallback
to the CPU that the caller did not ask for. Mode-averaged propagation with a response
that is not device-capable (plasma, Raman) gets the analogous
`_check_responses_device_capable!` for an *explicit* incompatible request.

Found on Metal hardware (see "What Metal caught"): `Stats.default`'s `Eω` argument is
used only to size its buffers at construction, but `Stats.jl`'s `EnvGrid` code path plans
an FFTW transform directly on a copy of it, which fails when `Eω` is a device array.
`prop_capillary_args` now passes `Luna.tohost(Eω)` in place of a device `Eω` for that one
call; `Stats.jl` itself is untouched.

`boundary_kwargs` is unchanged: it already filters the keywords it forwards to
`Luna.run` down to `(:boundary, :boundary_N, :boundary_length, :tcollar)`, so the new
keywords being present in the same captured `kwargs` (in `prop_capillary`/`prop_gnlse`,
which still capture everything not explicitly named) is harmless.

`saveargs` records the three new keywords as the caller passed them, like every other
keyword there.

## Unit scaling

Unchanged from `gpu/10-device-model`; this branch is the "unscaling happens in one place
only, the output wrapper" that PR was written in anticipation of. `ScaledOutput` is that
wrapper.

## Why the default CPU path did not move

- `Boundaries.RateAbsorber`'s temporal collar was rewritten from a loop restricted to
  `tidcs` to a full-array broadcast; the profile is exactly `1.0` (the rate exactly `0`)
  outside the collar by construction, so multiplying the whole array does the identical
  arithmetic as the loop did in the interior (`x * 1.0 == x` bit for bit) and the
  original arithmetic in the collar. `LegacyAbsorber`'s `ωwin`/`twin` mirrors are the
  identity on the host (`upload_like` returns the grid's own vectors unchanged).
- Replacing `ldiv!(Et, FT, Eω)` with `mul!(Et, IFT, Eω)`, `IFT = Utils.plan_ift(FT)`, is
  the same FFTW inverse plan applied the same way; `ldiv!` on an FFTW plan is defined in
  terms of exactly this.
- The removed-energy bookkeeping's new formula is not part of `Eω` or any recorded
  statistic, so a rounding-level (or larger) change there is invisible to the gate.
- `ScaledOutput` is never constructed on the default path (`isdevice(Eω)` is false and
  `isunity(scaling)` is true), so `output` is the same object it always was.
- `runscaling`'s fallback (`UNIT_SCALING`) and `Stats.default`'s host-template fix apply
  only when `Eω` is a device array.

## Tests

`julia --project=<worktree> -t 1`, with `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)`. Apple M1 Pro, Julia 1.13.0.

### The regression gate

```
LUNA_REGRESSION_BASE=c1ffc410 julia --project=$PWD -t 1 test/test_regression.jl
```

against a baseline generated from the base commit `c1ffc410`
(`test/regression/generate.jl c1ffc410`), and again against `gpu/int-A` (`782f55d1`) for
the cumulative view:

**460 pass, 0 fail, both baselines. Every case, both modes, both classes: `0.000e+00`.**

| case | `:fixed` Eω | `:fixed` stats | `:adaptive` Eω | `:adaptive` stats |
|---|---:|---:|---:|---:|
| every one of the 21 cases | 0 | 0 | 0 | 0 |

No case moved against either baseline. Step counts unchanged. The two cases this
branch's collar rewrite could plausibly have moved -- `modeavg_field_legacy` (`:legacy`
absorber) and any `:rate` case -- are exactly `0.000e+00`, not merely within the gate's
tolerance floor.

### New/changed test files

`test/test_boundaries.jl` (CPU gate unchanged; **all pre-existing testsets pass with the
same counts**): new testset **"RateAbsorber and LegacyAbsorber on JLArray"** (gated on
`JLArrays` being available), constructing both absorbers directly on host and `JLArray`
arrays with a shared host-backed `AbstractFFTs` plan shim and comparing one step's
output. 2 assertions, both at exactly the host result (0 relative difference).

`test/test_device.jl`: extended with **"boundaries and default statistics on JLArray"**
(mode-averaged Kerr, `boundary=:rate`, `Stats.default`, field-resolved and envelope) and
updated `kerrcase`/`gradientcase` (drop the now-unnecessary `ToHost` wrapper; `boundary`/
`stats` keywords). "Float32 on the CPU" no longer multiplies by `Eref` manually (the
output is unscaled automatically now) and checks `eltype === ComplexF32`. "a device run
refuses host-only machinery" no longer asserts `boundary=:rate` throws (it does not any
more); the normalisation/rescale refusals are unchanged.

| testset | assertions |
| --- | ---: |
| (13 unchanged testsets from `gpu/10`) | 157 |
| boundaries and default statistics on JLArray | 12 |
| Float32 on the CPU (updated) | 12 |

**169 pass, 0 fail** (15 testsets).

`test/test_output.jl`: four new testsets -- **willsave** (9), **PeriodicStats** (4),
**Float32 output eltype** (5), **HDF5 resume of a Float32 run** (13). All existing
testsets unchanged and pass. **102 pass, 0 fail** (the whole file, 23.3 s).

`test/test_interface.jl`: new testset **"device, precision and stats_period keywords"**
(16 assertions) covering an untouched call, `precision=Float32`, an explicit
`DeviceSpec(Array, Float64)`, `stats_period`, a pressure gradient, the multimode refusal,
and `prop_gnlse`'s refusal of anything but the default. **317 pass, 0 fail** (up from
301), full file run ≈4 minutes.

`test/test_metal.jl` (not part of the suite; run from an environment with Metal): see
"Metal hardware results" below.

### What Metal caught

Running on real hardware (not `JLArrays`) found two genuine bugs, neither of which any
CPU test or `JLArrays` test caught.

**1. `Stats.default` planning an FFTW transform on a device array.** In scope for this
branch to fix at the call site even though `Stats.jl` itself is out of scope:
`Stats.default`/`collect_stats`'s one-time construction call built its internal FFTW
plan from `similar(Eω)`/`copy(Eω)` where `Eω` was the actual (possibly device) state.
For a `RealGrid` this happens to work anyway (`plan_analytic` always allocates a host
`Array` there regardless of `Eω`'s type), but for an `EnvGrid` it inherits `Eω`'s array
type and tries to `FFTW.plan_ifft` directly on a copy of a device array, which fails
("Cannot access the contents of a private buffer" on Metal's private-storage memory).
`JLArrays` does not catch this: it silently presents host-backed memory to FFTW, so the
plan is built (if slowly) rather than refused. Fixed in `Interface.jl` (see above) by
passing a host-shaped template for that one construction call; `Stats.jl` untouched.

**2. `device`'s default resolving through `:auto` for paths that can never honour it.**
The first version of the `_cpu_only!` design (§"Interface.jl" above) gave the multimode/
radial `setup` methods and `prop_gnlse_args` a default of `Luna.device_request()` -- the
same default as the mode-averaged path. That is wrong: once a GPU package is loaded and
`Luna.settings["device"]` becomes `:auto`, `Luna.device_request()` resolves to the GPU,
so an *untouched* multimode or `prop_gnlse` call -- one that never mentioned `device` at
all -- started erroring the moment `using Metal` was run anywhere in the process, which
is a regression from `gpu/10`'s "always the CPU regardless of settings" behaviour and
plainly wrong (a working script should not break because an unrelated `using Metal` was
added). Caught by `test_metal.jl`'s "`Luna.set_device(:cpu)` opts out" testset, which
actually runs `prop_gnlse` and multimode `prop_capillary` under `Luna.set_device(:auto)`
with a real GPU registered -- a scenario no CPU-only test session can produce. Fixed by
giving `prop_capillary_args`'s own `device` keyword a sentinel default (`nothing` =
"not specified") and resolving it per path: `Luna.device_request()` only for a single
mode, `Luna.HostSpec()` (unconditionally) for everything else; `prop_gnlse_args`'s and
the multimode `setup` methods' own defaults are `Luna.HostSpec()` outright, since they
have only one path. An *explicit* `device`/`precision` argument still reaches
`_cpu_only!` and is validated as before.

## Metal hardware results

**Superseded by review round 1, finding 1 -- see "Changes after review round 1" below for
the corrected numbers.** The table below was measured with a testset bug (`href` had no
explicit `device`, so it resolved through the sentinel to `:auto` = Metal in this
Metal-loaded process, making the "exit criterion" comparison Metal against itself) and an
incorrect physical explanation for the resulting `0.0`; kept here, struck through in
spirit, only so the history of what changed is visible.

Apple M1 Pro, Metal.jl v1.11.1, environment built per `docs/src/gpu.md`'s "Running the
hardware tests".

**102 pass, 0 fail** (10 testsets):

| testset | assertions |
| --- | ---: |
| Metal registration | 8 |
| allocation and transfer on Metal | 7 |
| no stray Float64 in the kernels (extended with the boundary kernels) | 18 |
| mode-averaged Kerr on Metal | 20 |
| Metal against the Float64 CPU path | 4 |
| pressure gradient on Metal | 7 |
| **boundaries and default statistics on Metal** (new) | 14 |
| **prop_capillary on Metal** (new, the exit criterion) | 14 |
| **`Luna.set_device(:cpu)` opts out** (new) | 8 |
| Metal refuses what it cannot run | 2 |

| comparison (superseded numbers, see below) | `Eω` | energy |
| --- | ---: | ---: |
| mode-averaged Kerr, field-resolved, Metal vs CPU `Float32` (`boundary=:none`) | 2.7e-7 | -- |
| mode-averaged Kerr, envelope, Metal vs CPU `Float32` (`boundary=:none`) | 5.2e-8 | -- |
| Metal vs CPU `Float64`, He at 0.3 bar (`boundary=:none`) | 2.4e-7 | -- |
| field-resolved, `boundary=:rate` + default statistics | 3.9e-6 | 7.7e-7 |
| envelope, `boundary=:rate` + default statistics | 4.4e-6 | 5.3e-7 |
| `prop_capillary`, constant pressure ("exit criterion", `href` bug) | 0.0 (invalid) | 0.0 (invalid) |
| `prop_capillary`, gradient pressure ("exit criterion", `href` bug) | 0.0 (invalid) | 0.0 (invalid) |

The `boundary=:none` rows above (unaffected by the `href` bug -- `metalcase` always took
an explicit `spec`) still stand and reproduce `gpu/10`'s own numbers.

`Luna.set_device(:cpu)` opts out of a globally-registered GPU (exit criterion 3): verified
with `prop_capillary`, `prop_gnlse` and multimode `prop_capillary` all under
`Luna.set_device(:auto)` with Metal registered and functional. Unaffected by finding 1.

## Benchmarks

`benchmark/boundaries.jl`: the mode-averaged Kerr propagation (He, 1 bar, 1 cm, 20 fixed
steps), timing a whole `Luna.run` under `boundary=:none`/`Output.nostats`,
`boundary=:rate`/`Output.nostats` and `boundary=:rate`/`Stats.default`, on CPU
`Float64`/`Float32` (Metal numbers below, from a separate environment):

```
device           trange      state   none/nostats   rate/nostats     rate/stats
------------------------------------------------------------------------------------------
CPU Float64       400 fs       1025       4.564 ms       4.921 ms       7.742 ms
CPU Float32       400 fs       1025       3.815 ms       4.130 ms       7.054 ms
CPU Float64      1600 fs       4097      19.388 ms      20.979 ms      26.713 ms
CPU Float32      1600 fs       4097      15.465 ms      16.862 ms      22.629 ms
CPU Float64      6400 fs      16385     105.444 ms     113.224 ms     134.997 ms
CPU Float32      6400 fs      16385      70.767 ms      76.276 ms      97.553 ms
```

`boundary=:rate` (`RateAbsorber`'s broadcasts, no statistics) adds 6-9% over
`boundary=:none` at every size -- a full-array broadcast and a fused reduction, once per
accepted step, over the collar rewrite this branch made. The default statistics
(`Stats.default`, unchanged host code, `gpu/24`'s to speed up) dominate the difference
between the second and third columns: 30-60% on top of `boundary=:rate` alone, because
they do a full inverse FFT and several host reductions every accepted step regardless of
device. This is exactly the cost `stats_period` exists to amortise.

Including Metal (`julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/boundaries.jl")'`):

```
device           trange      state   none/nostats   rate/nostats     rate/stats
------------------------------------------------------------------------------------------
CPU Float64       400 fs       1025       4.544 ms       5.055 ms       7.858 ms
CPU Float32       400 fs       1025       3.833 ms       4.143 ms       7.050 ms
metal Float32     400 fs       1025      52.354 ms      64.722 ms      76.668 ms
CPU Float64      1600 fs       4097      19.861 ms      21.380 ms      27.471 ms
CPU Float32      1600 fs       4097      16.151 ms      17.604 ms      23.414 ms
metal Float32    1600 fs       4097      54.683 ms      63.592 ms      75.086 ms
CPU Float64      6400 fs      16385     103.856 ms     111.189 ms     133.085 ms
CPU Float32      6400 fs      16385      70.109 ms      75.234 ms      96.090 ms
metal Float32    6400 fs      16385      57.164 ms      65.514 ms      99.794 ms
```

Consistent with `gpu/10-device-model`'s own finding: Metal is launch-bound at these
sizes (a mode-averaged RHS is only a few thousand elements), so it is 10-14x slower than
either CPU path at 1025 samples and only catches up to CPU `Float64` at 16385 -- still
behind CPU `Float32`. The *relative* cost of `boundary=:rate` and the statistics on Metal
(24% and 46% over `none/nostats` at the smallest size, narrowing to 15% and 75% at the
largest) is the same story as the CPU numbers: a fixed number of extra kernel launches
per step, more expensive relative to the RHS when the RHS itself is cheap. Nothing here
changes GPU_PLAN.md §8's conclusion that a GPU is not the right choice for a
single-column run; the win is in the multi-column transforms, later branches' work.

## Known gaps

- **Only the mode-averaged transform is device-capable.** Radial, free-space and
  multimode propagation, and every response but Kerr, are unchanged host code, refused
  rather than run wrongly on a device or in reduced precision.
- **`RadialCollar`/`CartesianCollar` stay host-only.** `TransRadial`/`TransFree*` have no
  device path (`gpu/21`); the collar rewrite here is about kernel-discipline consistency
  and Float32 correctness on the CPU, not about running them on a device yet.
- **`prop_gnlse` is not device- or reduced-precision-capable.** It builds its own
  normalisation before `Luna.setup` knows the unit scaling, so it cannot yet be given a
  `spec`/`scaling` matching a real device or `Float32` the way the mode-averaged path's
  own default normalisation can. `_cpu_only!` refuses anything but the default rather
  than silently narrowing the request or producing a normalisation wrong by a factor of
  `E_ref`/`P_ref`.
- **Per-step statistics are host code.** `Stats.jl` is unchanged; a device run copies the
  field to the host every accepted step to compute the default statistics and warns
  once. `stats_period` (`Output.PeriodicStats`) reduces how often this happens.
  Device-capable statistics are `gpu/24`'s.
- **Tapers and gradients still upload `β` per stage from the host** (`gpu/10`'s layer 1
  of GPU_PLAN.md §4.5, unchanged here).
- **`ridcs`/`tidcs` became plain full-array broadcasts, not masks** (see "Deviations").

## Deviations from GPU_PLAN.md and the brief

1. **`tidcs`/`ridcs` index vectors removed rather than converted to masks.** The brief's
   wording ("the `ridcs`/`tidcs` index vectors become masks") is satisfied in spirit --
   the per-step application is a broadcast over the whole array, not a loop over
   specific indices -- but there is no explicit `Bool` mask array anywhere. Reasoning:
   the absorber profile is exactly `1.0` outside the collar (`α` exactly `0`), so the
   `removed`-energy formula (`abs2(e)*(1-f²)`) evaluates to an exact `0.0` there without
   needing to be told which elements those are, and the mutating multiply is `x *= 1.0`,
   an exact no-op. A mask would only save redundant arithmetic (computing `exp(0)` and
   multiplying by `1.0` for the whole non-collar region every step), not change the
   answer, and the collar is a small enough fraction of the grid that this was not
   worth the extra array and the masked-`ifelse` indirection. If a reviewer wants the
   explicit mask for performance or documentation clarity, it is a small, low-risk
   follow-up.
2. **`boundary_kwargs` unchanged.** The brief says to "add to the whitelist"; on reading
   it, the new keywords (`device`, `precision`, `stats_period`) are fully consumed by
   `prop_capillary_args`/`prop_gnlse_args` themselves (named keyword parameters there,
   not caught by a trailing `kwargs...`), so there is nothing left for
   `boundary_kwargs` -- which only forwards a filtered subset of `prop_capillary`'s own
   catch-all `kwargs` to `Luna.run` -- to do with them. `Luna.run` has no `device`/
   `precision`/`stats_period` keyword of its own to receive. No code change was needed
   or made; noted here in case the brief meant something this reading missed.
3. **`Stats.default`'s call site in `Interface.jl` was changed**, which is not named in
   the brief, to fix the Metal-only bug described above. `Stats.jl` itself is untouched.
4. **`prop_capillary_args`'s `device` default is `nothing` (a sentinel), not
   `Luna.device_request()` directly.** A first pass gave it `Luna.device_request()`
   outright, matching the brief's "replacing the explicit `device=Luna.HostSpec()`"
   literally; that is wrong (see "What Metal caught", finding 2) because the same
   variable is what the multimode/radial branch would see too, and that branch must stay
   on the CPU by default regardless of `Luna.settings["device"]`. The sentinel lets
   `prop_capillary_args` resolve the default differently per branch while keeping a
   single `device` keyword and a single `_cpu_only!` check for an explicit request.

### Documentation build

`include("docs/make.jl")` reports the same 9 pre-existing unresolved `@ref`s the base
branch (`c1ffc410`) does (`LinearOps.βz` x3, `loadFFTwisdom`, `saveFFTwisdom`,
`AbstractOutput`, `Luna.PhysData.crystal_internal_angle`, `norm_free`,
`LinearOps.make_const_linop`), and none from this branch. `Luna.ScaledOutput`,
`Luna.needs_host_y`, `Luna.needs_host_cache` and `Luna.runscaling` were added to
`docs/src/modules/Luna.md`'s `@docs` block so their `@ref`s resolve.

## Open questions for review

1. Is deviation 1 (no explicit mask) acceptable, or should `RateAbsorber`/`RadialCollar`/
   `CartesianCollar` carry an explicit `Bool` mask for the reduction regardless?
2. `_cpu_only!`'s error message and placement (`Interface.jl`, used by the multimode/
   radial `setup` methods and by `prop_gnlse_args`) -- is refusing at the `Interface`
   layer the right place, or should `Luna.setup` itself refuse for `TransModal`/
   `TransRadial` when given a non-host `device`?
3. `Output.willsave`'s restriction to `save_cond isa GridCondition` (a stateful save
   condition like `every_nth` is not speculatively evaluated, and conservatively returns
   `true`) -- acceptable, or should `willsave` instead take a `mutates::Bool` trait so a
   stateless custom save condition could also benefit?

## Changes after review round 1

The review's verdict was "request changes" on three majors, plus six minors/nits. All
nine are addressed here.

### 1 (major) — the exit-criterion testset compared Metal with Metal

`test_metal.jl`'s "prop_capillary on Metal" built its CPU reference with
`precision=Float32` and no `device`, which resolves through the sentinel to
`Luna.device_request()` -- `:auto` throughout `test_metal.jl` (Metal is loaded and
nothing in that testset had overridden the global setting), which resolves to Metal. Both
sides of the comparison were therefore the same run, and the `0.0` it measured said
nothing about the device path. The physical explanation offered for it (the Kerr phase
below `Float32`'s precision floor) was consequently wrong, since it was never being
tested; the review's own measurement with a genuine CPU Float32 reference gave 3.7e-6 at
the same parameters, not `0.0`.

Fixed:
- Every comparison in `test_metal.jl` now uses an explicit `device` for its host
  reference(s) -- `DeviceSpec(Array, Float32)` and, newly, `HostSpec()` (`Float64`) as
  well, so Metal is checked against both precisions on the host.
- The fibre length is the brief's own (`0.1` m), not the `1e-2` the file had (finding 9).
- `metalcase`/`metalgradientcase` gained a `fixed::Bool` keyword (`min_dz == max_dz`,
  bypassing the step-size controller) so a comparison's tolerance reflects arithmetic
  only, not a possibly different step sequence -- the same principle the regression
  gate's own `:fixed` mode uses. `metalgradientcase` already did this by default.
- "boundaries and default statistics on Metal" now runs the brief's own (weakly
  nonlinear) parameters *and* a visibly nonlinear one (He at 5 bar, 300 µJ, about ×2
  spectral broadening -- checked with a `broadening` helper, the ratio of the last save's
  rms spectral width to the first's), fixed steps, against both host references, for both
  constant and gradient pressure. One adaptive case (the strongly nonlinear gradient) is
  kept with a documented, looser tolerance (see finding 1's own note on this below).
- `precision`'s docstring now says explicitly that it does not select the CPU on its own
  when a GPU package is loaded.
- The wrong `0.0`/precision-floor claim is withdrawn from this document ("Metal hardware
  results", struck through above) and from `docs/src/gpu.md`, replaced with the measured
  numbers below.

One further problem surfaced while fixing this: 600 µJ at 5 bar over the full 10 cm at
only 20 fixed steps is numerically under-resolved (both host and device runs produced
`NaN`) -- a step-count problem with that combination of energy and grid, not a boundary
or device defect (300 µJ at 5 bar, used for the non-gradient visible-Kerr case, is well
resolved at the same step count). The gradient's visibly-nonlinear case uses 300 µJ
instead, matching the non-gradient one.

### 2 (major) — the default `prop_capillary` call failed once a GPU package was loaded

See "What Metal caught", finding 2, above (written before the review; the review found
the same bug independently on hardware). Fixed with `Nonlinear.device_capable` and the
sentinel now checking `all(Nonlinear.device_capable, resp)` as well as the mode count;
an explicit incompatible request now errors before `Luna.setup`, naming `device=:cpu`/
`Luna.set_device(:cpu)`. Verified interactively with a fake registered device forcing
`:auto`: a default (plasma on) call stays on the CPU in `Float64`; a Kerr-only call
follows `:auto`; an explicit device request with plasma raises the new message. New
`test_metal.jl` testset "Luna.set_device(:cpu) opts out" runs `prop_gnlse` and multimode
`prop_capillary` under a real `:auto` with Metal registered.

### 3 (major) — `stats_period` did not skip the device-to-host copy

`needs_host_y` was decided once, at `ScaledOutput` construction, from whether the wrapped
statistics function was `Output.nostats` -- a `PeriodicStats`-wrapped function never is,
so the copy ran on every accepted step regardless of the period; only the statistics
arithmetic was actually skipped. Fixed: `needs_host_y(o, t)` is now evaluated fresh every
step, and a new `Output.willfire(p::PeriodicStats, t)` predicts whether the *next* call to
`p` would run its wrapped function, without running it or mutating `p`. Verified directly
(not by timing, which is noisy under contention on this machine, but by construction): a
`ScaledOutput` built around a `PeriodicStats(period=3)` and called with different `y`
values each step shows its `ybuf` bit-for-bit unchanged on the two steps between fires and
updated only on a fire (`test_device.jl`, "stats_period skips the device-to-host copy",
JLArray, 10 assertions); the one-time host-statistics warning still fires exactly once.

### 4 (minor) — `stats_period` unvalidated; the "every Δz" form missing

Both `Interface.jl` call sites now go through a new `Output.maybe_periodic(statsfun,
period)`, which validates unconditionally (an integer must be `>= 1`, a non-integer
`> 0`) and only skips wrapping for the literal trivial value (integer `1`). `PeriodicStats`
gained the distance form: a non-integer `Real` period now means "every `period` metres of
propagation", firing when the propagation coordinate has advanced by at least that much
since the last fire (`Output.willfire` and the call operator both switch on it).
Documented in `prop_capillary`'s/`prop_gnlse`'s `stats_period` docstrings.
`test_output.jl`'s "PeriodicStats" grew from 4 to 19 assertions: `willfire`, both invalid
classes in both modes, the distance form, and `maybe_periodic`.

### 5 (minor) — `RadialCollar`'s precision parameter was inert

`spatialcollar`'s `RadialCollar` branch built its buffer with a hardcoded
`zeros(ComplexF64, ...)`; now `zeros(Complex{real(eltype(Et))}, ...)`, so `RadialCollar`'s
matrices and vectors actually follow the state's precision, as the PR already claimed.
Nothing is reachable at a different precision today (`TransRadial` has no device path),
so this is forward-looking, for `gpu/21`.

### 6 (minor) — no residency assertion on the new absorber structs

`RateAbsorber`, `LegacyAbsorber`, `RadialCollar` and `CartesianCollar` now call
`assert_resident` on their own buffers/mirrors at construction (a local `_specof` helper
derives the `DeviceSpec` from the reference array they are handed, since none of these
constructors otherwise has a reason to carry one). Self-consistent by construction today,
like `NonlinearRHS.TransModeAvg`'s own assertion; it catches a future mismatched
low-level construction that `JLArrays` would not.

### 7 (minor) — the collars allocated a profile array every accepted step

`RadialCollar` and `CartesianCollar` now hold a reusable `fac` scratch buffer, matching
`RateAbsorber`'s `tfac`, instead of allocating a fresh array every step.

### 8 (nit) — `test_boundaries.jl` did not restore `allowscalar`

Fixed with the same save/restore of the task-local `:ScalarIndexing` key
`test_device.jl` already does.

### 9 (nit) — the tested exit criterion was a tenth of the brief's

Covered by finding 1: `test_metal.jl` now uses `flength=0.1`, the brief's own length.

### Self-caught: `prop_capillary_args`'s docstring was silently orphaned

While filling in this section's numbers, rebuilding the docs (`include("docs/make.jl")`)
turned up a tenth `@ref` failure beyond the 9 pre-existing ones: `Interface.prop_capillary_args`,
plus two "no docstring found" warnings for the same binding. Cause: finding 2's
`_check_responses_device_capable!` helper was inserted between `prop_capillary_args`'s
docstring and the `function prop_capillary_args(...)` it documents, so Julia silently
attached the docstring to `_check_responses_device_capable!` instead -- a plain reordering
bug, not a review finding, introduced while fixing finding 2 and only caught by actually
rebuilding the docs afterwards. Fixed by moving the docstring immediately above `function
prop_capillary_args` again (`_check_responses_device_capable!` now precedes it). Confirmed
fixed: `include("docs/make.jl")` now reports exactly the same 9 pre-existing unresolved
`@ref`s as the base branch, none from this branch. `test_interface.jl` re-run after the
fix: still 317 pass, 0 fail -- a pure reordering of two independent top-level definitions.

### Re-run after the changes

| | |
|---|---|
| regression gate, `LUNA_REGRESSION_BASE=c1ffc410` | 460 pass, 0 fail, every case `0.000e+00` |
| regression gate, `LUNA_REGRESSION_BRANCH=gpu/int-A` | 460 pass, 0 fail, every case `0.000e+00` |
| `test_boundaries.jl` (worktree env) | 182 pass, 0 fail |
| `test_boundaries.jl` (JLArrays env) | 184 pass, 0 fail |
| `test_output.jl` | 117 pass, 0 fail |
| `test_device.jl` (JLArrays env) | 179 pass, 0 fail |
| `test_interface.jl` | 317 pass, 0 fail |
| `test_metal.jl` (M1 Pro, hardware) | 175 pass, 0 fail (10 testsets) |

Metal hardware numbers, corrected (finding 1) -- largest relative difference in `Eω` over
the saves, measured directly (a small standalone script reproducing each testset's own
parameters, not inferred from the pass/fail assertions):

| comparison | `Eω` | energy |
| --- | ---: | ---: |
| exit criterion (100 nJ, 10 cm, He 1 bar), Metal vs CPU `Float32` | 4.4e-6 | 4.2e-6 |
| exit criterion, Metal vs CPU `Float64` | 4.2e-6 | -- |
| visible Kerr (300 µJ, 10 cm, He 5 bar, ~×2.5 broadening), Metal vs CPU `Float32`, fixed steps | 4.3e-6 | 1.5e-6 |
| gradient, brief's parameters (100 nJ, Ar 1 bar→0), Metal vs CPU `Float32`, fixed steps | 3.1e-6 | 1.6e-6 |
| gradient, visible Kerr (300 µJ, Ar 5 bar→0), Metal vs CPU `Float32`, **adaptive** (looser tolerance, see finding 1) | 4.6e-4 | -- |
| `prop_capillary`, constant pressure, Metal vs CPU `Float32` / `Float64` | 4.4e-6 / 4.2e-6 | 4.2e-6 |
| `prop_capillary`, gradient pressure, Metal vs CPU `Float32` / `Float64` | 5.9e-6 / 2.0e-4 | -- |

The `prop_capillary` gradient row's `Float64` column (2.0e-4) is the same Float32-rounding-
perturbs-the-adaptive-controller effect as the low-level adaptive gradient row above it
(4.6e-4 vs `Float32`, within the file's own documented 5e-4 tolerance for that comparison)
-- both are well above the 1e-4 used for every fixed-step comparison in the file, and
`prop_capillary`'s own two hardware assertions for this case use 1e-4 (`Float32` reference)
and a separately documented 3e-4 (`Float64` reference) accordingly. Not a device defect:
every fixed-step comparison in `test_metal.jl`, including this same gradient at the same
energy and pressure, agrees to 1e-4 or better.

