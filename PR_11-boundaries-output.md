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

`prop_capillary`/`prop_gnlse` gain `device` (default `Luna.device_request()`),
`precision` (default `nothing`) and `stats_period` (default `1`) keywords, replacing the
`device=Luna.HostSpec()` `gpu/10` hardcoded at every call site. An untouched call
reproduces today's behaviour (the default resolves to the CPU unless the user has set
something); once a GPU package is loaded and `Luna.settings["device"]` becomes `:auto`,
a plain `prop_capillary(...)` call now actually runs on the GPU for the cases it can.

Multimode and radial propagation (`modes` a collection) and `prop_gnlse` are not device-
or reduced-precision-capable yet: their `setup`/`prop_gnlse_args` methods accept and
validate `device`/`precision` with a new `_cpu_only!` helper, erroring with a message
naming the actual limitation rather than an unrelated `MethodError` or a silent fallback
to the CPU that the caller did not ask for.

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

Largest relative difference in `Eω` over the saves, and in the energy statistic where
computed:

| comparison | `Eω` | energy |
| --- | ---: | ---: |
| mode-averaged Kerr, field-resolved, Metal vs CPU `Float32` (`boundary=:none`) | 2.7e-7 | -- |
| mode-averaged Kerr, envelope, Metal vs CPU `Float32` (`boundary=:none`) | 5.2e-8 | -- |
| Metal vs CPU `Float64`, He at 0.3 bar (`boundary=:none`) | 2.4e-7 | -- |
| field-resolved, `boundary=:rate` + default statistics | 3.9e-6 | 7.7e-7 |
| envelope, `boundary=:rate` + default statistics | 4.4e-6 | 5.3e-7 |
| `prop_capillary`, constant pressure (exit criterion) | 0.0 | 0.0 |
| `prop_capillary`, gradient pressure (exit criterion) | 0.0 | 0.0 |

All well inside the `1e-4` gate `test_metal.jl` asserts. `boundary=:none` reproduces the
`gpu/10` numbers to within a factor of 2 (that branch's own table: 2.7e-7/5.2e-8/2.4e-7 for
the same three comparisons -- unchanged, since `gpu/11` did not touch the mode-averaged
RHS). Adding `boundary=:rate` and the default statistics moves the field-resolved and
envelope cases from ~1e-7 to ~1e-6 -- the collar broadcasts and the extra host round trip
for statistics both add their own rounding, an order of magnitude the plan's own §8
estimate (~1e-4 for Float32 phase accumulation) has ample room for.

The exit-criterion `prop_capillary` comparisons (100 nJ, 1 cm, `plasma=false`,
`raman=false`) come back at *exactly* `0.0`, not merely small: at that pulse energy the
Kerr phase accumulated over 1 cm is far below `Float32`'s precision floor, so the
nonlinear right-hand side rounds to exact zero on both backends and the propagation
reduces to the linear operator's elementwise `exp` -- an operation with no
backend-dependent summation order, hence bit-identical. This is a property of the chosen
exit-criterion parameters, not a general claim; the field-resolved/envelope cases above,
at higher energy (1 µJ) where the nonlinearity is not negligible, show the ~1e-6 that
actually reflects the device path's rounding.

`Luna.set_device(:cpu)` opts out of a globally-registered GPU (exit criterion 3): verified
with `prop_capillary`, `prop_gnlse` and multimode `prop_capillary` all under
`Luna.set_device(:auto)` with Metal registered and functional.

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
