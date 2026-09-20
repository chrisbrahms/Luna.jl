# The plasma response on a device: three scans, one `ifelse`, and rates in `Float32`

Branch `gpu/13-plasma`, base `gpu/12-response-traits` (`7d72431c`). GPU_PLAN.md §2, §4.1,
§4.2, §4.3, §4.9, §6 Group D, §11.

Group D, parallel with `gpu/14-raman` and `gpu/15-chi2`. It makes
`Nonlinear.PlasmaCumtrapz` the first nonlinear response with a device kernel of its own
that is not a single fused broadcast, and makes the two ionisation rates a kernel can
evaluate — which together are the exit criterion: **the physics `prop_capillary` builds
by default for a field-resolved run in a non-Raman gas, Kerr and plasma, now runs on
Metal through the simple interface with `device=:auto`.**

**The regression gate is exactly `0.000e+00` for all 21 cases, in both modes and both
classes, against the base and against `gpu/int-A`.** The plasma cases did not move, but
not for a flattering reason: they carry no plasma (see "The plasma regression cases do
not ionise"). Measured separately, on propagations which do ionise, the branch moves `Eω`
by 0.9e-15 to 4.3e-15 relative.

## Motivation

`PlasmaCumtrapz` was `Columnwise`, and everything about it resisted a device:

- **three serial cumulative integrals.** `Maths.cumtrapz!` is
  `out[i] = out[i-1] + δt(y[i-1] + y[i])/2`: a dependency chain along the time axis,
  which is not a broadcast on any backend.
- **a branch per element** for the ionisation-loss term.
- **a rate function nothing could compile.** `IonRateADK` held nine `Float64` fields,
  `IonRatePPTAccel` a `CSpline` whose knots, values, coefficients and index closure were
  all `Float64`, and whose evaluation throws a `DomainError` built by string
  interpolation. Metal's kernel compiler rejects a `double` which survives optimisation
  and cannot format a message.
- **buffers sized for one column**, where a batched response is handed the whole block.

On a device it therefore went through `Nonlinear.HostResponse` — two host copies and a
host evaluation per right-hand side — which is correct and much slower than the CPU.

## What changed

### `src/Maths.jl`

- **`cumtrapz_scan!(out, y, δt)`**: the trapezoid integral along dimension 1 as one
  prefix scan and one broadcast, `δt*(cumsum(y) - (y + y[1])/2)`. `accumulate!` exists on
  `MtlArray` and `CuArray`; `JLArrays` has no scan, so `test_device.jl` supplies a
  host-backed one, exactly as it does for FFT plans. `Maths.cumtrapz!` is untouched and
  is still what `Stats` and `Ionisation.ionfrac!` use.
- **`spline_eval(c, x0)`**: `CSpline`'s evaluation without the optional bounds check.
  `(c::CSpline)(x0)` is now the check plus this, so there is one body and the two agree
  bit for bit.
- **`UniformIndex{T}`** replaces the closure `make_spline_ifun` returned for a uniform
  axis. The closure captured `xmin`, `xmax` and `N` as `Float64`s, which is exactly the
  stray-`Float64` case: it would have been read inside the kernel. Same arithmetic.
- **`todevice_spline(spec, c)`** and an `Adapt.adapt_structure` rule for `CSpline`. The
  first moves the knots at setup (precision included, which `Adapt` does not do); the
  second is for the kernel adaptor, which has to turn the arrays of a spline a broadcast
  closure captured into device pointers.

### `src/Ionisation.jl`

- **`IonRateADK{T}`** and **`IonRatePPTAccel{ST, T}`** are parametric in the real type of
  their constants. The ADK kernel's one `Float64` literal, `-4/3`, becomes
  `-(4*one(aE)/3)`; every other operand was already an `Int` times a field. The `Float64`
  values are unchanged.
- **`ionrate!(out, ir, E, Eref=1)`**: the array-level evaluation, one broadcast of
  `ratekernel(ir, Eref)` on any array type. The rate is not polynomial in the field, so
  the kernel reconstructs the physical field `Eref*e` rather than carrying a power of
  `E_ref` in a coefficient (GPU_PLAN.md §4.1). `Eref == 1` for every `Float64` run, where
  `1*e` is exact. Both rates' array call operators now go through it, so there is one
  body per rate.
- **`_pptaccel`** is the cached PPT kernel: no error path, and above `Emax` it saturates
  at the table's last value (`min(aE, Emax)` lands exactly on the last knot, where the
  spline returns `y[end]`). The host keeps today's error, raised once per call from a
  `maximum(abs, E)` check rather than once per element — GPU_PLAN.md §4.3's arrangement.
  This is the only place the branch trait decides anything, and it selects an *error
  report*, not a kernel.
- **`device_capable`**, **`device_rate`**, **`resident_arrays`** and
  `Adapt.adapt_structure` rules. A direct `IonRatePPT`, a table which ended up on a
  `FastFinder`, or a user's callable is refused with a message naming the alternatives
  and `device=:cpu`.

### `src/Nonlinear.jl`

`PlasmaCumtrapz` is `Batched()`. `_plasma_block!` is the whole response, for a block of
any shape on any array type: one broadcast for the rate, three `cumtrapz_scan!`s, one
`ifelse` broadcast for the loss term, one to accumulate into the output.

- `rate` and `fraction` (and `Em` for a two-component field) have a singleton
  polarisation axis and broadcast against both components, so no reshaping is needed for
  a block of any rank.
- The struct keeps the rate as given in `ratefunc` — the host, physical-unit object
  `Stats.electrondensity` evaluates — and the converted one in `ratedev`.
- `phase` is gone: `P` holds the phase-modulation term until the last scan and the
  polarisation after it, which is one block-sized buffer fewer than before.
- `coefficients` returns four scalars rather than one, because the response has no single
  polynomial degree: `(Eref, e_ratio*Eref, ionpot/Eref, ρ/(Pref*Eref))`. Each reduces to
  the physical constant at `Eref == Pref == 1`.
- `rescale(p, spec, scaling, Et)` sizes the buffers to the block and converts the rate.
- **CPU threading** (§4.9): above `PLASMA_THREAD_MINLEN` (2^14) elements and more than one
  column, `Threads.@threads :dynamic` over the columns, each task taking one column of
  every buffer. `test_device.jl` asserts a multi-column block equals the same columns
  passed one at a time, as **exact equality**.

`rescale_responses` and `resident_arrays_all` are new: the four-argument `rescale` and
`resident_arrays` mapped over a transform's response collection, through the
tuple-of-tuples a gas mixture is.

### `src/NonlinearRHS.jl`

`TransModal`, `TransRadial`, `TransFree` and `TransFree2D` call `rescale_responses` at
construction with their own block prototype. They did not call `rescale` at all before —
only `TransModeAvg` did — and a batched response cannot work without it, since its
buffers have to match the block, which for a radial or free-space transform is the whole
transverse grid rather than the single column the response's constructor is given. On
the default CPU path every fallback returns the response unchanged. `TransModeAvg` now
uses the same two helpers, which as a side effect gives a gas mixture's inner responses
the treatment they were missing.

### Tests and documentation

`test/test_device.jl` gains "the trapezoid scan", "ionisation rates in the run's
precision", "the plasma response" (all host), "plasma on JLArray" and "Kerr and plasma
propagation on JLArray", plus a host-backed `accumulate!` shim for `JLArray`.
`test/test_metal.jl` gains "plasma on Metal" including the two dynamic-range cases, the
rate structs in the stray-`Float64` smoke test, and a Kerr+plasma propagation in the
`prop_capillary` testset. `test/test_interface.jl`'s refusal test moves from plasma to
Raman, since plasma is device-capable now, and gains the `precision=Float32`-with-plasma
path which now runs.

`docs/src/developer/device_model.md` gains "The plasma response" with the dynamic-range
audit; `docs/src/gpu.md` gains "Ionisation rates on a device" and its "what runs where"
is corrected. `benchmark/device.jl` runs its sweep twice, Kerr and Kerr+plasma, and times
the plasma response alone over a block of columns.

## The plasma regression cases do not ionise

`modeavg_field_plasma`, `modeavg_field_adk`, `modeavg_field_vector` and
`multimode_field_plasma` are helium at 1 bar, 800 nJ in 10 fs through a 125 µm core:
about 3e11 W/cm², where the PPT rate is ~1e-60 1/s. The recorded `electrondensity`
statistic is **exactly** `0.0` at every step, and running the same case with
`plasma=false` gives a **bit-identical** `Eω`. The plasma response is called at every
right-hand side and contributes nothing that `Float64` can represent.

So the gate's `0.000e+00` on those four cases is real but says nothing about this branch.
The number that does is this, measured by running the same propagation on this branch and
on the base commit in the same environment, 20 fixed steps over 2 cm, 125 µm core, 10 fs:

| case | peak electron density [m⁻³] | `Eω`, global | `Eω`, per save |
| --- | ---: | ---: | ---: |
| Ar 1 bar, 300 µJ, PPT | 9.263792418965712e22 | 2.56e-15 | 3.03e-15 |
| Ar 1 bar, 300 µJ, ADK | 3.962828910370468e22 | 1.31e-15 | 1.50e-15 |
| Ar 1 bar, 300 µJ, PPT, elliptical (vector plasma) | 7.578097369373387e23 | 2.81e-15 | — |
| He 1 bar, 800 µJ, PPT | 2.874216387522913e21 | 7.68e-16 | 7.74e-16 |

Measured with **`set_fftw_wisdom(false)`**, as the gate is (see "review round 1, finding
11" below: the first version of this table was measured with the shared wisdom file
enabled, which changes the plans and moves these numbers by up to 1.7x).

0.1 to 3 per cent of the gas is ionised in those, and the difference is the scan's
summation order, three orders of magnitude inside the gate's 1e-12 tolerance. The
`electrondensity` statistic is bit-identical in the three scalar cases, because `Stats`
computes it with `Maths.cumtrapz!`, which this branch does not touch.

**Recommendation for a later branch** (not done here, since the gate matrix is
`gpu/00-harness`'s): raise the energy of the plasma cases, or switch them to argon, so
that they exercise the response they are named after.

## What the Metal tests found

The vector plasma response produced NaN on Metal and nowhere else. The ionisation-loss
term for a two-component field divides by `Em²`, and the guard was on `Em`. In the wings
of a pulse `Em` is 1e-22 of its peak (`exp(-50)`), so `Em²` is 1e-44 — subnormal in
`Float32`, which a device flushes to zero while `Em` itself is still a normal number. The
rate there is zero too, so the term became `0/0`, and one NaN travels through the three
scans which follow and destroys the whole column.

On the CPU, in `Float32` as much as in `Float64`, the subnormal survives and `0/x` is
`0`, so neither the `JLArray` tests nor the `Float32` host tests could see it. This is
GPU_PLAN.md §4.2 rule 7 — "every branch which adds or changes a kernel runs the Metal
hardware tests ... the only reliable detector" — earning its place.

The guard is now on the denominator, `Em² > 0`. In `Float64` the two conditions differ
only below 1e-162 V/m, and the regression gate is unchanged at `0.000e+00`, including
`modeavg_field_vector` and `multimode_field_plasma`.

## Deviations from GPU_PLAN.md and the brief

1. **No offset on the stored PPT table, and no format-version bump.** The brief asks for
   the table "held in `log` form with an offset so that `Float32` covers the range", and
   for a bump of the cache key's format version if the contents or format change. The
   table is already in `log` form, and the range that puts in it (−708 to +37) is covered
   by `Float32` with a resolution of 6.1e-5 at the bottom and 3.8e-6 at the top — the
   first at a rate of 1e-308 1/s, which nothing uses. An offset would also change the
   `Float64` path, which is otherwise bit-identical. **Nothing about the cached file
   changed, so the cache key is untouched** and every other worktree's cache stays valid.
   The measurements are in `device_model.md`.
2. **The four other transforms call `rescale`.** The brief says not to touch "transforms
   other than what the batched call needs". A batched response with block-sized buffers
   needs exactly this, and without it a radial or free-space plasma run — which works
   today — would break. The change is one line each and is a no-op on the host path.
3. **`rescale_responses`/`resident_arrays_all` are new in `Nonlinear.jl`.** Group C owns
   that file's contracts. These add nothing to the protocol and change no existing
   contract: they are the four-argument `rescale` and `resident_arrays` mapped over a
   collection, which every transform now needs. Flagged rather than negotiated because
   the alternative was four copies of the same recursion in `NonlinearRHS.jl`.
4. **`device_capable` is not overridden for `PlasmaCumtrapz`.** It is derived from `kind`
   (`gpu/12`), so a plasma response whose *rate* has no kernel — a direct `IonRatePPT`, or
   a user's callable — still reports `true`, and `prop_capillary` lets the call through to
   fail in `Nonlinear.rescale` with `Ionisation.device_rate`'s message instead of at the
   interface check. Overriding it would reintroduce the per-type declaration `gpu/12`
   removed. The message names the fix; the default PPT and ADK rates are both capable, so
   this is an edge case. Worth a decision in `gpu/int-D`.
5. **`PlasmaScalar!`/`PlasmaVector!` are gone**, replaced by `_plasma_block!`. They were
   internal helpers taking the response and one column; nothing outside `Nonlinear.jl`
   used them. `PlasmaCumtrapz`'s `phase` field is gone with them.
6. **Attribution trailer.** The commits carry `COMMON.md`'s `Claude Fable 5.1` line, as
   `gpu/12`'s round-2 commits do.

## Measurements

Apple M1 Pro, Julia 1.13.0, `set_fftw_mode(:estimate)`, one FFTW thread, one BLAS thread,
with other agents' jobs on the same machine.

### The regression gate

Baseline generated from `7d72431c` (the base commit) with `test/regression/generate.jl`,
and against `gpu/int-A` (`782f55d1`) for the cumulative view.

**460 pass, 0 fail against both. Every case, both modes, both classes: `0.000e+00`.**
Step counts unchanged.

### `test/test_device.jl`

`julia --project=<env with JLArrays> -t 1`, and again with `-t 4`.

| testset | assertions |
| --- | ---: |
| **the trapezoid scan** (new) | 7 |
| **ionisation rates in the run's precision** (new) | 35 |
| **the plasma response** (new) | 107 |
| **plasma on JLArray** (new) | 21 |
| **Kerr and plasma propagation on JLArray** (new) | 8 |
| the 18 testsets of the base branch | 263 |

**441 pass, 0 fail** (23 testsets) with one thread and with four, up from 263 on the base
branch. The threaded path is only taken with more than one thread; the same assertions
cover both, and the multi-column block is over `PLASMA_THREAD_MINLEN` so that the
four-thread run really is the threaded path against the serial one.

| comparison | max relative difference |
| --- | ---: |
| batched plasma vs the serial implementation it replaces (ADK and tabulated, scalar and vector) | < 1e-11 |
| `JLArray` vs host, plasma response, ADK and tabulated, scalar and vector | < 1e-10 |
| `JLArray` vs host, Kerr + plasma propagation, fixed steps | < 1e-10 |
| threaded columns vs one column at a time | **0** (exact) |
| scaled vs unscaled plasma, `Float64`, `Eref = 2^36` | < 1e-14 |
| `Float32` vs `Float64`, scaled plasma response | < 1e-4 |

### `test/test_metal.jl`

Apple M1 Pro, Metal.jl v1.11.1, in a separate environment as the file's header describes.

| testset | assertions |
| --- | ---: |
| no stray Float64 in the kernels (the two rate structs added) | 33 |
| **plasma on Metal** (new, including the two dynamic-range cases) | 35 |
| prop_capillary on Metal (Kerr+plasma added) | 45 |
| the other ten testsets | 147 |

**260 pass, 0 fail** (13 testsets), up from 194 on the base branch.

Measured separately, by a script reproducing the same parameters:

| | Metal vs CPU `Float32` | Metal vs CPU `Float64` | CPU `Float32` vs `Float64` |
| --- | ---: | ---: | ---: |
| plasma response, ADK, scalar | 5.4e-05 | 3.3e-05 | 2.1e-05 |
| plasma response, ADK, vector | 1.1e-05 | 1.4e-05 | 2.4e-05 |
| plasma response, tabulated, scalar | 9.5e-06 | 2.5e-05 | 1.6e-05 |
| plasma response, tabulated, vector | 3.3e-05 | 4.6e-05 | 1.8e-05 |
| Ar at the barrier-suppression field (rate 1.3e14 1/s, 18 % ionised) | 1.2e-05 | 4.3e-06 | 1.6e-05 |
| Kerr + plasma propagation, Ar 1 bar 150 µJ 1 cm, 20 fixed steps | **3.1e-07** | **7.8e-07** | 8.5e-07 |

The response-level differences are larger than the pure-Kerr ones `gpu/12` measured
(which were exactly zero, a single elementwise broadcast doing identical IEEE operations
in identical order): the scans are a parallel prefix sum on the device and a sequential
one on the host, so the summation order differs. They are the same size as the
`Float32`/`Float64` gap in the last column, i.e. the device adds nothing beyond what
single precision already costs. The propagation — which is what a user sees — agrees with
the `Float32` CPU run to 3e-7.

The helium dynamic-range case (0.3 bar, 1e9 V/m, far below threshold) gives a peak rate
of exactly zero and an output of exact zeros on Metal, finite throughout: the underflow
produces nothing, which is the point of testing it.

### Existing CPU test files

All CPU, all pass: `test_maths.jl` 138, `test_ionisation.jl` 25, `test_vectorplasma.jl`
2, `test_polarisation_field.jl` 8, `test_interface.jl` 327 (11 testsets, up from 323),
`test_multimode.jl` 6, `test_freespace.jl` 77.

### Benchmarks

`benchmark/device.jl`, M1 Pro, one FFTW and one BLAS thread, `:estimate`, no wisdom. The
sweep now runs twice, Kerr and Kerr+plasma. One right-hand side, one step, and a 20-step
propagation:

| device | trange | state | rhs (Kerr) | rhs (Kerr+plasma) | prop (Kerr) | prop (Kerr+plasma) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CPU `Float64` | 400 fs | 1025 | 17.1 µs | 63.1 µs | 4.55 ms | 11.9 ms |
| CPU `Float32` | 400 fs | 1025 | 11.2 µs | 54.2 µs | 3.83 ms | 10.8 ms |
| Metal `Float32` | 400 fs | 1025 | 289 µs | 701 µs | 52.3 ms | 89.1 ms |
| CPU `Float64` | 6400 fs | 16385 | 538 µs | 1.29 ms | 103 ms | 195 ms |
| CPU `Float32` | 6400 fs | 16385 | 266 µs | 957 µs | 69.3 ms | 154 ms |
| Metal `Float32` | 6400 fs | 16385 | 318 µs | 630 µs | 63.3 ms | 106 ms |

The plasma response roughly triples the cost of a mode-averaged right-hand side on the
CPU, which is the three scans and the rate. A single-column run is launch-bound on the
GPU, as it was for Kerr; Metal only wins at the largest grid, and by less than a factor
of two.

Where it does win is columns. `batched!` alone, on a 2048-sample block:

| columns | CPU `Float64`, 1 thread | CPU `Float64`, 8 threads | Metal `Float32` |
| ---: | ---: | ---: | ---: |
| 1 | 36.0 µs | 35.8 µs | 587 µs |
| 16 | 538 µs | 94.6 µs | 653 µs |
| 128 | 4.51 ms | 580 µs | 620 µs |

The column threading is close to linear (7.8x on eight threads at 128 columns) and the
GPU is flat in the number of columns, which is what the batched trait is for. Neither
matters for a mode-averaged run, which has one column; both will for the radial and
free-space transforms in Group E.

### Documentation build

`include("docs/make.jl")` reports exactly the pre-existing unresolved `@ref`s the base
branch does (`LinearOps.βz` x3, `loadFFTwisdom`, `saveFFTwisdom`, `AbstractOutput`,
`crystal_internal_angle`, `norm_free`, `LinearOps.make_const_linop`) and none from this
branch. There is no `docs/src/modules/Maths.md`, so the new `Maths` symbols are referred
to in code spans rather than with `@ref`.

## Known gaps

- **The plasma cases of the regression matrix do not ionise** (above). The branch is
  gated by them and by the four propagations in the table, but the matrix itself should
  be fixed.
- **`Stats.electrondensity` still integrates with `Maths.cumtrapz!` on the host**, so on
  a device run the per-step statistics still copy the field down and redo the rate on the
  host. That is `gpu/24`'s (device statistics), and it is why `ratefunc` is kept.
- **A gas mixture with plasma on a device** is not tested. `rescale_responses` now
  reaches the inner responses of a mixture, which it did not before, but the mixture's
  own `Et_to_Pt!` path on a device has no test.
- **Memory.** All four buffers are now sized for the block where they used to be sized
  for one column: `J` and `P` at the full block size, and `rate` and `fraction` at the
  block size with a singleton polarisation axis (so the full block for a scalar field,
  half of it for a two-component one). For a `(nto, 1, Nx, Ny)` free-space run that is
  four more block-sized `Float64` arrays than before, not two — the figure the Group E
  chunking decision should be taken against. Chunking the host path would bound it; it
  is not done here. (The count *per column* is one lower than before, since `phase` is
  gone.)
- **`preionfrac`** is carried through unchanged and is exercised only by
  `test_ionisation.jl`'s existing case.


## Changes after review round 1

The review's verdict was "request changes", on one blocker and two majors. All eleven
findings are addressed. **The regression gate is still exactly `0.000e+00` on all 21
cases**, in both modes and both classes, against `7d72431c`.

### 1 (blocker) — `ionrate!` dispatched on the type, not on capability

`ionrate!(out, ir::AbstractIonRate, E, Eref)` called `ratekernel`, which only
`IonRateADK` and `IonRatePPTAccel` have. Every other `AbstractIonRate` — including
Luna's own exported, documented `Ionisation.IonRatePPT`, and anything a user writes —
died with a `MethodError` inside `PlasmaCumtrapz` **on the plain host `Float64` path**,
which is a behaviour change with default settings and the exact opposite of what this
branch's own user page promises. Nothing exercised it.

`ionrate!` now dispatches on `device_capable(ir)`, as the review sketched: a rate with a
kernel takes the broadcast, anything else is called as `ir(out, E)` when `Eref == 1` and
`E` is a host array and refused, by name, otherwise. `IonRateADK`'s and
`IonRatePPTAccel`'s array call operators are back to `out .= ir.(E)` rather than
delegating to `ionrate!`, because the fallback calls them and a delegating one would
recurse — which is exactly what a cached rate on a non-uniform table would have done.

Tested in `test_device.jl`, "a rate with no device kernel": a direct `IonRatePPT`, a
`UserRate <: AbstractIonRate` written the way the documentation asks, and an
`IonRatePPTAccel` on a non-uniform table, each through `ionrate!` and inside a
`PlasmaCumtrapz` against the serial reference, plus the named refusal for a scaled run
and for `device_rate`.

### 2 (major) — a non-`Tuple` collection with a batched response

Fixed at both ends, as the review suggested:

- `TransModal`, `TransRadial`, `TransFree` and `TransFree2D` now `Tuple(...)` their
  responses at the `rescale_responses` call, as `TransModeAvg` already did;
- `Et_to_Pt!`'s legacy path refuses a `Batched` response (`_refuse_batched_legacy`,
  which also walks a mixture's inner tuples) with a message that says a batched response
  needs a `Tuple` collection, instead of letting it fail later as a shape mismatch. The
  response's own size-mismatch error also now names the tuple requirement.

Tested with a `Vector` of responses containing the plasma response through
`TransRadial`: the transform converts, the buffers come out block-sized, and the
polarisation is identical to the same responses passed as a tuple; the bare `Et_to_Pt!`
with a `Vector` and `idcs` raises the new message.

### 3 (major) — the interface refusal test pointed at Raman

`gpu/14-raman` makes Raman device-capable, after which no response `prop_capillary` can
build is columnwise and there is nothing left for a call-level test to point at. The test
now calls `Interface._check_responses_device_capable!` directly with a tuple containing a
user closure — columnwise by definition and staying that way — for `precision` alone,
`device` alone and both, asserting the message names `device=:cpu`, plus the two
not-refused cases. **`gpu/int-D` note:** there is no longer a `prop_capillary`-level
refusal test, because there is no longer a `prop_capillary` call which should be refused.

### 4 (minor) — the four `rescale_responses` call sites hard-coded the host

`TransModal`, `TransRadial`, `TransFree` and `TransFree2D` take `spec=HostSpec()` and
`scaling=UNIT_SCALING` keyword arguments, used at that call, with a comment saying Group
E changes the caller and not the line. The rest of each transform is still host-only and
`Luna.setup` still refuses a device for them.

### 5 (minor) — the memory note undercounted

Corrected to four block-sized buffers (`J` and `P` at full size, `rate` and `fraction`
with a singleton polarisation axis), from one per column before.

### 6 (minor) — the `Emax` error came out of a `@threads` loop

`_check_field_range` is now `Ionisation.check_field_range`, documented and public, and
`_plasma_run!` calls it once on the whole block — together with the field magnitude for a
two-component field, which it now also computes once for the block — before the threading
decision. `ionrate!` takes `check=false` for a caller which has already made the check.
The error keeps its type in the threaded case.

### 7 (minor) — the `Float32` limitation was not on the user page

`docs/src/gpu.md`, "Ionisation rates on a device", now says that a weak plasma can vanish
in single precision, names helium at ~1e14 W/cm² and the 2.6e-7 ratio, and says to use
`Float64` if a contribution that small matters.

### 8 (minor) — the Metal `prop_capillary` plasma test did not ionise

It now runs argon at 0.1 bar and 300 µJ through `prop_capillary_args` with fixed steps
(which `prop_capillary` has no keyword for), and asserts that the electron density is
more than 0.1 % of the gas — so the real cached PPT spline kernel is exercised at a rate
that is not zero — as well as that the device response holds an `IonRatePPTAccel` with
its knots on the GPU.

### 9-10 (nits)

`PLASMA_THREAD_MINLEN` has a docstring; `cumtrapz_scan!` checks `out !== y` rather than
only documenting it.

### 11 (nit) — the Ar PPT row did not reproduce

It was measured with the shared FFTW wisdom file enabled, which changes the plans and so
the rounding. With `set_fftw_wisdom(false)` — the gate's configuration — the row is
**2.56e-15 global, 3.03e-15 per save**, which is the review's 2.6e-15 / 3.0e-15 exactly;
the elliptical row is 2.81e-15 either way, which is why only one row disagreed. The table
above is re-measured with wisdom disabled and says so. With wisdom enabled the same
script reproduces the original 4.35e-15 bit for bit across runs, so neither number was
noise in the measurement — they are two different sets of FFT plans.

### Re-run after the changes

| | |
|---|---|
| regression gate, `LUNA_REGRESSION_BASE=7d72431c` | **460 pass, 0 fail**, every case `0.000e+00` (1m36) |
| `test_device.jl`, `-t 1` and `-t 4` | **471 pass, 0 fail** (25 testsets), up from 441 |
| `test_metal.jl` (M1 Pro, hardware) | **261 pass, 0 fail** (13 testsets) |
| `test_ionisation.jl` | **25 pass, 0 fail** |
| `test_interface.jl` | **329 pass, 0 fail** (11 testsets) |
