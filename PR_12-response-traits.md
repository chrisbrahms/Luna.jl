# The response protocol: kinds, one fused broadcast, and the host fallback

Branch `gpu/12-response-traits`, base `gpu/11-boundaries-output` (`0cff74b6`). GPU_PLAN.md
§4.1, §4.2, §4.3, §4.11, §6 Groups C and D, §11.

This is Group C: one small branch whose only job is to fix the contracts the three Group D
branches (plasma, Raman, χ⁽²⁾) will write responses against. It gives every nonlinear
response a **kind**, makes `NonlinearRHS.Et_to_Pt!` dispatch on it, fuses all the pointwise
responses of a transform into one broadcast per right-hand side, and adds
`Nonlinear.HostResponse`, the fallback which lets a response Luna has not made
device-capable — including one the user wrote as a closure — run on a GPU at all.

**The default CPU path is unchanged: the regression gate is exactly `0.000e+00` for all 21
cases, in both modes and both classes, against both `0cff74b6` and `gpu/int-A`.** No case
moved, including the two the fused broadcast could have moved (`modeavg_field_mixture`,
two pointwise terms; `modeavg_field_vector`, the vector forms). See "Why the default CPU
path did not move".

## Motivation

Two things stood in the way of Group D.

**There was no protocol.** `Et_to_Pt!` called every response the same way — once per
column, on the host, in physical SI units — and `gpu/10-device-model` made the Kerr
responses work on a device by making their call operators generic, which works for one
column but leaves the per-column loop, and no way for a response to say "call me with the
whole block" or "I am elementwise, fuse me". Each of the three Group D branches would have
invented its own answer, and they cannot be reconciled afterwards without diverging three
branches.

**There was no fallback.** `Nonlinear.rescale` errored for any response which was not one
of the three Kerr structs, so a device run was possible only with Kerr. GPU_PLAN.md §3
requires that "ad hoc responses and linear operators still work, falling back to the CPU on
a GPU run"; nothing implemented that.

## What changed

### `src/Nonlinear.jl` — the protocol

`Nonlinear.kind(r)` returns one of four singletons, `Columnwise()` by default:

| kind | what the dispatcher does | what the response supplies |
| --- | --- | --- |
| `Pointwise()` | fuses it into one broadcast over the whole block with the other pointwise responses, no buffer | `pointwise_kernel` (a `T -> T`), or `pointwise_expr` if it carries per-sample arrays |
| `VectorPointwise()` | two broadcasts, one per polarisation component, each fused across the group | `vector_kernel` (an `(ex, ey) -> SVector{2}`), or `vector_expr` |
| `Batched()` | calls it once with the whole `(nt, npol, ncols)` block | the call operator, plus its own full-size buffers in the run's array type |
| `Columnwise()` | calls it once per column, on the host; refused on a device array | nothing — the default, so any callable `(out, E, ρ)` works |

`kind(r, npol)` (with `npol` a `Val`) is the kind for a block with that many polarisation
components, defaulting to `kind(r)`. That is how `KerrField` is `Pointwise` for a scalar
field and `VectorPointwise` for a two-component one. The dispatcher resolves `npol` to a
compile-time constant before it groups the responses, so a response is free to report
different kinds for the two cases.

Three functions carry the units and the precision:

- **`coefficients(r, ρ, scaling)`** is now the one place a response's physical constants,
  the density and the powers of `E_ref`/`P_ref` meet. It runs on the host, in `Float64`.
  `Luna.polscale(scaling, n)` (new, in `Device.jl`) gives the factor a response of
  polynomial degree `n` in the field carries, `Eref^(n-1)/Pref`: cubic for the Kerr
  responses, quadratic for a χ⁽²⁾ one. It is exactly `1` for every `Float64` run.
- **the kernel** converts that scalar once, with `Luna.scalar(E, c)`, so nothing reachable
  from a device kernel is a `Float64`.
- **`rescale(r, spec, scaling)`** now converts only the *arrays* a response carries, and
  wraps a columnwise response for a device run. A response whose coefficients are all
  scalars needs no method: the fallback passes it through when `resident_arrays(r)` is
  empty, and errors only when there is an array and no rule for moving it.

The consequence is that **a response struct never enters a kernel** — only the scalars its
kernel captured and the arrays `rescale` converted. `KerrField` therefore keeps its
physical `Float64` `γ3` on a Metal run, and `test_metal.jl` checks the element type of the
*kernel's* coefficient rather than of the struct's field.

`device_capable(r)` is now derived (`!(kind(r) isa Columnwise)`) rather than declared per
type, so the two sets cannot drift apart. Its meaning shifts slightly: a columnwise
response *can* now run on a device, through the fallback, so `device_capable` is about
whether it has a kernel of its own, i.e. about speed.

### `src/Nonlinear.jl` — `HostResponse`

`rescale`'s fallback wraps a `Columnwise` response in a `Nonlinear.HostResponse` as soon as
the run is not host `Float64` in physical units, and logs one `@info` line naming the
response. The wrapper is `Batched`, so it receives the whole block, and per right-hand side
it

1. copies the block to a host buffer in the run's element type and multiplies by `E_ref`
   into a `Float64`/`ComplexF64` buffer — physical units, double precision, which is what
   the response was written for;
2. zeroes its own polarisation buffer and calls the response on it column by column,
   exactly as `Et_to_Pt!`'s `idcs` loop does;
3. multiplies by `1/(P_ref E_ref)` into the staging buffer, copies that back up and adds it
   to the output.

Buffers are allocated on the first call, because `rescale` sees the spec and the scaling
but not the block shape, and are read back through a function barrier so the per-element
code is still compiled for concrete types. The separate staging buffer in the run's element
type exists because `copyto!` between a host array and a device array does not convert the
precision (`copyto!(::MtlArray{Float32}, ::Array{Float64})` is not a conversion).

It is correct and slow, deliberately: two host copies and a host evaluation per right-hand
side, which on a GPU also serialises the step.

### `src/NonlinearRHS.jl` — the dispatch

`Et_to_Pt!` gains a `scaling` keyword (default `UNIT_SCALING`; only `TransModeAvg`, the
only scaled transform, passes anything else) and is rewritten around a list of
`(response, density)` pairs. Pairing each response with the density it sees turns the flat
case and the gas-mixture case (a tuple of tuples with a vector of densities) into the same
list, which is how the mixture's per-gas Kerr terms end up summed inside one broadcast
rather than applied one gas at a time.

The dispatcher then walks the list in order, takes the longest run of pointwise responses
at a time and materialises it as one broadcast, and applies everything else one response at
a time. The pointwise expressions are built as nested `Base.broadcasted` objects and summed
left-associated in tuple order, which Julia fuses into a single kernel on every backend.
The first group written *assigns* into `Pt` rather than zero-filling it and accumulating,
which is one pass over the block fewer.

A `responses` collection which is not a tuple keeps the historical per-response loop, and
is refused in a scaled run (it would apply physical-unit coefficients to a scaled state;
nothing in Luna reaches it, since `TransModeAvg` always holds a tuple, but a low-level
caller could).

`Et_to_Pt!(Pt, Et, resp, ρ, idcs)` now applies pointwise and batched responses to the whole
`(nt, npol, ncols)` block and loops over `idcs` only for columnwise ones. `idcs` covers
every column in all three transforms which use it, so hoisting the zero-fill out of the
loop is equivalent.

### `src/Interface.jl`

The simple interface deliberately does not use the fallback. An explicit `device` or
`precision` request to `prop_capillary` whose mode-averaged responses are not all
`device_capable` errors naming `device=:cpu`, rather than silently producing a run slower
than the CPU run the caller almost certainly wanted.

`precision=Float32` on its own now goes through that check too. That is the leftover nit
from review round 2 of `gpu/11-boundaries-output`: with `device` unspecified and plasma on,
`precision=Float32` used to fail deep inside `Nonlinear.rescale` with a message that did not
name the fix. It matters more now, because without the check it would no longer fail at all
— it would run, slowly, through `HostResponse`. A call which mentions neither keyword is
unaffected and still resolves exactly as it did.

The error message and the `device`/`precision` docstrings are reworded: "a response with no
`Nonlinear.rescale` method" is no longer accurate, and the message now says what the
low-level interface will do instead.

### Documentation

`docs/src/developer/device_model.md` gains "The response protocol" (the table filled in,
the three functions which carry the units, a worked example of a new pointwise response,
and the deviation below) and "The host fallback". `docs/src/gpu.md`'s "What runs where" is
corrected — responses are no longer refused on a device, transforms still are — and gains
"An ad hoc response on a device", which says plainly that the fallback is correct and slow
and that the simple interface refuses rather than using it.

## Deviations from GPU_PLAN.md and the brief

1. **The vector-pointwise form is two broadcasts, not one `reinterpret`ed `SVector{2}`
   write.** GPU_PLAN.md §4.3 (and the brief, item 1) describe the vector broadcast as
   returning an `SVector{2}` "written into a `reinterpret`ed `(2, nt, ncols)` view of
   `out`". Luna's buffers are `(nt, npol, ncols)`: the polarisation index is the **slow**
   axis, so there is no contiguous leading axis of length 2 to reinterpret, and producing
   one would mean a transpose and a field-sized buffer — exactly what "pointwise, no
   buffer" is supposed to avoid. The *kernel* contract is kept as written (the response
   returns an `SVector{2}` of the two lab-frame components, which is what Group D's χ⁽²⁾
   branch needs); the dispatcher materialises component 1 and component 2 into the two
   output views with one broadcast each, both fused across the whole pointwise group. For
   the two-component Kerr forms this is the same arithmetic and the same number of passes
   as the pair of broadcasts they were already written as. **This is the one contract Group
   D has to know about that differs from the plan's wording.**
2. **`coefficients(r, ρ, scaling)` takes the scaling, and `rescale` no longer folds it in.**
   The brief says `coefficients` should be "the one place where physical constants and
   `E_ref`/`P_ref` powers are combined ... (gpu/10's `rescale` may be renamed/folded into
   this; keep one path)". Taken literally: `rescale` now converts arrays only, the response
   keeps its physical constants, and `coefficients` combines everything per call. The cost
   is that `Et_to_Pt!` needs the scaling, which it takes as a keyword from `TransModeAvg`;
   no other transform passes one. The benefit beyond tidiness is that a `Float32` run never
   has a small `Float32` intermediate: `γ3·Eref²/Pref` used to be converted to `Float32`
   and then multiplied by `ρ·ε₀` in `Float64`, whereas now the whole product is formed in
   `Float64` and converted once.
3. **"Make the mixture Kerr a pointwise response summing per-gas coefficients" is
   implemented as one broadcast whose body sums the per-gas terms**, each with its own
   density and its own coefficient — not as a single response object holding the summed
   coefficient. Summing the coefficients first (`(ρ₁γ₁ + ρ₂γ₂)E³` instead of
   `ρ₁γ₁E³ + ρ₂γ₂E³`) would change the arithmetic and move `modeavg_field_mixture` at
   rounding level, and would only work when every gas has the same response type. The
   effect asked for — the mixture fuses like everything else, one broadcast for the whole
   mixture — is what happens.
4. **Commits.** The brief suggests five (traits + dispatch; `HostResponse`; Kerr/mixture as
   pointwise; tests; docs + PR). There are four: the traits, the dispatch and the Kerr
   kernels are one commit, because the dispatch is untestable without at least one
   non-columnwise response and the Kerr kernels are where `coefficients` is defined.
   Splitting them further would have meant an intermediate commit which compiles but whose
   new code is never reached.
5. **Attribution trailer.** The commits carry `Co-Authored-By: Claude Opus 5 (1M context)`,
   which is the model that wrote them and what this session's harness specifies, not the
   `Claude Fable 5.1` line in `briefs/COMMON.md`. Everything else about the trailers,
   including the session URL, is as COMMON.md says.

## Why the default CPU path did not move

- **The fused broadcast is the same arithmetic.** Terms are summed left-associated in tuple
  order, which is the order the `for resp! in responses` loop accumulated them in. The only
  difference is the leading `0 +` the zero-fill used to contribute, which is exact for
  every value; the one case where it is observable is the *sign* of an exact zero
  (`0.0 + (-0.0)` is `+0.0`, the assignment keeps `-0.0`), and a difference between `+0.0`
  and `-0.0` has magnitude zero, which is what the gate measures. In practice the gate is
  `0.000e+00`, so even that did not arise.
- **`polscale` is exactly `1` on the `Float64` path** (`Eref == Pref == 1`), and
  multiplying by `1.0` and dividing by `1.0` are exact, so moving the scaling out of
  `rescale` and into `coefficients` changes no value.
- **The Kerr kernel bodies are the expressions they were**, in the same association:
  `fac*(ex^2 + ey^2)` then `*ex`, `3/4` folded into the scalar factor, and so on. The
  columnwise call operators (`KerrScalar!`, `KerrVector!`, `KerrScalarEnv!`,
  `KerrVectorEnv!`, still present and still public) now call the same kernel functions, so
  the two paths cannot drift apart — one kernel body per response, as GPU_PLAN.md §4.2
  rule 1 requires.
- **Nothing on the default path constructs a `HostResponse`**: `rescale` returns the
  response unchanged for host `Float64` in physical units.
- **The `idcs` loop's zero-fill was hoisted**, which is equivalent because `idcs` covers
  every column of `Pt` in all three transforms which pass one.

## Tests

`julia --project=<worktree> -t 1`, `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)`, `BLAS.set_num_threads(1)`. Apple M1 Pro, Julia 1.13.0, with
other agents' jobs on the same machine.

### The regression gate

Baseline generated from the base commit with `test/regression/generate.jl 0cff74b6`, and
again against `gpu/int-A` (`782f55d1`, the merge base) for the cumulative view:

**460 pass, 0 fail against both baselines. Every case, both modes, both classes:
`0.000e+00`.**

| case | `:fixed` Eω | `:fixed` stats | `:adaptive` Eω | `:adaptive` stats |
|---|---:|---:|---:|---:|
| every one of the 21 cases | 0 | 0 | 0 | 0 |

**Largest difference over all cases and modes: `0.000e+00`,** against `0cff74b6` and
against `782f55d1` (`gpu/int-A`). Step counts unchanged (the gate checks them separately
and fails hard if they move).

The two cases which could have moved are the ones to look at:

| case | what this branch changed about it | `:fixed` Eω | `:adaptive` Eω |
|---|---|---:|---:|
| `modeavg_field_mixture` | two pointwise Kerr terms, one per gas, now summed inside one broadcast instead of applied one gas at a time | 0 | 0 |
| `modeavg_field_vector` | the vector Kerr forms now go through `vector_kernel` and the component broadcasts of the fused group | 0 | 0 |
| `modeavg_env_thg` | `KerrEnvTHG` now supplies a `pointwise_expr` with `C` as a broadcast argument | 0 | 0 |
| `free2d_field_chi2`, `free2d_env_chi2` | two-component columnwise responses through the new `Val(2)` dispatch path | 0 | 0 |
| `modeavg_field_plasma`, `modeavg_field_raman`, `modeavg_field_nothg` | Kerr (fused, assigning) followed by a columnwise response (accumulating) | 0 | 0 |
| `radial_*`, `free3d_env_kerr`, `multimode_field_plasma` | the `idcs` path, with the zero-fill hoisted out of the column loop | 0 | 0 |

### `test/test_device.jl`

Run twice: in the worktree environment (the `JLArrays` half skips itself) and in an
environment with `JLArrays` added, which is what `Pkg.test()` gives it.

| testset | assertions |
| --- | ---: |
| backend trait | 16 |
| device spec and settings | 23 |
| allocation and transfer | 18 |
| residency assertions | 6 |
| unit scaling | 8 |
| **response kinds and the fused broadcast** (new) | 29 |
| FFT planner dispatch | 9 |
| constβ is checked, not trusted | 2 |
| JLArray basics | 14 |
| RK45 kernels on JLArray | 11 |
| mode-averaged Kerr on JLArray | 24 |
| pressure gradient on JLArray | 10 |
| **pointwise responses on JLArray** (new) | 14 |
| **a user closure response through HostResponse on JLArray** (new) | 9 |
| boundaries and default statistics on JLArray | 12 |
| stats_period skips the device-to-host copy | 10 |
| a device run refuses host-only machinery (extended) | 10 |
| Float32 on the CPU (updated) | 15 |

**240 pass, 0 fail** (18 testsets, 97 s), up from 179 on the base branch. Without
`JLArrays`: 126 pass, 0 fail.

Every comparison in "response kinds and the fused broadcast" is **exact equality** with
the per-response loop, not a tolerance: the fused broadcast, the vector forms, the
carrier-array response, the gas mixture, a columnwise response before and after a
pointwise one, and the multi-column `idcs` path.

| comparison | max relative difference in `Eω` |
| --- | ---: |
| `JLArray` vs host, Kerr + a user closure through `HostResponse` | 0 |
| `JLArray` vs host, pointwise `Et_to_Pt!` (Kerr field/env/THG, scalar and vector) | 0 |
| scaled vs unscaled `Et_to_Pt!` (`Eref = 1024`, `Pref = ε₀`), two responses of different degree | 1.1e-16 |

### `test/test_metal.jl`

Apple M1 Pro, Metal.jl v1.11.1, environment built per `docs/src/gpu.md`.

| testset | assertions |
| --- | ---: |
| Metal registration | 8 |
| allocation and transfer on Metal | 7 |
| no stray Float64 in the kernels (Kerr part rewritten to go through `Et_to_Pt!`) | 22 |
| mode-averaged Kerr on Metal | 24 |
| Metal against the Float64 CPU path | 6 |
| pressure gradient on Metal | 9 |
| boundaries and default statistics on Metal | 68 |
| **vector pointwise responses on Metal** (new) | 2 |
| **a user closure response through HostResponse on Metal** (new) | 9 |
| prop_capillary on Metal | 25 |
| `Luna.set_device(:cpu)` opts out | 8 |
| Metal refuses what it cannot run (updated) | 3 |

**191 pass, 0 fail** (12 testsets), up from 175 on the base branch.

Measured separately (a standalone script reproducing the two new testsets' own
parameters, not inferred from the assertions):

| comparison | relative difference |
| --- | ---: |
| `Et_to_Pt!`, vector `KerrField`, Metal vs CPU `Float32` | **0** |
| `Et_to_Pt!`, vector `KerrEnv`, Metal vs CPU `Float32` | **0** |
| Kerr + user closure through `HostResponse`, Metal vs CPU `Float32` (`Eω`, fixed steps) | 2.0e-7 |
| same, Metal vs CPU `Float64` | 7.6e-7 |
| same, CPU `Float32` vs CPU `Float64` | 8.1e-7 |
| the closure's own contribution (CPU `Float64`, with vs without it) | 9.0e-6 |

The two zeros are genuine and expected, not a comparison of Metal with itself: a
vector-pointwise `Et_to_Pt!` is a pure elementwise `Float32` broadcast with no reduction
and no FFT, so Metal and the CPU do the identical IEEE operations in the identical order.
The host side is an explicit `DeviceSpec(Array, Float32)` and the device side an
`MtlArray`; "no stray Float64 in the kernels" separately asserts the result is not all
zero, so this is not two empty arrays agreeing. The `HostResponse` rows are the ordinary
`Float32` phase-accumulation difference — the device path adds nothing beyond what
`Float32` on the CPU already costs (2.0e-7 against the `Float32` CPU run, against 8.1e-7
between the two CPU precisions) — and the last row shows the closure is a real
contribution rather than a rounding-level one.

### Existing CPU test files

All CPU, all pass, one process, in order:

| file | result | time |
| --- | ---: | ---: |
| `test_kerr.jl` | pass (2 bare `@test`s) | 1.4 s |
| `test_polarisation.jl` | 15 pass | 0.6 s |
| `test_polarisation_env.jl` | 4 pass | 7.6 s |
| `test_polarisation_field.jl` | 8 pass | 76.3 s |
| `test_vectorplasma.jl` | 2 pass | 33.9 s |
| `test_mixtures.jl` | 2049 pass | 5.6 s |
| `test_chi2.jl` | 16 pass | 2.2 s |
| `test_raman.jl` | pass (7 bare `@test`s) | 1.8 s |
| `test_gnlse.jl` | 4 pass | 21.6 s |
| `test_multimode.jl` | 6 pass | 208.9 s |
| `test_freespace.jl` | 77 pass | 304.3 s |
| `test_interface.jl` | 317 pass | 245.7 s |

`test_kerr.jl` and `test_raman.jl` use bare `@test`s outside a `@testset`, so they print
no summary; they throw on failure and did not. `test_multimode.jl`, `test_freespace.jl`
and `test_chi2.jl` are the ones which exercise the `idcs` and two-component paths of the
new dispatch; `test_interface.jl` covers the `device`/`precision` keyword change and is
unchanged at 317.

## Measurements

**Per-call allocation of `Et_to_Pt!`** on a 1024-sample scalar block, measured with
`@allocated` after warm-up:

| responses | bytes per call |
| --- | ---: |
| 1 pointwise | 16 |
| 2 pointwise | 32 |
| 3 pointwise | 32 |

Constant, not proportional to the block size: it is the boxed closure each pointwise
kernel hands to `Broadcast.materialize!`, of the same order as `RK45`'s `_zipreduce`
(16 B per call, accepted in `gpu/10-device-model`). Nothing field-sized is allocated.

No propagation benchmarks are reported. The branch removes one pass over the oversampled
block per right-hand side (the zero-fill, when the first group is pointwise) and merges
the per-response passes of a multi-response pointwise group into one, so it should be
slightly faster; but the only case in Luna today with more than one pointwise response is
the gas mixture, and with other agents' jobs on this machine the difference is inside the
noise. The two regression-gate runs (21 cases x 2 modes each) took 1m37 and 1m36.

### Documentation build

`include("docs/make.jl")` reports exactly the 9 pre-existing unresolved `@ref`s the base
branch does (`LinearOps.βz` x3, `loadFFTwisdom`, `saveFFTwisdom`, `AbstractOutput`,
`Luna.PhysData.crystal_internal_angle`, `norm_free`, `LinearOps.make_const_linop`) and
none from this branch. `Luna.polscale` was added to `docs/src/modules/Luna.md`'s `@docs`
block; everything new in `Nonlinear.jl` is picked up by that module's `@autodocs`.

## Known gaps

- **Only the mode-averaged transform passes a `scaling` to `Et_to_Pt!`.** The radial,
  free-space and multimode transforms are still host `Float64` only, so they pass the
  identity. They will need the keyword when they gain a device path (Group E).
- **`HostResponse` is not fast and is not meant to be.** It holds five host-sized buffers
  and does two copies and a host evaluation per right-hand side. Group D replaces the
  responses which matter (plasma, Raman, χ⁽²⁾) with real kernels; the fallback stays for
  whatever a user writes.
- **`HostResponse` allocates its buffers on the first call**, from untyped fields read
  through a function barrier. `rescale` does not know the block shape, and adding it to the
  signature would change every response's `rescale` method for the benefit of one type.
- **No batched response other than `HostResponse` exists yet**, so the `Batched` contract is
  exercised only through it. Group D is what tests it properly.
- **`Chi2Field`/`Chi2Env`, `PlasmaCumtrapz` and the Raman responses are untouched** and
  remain `Columnwise()`, which is Group D's scope. On a device they now go through
  `HostResponse` instead of erroring.
- **The `device_capable` gate in `Interface.jl` is stricter than the low-level path.** A
  `prop_capillary` call which explicitly asks for a device with plasma on errors, while the
  low-level equivalent runs (slowly). That is deliberate and documented, but it is an
  inconsistency a future branch may want to turn into a warning once the fallback is less
  pathological.

## Open questions for review

1. Deviation 1 (two broadcasts instead of the plan's `reinterpret`ed `SVector{2}` write).
   The layout makes the plan's version impossible without a transpose; is the kernel
   contract as implemented the right one for Group D's χ⁽²⁾ branch?
2. Deviation 2 (`scaling` as a keyword of `Et_to_Pt!`). The alternative is to keep
   `gpu/10`'s arrangement, where `rescale` folds the scaling into the response's constants
   and `coefficients` only adds the density — one fewer argument to thread, one more place
   where a power of `E_ref` appears.
3. `rescale`'s new default for a device kind with no arrays (pass through when
   `resident_arrays` is empty). It removes boilerplate for every scalar-coefficient
   response, but it means a response which carries an array and forgets to list it in
   `resident_arrays` is passed through silently rather than refused. The transform's
   residency assertion does not catch that either, for the same reason.
4. The `@info` line per wrapped response fires at transform construction. For a scan that
   is once per point; `Logging.@info` with a `maxlog` was not used because the wrapping is
   genuinely per-run information. Is once per `Luna.setup` the right frequency?
