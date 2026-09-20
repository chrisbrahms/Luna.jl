# The device and precision model: mode-averaged propagation on any array type

Branch `gpu/10-device-model`, base `gpu/int-A` (`782f55d1`). GPU_PLAN.md §2, §3, §4.1,
§4.2, §4.4, §4.5 layer 1, §4.6, §4.11, §6 Group B, §7.

This is the first branch of Group B and the one every later branch builds on. It makes the
mode-averaged propagation — `TransModeAvg`, the Kerr responses, constant and z-dependent
linear operators, and the whole of `RK45` — run unchanged on `Array{Float64}`, on
`Array{Float32}` (scaled), on `JLArray` (test) and on `MtlArray{Float32}` (Metal, tested
on hardware). The CUDA extension loads and registers but is untested: there is no CUDA
device here.

**The default CPU path is unchanged — bit-for-bit wherever the transform lengths are
powers of two, which is every case Luna runs today; see "Why the default CPU path did not
move" for the one construction that is only rounding-level. The regression gate is exactly
`0.000e+00` for all 21 cases, in both modes and both classes.**

Review round 1 requested changes; they are all in, and what changed is in
"Changes after review round 1" at the end.

## Motivation

Luna's per-step work was written against concrete `Array{Float64}`/`Array{ComplexF64}`
buffers, scalar `for` loops, FFTW plans on host arrays and closures over host vectors. This
branch replaces that, for the mode-averaged path, with code which is generic over an array
type and a real element type, so that the same implementation runs on the host and on a
GPU. There are no CPU-only copies of any kernel.

Two constraints shape it (GPU_PLAN.md §2): Metal refuses `Float64` arrays and its kernel
compiler rejects any `double` which survives optimisation, so nothing reachable from a
kernel argument may hold one; and SI units do not fit in `Float32` — `ρ ε₀ γ₃` for helium
at 0.3 bar is 3.3e-39, below the smallest normal `Float32`, and GPUs flush subnormals to
zero.

## What changed

### `src/Device.jl` (new, included after `Utils`)

Everything which decides where a propagation runs and in what units.

- `DeviceSpec{A, T}`: a singleton carrying the array type and the real element type, so
  everything derived from it is resolved at compile time. `Luna.device()` /
  `Luna.set_device(x)` on `settings["device"]`, which accepts `:cpu` (the meaning of an
  absent key), `:auto`, `:metal`, `:cuda` or a spec.
- A registry of vendor hooks, `register_device!(name, spec; functional, synchronize,
  reclaim, memory_status)`, filled by the extensions from their `__init__` and called
  through `invokelatest` at setup or teardown only. Registering sets
  `settings["device"] = :auto` **only if the key is absent**. The fork's lazy
  `resolve_arraytype(:cuda)` is not adopted (world-age hazard).
- `alloc`, `todevice`, `tohost`, `scalar`, `upload_like`, `mask_like`, `assert_resident`.
- `UnitScaling`/`unitscaling`: the reference field and polarisation a reduced-precision run
  is expressed in. The identity for every `Float64` run.
- `GridVectors`/`gridvectors`: the mirror of the grid vectors a transform broadcasts
  against, aliasing the grid's own vectors on the host.
- `HostMirror`/`upload!`: a `Float64` host buffer plus its device copy, for quantities
  still produced by host scalar code on every right-hand side.
- `log_device`: the one-line "Propagating on … in … precision" at setup, and the
  "`:auto` requested but no GPU package is loaded in this process" message for a `Scans`
  worker.

### `src/Utils.jl`

The backend trait (`Backend`, `CPUBackend`, `DeviceBackend`, `backend`, `isdevice`), used
only for FFT planning, host-copy fallbacks and residency checks — never to select a
kernel — and the planners `plan_ft`/`plan_ift` with `iplan`/`iscale`, which split an
inverse plan into its unnormalised plan and its `1/N`.

### `ext/LunaMetalExt.jl`, `ext/LunaCUDAExt.jl` (new)

Register `DeviceSpec(MtlArray, Float32)` and `DeviceSpec(CuArray, Float64)` and the four
vendor operations. Declared under `[weakdeps]`/`[extensions]` with compat entries; neither
package is ever installed with Luna. `AbstractFFTs` becomes a direct dependency (it was
already indirect through FFTW) for the generic device planners.

### `src/RK45.jl`

Every per-step operation becomes a broadcast or a reduction, without changing any
arithmetic on the CPU.

- `combine!` accumulates `y + Σ dt bⱼ kⱼ` in one fused broadcast instead of a sequence of
  `.+=` passes; `errorestimate!` does the same for the embedded error estimate;
  `interpolant!` does it for the dense output, into the stepper's existing `yi` buffer.
  Terms are left-associated in tableau order and the zero weights are skipped exactly as
  the loops did, so all three are bit-identical. An assertion on the tableau records which
  weights the skips assume are zero.
- The four error norms become single reductions through `_zipreduce`.
- `make_fbar!` propagates in place, one field-sized buffer fewer. The first RHS evaluation
  in the `PreconStepper` constructor is made on a copy, so the caller's initial condition
  is not touched.
- `make_prop!` evaluates a closure operator into a host buffer and copies it up when the
  state is on a device, and converts the step to the state's real element type. `lastt2` is
  a `Ref`.
- `solve(output=true)` copies each save to the host.

`interpolate` now returns the stepper's own buffer rather than a fresh array. That is the
documented change for low-level users who kept the returned array; every output handler in
Luna copies it immediately, as does `solve(output=true)`.

### `src/NonlinearRHS.jl`

- `to_time!`, `to_freq!`, `copy_scale!` and `copy_scale_both!` dispatch on the element type
  rather than on `::Array`, and the copies are broadcasts over views of the first and last
  `N` samples. `to_time!` takes an explicit inverse plan and folds its `1/N` into the scale
  factor of the oversampling copy; every transform now holds that plan as `IFT`.
- `TransModeAvg` is parametric in its buffer array type, allocates through `alloc`, holds
  a `GridVectors` mirror, converts its responses with `Nonlinear.rescale`, and asserts
  residency of everything it holds.
- `norm_mode_average` and `norm_mode_average_gnlse` return callable structs (`NormModeAvg`,
  `NormModeAvgGNLSE`) whose body is one masked broadcast. `constβ=true` folds the
  propagation constant in at setup; otherwise `βfun!` is evaluated on the host per call and
  uploaded through a `HostMirror`.
- `check_norm` refuses a caller-supplied normalisation which was not built for the run's
  spec and scaling.

### `src/Nonlinear.jl`

`KerrField`, `KerrEnv` and `KerrEnvTHG` are structs parametric in the real element type,
with `rescale` methods and, where they carry an array, an `Adapt` rule. `Kerr_field(γ3)`
and the other constructors are unchanged and return the structs. The vector forms are two
broadcasts over column views. `rescale`'s fallback passes any other response through
unchanged for an unscaled `Float64` run and errors otherwise, so an ad hoc closure response
keeps working on the default CPU path. The response traits (`pointwise`/`batched`) are
`gpu/12`'s.

### `src/Luna.jl`

`setup` for mode-averaged propagation takes `device`, `precision` and `constβ`; it plans on
the chosen array type, builds the input on the host and converts it, computes the unit
scaling and logs the device and precision once. The `RealGrid` and `EnvGrid` methods
delegate to one implementation.

`run` band-limits with a masked broadcast, uploads a constant operator **after**
`Boundaries.setup` (which wraps it), re-uploads a resumed cached field, skips the host-side
`Et` on a device, and refuses a device run with absorbing boundaries.

### `src/Interface.jl`

One line: the constant mode-averaged branch passes `constβ=true`. `prop_capillary` and
`prop_gnlse` are otherwise untouched; the device keywords are `gpu/11`'s job.

## Unit scaling

The state is `e = E/E_ref` and the polarisation buffer holds `p = P/(P_ref·E_ref)`.
`E_ref = P_ref = 1` for every `Float64` run. For `Float32`, `E_ref` is the peak of the
time-domain input at `z0` rounded to a power of two and `P_ref = ε₀`.

**Deviation from GPU_PLAN.md §4.1.** The plan fixes the polarisation unit as `P/(ε₀E_ref)`
unconditionally. Here `P_ref` is a parameter of the scaling in the same way `E_ref` already
is, and is 1 rather than `ε₀` for a `Float64` run. The reason is the same as the plan's own
reason for `E_ref = 1` there: with `P_ref = ε₀` the `Float64` Kerr coefficient would be
`ρ·γ₃` and the normalisation would carry an extra `ε₀`, which is a rounding-level change to
the default CPU path for no benefit. With `P_ref = 1` on that path the coefficients are the
expressions they always were, and the gate is exactly zero. There is still one coefficient
path and one kernel body per response; the units are a property of the constants it holds.

Dynamic-range audit, helium (the worst case in Luna's range, `γ₃(He)` being the smallest):

| quantity | 1 bar | 0.3 bar |
| --- | ---: | ---: |
| `ρ ε₀ γ₃` unscaled | 1.1e-38 | 3.3e-39 |
| smallest normal `Float32` | 1.2e-38 | 1.2e-38 |
| scaled `γ₃` (`E_ref = 8192`, `P_ref = ε₀`) | — | 3.9e-34 |
| `ρ ε₀ γ₃` as the kernel sees it | — | 2.5e-20 |

At 0.3 bar the unscaled coefficient does not exist in `Float32` at all; scaled, it is
eighteen orders of magnitude above the subnormal threshold. The audit is asserted in
`test_device.jl` and `test_metal.jl` rather than only written down.

Unscaling happens in one place only, the output wrapper, which is `gpu/11`'s. Until then a
`Float32` run saves the scaled state, and the exit tests multiply by `E_ref` themselves.

## Why the default CPU path did not move

Every change to the CPU arithmetic was chosen to be exactly equivalent, and the gate
confirms it:

- The fused stage combines, error estimate and dense output are left-associated in the
  original order with the same scalar coefficients, and skip the same zero-weight terms, so
  they perform the identical sequence of operations.
- The norms are folded over a lazy `Broadcast.Broadcasted`. With an explicit `init`, Base
  reduces that with `mapfoldl` — a serial fold in index order — so it is the same
  arithmetic as the scalar loops, and it allocates only the boxed tuple result (16 bytes
  per call), not anything field-sized. (Base's *multi-array* `mapreduce`, which
  GPU_PLAN.md §4.6 assumed was a zipped fold, was measured to materialise `map(f, As...)`
  first: 2.4 MB per call at n = 1e5. It is not used. `GPUArrays` has a `mapreduce` method
  for a `Broadcasted` of its own style, so the device path is still its tree reduction.)
- Folding the inverse plan's `1/N` into the oversampling copy is exact **when `1/N` is a
  power of two**: scaling by a power of two is exact, and an exactly scaled FFT input gives
  an exactly scaled output. That covers every transform over the time axis alone, Luna's
  time grids being powers of two by construction (`Grid.jl:57, 78, 178, 200`). A
  multi-axis free-space transform normalises by `1/(Nt·Nx·Ny)` instead, and
  `Grid.FreeGrid`/`Grid.Free2DGrid` accept any `Nx`, `Ny`; for a transverse length which is
  not a power of two the folded and unfolded routes differ at rounding level. Measured
  directly against `copy_scale!` + `ldiv!`: 0 for length 16 and 32, 5.6e-17 for 24, 6.9e-17
  for 48. Every free-space grid in the gate, the tests and the examples happens to use a
  power of two, so the gate cannot see this; it is within what GPU_PLAN.md §3 allows, and
  it is stated here, in the `to_time!` docstring and in the developer guide rather than
  left to the reader. (Found by review round 1, finding 3.)
- `pre/β·√aeff` precombined at setup under `constβ` uses the same operands in the same
  order as the per-element expression it replaces.
- The masked broadcasts (`norm!`, the band limit) write an exact zero where the indexed
  assignment did and leave the in-band elements untouched.
- `E_ref = P_ref = 1` makes the scaling arithmetically invisible.

## Tests

`julia --project=<worktree> -t 1`, with `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)`, `Luna.set_fftw_wisdom(false)` and
`LinearAlgebra.BLAS.set_num_threads(1)`. Apple M1 Pro, Julia 1.13, with other agents' jobs
on the same machine, so the times are upper bounds.

### The regression gate

```
LUNA_REGRESSION_BRANCH=gpu/int-A julia --project=$PWD -t 1 test/test_regression.jl
```

against a baseline generated from the merge-base `782f55d1`
(`test/regression/generate.jl 782f55d1`).

**460 pass, 0 fail. Every case, both modes, both classes: `0.000e+00`.**

| case | `:fixed` Eω / tol | `:fixed` stats / tol | `:adaptive` Eω / tol | `:adaptive` stats / tol |
|---|---|---|---|---|
| `modeavg_field_kerr` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 5.9e-07 |
| `modeavg_field_nothg` | 0 / 1.0e-12 | 0 / 1.0e-12 | 0 / 2.3e-12 | 0 / 2.1e-04 |
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

**Largest difference over all cases and modes: 0.000e+00.** No case moved. Step counts are
unchanged (the gate checks them separately and fails hard if they move).

### New test files

`test/test_device.jl` (added to `runtests.jl`; the `JLArray` half needs `JLArrays`, which
is in `[extras]`/`[targets]`, and skips itself otherwise):

| testset | assertions |
| --- | ---: |
| backend trait | 16 |
| device spec and settings | 24 |
| allocation and transfer | 18 |
| residency assertions | 6 |
| unit scaling | 8 |
| FFT planner dispatch | 9 |
| constβ is checked, not trusted | 2 |
| JLArray basics | 14 |
| RK45 kernels on JLArray | 11 |
| mode-averaged Kerr on JLArray | 24 |
| pressure gradient on JLArray | 12 |
| a device run refuses host-only machinery | 5 |
| Float32 on the CPU | 11 |

**157 pass, 0 fail.** The propagation comparisons:

| comparison | max relative difference in `Eω` |
| --- | ---: |
| `JLArray` vs host, field-resolved Kerr | 0 |
| `JLArray` vs host, envelope Kerr | 0 |
| `JLArray` vs host, pressure gradient (z-dependent operator) | 0 |
| CPU `Float32` (scaled) vs CPU `Float64`, He at 0.3 bar | 3.7e-7 |

`test/test_metal.jl` (not part of the suite; run from an environment with Metal):

| testset | assertions |
| --- | ---: |
| Metal registration | 8 |
| allocation and transfer on Metal | 7 |
| no stray Float64 in the kernels | 16 |
| mode-averaged Kerr on Metal | 20 |
| Metal against the Float64 CPU path | 4 |
| pressure gradient on Metal | 7 |
| the simple interface stays on the CPU | 5 |
| Metal refuses what it cannot run | 1 |

**68 pass, 0 fail**, on an Apple M1 Pro with Metal.jl v1.11.

| comparison | field-resolved | envelope |
| --- | ---: | ---: |
| Metal `Float32` vs CPU `Float32` | 2.7e-7 | 5.2e-8 |
| Metal `Float32` vs CPU `Float64` | 3.0e-7 | 7.3e-8 |
| CPU `Float32` vs CPU `Float64` | 3.5e-7 | 7.3e-8 |
| Metal vs CPU `Float64`, He at 0.3 bar | 2.4e-7 | — |

That is two to three orders of magnitude better than the ~1e-4 GPU_PLAN.md §8 expected for
`Float32` phase accumulation, on these short propagations. It will get worse with distance;
the per-case numbers are what the later branches report.

### Existing test files

All CPU, all pass:

| file | result | time |
| --- | ---: | ---: |
| `test_rk45.jl` | 48 pass | 6.9 s |
| `test_kerr.jl` | 2 pass | 1.2 s |
| `test_linops.jl` | 196 pass | 12.4 s |
| `test_gradient.jl` | 7 pass | 10.4 s |
| `test_tapers.jl` | 2 pass | 16.5 s |
| `test_output.jl` | 71 pass | 13.7 s |
| `test_gnlse.jl` | 4 pass | 21.4 s |
| `test_mixtures.jl` | 2049 pass | 5.0 s |
| `test_stats.jl` | 1 pass | 3.2 s |
| `test_noise.jl` | 49 pass | 15.6 s |
| `test_utils.jl` | 33 pass | 3.8 s |
| `test_interface.jl` | 301 pass | 4m02 |
| `test_multimode.jl` | 6 pass | 3m38 |
| `test_freespace.jl` | 77 pass | 5m12 |
| `test_boundaries.jl` | 182 pass | 15.0 s |
| `test_fields.jl` | 179 pass | 50.3 s |
| `test_chi2.jl` | 16 pass | 1.5 s |
| `test_raman.jl` | 7 pass | 1.1 s |
| `test_modes.jl` | 587 pass | 24.3 s |
| `test_radialgrid.jl` | 159 pass | 4.3 s |
| `test_grid.jl` | 88 pass | 5.5 s |
| `test_polarisation.jl` | 15 pass | 0.2 s |
| `test_polarisation_env.jl` | 4 pass | 5.7 s |
| `test_vectorplasma.jl` | 2 pass | 34.3 s |
| `test_linearprop.jl` | 2 pass | 1.5 s |
| `test_polarisation_field.jl` | 8 pass | 1m17 |

`using Luna` without Metal or CUDA installed was checked in a fresh process: neither is in
`[deps]` (`Pkg.status` shows only `AbstractFFTs`, `Adapt` and `GPUArraysCore` of the three
new-ish dependencies), `Luna.devicenames()` is empty, `Luna.device()` is
`DeviceSpec(Array, Float64)`, and `prop_capillary` runs. With `using Metal` in the same
process the extension registers, `settings["device"]` becomes `:auto` and `Luna.device()`
is `DeviceSpec(MtlArray, Float32)`; `Luna.set_device(:cpu)` opts back out.

The documentation build (`include("docs/make.jl")`) reports the same six pre-existing
unresolved `@ref`s it did on the base branch, and none from this branch.

Two existing tests called `NonlinearRHS.to_time!` with the forward plan
(`test_interface.jl:426`, `test_freespace.jl:404`) and now pass the transform's `IFT`.
That is the only change to an existing test.

### CI

`.github/workflows/run_tests.yml` gains a macOS job which installs Metal into a separate
environment stacked on the package and runs `test_metal.jl`. The Linux and Windows jobs do
not install it, and Metal never enters `Project.toml`.

## Benchmarks

`benchmark/run.jl`, the same cases as `PR_00-harness.md` measured on `evanescent`:

| Case | state | rhs (base → here) | step (base → here) | prop (base → here) |
| --- | ---: | ---: | ---: | ---: |
| `modeavg_field_kerr` | 2049 | 41.2 → 37.9 µs | 495 → 463 µs | 36.7 → 35.6 ms |
| `modeavg_field_nothg` | 2049 | 100 → 98.4 µs | 852 → 832 µs | 45.2 → 49.0 ms |
| `modeavg_field_raman` | 2049 | 251 → 242 µs | 1.76 → 1.72 ms | 72.9 → 68.0 ms |
| `modeavg_field_vector` | 4098 | 5.03 → 4.94 ms | 31.2 → 30.5 ms | 6.09 → 6.47 s |
| `modeavg_env_kerr` | 2048 | 21.1 → 18.4 µs | 366 → 335 µs | 31.1 → 35.0 ms |
| `modeavg_env_thg` | 2048 | 40.3 → 38.5 µs | 467 → 434 µs | 35.6 → 38.0 ms |
| `gnlse_sech` | 2048 | 20.7 → 18.5 µs | 516 → 487 µs | 31.4 → 38.2 ms |
| `gradient_field_kerr` | 2049 | 80.8 → 77.4 µs | 1.16 → 1.12 ms | 235 → 226 ms |
| `taper_field_kerr` | 2049 | 309 → 300 µs | 3.65 → 3.58 ms | 101 → 107 ms |

The RHS and step times improve by 3–13 %: one pass fewer over the oversampled array (the
folded inverse normalisation), one buffer fewer and one copy fewer per stage (the in-place
`make_fbar!`), and the fused stage combines. `prop` is dominated by setup at these sizes
and is noisy under contention; nothing in it regressed beyond that noise.

`benchmark/device.jl`, the same mode-averaged Kerr propagation on each device and
precision, sweeping the time-grid size:

```
Mode-averaged Kerr, He at 1 bar, 0.01 m, 20 fixed steps, boundary=:none
1 Julia thread, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom

device           trange      state          rhs         step         prop
------------------------------------------------------------------------
CPU Float64       400 fs       1025    17.792 µs   223.042 µs     4.926 ms
CPU Float32       400 fs       1025    11.666 µs   183.666 µs     3.940 ms
metal Float32     400 fs       1025   386.167 µs     2.376 ms    52.921 ms
CPU Float64      1600 fs       4097    82.917 µs   940.625 µs    19.635 ms
CPU Float32      1600 fs       4097    51.375 µs   768.666 µs    15.784 ms
metal Float32    1600 fs       4097   405.208 µs     2.618 ms    63.552 ms
CPU Float64      6400 fs      16385   551.417 µs     5.162 ms   105.844 ms
CPU Float32      6400 fs      16385   271.833 µs     3.475 ms    71.100 ms
metal Float32    6400 fs      16385   517.750 µs     3.304 ms    77.051 ms
```

Two things to read out of this.

**`Float32` on the CPU is 1.5-2x faster than `Float64`**, on RHS and on step, at every
size. That is worth having on its own, and it is the same scaling layer the GPU needs.

**Metal is launch-bound on a mode-averaged run, exactly as GPU_PLAN.md §8 predicted.**
The RHS costs ~390 µs however small the problem is: a single-column mode-averaged
right-hand side is about ten kernel launches and two FFTs over a few thousand elements,
and the launches dominate. Only at the largest size does it come within 10 % of the CPU,
and it is still slower. The win is in the multi-column transforms (radial, free-space,
fixed-rule multimode), which are later branches; this branch's job is that the mode-averaged
path *runs* correctly on a device, which it does. A GPU should not be used for a
single-column run, and the user page says so.


## Known gaps

Stated as such; each is another branch's scope in GPU_PLAN.md.

- **The simple interface is not device-capable and deliberately stays on the CPU.**
  `Interface` passes `device=Luna.HostSpec()` explicitly, so `prop_capillary` and
  `prop_gnlse` give exactly the result they always did whatever `settings["device"]` says.
  They gain the `device`/`precision` keywords in `gpu/11`, together with the boundaries
  and the output wrapper. Tested on hardware (`test_metal.jl`, "the simple interface stays
  on the CPU").
- **Only mode-averaged Kerr runs on a device.** The radial, free-space and multimode
  transforms, and the plasma, Raman and χ⁽²⁾ responses, are unchanged host code. They keep
  compiling and passing their CPU tests, and they are refused rather than run wrongly on a
  device.
- **`boundary=:none` only on a device**, and no per-step statistics: the absorbers and the
  statistics are host scalar code. `Luna.run` errors with a message naming the keyword.
  `gpu/11`.
- **`Output` is device-unaware.** A device run needs an output which copies to the host;
  both test files define a five-line one, and the real `ScaledOutput` (which also unscales)
  is `gpu/11`'s. `MemoryOutput`/`HDF5Output` still allocate `ComplexF64`, so a `Float32`
  run is widened on the way into the output. Also `gpu/11`.
- **A `Float32` run saves the scaled state.** Unscaling is in the output wrapper, so until
  `gpu/11` the caller multiplies by `transform.scaling.Eref`. The tests do.
- **Tapers and gradients upload `β` per stage from the host**, and a closure `linop!` is
  evaluated into a host buffer per stage. That is layer 1 of GPU_PLAN.md §4.5 and the
  interim state until `gpu/23` tabulates them. It works but is not fast on a device.
- **The CUDA extension is untested on hardware.** It loads, registers and has the same
  shape as the Metal one, but there is no CUDA device on this machine. `test/test_cuda.jl`
  is a later branch's.
- **`interpolate` returns the stepper's buffer.** A low-level user who kept the returned
  array across steps now has to copy it. Documented on the function.
- **`Nonlinear.rescale` has no fallback for an arbitrary response.** On a device or in
  single precision, a response which is not one of the Kerr structs errors with a message.
  The `HostResponse` fallback which copies a column block to the host is `gpu/12`.

## Open questions for review

1. `P_ref` as a parameter of the scaling rather than a fixed `ε₀` (see **Unit scaling**).
   It is a deliberate deviation from GPU_PLAN.md §4.1, taken so that the `Float64` gate
   stays exactly zero. If the reviewer prefers the plan as written, the cost is a
   rounding-level move of every mode-averaged case.
2. `constβ` is a new keyword of `Luna.setup` and of `norm_mode_average`, set from
   `Interface` in the constant branch. GPU_PLAN.md §4.4 says `βfun!` should be uploaded
   "only when the operator is z-dependent" but does not say how the transform learns which
   it is; an explicit keyword seemed better than probing `βfun!` at two values of `z`.
3. `AbstractFFTs` as a direct dependency. It was already indirect through FFTW, and
   GPU_PLAN.md §3 names only `Adapt` and `GPUArraysCore`, so it is one more than the plan
   says.
4. `_adapt` moves an array with the array type's own constructor rather than
   `Adapt.adapt`. `JLArrays` defines no `adapt_storage` rule for its bare array type and
   `Adapt.adapt` then silently returns the host array, which would put a host array into a
   device kernel. `Adapt` is still used for the `adapt_structure` rules on `GridVectors`
   and `KerrEnvTHG`.

## Changes after review round 1

The review's verdict was "request changes" on two majors, plus six minors and six nits.
All fourteen are addressed. Nothing in the review's reproduction of the numbers changed:
the regression gate is still exactly `0.000e+00` on all 21 cases, in both modes and both
classes.

### 1 (major) — loading Metal or CUDA broke `prop_capillary`/`prop_gnlse`

The extension sets `settings["device"] = :auto` for the whole process and
`Interface.setup` inherited it through `Luna.setup`'s default keyword, so the simple
interface failed (or, for `modes=4`, silently ran on the CPU) as soon as a GPU package was
loaded — for whatever reason the user loaded it.

`Interface` now passes `device=Luna.HostSpec()` explicitly at all three call sites: both
mode-averaged branches of `prop_capillary` (`Interface.jl:896, 909`) and `prop_gnlse`
(`Interface.jl:1125`). Until `gpu/11` plumbs `device` and `precision` through the simple
interface, it stays on the CPU in double precision whatever is loaded, which is what
"the simple interface is unchanged" requires — the boundaries, the statistics and every
response but Kerr are still host code.

`test_metal.jl` gains "the simple interface stays on the CPU": with Metal loaded, it runs
`prop_capillary` (Kerr only, `plasma=false`) and `prop_gnlse` under `set_device(:cpu)` and
again under `set_device(:auto)` — checking first that `:auto` really does resolve to
`MtlArray` — and asserts that both give `ComplexF64` output **equal** to the CPU run.
`docs/src/gpu.md` now opens with a warning saying so, and it is in "Known gaps".

### 2 (major) — the "`:auto` requested but no GPU package is loaded" message was dead code

`setup`'s `device` keyword defaulted to `Luna.device()`, which had already resolved the
setting, so `log_device`'s `request === :auto` test could never be true from any Luna code
path. GPU_PLAN.md §3 requires the message and §8 relies on it.

The default is now the *unresolved* request, `Luna.device_request()` (a new one-line
accessor: `get(settings, "device", :cpu)`), and `setup_mode_average` resolves it itself.
`test_device.jl`'s "device spec and settings" now sets `settings["device"] = :auto` in a
process with no GPU package and `@test_logs` for the message, and for the ordinary
"Propagating on …" line in the `:cpu` case.

### 3 (minor) — the bit-identity claim was unconditional

Qualified in three places (the headline, "Why the default CPU path did not move", the
`to_time!` docstring and `device_model.md` rule 5), with the reviewer's measurement: exact
for a power-of-two transform length, 5.6e-17 for 24 and 6.9e-17 for 48. Every grid Luna
runs today uses a power of two.

### 4 (minor) — `check_norm` ignored the scaling

`NormModeAvg` and `NormModeAvgGNLSE` now carry the `UnitScaling` they were built for, and
`check_norm` compares it as well as checking residency. Both failure modes now go through
one message, the same one the generic fallback gives, so a normalisation built with the
right array type but the wrong units is refused instead of producing a polarisation wrong
by `Pref`. `Luna.all_resident` (the predicate form of `assert_resident`) was added for
this. Tested in "a device run refuses host-only machinery".

### 5 (minor) — no device test of the z-dependent path

`test_device.jl` and `test_metal.jl` gain a pressure-gradient case
(`Capillary.gradient` + `LinearOps.make_linop`, `constβ=false`, fixed steps,
`boundary=:none`), which is the only thing that exercises `NormModeAvg`'s `HostMirror`
branch and `RK45.make_prop!`'s host-buffer branch — the two pieces of §4.5 layer 1 the
exit criteria name. Both testsets also assert the shape of the mirror (host `Float64`
buffer, staging buffer in the device precision, device array; and, on the host, that the
device array *is* the host buffer and `upload!` does nothing).

`JLArray` vs host, fixed steps: max relative difference in `Eω` **0**. Metal vs the
`Float32` CPU path: within 1e-4, as for the constant-operator cases.

### 6 (minor) — the documented way to run the test files did not work

`docs/src/gpu.md` and the headers of both test files now give the same package list the CI
job installs (`Test`, `GPUArraysCore`, `Adapt` for `test_metal.jl`; those plus
`AbstractFFTs`, `FFTW` and `JLArrays` for `test_device.jl`), as a copy-pastable
`Pkg.develop`/`Pkg.add` pair.

### 7 (minor) — the residency assertion had holes

`TransModeAvg` now asserts `gv.twin` as well, and the arrays its rescaled responses carry,
through a new `Nonlinear.resident_arrays(response)` hook which returns a tuple (empty by
default, `(k.C,)` for `KerrEnvTHG`).

### 8 (minor) — `constβ=true` was unguarded

`NonlinearRHS.check_constβ` evaluates `βfun!` a second time at
`NonlinearRHS.CONSTβ_PROBE_Z = 1e-3` m and errors if the result differs from the one at
`z = 0`, naming `constβ=false` as the fix; a `βfun!` which throws at that `z` gets its own
message. A silent physics error (β frozen at the entrance) becomes an error at setup.
Tested both ways in "constβ is checked, not trusted".

### 9, 11, 13, 14 (nits)

- The `_zipreduce` docstring and the PR now say 16 bytes per call (the boxed tuple), not
  "allocates nothing".
- `_adapt(::Type{Array}, x)` copies a device array to the host rather than returning it,
  so `todevice(HostSpec(), x)` means what its name says. (Still unreachable from Luna.)
- `test_device.jl` restores the task-local `:ScalarIndexing` key at the end of the file, so
  one file no longer changes what the rest of `Pkg.test()` sees. The process-global default
  is `ScalarDisallowed` in a non-interactive run anyway, which is what `allowscalar(false)`
  sets it to; that is recorded in a comment.
- `Project.toml`'s `CUDA` and `Metal` compat entries are in alphabetical order.

### 10, 12 (nits)

- `Utils.iplan`/`Utils.iscale` now error for anything that is not a `ScaledPlan` instead of
  passing a forward plan through with a factor of 1, so `to_time!` called with `FT` says so
  rather than transforming the wrong way. Every inverse plan Luna makes is a `ScaledPlan`
  (FFTW, Metal, CUDA and the JLArray shim all return one), dispatch decides, and the check
  costs nothing per step. Tested in "FFT planner dispatch".
- `errorestimate!` carries a note saying that its skip of the second term encodes DOPRI5's
  zero pattern, that the `@assert`s above hold the tableau to it, and what would have to
  change if the tableau ever became a runtime choice.

### Re-run after the changes

| | |
|---|---|
| `test_regression.jl` vs `782f55d1` | **460 pass, 0 fail**, every case `0.000e+00` |
| `test_device.jl` | **157 pass, 0 fail** (13 testsets) |
| `test_metal.jl` (M1 Pro) | **68 pass, 0 fail** (8 testsets) |
| `test_rk45.jl` | 48 pass |
| `test_interface.jl` | 301 pass |
| `test_gnlse.jl`, `test_gradient.jl`, `test_tapers.jl`, `test_linops.jl` | pass |
