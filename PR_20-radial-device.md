# `gpu/20-radial-device`: the radial transform, the free-space normalisation and the transverse collar on a device

Base: `gpu/int-D` at `90826dc4`.

Radially symmetric free-space propagation — field-resolved and envelope, Kerr, plasma and
Raman, with `boundary=:rate` — now runs on Metal and on a `JLArray`, through the low-level
interface. The default CPU path is bit-for-bit what it was: the regression gate against
the base commit is exactly `0.000e+00` on all 460 rows.

The three pieces this branch also has to make general, because `gpu/21-free-device` is cut
from it, are `NonlinearRHS.FreeSpaceNorm` (shared by all three free-space transforms),
`Boundaries.spatialcollar` and `Grid.radial_matmul!`. `TransFree`/`TransFree2D` and
`Boundaries.CartesianCollar` are untouched apart from the shared normalisation.

## Commits

| commit | what |
| --- | --- |
| `5e8a021e` | `TransRadial` and `FreeSpaceNorm` on the device; `Luna.setup`'s radial `device`/`precision` |
| `7656fe6a` | `Boundaries.RadialCollar` on the device |
| `ba57b1b4` | `test_device.jl`: the radial transform on `JLArray` |
| `a1df38ea` | `test_metal.jl`: the radial transform on hardware |
| `ec4ae61e` | Regression gate: a multi-column Raman case, and case filters |
| `05655e24` | `benchmark/radial.jl` |
| `a331141d` | `docs/src/developer/device_model.md`, `docs/src/gpu.md` |

## What changed

### `NonlinearRHS.TransRadial` (`src/NonlinearRHS.jl`)

Parametric in its buffer array type, like `TransModeAvg`: every buffer is allocated with
`Luna.alloc`, the grid vectors are mirrored (`Luna.gridvectors`), the inverse plan is held
explicitly, the responses are given the run's `spec` and `scaling` through
`Nonlinear.rescale_responses`, the noise field is uploaded and divided by `Eref`, and the
constructor asserts residency of every buffer, mirror, transform matrix and response array.
The transform carries its `UnitScaling`, and `Luna.runscaling` knows about it, so a
`Float32` radial run gets a `ScaledOutput` like a mode-averaged one.

**The Hankel step is one GEMM per direction.** It was one `mul!` per polarisation component
on `view(Eto, :, ip, :)`; it is now one `Grid.radial_matmul!` on the block reshaped to
`(nto·npol, nr)`. On Metal that is the difference between MPSGraph's accelerated matmul and
a fallback: the accelerated path needs plain zero-offset operands of equal element type, and
a `view` is neither. The transform therefore holds its own copies of the grid's `Tfwd`/`Tbwd`
**in the time-domain element type** (real on a `RealGrid`, complex on an `EnvGrid`) rather
than the grid's `Float64` ones. `Grid.radial_matmul!` already existed (`gpu/02-radialgrid`
wrote it) and needed no change; it is now what the per-step path uses.

**The frequency-domain normalisation** was
`nl .*= ωwin .* (-im.*ω) ./ (2 .* normfun(z))`, which rebuilt two host vectors on every
right-hand side. It is now one fused broadcast over a precombined `prefac = ωwin·(-iω)·Pref`
(built once, mirrored) and the normalisation array, written in the same association, with
the `2` converted by `Luna.scalar`.

### `NonlinearRHS.FreeSpaceNorm` (shared with `gpu/21`)

Takes a `spec`; its output array and the vectors its kernel broadcasts against live there.
The isotropic fill is one broadcast over

- `n(ω; z)`, evaluated by host scalar code into a `(Nω, Npol)` `Luna.HostMirror` and
  uploaded once per call (it runs once per `z`, not once per element, so it does not need a
  kernel);
- mirrors of `ω` and `grid.sidx`, and of `kperp2` and the k-space window reshaped to
  `(1, 1, Nk...)`.

The kernel is the old `normfactor` with every constant converted to the element type — `c`,
`μ₀`, `κmax`, `ℓ` — and `βz` written as
`βsq < 0 ? complex(0, -√-βsq) : complex(√βsq, 0)` rather than with an `im` literal, which is
a `Complex{Int}` and would put a 64-bit integer into a kernel. The value is identical in
`Float64`.

**The k-window mirror is rebuilt lazily.** `Boundaries.setup` calls `reflength!` *after* the
transform is constructed, to set `ℓ`, `κmax` and the k-space absorber profile; `reflength!`
clears a `mirrored` flag and the next `fillnorm!` refills the mirror. Nothing else the
kernel broadcasts against changes after construction.

The **crystal-optics** variant (root-finding per `(ω, kx)`) stays host scalar code and
copies its result up through a staging buffer. It is reached only from `TransFree`/
`TransFree2D`, which are host-only until `gpu/21`.

`norm_radial`/`norm_free`/`norm_free2D` and the `const_` variants take a `spec` keyword.
Because the normalisation is a *positional* argument of `Luna.setup` — every low-level
radial script builds one before it knows what device the run will use — `Luna.setup` calls
the new `NonlinearRHS.retarget`, which rebuilds a `FreeSpaceNorm` for the run's spec,
carrying over anything `reflength!` has already set, and returns anything else unchanged on
the default host `Float64` path (and errors with a message naming the constructors
otherwise). `check_norm` gains a `FreeSpaceNorm` method, which checks residency only: the
unit scaling is not the normalisation's, since `βz/(μ₀ω)` is physics and `Pref` is folded
into the transform's `prefac`.

### `Luna.setup` (`src/Luna.jl`)

The four radial methods (two grid types × `QDHT`-or-`RadialGrid`) become two forwarding
methods and one `setup_radial`, which takes `device` and `precision`, logs the choice,
retargets the normalisation, plans the oversampled and state transforms on the run's array
type through `Utils.plan_ft`, keeps host plans for building the input fields (`Fields` is
host scalar code), computes the unit scaling from the peak of the input taken back to
`(t, pol, r)`, and uploads the initial field. Behaviour on the default CPU path is
unchanged, including which FFTW plans are made.

### `Boundaries.RadialCollar` (`src/Boundaries.jl`)

Its transform matrices, absorption rate, integration weights and scratch factor become
parametric in their array type and are moved with `Luna.todevice` after being converted to
the state's precision on the host (a real-to-complex conversion is not one `todevice`
performs). `spatialcollar` builds the collar's buffer with `similar(Et, ...)` rather than a
host `ComplexF64` array. The per-step code is unchanged — it was already one
`radial_matmul!` per direction, broadcasts, and `mapreduce` folds over lazy `Broadcasted`
arguments.

### `Interface.jl`

Comment only. `prop_capillary`/`prop_gnlse` never build a radial run, so the `_cpu_only!`
guard is now about multimode propagation alone.

## Rounding

The brief expected the radial gate cases to move at rounding level from the single-GEMM
formulation. **They did not move at all.** For one polarisation component — which is what
both radial gate cases have — `reshape(Eto, nto*1, nr)` is the same
matrix, with the same leading dimension, as `view(Eto, :, 1, :)` was, so BLAS sees an
identical problem and returns identical bits. The gate against the base commit is
`0.000e+00` on every row, including the two radial cases.

For two components the reshape changes the GEMM's `m` from `nto` to `2nto`, which may change
BLAS's blocking and hence the summation order. `test_device.jl` measures it directly
(reshaped GEMM against the per-view loop, `Float64` and `ComplexF64`): **exactly zero** for
one component and **below 1e-14 relative** for two. No gate case is a two-component radial
run, so nothing in the matrix sees even that. Two-component radial propagations are run by
`test_freespace.jl` (`pol = true`) and by two of the examples, and they pass.

The two other candidate sources of movement are also exactly zero in the gate: the fused
normalisation broadcast (same association, same values — `Pref` is `1.0` at `Float64`) and
the `FreeSpaceNorm` broadcast (same expression, same constants at `Float64`).

## Regression gate

M1 Pro, Julia 1.13.0, `-t 1`, `set_fftw_mode(:estimate)`, one FFTW thread, one BLAS thread,
wisdom off, `shotnoise=false`, through the timing wrapper.

| baseline | cases | result | largest difference |
| --- | --- | --- | --- |
| `fa556e6f` (source-identical to `90826dc4`, this branch's base) | 21 | **460 pass, 0 fail** | `0.000e+00` on every row |
| `fdf8dbe3` (`evanescent`) | 21 | **460 pass, 0 fail** | 7.240e-05 (`multimode_field_plasma` adaptive statistics) |
| `fa556e6f`, `radial_field_raman` only (own baseline) | 1 | **6 pass, 0 fail** | `0.000e+00` |

Against the base, **every row of every case in both modes is exactly `0.000e+00`**. No case
moved, no step count moved.

Against `evanescent` the record is the `gpu/int-D` one reproduced to every digit, with no
new rows — the same twelve non-zero rows, the same values:

| case | mode | `Eω` | `stats` |
| --- | --- | ---: | ---: |
| `modeavg_field_plasma` | fixed | 1.131e-15 | 4.769e-15 |
| `modeavg_field_plasma` | adaptive | 5.616e-09 | 7.103e-06 |
| `modeavg_field_adk` | fixed | 1.245e-15 | 5.202e-15 |
| `modeavg_field_adk` | adaptive | 2.283e-13 | 8.550e-08 |
| `modeavg_field_vector` | fixed | 1.053e-15 | 6.300e-13 |
| `modeavg_field_vector` | adaptive | 1.962e-10 | 3.294e-05 |
| `multimode_field_plasma` | fixed | 3.785e-14 | 7.340e-13 |
| `multimode_field_plasma` | adaptive | 7.636e-08 | 7.240e-05 |
| `radial_field_kerr` | fixed | 1.901e-15 | 0 |
| `radial_field_kerr` | adaptive | 6.713e-11 | 0 |
| `radial_env_kerr` | fixed | 2.866e-15 | 0 |
| `radial_env_kerr` | adaptive | 4.017e-11 | 0 |
| `free2d_field_chi2` | fixed | 7.683e-16 | 0 |
| `free2d_field_chi2` | adaptive | 4.841e-14 | 0 |
| `free2d_env_chi2` | fixed | 8.348e-16 | 0 |
| `free2d_env_chi2` | adaptive | 3.151e-14 | 0 |

The radial rows are `gpu/02-radialgrid`'s, unchanged; this branch adds nothing to them.

### The new case: `radial_field_raman`

The matrix had no multi-column Raman case: every Raman case in it (mode-averaged field,
mode-averaged envelope, GNLSE) has a single transverse column, so none of them exercised
`RamanPolarField`'s block-sized buffers or its one-pair-of-FFTs-per-right-hand-side path.
`gpu/int-D` flagged that as work for Group E.

`radial_field_raman` is radial Kerr + Raman in nitrogen at 1 bar, on the same grid as the
two radial Kerr cases, at 50 µJ. The energy is measured, not guessed — against the same case
with the Raman response removed, at the end of the propagation:

| energy | Raman contribution | adaptive steps (unperturbed / +1 ulp) |
| ---: | ---: | --- |
| 5 µJ | 8.6e-04 | 23 / 23 |
| 20 µJ | 3.5e-03 | 23 / 23 |
| **50 µJ** | **9.3e-03** | **23 / 23** |
| 100 µJ | 2.1e-02 | 23 / 23 |

50 µJ is 1.6 GW, comfortably below the critical power for self-focusing in nitrogen. The
case costs about 0.2 s per mode.

Tolerances, from `test/regression/sensitivity.jl` (100× the one-ulp sensitivity, floor
1e-12):

| case | mode | `:Eω` | `:stats` | measured sensitivity |
| --- | --- | ---: | ---: | ---: |
| `radial_field_raman` | fixed | 1.0e-12 | 1.0e-12 | 2.4e-15 / 0 |
| `radial_field_raman` | adaptive | 1.1e-09 | 1.0e-12 | 1.1e-11 / 0 |

In line with the two radial Kerr cases.

**Its baseline has to be generated by the integrator.** No pre-`gpu/20` baseline directory
contains it, and this branch deliberately did not add it to the shared ones. It was
generated into a scratch directory instead (`generate.jl fa556e6f <dir> radial_field_raman`,
which reuses the existing baseline worktree and writes only that case's file), and the
21-case runs above were done with `LUNA_REGRESSION_SKIP=radial_field_raman`.

To make that possible the gate gains `LUNA_REGRESSION_ONLY` and `LUNA_REGRESSION_SKIP`
(comma-separated case names, `ONLY` applied first). The header prints the selection whenever
it is not the whole matrix, an unknown name is an error, and neither variable can make a
failing case pass. `test/regression/README.md` documents the pattern.

## Tests

All on M1 Pro / Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread, wisdom
off, through the timing wrapper.

| what | result | before |
| --- | --- | --- |
| `test/test_regression.jl` × 3 | 460/0, 460/0, 6/0 — the record above | |
| `test/test_device.jl` (JLArrays) | **675 pass, 0 fail**, 36 testsets | 611 / 31 |
| `test/test_metal.jl` (Metal 1.11, hardware) | **452 pass, 0 fail**, 20 testsets | 386 / 16 |
| `test/test_freespace.jl` | **77 pass, 0 fail** | 77 pass |
| `test/test_boundaries.jl` | **182 pass, 0 fail** | 182 pass |
| `include("docs/make.jl")` | the 9 pre-existing unresolved `@ref`s, none new | |
| the five radial examples | all run | |

### What was added

`test/test_device.jl`

- `the Hankel step as one GEMM` (host): the reshaped GEMM against the per-polarisation
  `mul!` on a view, `Float64` and `ComplexF64`, one and two components. Exactly equal for
  one, `< 1e-14` for two.
- `the Hankel GEMM on JLArray`: the same product on a device array under
  `allowscalar(false)`, including the `out === A` case, and an assertion that the reshape
  stays a `JLArray` matrix.
- `radial propagation on JLArray`, `radial plasma and Raman on JLArray`,
  `radial envelope Raman on JLArray`: a `radialcase` helper (field-resolved and envelope,
  Kerr / Kerr+plasma / Kerr+Raman, `boundary=:rate`, fixed steps) run on the host and on
  `JLArray`, plus residency assertions on the transform, its mirrors, its Hankel matrices
  and the retargeted normalisation.

Measured differences (per save, normalised by the largest `|Eω|` in that save — the gate's
metric; an elementwise relative difference is meaningless in the k-channels the evanescent
taper has emptied):

| case | JLArray vs host |
| --- | ---: |
| radial field Kerr | 8.589e-16 |
| radial envelope Kerr | 3.704e-16 |
| radial field Kerr+plasma | 9.844e-16 |
| radial field Kerr+Raman | 6.223e-16 |
| radial envelope Kerr+Raman | 8.383e-16 |

The plasma case runs on a tighter grid at a smaller waist than the Kerr ones so that argon
actually ionises, and the test asserts the ionised fraction is non-zero — otherwise it would
silently become a comparison of two Kerr-only runs, which is what the `gpu/13-plasma` review
found in the gate's original plasma cases.

`test/test_metal.jl`

- `the Hankel GEMM on Metal`: `Float32` and `ComplexF32`, one and two components, against
  the host product, with assertions that what the GEMM sees is a plain `MtlMatrix` of the
  same element type as the transform matrix. `ComplexF32 × ComplexF32` through MPSGraph is
  the least-exercised of the paths Luna uses (GPU_PLAN.md §8); it works, and agrees with the
  host to better than 1e-5.
- `no stray Float64 in the radial kernels`: the element type of every array the free-space
  normalisation and the transverse collar hold, the normalisation kernel compiled and run on
  the GPU against the host at the same parameters, on a grid fine enough
  (`R = 100 µm`, 64 points, 400–4000 nm) that the evanescent branch — the one which carries
  the taper — is reached; the lazy k-window mirror rebuilt after `reflength!`; and
  `apply_kspace!` run on a device state.
- `radial propagation on Metal` and `radial Kerr on Metal against the Float64 CPU path`.

Measured:

| case | Metal vs CPU `Float32` | Metal vs CPU `Float64` | CPU `Float32` vs `Float64` |
| --- | ---: | ---: | ---: |
| radial field Kerr | 1.147e-06 | 7.915e-07 | 9.631e-07 |
| radial envelope Kerr | 1.992e-06 | 1.861e-06 | 5.959e-07 |

Two orders inside the 1e-4 the tests assert, and the Metal-vs-`Float64` difference is the
same size as the `Float32`-vs-`Float64` one, i.e. it is single precision and not the device.

## Benchmark

`benchmark/radial.jl`. M1 Pro, Julia 1.13.0, `-t 1`, one FFTW thread, one BLAS thread,
`:estimate`, no wisdom; 20 fixed steps over 1 cm of argon at 1 bar, a 100 fs / 400–2000 nm
`RealGrid` (129 frequency samples per column, so the state is `129 x 1 x N`),
`R = 1 mm`, `w0 = 200 µm`, `boundary=:none`.

| N | | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | --- | ---: | ---: | ---: |
| 64 | Hankel GEMM | 95.7 µs | 48.2 µs | 170.7 µs |
| | right-hand side | 365.9 µs | 208.8 µs | 494.0 µs |
| | one step | 3.460 ms | 2.499 ms | 3.545 ms |
| | propagation | 72.1 ms | 51.6 ms | 93.6 ms |
| 256 | Hankel GEMM | 1.395 ms | 703.6 µs | 287.1 µs |
| | right-hand side | 3.480 ms | 1.858 ms | 773.3 µs |
| | one step | 26.96 ms | 16.79 ms | 5.439 ms |
| | propagation | 556.0 ms | 341.8 ms | 114.6 ms |
| 1024 | Hankel GEMM | 22.23 ms | 11.12 ms | 592.6 µs |
| | right-hand side | 47.68 ms | 24.09 ms | 1.215 ms |
| | one step | 313.8 ms | 169.4 ms | 8.404 ms |
| | propagation | 6.51 s | 3.47 s | 210.6 ms |

Crossover, from the narrowing sweep (`LUNA_BENCH_NRADIAL=96,128,192`), propagation wall
time:

| N | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | ---: | ---: | ---: |
| 64 | 72.1 ms | 51.6 ms | 93.6 ms |
| 96 | 125.9 ms | 86.4 ms | 83.5 ms |
| 128 | 192.0 ms | 126.8 ms | 102.8 ms |
| 192 | 348.3 ms | 221.4 ms | 89.7 ms |
| 256 | 556.0 ms | 341.8 ms | 114.6 ms |
| 1024 | 6.51 s | 3.47 s | 210.6 ms |

**Metal passes the `Float64` host between 64 and 96 radial points, and the `Float32` host
between 96 and 128.** At 1024 points it is 31× the `Float64` host and 16× the `Float32` one,
and the Hankel GEMM alone is 38× faster. This is the first geometry in the project where a
GPU is worth using: the mode-averaged transform has one transverse column and is
launch-bound (`benchmark/device.jl`), and a radial run has one per radial point.

The Metal propagation time is not monotonic in `N` between 96 and 256 (83.5 ms at 96,
102.8 ms at 128, 89.7 ms at 192, 114.6 ms at 256). In that range the run is dominated by
per-step launch and synchronisation overhead, which does not grow with `N`, so the variation
is noise on a roughly constant floor; the trend is clean from 256 up, where the kernels
dominate.

## Known gaps and open questions

- **`Boundaries.clampdecay`, `addloss` and `addloss_k` are host `Float64` code.** They are
  applied to the linear operator by `Boundaries.setup`, which `Luna.run` calls *before*
  uploading the operator, so nothing in Luna's own flow reaches them with a device or
  `Float32` array and the precision is right. A low-level caller who uploads the operator
  first and then calls `Boundaries.setup` by hand gets a `Float64` promotion on the host
  and, on Metal, an error from `similar` (`clampdecay`'s
  `complex(max(real(linop), -ratemax), imag(linop))` promotes `ComplexF32` to `ComplexF64`
  because `ratemax` is a `Float64`). Found while writing `benchmark/radial.jl`, which did
  exactly that (copied from `benchmark/device.jl`, where it is harmless because a modal
  transform has no evanescent channels and `evanescent` returns the operator unchanged). The
  benchmark now keeps both forms of the operator. Not fixed here — it is a one-line
  conversion in three functions, but it is outside this brief and `gpu/23-tabulated-linop`
  is rewriting how the operator reaches the device anyway.
- **The refractive index is uploaded on every `fillnorm!` call.** For a constant index
  (`const_norm_radial`, which is what every radial example and the gate cases use) that is
  once per propagation. For a z-dependent one it is one `(Nω, Npol)` host evaluation and one
  upload per call — the same interim arrangement as `NormModeAvg`'s `β` mirror, and the same
  thing `gpu/23` removes.
- **`TransRadial` does not alias `Pωo` with `Eωo`.** GPU_PLAN.md §4.4 adopts the fork's
  `Pωo === Eωo` aliasing for the free-space transforms, which removes one field-sized buffer.
  It is written there as part of the `TransFree` work; `gpu/21` should do it for all three at
  once, since the argument (the inverse transform consumes `Eωo` and nothing reads it again)
  has to be checked against each transform's own call sequence.
- **No two-component radial gate case.** The reshaped GEMM's only possible rounding
  difference is for two polarisation components, and it is measured directly in
  `test_device.jl` (< 1e-14) rather than end to end. `test_freespace.jl` does run
  two-component radial propagations, and they pass.
- **The `Float32` radial path is not exercised on the regression gate**, only in the tests:
  the gate is a `Float64` contract by construction.
- **Plasma and Raman are not in the Metal radial testset**, only in the `JLArray` one. Their
  kernels are already covered on Metal by `gpu/13`/`gpu/14`'s own testsets (on device blocks
  and through mode-averaged propagations), and what this branch adds is the transform around
  them; a radial Metal plasma run takes about 40 s, which is more than the file's budget
  justifies for a third copy of the same kernel coverage. `metalradialcase` takes `plasma`
  and `raman` keywords, so adding them is a two-line change if a reviewer wants them.

## Deviations from GPU_PLAN.md

- §4.4 says the radial gate cases will move at rounding level. They do not move at all; the
  reasoning and the direct measurement are under "Rounding" above.
- The normalisation is retargeted by `Luna.setup` (`NonlinearRHS.retarget`) rather than
  built for the device by the caller. The plan does not say how a `FreeSpaceNorm`, which is
  a positional argument built before `Luna.setup` is called, learns what device the run uses;
  retargeting keeps every existing low-level script, example and gate case working unchanged
  while still allowing `spec` to be passed to `norm_radial` directly.
- `Pωo === Eωo` aliasing is left to `gpu/21` (above).
