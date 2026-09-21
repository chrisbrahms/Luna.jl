# `gpu/21-free-device`: the Cartesian free-space transforms and collar on a device

Base: `gpu/20-radial-device` at `73e8dba6`.

Cartesian free-space propagation — 2-D (`x`-`z`) and full 3-D, field-resolved and
envelope, Kerr and χ⁽²⁾, with `boundary=:rate` — now runs on Metal and on a `JLArray`,
through the low-level interface. With `gpu/20`'s radial transform, that is every
geometry in Luna except the multimode one.

The default CPU path is bit-for-bit what it was: the regression gate against the base
commit is exactly `0.000e+00` on all 460 rows, including the three free-space cases,
which the brief expected to move.

## Commits

| commit | what |
| --- | --- |
| `70e83f0e` | `test_metal.jl`: the multi-axis FFT plans against FFTW |
| `be2bbbf5` | `TransFree` and `TransFree2D` on the device |
| `e2f3baa2` | The Cartesian collar on the device, and `test_device.jl` |
| `c8d360c9` | `test_metal.jl`: the Cartesian free-space transforms on hardware |
| `9b09ac45` | `benchmark/free.jl`, memory accounting, docs and this file |
| `73b6365c` | PR fix: `test_boundaries.jl` row and the last commit hash |
| `a39604f4` | Review round 1 fixes (this row corrected in the commit after it) |

## The plan test, first

GPU_PLAN.md §2 records that Metal.jl's own tests cover the multi-axis real regions
`(1, 3)` and `(1, 4)` but **not** `(1, 3, 4)`, which is what `TransFree` needs, and the
brief asked for that to be settled before anything was built on it.

`test_metal.jl`, "multi-axis FFT plans on Metal", plans `(1, 3)` and `(1, 3, 4)` real and
complex transforms on an `MtlArray` and compares them with FFTW on the same random data:
forward and inverse, one and two polarisation components, transverse axes which are
powers of two and not, and an odd time axis (where a real inverse plan's output length is
not its input's).

**All four plan kinds are supported and correct.** Agreement with the `Float64` host
transform is 2.0e-7 to 5.3e-7 relative — `Float32` roundoff — and `inv(plan).scale` is
`1/N` over the *output* lengths, which is what `to_time!` assumes when it folds the
normalisation into the oversampling copy. 48 assertions. **No fallback to two plans is
needed**, and none is in the branch.

## What changed

### `NonlinearRHS.TransFree` and `TransFree2D` (`src/NonlinearRHS.jl`)

Both become parametric in their buffer array type, like `TransModeAvg` and `TransRadial`:
every buffer allocated with `Luna.alloc`, the grid vectors mirrored (`Luna.gridvectors`),
the inverse plan held explicitly, the responses given the run's `spec` and `scaling`
through `Nonlinear.rescale_responses`, the noise field uploaded and divided by `Eref`, and
a constructor which asserts residency of every buffer, mirror and response array. Both
carry a `UnitScaling`, and `Luna.runscaling` knows about them, so a `Float32` free-space
run gets a `ScaledOutput` like a mode-averaged one.

**One field-sized buffer fewer: `Pωo === Eωo`.** The oversampled frequency-domain buffer
does double duty. `to_time!` writes the field into it, applies the inverse plan, and
nothing reads it again; `to_freq!` then writes the nonlinear polarisation into the same
array. The two never appear as the input and the output of the same FFT call, which a
device plan would reject. Each Cartesian transform holds **three** field-sized arrays
(`Eto`, `Pto`, and the shared frequency-domain buffer) where it held four.
`test_device.jl` asserts the count by identity, and five when the modified shot-noise
model adds `Et_noise` and `Et_nl`.

**The frequency-domain normalisation** was
`nl .*= ωwin .* (-im.*ω) ./ (2 .* normfun(z))`, which rebuilt two host vectors on every
right-hand side. It is now the same fused `fsnorm!` broadcast over a precombined, mirrored
`prefac = ωwin·(-iω)·Pref` that `gpu/20` gave `TransRadial`.

**One body for both.** The per-step sequence is identical in 2-D and 3-D — they differ
only in how many transverse axes they have — so it is written once, as
`NonlinearRHS.freetransform!`. That is the `TODO: this can probably be combined with the
case for TransFree` the 2-D call operator carried. They stay separate types because
`Boundaries.spacegrid` and `Luna.setup` dispatch on them and because the FFT region
differs.

**The `scale` field and constructor argument are gone.** Nothing read them except the 3-D
noise setup, which now goes through `to_time!` and computes the same factor itself.

### `Luna.setup` (`src/Luna.jl`)

The four Cartesian free-space methods (two grid types × 2-D/3-D) become four forwarding
methods and one `setup_free`, which takes `device` and `precision`, logs the choice,
retargets the normalisation (`NonlinearRHS.retarget`, `gpu/20`), plans the oversampled and
state transforms on the run's array type through `Utils.plan_ft`, keeps host plans for
building the input fields (`Fields` is host scalar code), computes the unit scaling from
the peak of the input taken back to `(t, r⊥)` by the joint inverse transform, and uploads
the initial field. The transverse shape and the FFT region are the only things that differ
between the two geometries (`freeshape`, `freeregion`, `freetransform`).

Behaviour on the default CPU path is unchanged, including which FFTW plans are made, in
what order, with what flags.

### `Boundaries.CartesianCollar` (`src/Boundaries.jl`)

**No code change was needed.** `gpu/11` had already written it as one broadcast over a
rate mirrored with `Luna.upload_like`, one `sum(abs2, ·)` and one `mapreduce` over a lazy
`Broadcasted`, with a residency assertion at construction; its docstring said it was tested
on the host only because the Cartesian transforms had no device path. It is now the
transverse boundary of a device free-space run and says so. Unlike the radial collar it
does no transform: the Cartesian grids transform time and space together, so when
`RateAbsorber` applies the temporal collar the state is already in `(t, x[, y])`.

### `Interface.jl`

Comment only: `prop_capillary`/`prop_gnlse` never build a free-space run, so all three
free-space device paths are reached through the low-level interface. `_cpu_only!` is now
about multimode propagation alone.

## A bug found and fixed

**`TransFree`'s modified shot-noise path could not have worked.** Its buffers were built
without the polarisation axis —

```julia
Eωo_noise = zeros(ComplexF64, (length(grid.ωo), Nx, Ny))   # 3-D
Et_noise  = zeros(TT, (length(grid.to), Nx, Ny))           # 3-D
ldiv!(Et_noise, FT, Eωo_noise)                             # FT is the 4-D (1,3,4) plan
```

— so `ldiv!` against the 4-D plan and the later `@. t.Et_nl = t.Eto + t.Et_noise`
broadcast were both shape mismatches. Nothing in the test suite, the gate or the examples
passes `noise_field` to a Cartesian free-space setup, so it had never been run.

Since the constructor was being rewritten anyway, the path is now the `TransRadial`
arrangement: an `(nω, npol, nk...)` noise field taken to the oversampled real-space time
domain by the transform's own inverse plan, with an explicit shape check and a message
naming the expected shape. `TransFree2D`, which had no noise support at all, has the same.
Both are tested on `JLArray`.

## Rounding

The brief expected the free-space gate cases (`free3d_env_kerr`, `free2d_field_chi2`,
`free2d_env_chi2`) to move at rounding level from the buffer aliasing or from plan
batching. **They did not move at all**, and none of the three candidate sources can move
them:

- **The aliasing changes no arithmetic.** `to_time!` starts with `fill!(Aωo, 0)` and then
  `copy_scale!`, so what the inverse plan sees does not depend on what was in the buffer;
  and nothing reads the field's copy after `to_time!` returns.
- **The plans are the same plans.** `Utils.plan_ft` on a host array is
  `FFTW.plan_rfft`/`plan_fft` with `settings["fftw_flag"]` — the same call the four
  `setup` methods made — over the same regions, on buffers of the same shapes, created in
  the same order.
- **`fsnorm!` keeps the association** of the expression it replaces (`pre ./ (2 .* norm)`)
  and `Pref` is `1.0` at `Float64`, so `prefac` is the same vector the per-step expression
  used to rebuild.

This is the same argument `gpu/20` made for the radial cases, and the gate confirms it.

## Regression gate

M1 Pro, Julia 1.13.0, `-t 1`, `set_fftw_mode(:estimate)`, one FFTW thread, one BLAS
thread, wisdom off, `shotnoise=false`, through the timing wrapper. 21 cases with
`LUNA_REGRESSION_SKIP=radial_field_raman` (its baseline is the integrator's, as
`PR_20-radial-device.md` records).

| baseline | cases | result | largest difference |
| --- | --- | --- | --- |
| `fa556e6f` (source-identical to `73e8dba6`, this branch's base) | 21 | **460 pass, 0 fail** | `0.000e+00` on every row |
| `fdf8dbe3` (`evanescent`) | 21 | **460 pass, 0 fail** | 7.240e-05 (`multimode_field_plasma` adaptive statistics) |

Against the base, **every row of every case in both modes is exactly `0.000e+00`**. No
case moved, no step count moved.

Against `evanescent` the record is `gpu/20`'s reproduced to every digit, with no new rows
— the same sixteen non-zero rows (eight cases in two modes), the same values:

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

`free3d_env_kerr` is `0.000e+00` against both baselines, as it was on the base.

## Tests

All on M1 Pro / Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread,
wisdom off, through the timing wrapper. The `JLArrays` and `Metal` environments are
temporary stacked ones with `Hankel` pinned to the worktree's 0.5.9.

| what | result | before |
| --- | --- | --- |
| `test/test_regression.jl` × 2 | 460/0 and 460/0 — the record above | |
| `test/test_device.jl` (JLArrays) | **783 pass, 0 fail**, 44 testsets | 687 / 38 |
| `test/test_metal.jl` (Metal 1.11.1, hardware) | **567 pass, 0 fail**, 25 testsets | 452 / 20 |
| `test/test_freespace.jl` | **77 pass, 0 fail** | 77 pass |
| `test/test_chi2.jl` | **42 pass, 0 fail** | 42 pass |
| `test/test_boundaries.jl` (with JLArrays in the environment) | **184 pass, 0 fail** | 184 pass |
| `include("docs/make.jl")` | 8 unresolved `@ref`s, one **fewer** than the base's 9 | 9 |
| the free-space and radial examples | 12 of 13 run (below) | |

`docs/make.jl` was built on this HEAD and on `73e8dba6` in a temporary worktree. The base
has nine unresolved `@ref`s; this branch has the same eight of them and no new ones. The
one that goes away is `[`norm_free`](@ref)` in `docs/src/modules/NonlinearRHS.md`, which
came from the old `TransFree2D` docstring this branch rewrote.

The thirteenth example is `examples/low_level_interface/freespace/chi2_benchmarking.jl`,
which fails on `@profview` (a `ProfileView` macro that is not in the environment).
Pre-existing and unrelated; the other five 2-D/3-D free-space examples and all five radial
ones run.

### What was added

`test/test_device.jl` (all against the host on `JLArray`, `allowscalar(false)`):

- `the free-space transforms hold three field-sized buffers` — the buffer-count
  assertion, by `objectid`, for all three transform/grid combinations, plus five with a
  noise field.
- `2-D free-space χ⁽²⁾ on JLArray` — type I SHG in BBO on a `Grid.Free2DGrid`,
  field-resolved and envelope, `boundary=:rate`, fixed steps. The same case as
  `test_freespace.jl`'s "BBO SHG" testset and the `free2d_*_chi2` gate cases, on a smaller
  grid. Exercises the region-`(1, 3)` plans, the two-component χ⁽²⁾ response, the
  crystal-optics normalisation staged through the host, the k-space absorber, the
  evanescent source taper and `CartesianCollar`, with residency assertions on every buffer,
  mirror and the retargeted normalisation.
- `3-D free-space Kerr on JLArray` — the `free3d_env_kerr` geometry on a coarser grid;
  region-`(1, 3, 4)` plans and the isotropic `FreeSpaceNorm` broadcast (which needs no
  host staging buffer, unlike the crystal-optics one).
- `the Cartesian collar on JLArray` — both transverse shapes against the host, including
  the absorbed-energy bookkeeping.
- `the 2-D free-space noise field on JLArray` — the upload, the device inverse plan, the
  field+noise buffer, and that a noise field of the wrong shape is refused.

Measured, per save, normalised by the largest `|Eω|` in that save (the gate's metric):

| case | JLArray vs host |
| --- | ---: |
| 2-D free-space field χ⁽²⁾ (BBO) | **0.000e+00** |
| 2-D free-space envelope χ⁽²⁾ (BBO) | **0.000e+00** |
| 3-D free-space envelope Kerr | **0.000e+00** |
| 2-D free-space field χ⁽²⁾ with a shot-noise field | **0.000e+00** |
| 3-D free-space envelope Kerr with a shot-noise field | **0.000e+00** |

Exactly zero, where the radial cases are at 1e-16. There is no GEMM in a Cartesian
right-hand side, and the JLArray FFT shims wrap host plans, so every operation is the same
host arithmetic in the same order — the device path differs from the host one only in
where the arrays live.

`test/test_metal.jl`:

- `multi-axis FFT plans on Metal` — above.
- `no stray Float64 in the Cartesian free-space kernels` — the element type of every array
  the free-space normalisation and the Cartesian collar hold; the isotropic fill compiled
  and run on the GPU against the host on a grid fine enough (`R = 10 µm`, 32 × 16,
  400–4000 nm) that the evanescent branch — the one carrying the taper — is reached; the
  lazy k-window mirror rebuilt after `reflength!`; the crystal-optics fill with its host
  `ohost` staging buffer; and `apply_realspace!` on a device state in both transverse
  shapes.
- `2-D free-space χ⁽²⁾ on Metal`, `3-D free-space Kerr on Metal`, `free space on Metal
  against the Float64 CPU path`.

The Metal references are **explicit host specs** (`DeviceSpec(Array, Float32)` and
`HostSpec()`), not the `device` sentinel, which with Metal loaded resolves to the GPU and
would compare Metal with itself — the `gpu/11` review's finding 1. Fixed steps throughout,
so the comparison is of arithmetic and not of step sequences.

Measured:

| case | Metal vs CPU `Float32` | Metal vs CPU `Float64` | CPU `Float32` vs `Float64` |
| --- | ---: | ---: | ---: |
| 2-D free-space field χ⁽²⁾ (BBO) | 3.179e-06 | 2.117e-06 | 1.251e-06 |
| 2-D free-space envelope χ⁽²⁾ (BBO) | 3.277e-06 | 2.711e-06 | 8.442e-07 |
| 3-D free-space envelope Kerr | 1.975e-06 | 1.980e-06 | 3.021e-07 |

Between one and two orders inside the 1e-4 the tests assert, and the Metal-vs-`Float64`
difference is the same size as the `Float32`-vs-`Float64` one, i.e. it is single precision
and not the device.

One test-only fix: the k-window passed to `reflength!` is clamped at `exp(-MAX_αℓ/2)`
exactly as `Boundaries.setup` clamps it. The raw Planck taper `Boundaries.kprofile`
returns is exactly zero at the Nyquist wavevector of an FFT grid, and the normalisation
divides by it — `Boundaries.setup` has always clamped it for that reason; the first draft
of the test did not, and got `Inf`.

## Memory

Counted from the shapes (`scratchpad` script, reproduced in
`docs/src/developer/device_model.md`) for the 3-D example
`examples/low_level_interface/freespace/full3D.jl`: field-resolved, 400–2000 nm, 0.2 ps
(`nt = 512`, `nto = 1024`, `nω = 257`, `nωo = 513`), 128 × 128 transverse, one
polarisation, in `Float32`:

| item | Kerr | Kerr + plasma | Kerr + Raman |
| --- | ---: | ---: | ---: |
| transform `Eto`, `Pto` | 128 MB | 128 MB | 128 MB |
| transform `Eωo === Pωo` | 64 MB | 64 MB | 64 MB |
| normalisation `out` | 32 MB | 32 MB | 32 MB |
| stepper (`y`, `yn`, `yi`, `yerr`, 7 stages) | 353 MB | 353 MB | 353 MB |
| linear operator | 32 MB | 32 MB | 32 MB |
| absorber `Et` | 32 MB | 32 MB | 32 MB |
| response block buffers | — | 256 MB | 384 MB |
| **total** | **0.63 GB** | **0.88 GB** | **1.00 GB** |

**Field-sized buffers per transform: four before, three after** (`Eto`, `Pto`, `Eωo`,
`Pωo` → `Eto`, `Pto`, `Eωo === Pωo`); five with the modified shot-noise model, which
was six. The aliasing is worth 64 MB of 0.88 GB, about 7 %, on the run above, and a larger
fraction of a transform-dominated one.

**A 3-D plasma run at the examples' sizes fits a 16 GB Metal device eighteen times over**,
so the response block is **not** chunked. Two reasons beyond the headroom:

- The total scales as `Nx·Ny`: 0.88 GB at 128 × 128 is 3.5 GB at 256 × 256 and 14.0 GB at
  512 × 512, so 16 GB runs out somewhere around 512 × 512 — an enormous 3-D grid.
- Chunking would move that limit by much less than it looks. The batched response buffers
  are 29 % of the total and the stepper's eleven state-sized arrays, which cannot be
  chunked without changing `RK45`, are 39 %. Chunking the response block along the
  transverse axes would buy about 1.4× in grid area for a per-step gather/scatter on every
  backend, including the host.

A `Float64` host run is exactly twice these numbers (1.75 GB for the plasma case).

**Setup costs a little more than the propagation holds** (review 1, findings 4 and 5), all
of it collectable once `Luna.setup` returns: on the host, the `Float64` prototypes the
input-field plans are made against (`(nt, npol, Nk...)` and `(nt, 2, Nk...)`, 64 MB and
128 MB at this grid) and the initial state before it is uploaded (64 MB); on the device,
one state-shaped time-domain block for the state's own plan (32 MB). The *oversampled*
block is no longer among them: `setup_free` hands the array it planned `FTo` against to
the transform as its `Eto` rather than leaving it to the garbage collector, which is 64 MB
less peak memory at setup on this grid. `setup_radial` still allocates its planning
prototype and drops it — that is `TransRadial`, out of scope here.

## Benchmark

`benchmark/free.jl` (new, modelled on `benchmark/radial.jl`). M1 Pro, Julia 1.13.0,
`-t 1`, one FFTW thread, one BLAS thread, `:estimate`, no wisdom; 10 fixed steps over 1 cm
of argon at 1 bar, a 100 fs / 400–2000 nm **envelope** grid (`nω = 128`), `R = 1 mm`,
`w0 = 200 µm`, `boundary=:none`. `LUNA_BENCH_NFREE` overrides the sweep.

| transverse grid | | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | --- | ---: | ---: | ---: |
| 32 × 32 | joint inverse FFT | 1.148 ms | 916.5 µs | 313.4 µs |
| | right-hand side | 3.680 ms | 2.747 ms | 465.6 µs |
| | one step | 46.15 ms | 39.77 ms | 3.297 ms |
| | propagation | 477 ms | 412 ms | 56.0 ms |
| 64 × 64 | joint inverse FFT | 5.246 ms | 3.887 ms | 424.3 µs |
| | right-hand side | 16.21 ms | 11.52 ms | 851.8 µs |
| | one step | 195.0 ms | 164.1 ms | 7.413 ms |
| | propagation | 2.03 s | 1.69 s | 109 ms |
| 128 × 128 | joint inverse FFT | 36.05 ms | 17.52 ms | 1.057 ms |
| | right-hand side | 93.23 ms | 50.47 ms | 2.688 ms |
| | one step | 1.004 s | 700.3 ms | 25.93 ms |
| | propagation | 10.56 s | 7.27 s | 352 ms |

**There is no crossover to report.** The smallest grid in the sweep already has 1024
transverse columns — past where the radial transform crosses over.

The figure to quote is the **per-step** one, which `BenchmarkTools` repeats many times and
which reproduces to about 3 %: Metal is **14×** the `Float64` host at 32 × 32, 26× at
64 × 64 and **39×** at 128 × 128. End to end the propagation is 8.5×, 18.7× and 30.0×
faster — lower because a propagation also does its setup, its output and, on a device, the
host copy of each saved field, none of which the GPU helps with. The gap is the FFT: at
128 × 128 the joint inverse transform alone is 34× faster on the GPU. This is the geometry
the GPU work was for.

**On the scatter of the `prop` column** (review 1, finding 3). It is a single `@elapsed`
per sample, not a `BenchmarkTools` loop, so it pays for compilation on the first run of a
process, for the device's first-touch graph and kernel caching, and for whatever the
garbage collector does among the fresh field-sized buffers each sample allocates. Review 1
measured a factor of **2.2** between two runs of the 64 × 64 Metal row with the original
`samples=2` and no warm-up. `proptime` now discards a warm-up run and takes the minimum of
five (`LUNA_BENCH_PROPSAMPLES` overrides it). Re-measured: two full sweeps and an
independent `LUNA_BENCH_NFREE=64` run agree to within 10 % on every row — Metal at
64 × 64 is 108.8 ms, 108.8 ms and 121.8 ms; CPU `Float64` 2.03 s, 2.05 s and 2.23 s — so
the column is worth two figures and no more. The table above is the second full sweep;
the `fft`, `rhs` and `step` columns agree between all runs to a few per cent.

## Documentation

- `docs/src/developer/device_model.md`: a new "The Cartesian free-space transforms"
  section — the per-step sequence, the multi-axis regions and what is known about them on
  Metal, the `Pωo === Eωo` aliasing and why it is safe, the memory table and the argument
  against chunking, the transverse collar, and the benchmark. The testing list names
  `benchmark/free.jl`.
- `docs/src/gpu.md`: the "work in progress" admonition now says everything but multimode
  propagation is device-capable; "What runs where" says the same and shows how to build a
  Cartesian free-space run; "Performance" gains the free-space benchmark table and a
  "Memory" subsection.

## Known gaps and open questions

- **`TransRadial` still does not alias `Pωo` with `Eωo`.** `PR_20-radial-device.md` left
  the aliasing to this branch "for all three at once"; the brief for this branch lists
  `TransRadial` under "do not touch ... beyond the shared norm", so it is not done here.
  The argument holds for it unchanged — `to_time!` is the only reader of `Eωo` and
  `to_freq!` writes `Pωo` afterwards — and the change is one constructor line plus the
  struct field. It would save one `(nωo, npol, nr)` complex array per radial run and is a
  natural follow-up on `gpu/int-E`.
- **`Boundaries.clampdecay`, `addloss` and `addloss_k` are still host `Float64` code**,
  as `PR_20` records for `gpu/23-tabulated-linop`. Unreachable through `Luna.run`, which
  calls `Boundaries.setup` before uploading the operator; `benchmark/free.jl` builds its
  hand-made stepper from the clamped host operator for the same reason `benchmark/radial.jl`
  does.
- **The crystal-optics normalisation is host scalar code**, by design (GPU_PLAN.md §4.4):
  a root-find per `(ω, kx)` for the internal angle. It runs once per propagation for a
  `const_norm_free2D`/`const_norm_free`, and its result is uploaded through `ohost`. For a
  z-dependent crystal index it is one host fill and one upload per right-hand side, which
  would be worth batching if anyone runs one.
- **The refractive index is uploaded on every `fillnorm!` call** for the isotropic fill —
  `gpu/20`'s note, unchanged; once per propagation for a constant index.
- **Plasma and Raman are not in the Metal free-space testsets**, only χ⁽²⁾ and Kerr.
  Their kernels are covered on Metal by `gpu/13`/`gpu/14`, and on multi-column blocks by
  `gpu/20`'s radial JLArray testsets; what this branch adds is the transform around them.
  The `JLArray` free-space tests do not run them either, which is the gap — a 3-D plasma
  case would be the most expensive test in either file.
- **The `Float32` free-space path is not exercised by the regression gate**, only by the
  tests: the gate is a `Float64` contract by construction.
- **No `npol = 2` Cartesian gate case.** The two χ⁽²⁾ gate cases are two-component 2-D
  runs, so the shape is covered there; there is no two-component *3-D* case anywhere
  except `test_freespace.jl`'s host-only `pol = true` sweep.
- **The new `noise_field` keywords are tested on `JLArray` but not on Metal.** Both
  dimensionalities are covered there (2-D and, since review round 1, 3-D, which is the one
  the branch repaired); the gate runs with `shotnoise=false` by construction.

## Review round 1

`scratchpad/reviews/gpu-21-free-device-1.md`, verdict **approve with minor fixes**. Every
load-bearing claim was independently reproduced — both gates (460/0, `0.000e+00` on every
row against the base; the `evanescent` record to every digit), `test_device` 775,
`test_metal` 567, `test_freespace` 77, `test_chi2` 42, `test_boundaries` 184, the
Metal-vs-CPU agreements to every digit, and the memory arithmetic — **except the
benchmark's `prop` column**, which scattered by 2.2× between two of the reviewer's own
runs. The eight findings are addressed in one follow-up commit:

| # | finding | fix |
| --- | --- | --- |
| 1 | the 3-D shot-noise path, the one this branch repaired, had no test (the noise testsets and the wrong-shape check were all 2-D) | `free3dcase` gains `noise_field`; new testset `the 3-D free-space noise field on JLArray` — `Et_noise` resident, `(nto, 1, 8, 6)`, bit-identical to the host's, one device right-hand side finite, the propagation at 0.000e+00, and a noise field without the polarisation axis refused |
| 2 | `Luna.run`'s comment still named two of the four transforms that carry a scaling (stale twice over: review 20's finding 4 was the other half of the same sentence) | names all four, and `TransModal` as the exception |
| 3 | the `prop` column does not reproduce, and three speed-up figures in the PR and both docs pages rest on it | `proptime` discards a warm-up run and takes the minimum of five (`LUNA_BENCH_PROPSAMPLES`); re-measured over two sweeps and a single-size run, now within 10 %; the quoted ratios come from `step`, which reproduces to 3 %, with the end-to-end ones given alongside and a paragraph on the scatter in the PR, in `benchmark/free.jl`'s docstring and in `device_model.md` |
| 4 | `setup_free` allocated two field-sized blocks only to plan against them | the oversampled one is handed to the transform as its `Eto` (`freebuffers` takes it over, checks its shape and zero-fills it); the state-shaped one is documented |
| 5 | the memory table was device-side only | a paragraph on what setup allocates on the host and on the device, in the PR and both docs pages |
| 6 | `freeprefac` duplicated the expression `TransRadial` inlined | one `fsprefac`, documented, next to `fsnorm!`, used by all three free-space transforms |
| 7 | the `gpu.md` sentence introducing the radial table had been edited to describe the Cartesian one, and was over-width | rewritten to introduce the radial table; the two over-width lines this branch added are wrapped |
| 8 | `test_device.jl`'s "backend trait" testset errors if a GPUArrays-backed package is loaded before the file (pre-existing) | the device wrapper instances are constructed directly (`SubArray`, `Base.ReshapedArray`) instead of through `Base.view`/`Base.reshape`, which `GPUArrays` takes over for `AbstractGPUArray`; verified both ways round — 783/0 with `using JLArrays` first and with it loaded by the file |

Re-run after the fixes: gate against `fa556e6f` (21 cases) **460 pass, 0 fail,
`0.000e+00` on every row**; `test_device.jl` **783 pass, 0 fail, 44 testsets** (and the
same with `JLArrays` loaded first, which used to be 773 / 2 errored);
`test_freespace.jl` **77 pass, 0 fail**; `test_chi2.jl` **42 pass, 0 fail**;
`test_boundaries.jl` **184 pass, 0 fail**; `benchmark/free.jl` re-measured (the table
above is the new run).

The reviewer's two carried-forward notes stand: `TransRadial` still holds six field-sized
buffers and does not alias `Pωo` with `Eωo` (out of scope here, two lines on `gpu/int-E`),
and the `Boundaries.clampdecay`/`addloss`/`addloss_k` `Float64` promotion is still
`gpu/23`'s.

## Deviations from GPU_PLAN.md

- §4.4 and the brief expected the free-space gate cases to move at rounding level from the
  aliasing or from plan batching. They do not move at all; the reasoning and the gate are
  under "Rounding" above.
- §2 flags `(1, 3, 4)` on Metal as untested and the brief asked for a fallback to two
  plans `(1,)` then `(3, 4)` if it were wrong or unsupported. It is neither, so there is
  no fallback and no branch in the code; the measurement is the first testset in the
  free-space part of `test_metal.jl`.
- The brief asked for chunking of the response block "if a 3-D plasma run at the examples'
  sizes does not fit a 16 GB Metal device, or say why not". It fits with 18× headroom, so
  there is no chunking; the numbers and the argument are under "Memory".
- `TransFree`/`TransFree2D` lose their `scale` constructor argument and field (dead once
  the noise path uses `to_time!`), and `TransFree2D` gains a `noise_field` keyword it did
  not have. Both are changes to the low-level constructor signatures, which only
  `Luna.setup` uses.
