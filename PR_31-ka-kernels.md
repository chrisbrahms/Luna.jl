# Threaded pointwise kernels on the CPU

Base: `gpu/30-neff-grid` (5b38450f). Plan §4.9 and §6 Group F (`gpu/31-ka-kernels`).

## Motivation

The plan made this branch conditional on profiling showing that broadcasts dominate the
CPU step. They do. Sampling profile (`benchmark/threaded/profile.jl`, 1 thread, FFTW/BLAS 1 thread,
Kerr + PPT plasma in argon; share of samples):

| case | broadcasts | scans | FFT | GEMM | largest broadcast callers |
|---|---|---|---|---|---|
| mode-averaged field, 1 column | 52 % | 8 % | 11 % | – | propagator `exp` 16 %, ionisation rate 11 %, plasma block 11 % |
| 4-mode `:fixed`, 16 points | 69 % | 15 % | 2 % | 5 % | plasma block 25 %, ionisation rate 25 %, scan correction 9 % |
| radial, 256 points | 35 % | 4 % | 5 % | 46 % | propagator 9 %, plasma block 8 %, rate 5 % |

Every one of these broadcasts runs on one thread in Base. Only the plasma response's
column loop (multi-column blocks) was threaded before.

## Design

`Utils.threaded(dest; minlen)` wraps a broadcast destination:

```julia
@. $(Utils.threaded(y; minlen=Utils.THREAD_MINLEN_HEAVY)) = y * exp(linop*dt)
```

On a host array longer than `minlen`, with more than one Julia thread, it runs the same
`Broadcasted` object over slabs of `dest` on up to `Threads.nthreads()` tasks; otherwise it is
exactly the plain broadcast (on a device it always is). It does what Base's `copyto!` does
before its loop (axes fixed to `dest`, unalias/extrude via `Broadcast.preprocess`), so an
argument that aliases `dest` behaves as in Base. **Each element is computed by the same
scalar code whichever task runs it, so results are bit-identical to the serial broadcast
for any thread count** — the regression gate shows exactly 0 at 1 and at 4 threads.

A nesting guard (`Utils.serial_region`, a process-wide atomic depth counter) keeps
broadcasts serial inside an already-threaded region; the plasma column loop uses it.
`Luna.set_threaded_broadcasts(false)` turns the whole thing off (`settings`), for
benchmarking; it changes speed only.

### Why not KernelAbstractions

The plan named KernelAbstractions. Measured on this machine (`benchmark/threaded/executors.jl`), a
KA `@kernel` over the same `Broadcasted` on `CPU()` is no faster than `Threads.@spawn` over
chunks (complex `exp` at 2²⁰ elements: serial 12.4 ms, KA 3.31 ms, `@spawn` 3.40 ms at 4
threads; 2.81 / 2.61 ms at 8), and KA would be a new hard dependency of core Luna, which §3
excludes (only `Adapt`/`GPUArraysCore`). The kernel body is the broadcast expression, as
before; only the host executor is new, so §4.2 rule 2 (one implementation per operation)
holds.

### Thresholds

Task spawning costs ≈15–25 µs, more at 8 threads, where equal chunks also wait for the
M1's efficiency cores. The first version used 2¹⁶/2¹² (from the microbenchmark alone) and
slowed a short mode-averaged run to ×0.67 at 8 threads; the thresholds are now
`THREAD_MINLEN = 2¹⁷` (cheap arithmetic) and `THREAD_MINLEN_HEAVY = 2¹⁵` (complex `exp`,
spline evaluation), with at least `minlen ÷ 4` elements per task. Chunks are slabs along the
last non-singleton dimension iterated as `CartesianIndices` blocks; a linear index per
element (the first version) cost an integer division per dimension and made the
multi-column 4-mode case ≈7 % slower.

## Call sites

| site | threshold |
|---|---|
| `RK45.make_prop!(::AbstractArray)` and the four/two propagator broadcasts of `LinearOps._make_prop!` | heavy |
| `Ionisation.ionrate!` (device-capable rates) | heavy |
| plasma ionisation fraction `pf + 1 - exp(-fraction)` | heavy |
| plasma polarisation, output accumulate, loss term (`_plasma_loss!`, both forms) | default |
| `Maths.cumtrapz_scan!` correction broadcast | default |
| `RK45.combine!` (all 8 forms), `RK45.errorestimate!` | default |
| `NonlinearRHS._materialise!` (fused pointwise responses: Kerr, χ⁽²⁾, THG) | default |
| `NonlinearRHS.fsnorm!` (free-space normalisation) | default |

`y *= e` forms are rewritten as `y = y * e` (the same operation). The prefix scan itself
(`accumulate!`) stays serial per column; multi-column blocks already thread it through the
plasma column loop.

## Measurements

`benchmark/threaded/speed.jl`, M1 Pro (8 performance + 2 efficiency cores), best of two
runs per setting in one process, `Luna.set_threaded_broadcasts` off vs on; every output
bit-identical between the two. Cases (`benchmark/threaded/cases.jl`): `modeavg` 0.5 m
capillary, 1 ps window; `modal` 4 modes `:fixed`; `radial` 256 points; `modeavg_long` 8 ps
window (nt ≈ 2¹⁶); `radial_big` 1024 points. All Kerr + PPT plasma in argon.

| case | `-t 1` | `-t 4`, FFTW/BLAS 1 | `-t 8`, FFTW/BLAS 1 | `-t 8`, FFTW/BLAS 8 |
|---|---|---|---|---|
| modeavg | 1.00 | 1.00 (2.14 s) | 1.00 (2.15 s) | 1.00 (3.27 s) |
| modal | 1.00 | 0.99 (2.32 s) | 0.96 (1.75 s) | 0.98 (5.47 s) |
| radial | 1.00 | 1.09 (1.71 s) | 1.11 (1.65 s) | 1.40 (2.57 s) |
| modeavg_long | 0.99 | 1.26 (6.29 s) | 1.18 (6.69 s) | 1.24 (6.02 s) |
| radial_big | 1.00 | 1.05 (17.8 s) | 1.05 (17.6 s) | 1.23 (9.40 s) |

Speed-up = off/on (on-time in brackets). The gains are modest: the threaded broadcasts are
only part of the step, a mode-averaged state is mostly below the thresholds, and the
multimode case already threads its plasma columns. The one measured slowdown is `modal` at
`-t 8` (×0.96, one run; ×0.99 at `-t 4`). Log: `scratchpad/F31/speed_final.log`, one run
per configuration after the chunking change (each speed-up is the ratio of the best of two
timings per setting).

**Aside, pre-existing, not changed here:** FFTW's default thread count (`fftw_threads = 0`
means 4 × `Threads.nthreads()`) is harmful for these grid sizes: `modal` takes 5.4 s with
8 FFTW/BLAS threads against 1.7 s with 1, `modeavg` 3.3 s against 2.2 s. Only
`radial_big` gains from FFTW/BLAS threads (9.4 s against 17.6 s). Recorded for the user
page in `gpu/32-docs`.

## Tests

- `test/test_utils.jl`: new `threaded_checks.jl` (in-place with `dest` as an argument, a
  broadcast along other axes, a scalar right-hand side, a view destination, an aliasing
  shifted view, the nesting guard, the setting); run in process when Julia has several
  threads, otherwise in a child process with `-t 2`, so `Pkg.test()` exercises it.
  8/8 at 4 threads (re-run after the chunking change: 8/8), 1/1 (child) at 1 thread.
- CPU, 1 thread: `test_rk45` 51, `test_linops` 364, `test_ionisation` 37, `test_multimode`
  15, `test_vectorplasma` 2, `test_freespace` 77 — all pass.
- CPU, 4 threads: `test_rk45`, `test_ionisation`, `test_multimode`, `test_vectorplasma`,
  `test_freespace`, `test_kerr` — all pass (before the chunking change, and re-run after
  it: all pass; `test_utils.jl` 40/41 with the runner artifact below).
- `test_device.jl` at 4 threads: 1432/1432. `test_metal.jl`: all pass.
- Regression gate against ced6b299. Most of its 23 cases are mode-averaged or GNLSE with
  arrays below both thresholds, so at `-t 4` it mostly exercises the serial path; the
  stronger bit-identity evidence is the off/on comparison of `speed.jl` on the five larger
  cases above (identical outputs in every configuration) and `threaded_checks.jl`
  (`minlen=16`). Result: 506/506, largest difference 0.000e+00 at `-t 1` and at
  `-t 4` (FFTW/BLAS pinned to 1 thread by the gate); re-run at `-t 4` after the chunking
  change: 506/506, 0.000e+00 (177 s). `test_device.jl` (4 threads, 370 s) and `test_metal.jl`
  (394 s) re-run after it: all pass.
- `test_utils.jl`'s "FFTW wisdom switch" fails only under the implementer's runner, which
  disables wisdom before including the file; it passes under `Pkg.test()`.

## Known gaps

- The per-column prefix scan is serial for single-column (mode-averaged) blocks; a parallel
  scan would change the summation order and so the results.
- `TransModeAvg`'s own broadcasts, the boundaries and the statistics are not wrapped (each
  ≤ 2 % in the profile, and mode-averaged arrays are mostly below `THREAD_MINLEN`).
