# Thread and FFTW-planning defaults for CPU runs, from a measured study

Base: `gpu/32-docs` (7e8f7d05).

## Motivation

Luna's CPU defaults were FFTW threads = 4 × the Julia thread count for every plan,
BLAS at OpenBLAS's own default (8 here), and `:patient` FFTW planning. A study on an
Apple M1 Pro (`benchmark/threads/THREADS_REPORT.md`, three rounds of independent review)
found each of these costly: up to 3.5× on radial runs, 2.8× on multimode runs, 2× on
small grids, and up to 1.8× on 3-D grids from `:patient` plans.

## Changes

- **FFTW planning default `:measure`** (`Luna.settings["fftw_flag"]`). `:patient` plans
  were up to 1.45× slower for batched and multi-dimensional transforms on the test
  machine; `:measure` was within 5 % of the best mode in every whole run.
- **FFTW threads per plan** (`Utils.FFTWthreads(n)`, `Utils.FFTW_THREAD_MINLEN = 2¹⁷`):
  with `set_fftw_threads(0)` (the default) a plan of at least 2¹⁷ elements gets
  `Threads.nthreads()` threads, anything smaller one, and plans made inside
  `Utils.serial_fftw` (the radial and multimode `setup`s, whose steps are dominated by
  matrix products) one. `Utils.plan_ft` and `Utils.plan_ift` set the count per plan. An
  explicit `set_fftw_threads(n)` applies to every plan and is now honoured with one Julia
  thread too (FFTW then uses its own threads). The 4× multiplier is gone. The wisdom file
  is `FFTWcache_auto` for the automatic counts (FFTW keys its wisdom by thread count) and
  `FFTWcache_<n>threads` for an explicit one.
- **BLAS threads per run** (`Luna.blas_threads`, `Luna.set_blas_threads`,
  `Utils.with_BLAS_threads(f, nthreads)`): `Luna.run` wraps the propagation in
  `with_BLAS_threads`, which sets the count and restores the previous one in a `finally`,
  so also when the propagation throws. Automatic count (`set_blas_threads(0)`, default):
  unchanged on a device, with one Julia thread, and for transforms without matrix
  products; `B₀ ÷ 2` for `TransRadial`; `max(1, B₀ − J)` for `TransModal`/
  `TransModalFixed`, with `B₀` OpenBLAS's default recorded at `__init__`
  (`Utils.BLAS_DEFAULT`). The count used is logged once per run.
- **Regression harness and benchmark scripts** pin BLAS through `set_blas_threads(1)`
  (guarded, for older baseline commits); `benchmark/threads/runs.jl` accepts `0` (auto) in
  `FFTW_LIST`/`BLAS_LIST`, and a `BLAS_LIST` now applies to every case.
- **Docs**: `gpu.md` "On the CPU: threads" rewritten for the automatic defaults, with the
  advice to run `benchmark/threads` on one's own machine for speed-sensitive work; README
  sentence; `set_blas_threads`/`blas_threads` in the API page.
- `THREADS_REPORT.md` moved to `benchmark/threads/`, with an "Implemented" section.

## Tests

- Regression gate against `gpu/int-E2` ced6b299: 506/506, largest difference 0.000e+00 at
  `-t 1` and at `-t 4` (the gate pins `:estimate`, FFTW 1 and now BLAS 1 explicitly).
- CPU, `-t 1` and `-t 8`: `test_utils` (new "thread counts" testset: the FFTW rule and
  its override, `serial_fftw` and `with_BLAS_threads` restoring after a throw, the BLAS
  rule; `threaded_checks.jl` in a `-t 2` child process checks the multi-thread FFTW
  rule), `test_rk45`, `test_linops`, `test_multimode`, `test_freespace`,
  `test_interface`, `test_boundaries`: all pass (the runner's known "FFTW wisdom switch"
  artefact aside).
- `test_device.jl` (`-t 4`) and `test_metal.jl`: pass.
- Smoke, `-t 4`: radial and 4-mode runs log "BLAS threads for the propagation: 4
  (restored to 8 afterwards)" and leave BLAS at 8.

## Measurements

Previous against new defaults, measured side by side in one session
(`benchmark/threads/results/runs_olddefaults.csv`, `runs_defaults.csv`; best of 3;
repeatability ≈ ±5 %), speed-up old/new:

| case | `-t 1` | `-t 4` | `-t 8` |
|---|---:|---:|---:|
| mode-averaged 1 ps | 1.01 | 1.35 | 1.92 |
| mode-averaged 8 ps | 1.01 | 0.98 | 1.03 |
| envelope Raman | 1.01 | 1.50 | 2.09 |
| GNLSE | 1.00 | 1.54 | 2.30 |
| 4 modes `:fixed` | 0.96 | 1.86 | 2.37 |
| 4 modes `:adaptive` | 1.06 | 0.99 | 1.02 |
| radial 256 | 0.99 | 1.61 | 1.43 |
| radial 1024 | 1.04 | 1.35 | 1.13 |
| 2-D χ⁽²⁾ | 1.03 | 1.02 | 1.18 |
| 3-D 16 × 8 | 1.10 | 1.11 | 1.76 |
| 3-D 64 × 64 | 1.09 | 1.41 | 3.79 |

No case is slower beyond the repeatability (worst 0.96, `:fixed` multimode at `-t 1`,
where the thread counts are the same and only the planning mode differs). A first
measurement of the new defaults was discarded: the laptop slept during it (wall times of
test files of up to 398 min), and the benchmark script let a BLAS count set for one case
leak into the next (fixed: a `BLAS_LIST` now applies to every case).

## Not changed (recorded in the report for separate decisions)

- The multimode default (`modal_integral=:adaptive`); `:fixed` with `nr = 128` was as
  accurate and 2.3–2.6× faster at 8 threads in the one strongly nonlinear case measured.
- `Luna.tune_threads()` (per-machine thresholds): later.
- The fixed rule's statistic evaluating the right-hand side once per step even without a
  Kronrod estimate (pinned, next).
