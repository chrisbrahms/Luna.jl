# Julia, FFTW and BLAS threads on the CPU: measurements and a proposed default

Branch `gpu/33-threads` (base `gpu/32-docs` 7e8f7d05). Machine: Apple M1 Pro, 8 performance
+ 2 efficiency cores, Julia 1.13, FFTW.jl with `fftw_jll`, OpenBLAS (default 8 threads).
Scripts and raw CSVs in `benchmark/threads/`; `summarise.jl` prints the tables.

## How the three pools interact

- **Julia threads** (`julia -t J`): Luna's threaded broadcasts and the plasma column loop.
- **FFTW threads**: with `J > 1`, FFTW.jl runs FFTW's threads as Julia *tasks* on the same
  pool (`Threads.@spawn`), so an FFTW count above `J` only makes finer tasks; with `J = 1`
  FFTW uses its own pthreads, but Luna currently forces 1 there (`Utils.FFTWthreads`). The
  count is fixed when a plan is made. Luna's default is `4J`.
- **BLAS threads** (OpenBLAS's own pthreads, default 8 here): the radial Hankel GEMM and
  the multimode synthesis/projection. They compete with the Julia pool for cores.

Every configuration below gave **bit-identical** output to the single-threaded run of the
same case (`maxreldiff` column, all 0).

## Microbenchmarks (`fft.jl`, `blas.jl`)

FFT alone, `:estimate` plans, time with the best FFTW count / time with one:

| transform | elements | J=1 (pthreads) | J=4 | J=8 | J=10 |
|---|---:|---:|---:|---:|---:|
| real 1-D 8192 | 8k | 1.00 | 1.00 | 1.00 | 1.00 |
| complex 1-D 16384 | 16k | 0.78 | 0.46 | 0.55 | 0.93 |
| real 1-D 65536 | 64k | 0.69 | 0.60 | 0.81 | 0.82 |
| real 1-D 262144 | 256k | 0.33 | 0.43 | 0.40 | 0.42 |
| real 4096 × 32 batched | 128k | 0.39 | 0.33 | 0.25 | 0.53 |
| complex 1024 × 64 × 64 | 4M | 0.14 | 0.26 | 0.14 | 0.15 |

In isolation threading pays from ≈2¹⁵ elements, best near `J` FFTW threads (the `2J`/`4J`
multipliers never win by more than noise), and 10 Julia threads is worse than 8 (the two
efficiency cores hold back the slowest task). `:patient` planning time grows steeply with
the FFTW count: 0.44 s → 6.1 s → 28 s for a 2¹⁸ complex transform at 1, 8, 40 threads —
paid at every `setup` without wisdom, and once per thread count with it.

GEMM alone: 8 BLAS threads are 5–8× faster than 1 for every Luna shape, even the 4-mode
`8192×4×16` product (≈2×), independently of `J`.

## Whole propagations (`runs.jl`)

Wall time (s), best of ≥ 2, `:estimate`, wisdom off. "Old" is today's default (FFTW `4J`,
BLAS 8), "proposed" the rule below, "best" the fastest configuration measured.

| case (per-plan FFT elements) | J | old | proposed | best |
|---|---|---:|---:|---:|
| mode-averaged 1 ps (16k) | 1 / 4 / 8 | 2.17 / 2.64 / 4.25 | 2.17 / 2.14 / 2.16 | 2.17 / 2.14 / 2.16 |
| mode-averaged 8 ps (128k) | 1 / 4 / 8 | 7.94 / 5.12 / 6.30 | 7.54 / 5.14 / 6.00 | 7.41 / 5.12 / 6.00 |
| envelope Raman (2k), GNLSE (2k) | 1–8 | up to 2.2× slower at `J=8` | = best | |
| 4 modes `:fixed` | 1 / 4 / 8 | 6.32 / 4.84 / 4.86 | 6.26 / 2.34 / 1.72 | 6.26 / 2.13 / 1.72 |
| 4 modes `:adaptive` | 1 / 4 / 8 | 17.6 / 15.4 / 18.0 | 17.6 / 15.4 / 17.2 | 17.6 / 15.1 / 17.0 |
| radial 256 | 1 / 4 / 8 | 1.42 / 3.33 / 1.96 | 1.42 / 0.95 / 1.33 | 1.42 / 0.95 / 1.33 |
| radial 1024 | 1 / 4 / 8 | 7.07 / 11.8 / 7.92 | 7.07 / 6.43 / 6.02 | 7.07 / 6.22 / 6.02 |
| 2-D free χ⁽²⁾ (64k) | 1 / 4 / 8 | 0.43 / 0.48 / 0.53 | 0.43 / 0.48 / 0.47 | 0.43 / 0.48 / 0.47 |
| 3-D free 16 × 8 (16k) | 1 / 4 / 8 | 0.15 / 0.17 / 0.21 | 0.15 | 0.15 |
| 3-D free 64 × 64 (512k) | 1 / 4 / 8 | 5.40 / 2.07 / 1.67 | 4.47 / 2.06 / 1.67 | 4.31 / 2.06 / 1.67 |

The proposed rule is within 10 % of the best measured configuration in every case (worst:
4-mode `:fixed` at `J=4`, 2.34 vs 2.13 s) and never slower than the old default beyond
noise; the old default is up to 3.5× slower (radial 256 at `J=4`).

What drives it:
1. **BLAS against the Julia pool.** In the 4-mode `:fixed` run the GEMMs are small and
   interleaved with the threaded plasma loop; 8 BLAS threads make the whole run 2.8×
   slower at `J=8` although the GEMM alone is 2× faster. Large radial GEMMs gain 3× from
   BLAS threads even at `J=1`.
2. **FFTW threads help only large non-GEMM transforms** (8 ps window, 64 × 64 3-D grid) and
   hurt small ones (up to 2× at `J=8`). In the whole runs the threshold is higher than in
   isolation: the 2-D χ⁽²⁾ case (64k elements per plan) is 5–10 % slower with FFTW threads,
   the 8 ps case (128k) 16–25 % faster. With GEMMs in the same step (radial), FFTW threads
   plus BLAS threads oversubscribe (radial 1024 at `J=8`: 6.0 → 10.9 s).
3. **FFTW's own pthreads at `J=1`** give 4–20 % on the large non-GEMM cases (3-D 64 × 64:
   5.4 → 4.3 s with 8) and cost 30 % on the 1 ps case with 8, which the size threshold avoids.
4. **`J`**: `-t 4` is as fast as or faster than `-t 8` for mode-averaged, adaptive multimode
   and small radial runs; `-t 8` wins for `:fixed` multimode and the large free-space and
   radial grids. Nothing at `-t 10` was better than `-t 8`.

**Metal** (`runs_metal.csv`): host thread settings barely matter on a device run. FFTW
threads change nothing; BLAS 8 helps by 10–30 % through host work at setup (building the
radial matrices).

## Proposed rule ("auto", `fftw_threads = 0` and a new `blas_threads = 0`)

Per transform, at `Luna.setup`, with `J = Threads.nthreads()` and `N` the elements of each
FFT plan:

| state | FFTW threads per plan | BLAS threads |
|---|---|---|
| device (Metal/CUDA) | 1 | unchanged (OpenBLAS default) |
| CPU, transform without GEMMs (mode-averaged, GNLSE, Cartesian free space) | `N ≥ 2¹⁷` ? (`J > 1` ? `J` : 4) : 1 | unchanged |
| CPU, radial (large GEMM) | 1 | `J > 1` ? `J` : unchanged (OpenBLAS default) |
| CPU, multimode (small GEMMs) | 1 | 1 |

An explicit `set_fftw_threads(n)` / `BLAS.set_num_threads(n)` by the user wins. The 4×
multiplier goes. At `J = 1` Luna would, for the first time, let FFTW use pthreads (4) for
large plans.

## Decisions for the user

1. The rule and its thresholds as above (FFT threshold 2¹⁷ elements; `J=1` large plans get 4
   FFTW pthreads; multimode BLAS = 1 even though `J=4`/BLAS 4 was 9 % faster once).
2. Setting BLAS threads from Luna changes global state that other code in the session sees.
   Options: set it at `setup` and log once (proposed), or only recommend it in the docs.
3. `Luna.tune_threads()` (fit the two thresholds per machine, stored in the scratch cache):
   implement now, or leave for later.
4. The docs could recommend `-t 4` to `-t 8` (not `-t auto`, which includes efficiency cores)
   for laptops; Luna cannot choose `J` itself.
