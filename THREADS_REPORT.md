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

GEMM alone: 8 BLAS threads are 5–8× faster than 1 for the radial and 16-mode shapes and
1.7–5× for the 4-mode ones (`8192×4×16`: ≈2×), independently of `J`.

## Whole propagations (`runs.jl`)

Wall time (s), best of ≥ 2, `:estimate`, wisdom off. "Old" is today's default (FFTW `4J`,
BLAS 8), "proposed" the rule below, "best" the fastest configuration measured.

| case (per-plan FFT elements) | J | old | proposed | best |
|---|---|---:|---:|---:|
| mode-averaged 1 ps (16k) | 1 / 4 / 8 | 2.17 / 2.64 / 4.25 | 2.17 / 2.14 / 2.16 | 2.17 / 2.14 / 2.16 |
| mode-averaged 8 ps (128k) | 1 / 4 / 8 | 7.94 / 5.12 / 6.30 | 7.94 / 5.14 / 6.00 | 7.41 / 5.12 / 6.00 |
| envelope Raman (2k), GNLSE (2k) | 1–8 | up to 2.2× slower at `J=8` | = best | |
| 4 modes `:fixed` | 1 / 4 / 8 | 6.32 / 4.84 / 4.86 | 6.26 / 2.34 / 1.72 | 6.26 / 2.13 / 1.72 |
| 4 modes `:adaptive` | 1 / 4 / 8 | 17.6 / 15.4 / 18.0 | 17.6 / 15.4 / 17.2 | 17.6 / 15.1 / 17.0 |
| radial 256 | 1 / 4 / 8 | 1.42 / 3.33 / 1.96 | 1.42 / 0.95 / 1.33 | 1.42 / 0.95 / 1.33 |
| radial 1024 | 1 / 4 / 8 | 7.07 / 11.8 / 7.92 | 7.07 / 6.43 / 6.02 | 7.07 / 6.22 / 6.02 |
| 2-D free χ⁽²⁾ (64k) | 1 / 4 / 8 | 0.43 / 0.48 / 0.53 | 0.43 / 0.48 / 0.47 | 0.43 / 0.48 / 0.47 |
| 3-D free 16 × 8 (16k) | 1 / 4 / 8 | 0.15 / 0.17 / 0.21 | 0.15 | 0.15 |
| 3-D free 64 × 64 (512k) | 1 / 4 / 8 | 5.40 / 2.07 / 1.67 | 5.40 / 2.06 / 1.67 | 4.31 / 2.06 / 1.67 |

The proposed rule is within 10 % of the best measured configuration in every case at
`J > 1` (worst: 4-mode `:fixed` at `J=4`, 2.34 vs 2.13 s) and never slower than the old
default beyond noise; the old default is up to 3.5× slower (radial 256 at `J=4`). At `J = 1`
the two large non-GEMM cases are 7 % and 25 % above the best, which needed FFTW pthreads
(see below).

What drives it:
1. **BLAS against the Julia pool.** In the 4-mode `:fixed` run the GEMMs are small and
   interleaved with the threaded plasma loop; 8 BLAS threads make the whole run 2.8×
   slower at `J=8` although the GEMM alone is 2× faster. Large radial GEMMs gain 3× from
   BLAS threads even at `J=1`.
2. **FFTW threads help only large non-GEMM transforms** (8 ps window, 64 × 64 3-D grid) and
   hurt small ones (up to 2× at `J=8`). In the whole runs the threshold is higher than in
   isolation: the 2-D χ⁽²⁾ case (64k elements per plan) is 0–2 % slower with FFTW threads
   at `J=4` and 11–13 % at `J=8`, the 8 ps case (128k) 20 % faster at `J=4` and 9 % at
   `J=8`. One case on each side, of different geometry, brackets the threshold between 64k
   and 128k elements; the 128k case sits on it. With GEMMs in the same step (radial), FFTW threads
   plus BLAS threads oversubscribe (radial 1024 at `J=8`: 6.0 → 10.9 s).
3. **FFTW's own pthreads at `J=1`** give 7–20 % on the large non-GEMM cases (3-D 64 × 64:
   5.4 → 4.3 s with 8) and cost 30 % on the 1 ps case with 8, which the size threshold avoids.
4. **`J`**: `-t 4` is as fast as or faster than `-t 8` for mode-averaged, adaptive multimode
   and small radial runs; `-t 8` wins for `:fixed` multimode and the large free-space and
   radial grids. Nothing at `-t 10` was better than `-t 8`.

**Metal** (`runs_metal.csv`): host thread settings matter little on a device run. FFTW
threads change little (the exception: radial at `J=8`/BLAS 8, 0.48 → 0.78 s with 16 FFTW
threads, i.e. more threads only hurt); BLAS 8 helps by 10–30 % through host work at setup
(building the radial matrices).

## Multimode: adaptive against fixed transverse integral (`modal_rules.jl`)

The whole-run matrix above used a weakly nonlinear multimode case, and so did the README
example: there the adaptive rule stops at the same level at every step (31 points: the
33-point Clenshaw–Curtis level of `pcubature` with the two end points excluded; it starts
from 3) and HE₁₁ keeps > 99 % of the energy, so the transverse integral is easy and any
reasonable rule is accurate. For a meaningful comparison the peak power has to approach the critical
power for self-focusing, `Tools.Pcr` (Fibich–Gaeta, with `Tools.getN0n0n2`): 11.3 GW for
argon at 1 bar and 800 nm, 271 GW for helium (so the README case, 11 GW in helium, is at
0.04 P_cr). Probe (8 HE₁ₘ modes, 125 µm, Ar 1 bar, 30 fs, 30 cm, adaptive rule at the
default `radial_integral_rtol = 1e-3`; ratios with `Tools.Pcr`; 125 µm is the core
radius, the first argument of `prop_capillary`; at 0.95 P_cr the initial peak intensity
is ≈1e14 W/cm², the onset of argon ionisation):

| P/P_cr | energy | adaptive points median / 90 % / max | HE₁₁ energy share at the end |
|---|---:|---|---:|
| 0.51 | 184 µJ | 31 / 63 / 63 | 0.987 |
| 0.81 | 294 µJ | 63 / 127 / 255 | 0.730 (HE₁₂ 0.259) |
| 0.97 | 350 µJ | 127 / 127 / 255 | 0.634 (HE₁₂ 0.343) |

(The probe used P_cr = 11.5 GW, Marburger's constant; the timings below use
`Tools.Pcr` = 11.3 GW, so their 0.95 row is ≈2 % lower in energy than the 0.97 probe.)

Timing at 0.8 and 0.95 P_cr (energies from `Tools.Pcr`), FFTW 1 thread, one run each
(30–550 s, so run-to-run noise is small relative to the differences). Error: final `Eω`
against `:fixed` with `nr = 256`, largest difference over all modes / peak over all modes.

| P/P_cr | rule | error | `-t 1` | `-t 8`, BLAS 1 | `-t 8`, BLAS 8 |
|---|---|---:|---:|---:|---:|
| 0.8 | adaptive, rtol 1e-3 (default) | 6.6e-5 | 100.6 s | 58.0 s | 63.4 s |
| 0.8 | adaptive, rtol 1e-4 | 2.0e-5 | 369.8 s | 189.3 s | 214.8 s |
| 0.8 | fixed, nr 32 | 4.7e-4 | 29.4 s | 11.6 s | 24.2 s |
| 0.8 | fixed, nr 64 (default) | 5.4e-5 | 53.7 s | 16.3 s | 38.1 s |
| 0.8 | fixed, nr 128 | 2.3e-5 | 101.6 s | 26.7 s | 63.4 s |
| 0.95 | adaptive, rtol 1e-3 (default) | 1.5e-4 | 133.9 s | 73.9 s | 82.6 s |
| 0.95 | adaptive, rtol 1e-4 | 5.2e-5 | 543.8 s | 277.3 s | 343.3 s |
| 0.95 | fixed, nr 32 | 1.1e-3 | 33.9 s | 12.6 s | 25.9 s |
| 0.95 | fixed, nr 64 (default) | 1.6e-4 | 57.7 s | 17.5 s | 41.5 s |
| 0.95 | fixed, nr 128 | 5.1e-5 | 108.4 s | 28.5 s | 74.7 s |

(`-t 1` column: BLAS 1; BLAS 8 at `-t 1` is within 3 % of it in every row.) Both rules'
default statistics evaluate the right-hand side once more per accepted step (the
reconstruction error for adaptive, the transverse-integral statistic for fixed), so that
overhead is the same on both sides.

Readings:
- At matched accuracy the fixed rule at its default `nr = 64` (5.4e-5, 1.6e-4) is as
  accurate as the adaptive default (6.6e-5, 1.5e-4) and **1.9–2.3× faster at `-t 1`,
  3.6–4.2× faster at `-t 8`**. At its default tolerance the adaptive rule refines to
  63–255 points (1023 at rtol 1e-4) where the integrand is hardest, and each refinement
  round hands the threads only a few points.
- Threads: `nr = 64` goes 53.7 → 16.3 s from `-t 1` to `-t 8` (3.3×); adaptive 100.6 →
  58.0 s (1.7×).
- BLAS: 8 BLAS threads at `-t 8` make the fixed rule 2.1–2.6× slower and the adaptive
  rule 9–24 % slower; no effect at `-t 1`. This is the multimode BLAS = 1 rule.
- Error floor: adaptive at 1e-4 and fixed at `nr = 128` agree with each other and with
  the reference to 2–5e-5; below that the metric is limited by something common to all
  runs (the propagation's own step-size tolerance, or the reference itself), so the
  accuracy of the finer rules is not resolved here. The conclusions use only the rules
  above that floor.
- In the weak case (adaptive stays at 31 points) the two rules cost about the same (4-mode,
  `-t 1`: adaptive 3.9 s, fixed `nr = 32` 3.3 s).

**Found on the way (reviewer):** with the fixed rule's default `modal_kronrod = false`,
the default statistic `Stats.transverse_integral_error` still evaluates the whole
right-hand side once per accepted step (`Stats.jl`, `TransverseIntegralError`) and then
records NaN: about one evaluation in seven is wasted. Skipping it when there is no error
estimate would make the fixed rule faster still; not changed here.

This concerns the numerical method, not only threads, and is recorded for a separate
decision (narrowed pending the accuracy check below): on a multicore CPU, `modal_integral = :fixed` with the default `nr = 64` is the
faster multimode rule at the default accuracy in the strongly nonlinear regime measured
here.

## Conclusions: what matters for CPU speed

In order of effect, from the measurements above (this machine):
1. **The transverse rule of multimode runs** (up to 4.2× at matched peak-normalised
   accuracy, with threads; 1.9–2.3× single-threaded), for 8 HE₁ₘ modes at ≤ 0.95 P_cr.
   Not a thread setting, but the largest single factor measured; see the accuracy check
   for how far it holds.
2. **BLAS threads against Julia threads**: with `J > 1`, BLAS = 1 for multimode (up to
   2.8× in the whole-run matrix, 2.6× above); BLAS = `J` for radial (radial 256 at `-t 4`:
   0.95 s against 1.28 s with 8 and 1.68 s with 1); at `J = 1` OpenBLAS's default (8) is right
   (radial 1024: 7.1 s against 20.2 s).
3. **FFTW threads**: only for large plans without GEMMs (≥ 2¹⁷ elements); otherwise 1.
   The current default (`4J` for every plan) costs up to 2× on small grids and up to 3.5×
   on radial runs together with BLAS 8, and `:patient` planning time grows with the FFTW
   count (25–64× that of one thread at 40 threads, i.e. `4J` for `J = 10`).
4. **Julia threads**: large gains for fixed-rule multimode (3.3×) and large free-space
   grids (3-D 64 × 64: 5.4 → 1.7 s, 3.2×); moderate for radial (1.2–1.5×, the GEMM already
   uses BLAS threads at `-t 1`) and long mode-averaged windows (8 ps: 7.9 → 5.1 s); none
   for small mode-averaged, GNLSE or small free-space runs, which the thresholds leave
   single-threaded. `-t 4` to `-t 8` (at most the performance cores) is the useful range
   here: `-t 4` was as fast as `-t 8` or faster in half the cases, `-t 10` never faster.
5. **What does not matter**: FFTW threads above `J`; any thread setting for small grids
   beyond not oversubscribing; host thread settings on a Metal run (≤ 30 %, through
   setup).

## Proposed rule ("auto", `fftw_threads = 0` and a new `blas_threads = 0`)

Per transform, at `Luna.setup`, with `J = Threads.nthreads()` and `N` the elements of each
FFT plan:

| state | FFTW threads per plan | BLAS threads |
|---|---|---|
| device (Metal/CUDA) | 1 | unchanged (OpenBLAS default) |
| CPU, transform without GEMMs (mode-averaged, GNLSE, Cartesian free space) | `N ≥ 2¹⁷` ? `J` : 1 | unchanged |
| CPU, radial (large GEMM) | 1 | `J > 1` ? `J` : unchanged (OpenBLAS default) |
| CPU, multimode (small GEMMs) | 1 | 1 |

An explicit `set_fftw_threads(n)` / `BLAS.set_num_threads(n)` by the user wins. With
per-plan counts the FFTW wisdom cache (today one file per global count,
`FFTWcache_<n>threads`, `Utils.jl`) has to hold plans for two counts (1 and `J`). The 4×
multiplier goes. At `J = 1` FFTW stays single-threaded, as today (decided with the user:
FFTW's own pthreads would gain 7–20 % on the two large non-GEMM cases at `-t 1`, but a run
started without `-t` should stay single-threaded, and the evidence is four cases; the gain
remains available through an explicit `set_fftw_threads(n)`, which is to be documented;
today `Utils.FFTWthreads` returns 1 whenever `nthreads() == 1`, so the implementation has
to honour an explicit `n` there).

## Decisions for the user

1. The rule and its thresholds as above (FFT threshold 2¹⁷ elements; multimode BLAS = 1
   even though `J=4`/BLAS 4 was 9 % faster once). Decided: no FFTW threads at `J = 1`.
2. Setting BLAS threads from Luna changes global state that other code in the session sees.
   Options: set it at `setup` and log once (proposed), or only recommend it in the docs.
3. `Luna.tune_threads()` (fit the two thresholds per machine, stored in the scratch cache):
   implement now, or leave for later.
4. The docs could recommend `-t 4` to `-t 8` (not `-t auto`, which includes efficiency cores)
   for laptops; Luna cannot choose `J` itself.
