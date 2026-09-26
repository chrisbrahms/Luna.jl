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
BLAS 8), "proposed" the **first version** of the rule (FFTW threshold 2¹⁷, BLAS `J` for
radial, 1 for multimode), "best" the fastest configuration in `runs.csv`. The revised rule
and its comparison against all data are under "Proposed rule" below.

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
(30–550 s). The same configuration timed by `modal_accuracy.jl` in another session differs
by up to 22 % (fixed `nr = 128`: 28.5 vs 34.7 s; adaptive: 73.9 vs 80.9 s), so ratios
below are good to about that. Error: final `Eω`
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
- On the peak-normalised metric the fixed rule at its default `nr = 64` (5.4e-5, 1.6e-4)
  matches the adaptive default (6.6e-5, 1.5e-4) and is 1.9–2.3× faster at `-t 1`,
  3.6–4.2× at `-t 8`; **but that metric does not see weak modes or spectral wings** —
  see the accuracy check below, which changes this conclusion. At its default tolerance the adaptive rule refines to
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

### Accuracy check with mode- and spectrum-resolved metrics (`modal_accuracy.jl`)

0.95 P_cr, `-t 8`, BLAS 1. Two independent references at a 10× tighter propagation
tolerance (`rtol = 1e-7`): `:fixed nr = 256` (refF) and `:adaptive` with
`radial_integral_rtol = 1e-5` (refA, 892 s). They agree to 9e-6 (global), 2.5e-3
(per mode), 5.5e-4 (mode energy), 1.1e-3 dB (spectrum above −40 dB), which is the floor of
each metric. Reference mode energy shares: 0.64, 0.34, 1.6e-2, 3.7e-3, 2.8e-4, 3.2e-5,
5.4e-6, 1.9e-6. Error against refF / refA:

| rule | steps | time | global | per mode | mode energy | spectrum (dB) |
|---|---:|---:|---|---|---|---|
| adaptive 1e-3 (default) | 754 | 80.9 s | 1.3e-4 | 1.2–1.3e-2 | 1.3e-3 | 0.0065–0.0068 |
| adaptive 1e-4 | 741 | 274.6 s | 2.8–3.2e-5 | 2.5–3.8e-3 | 0.6–1.1e-3 | 0.0010–0.0015 |
| fixed nr 32 | 777 | 12.8 s | 1.1e-3 | 8.3e-2 | 1.3e-2 | 0.050 |
| fixed nr 64 (default) | 758 | 18.0 s | 1.6–1.7e-4 | 2.3e-2 | 3.6–4.1e-3 | 0.013 |
| fixed nr 128 | 748 | 34.7 s | 4.5–4.7e-5 | 0.9–1.0e-2 | 6.8–7.0e-4 | 0.0025–0.0028 |
| fixed nr 256 (default rtol) | 742 | 59.0 s | 2.1–2.4e-5 | 2.2–3.2e-3 | 1.6–4.3e-4 | 0.0007–0.0012 |
| fixed nr 64, rtol 1e-7 | 1448 | 37.4 s | 1.7–1.8e-4 | 2.4–2.5e-2 | 3.5e-3 | 0.012 |

- `nr = 64` is **2–3× less accurate than the adaptive default** on the per-mode and energy
  metrics (and 2× on the spectrum), so it is not "as accurate"; the peak-normalised
  metric hid that.
- `nr = 128` is **at least as accurate as the adaptive default on all four metrics
  (1.3–2.8× smaller error; the energy error, 7e-4, is close to the 5.5e-4 floor) and
  2.3–2.6× faster at `-t 8`** (34.7 vs 80.9 s here, 28.5 vs 73.9 s in `modal_rules.csv`);
  at `-t 1` the same comparison is 1.24× (108.4 vs 133.9 s, `modal_rules.csv`).
- The per-mode metric is set by the weakest modes (energy shares down to 1.9e-6), so a
  per-mode error of 1e-2 is a 1 % error in a mode carrying a millionth of the energy.
- Tightening the propagation tolerance does not change `nr = 64`'s error, so the errors
  above are the transverse rule's, not the step controller's. All rules take a similar
  number of steps (741–777), so the time differences are per right-hand side.
- In absolute terms every rule except `nr = 32` is within 0.013 dB on the spectrum and
  ≤ 0.4 % on mode energies; whether that matters is the user's call.

**Found on the way (reviewer):** with the fixed rule's default `modal_kronrod = false`,
the default statistic `Stats.transverse_integral_error` still evaluates the whole
right-hand side once per accepted step (`Stats.jl`, `TransverseIntegralError`) and then
records NaN: about one evaluation in seven is wasted. Skipping it when there is no error
estimate would make the fixed rule faster still; not changed here.

This concerns the numerical method, not only threads, and is recorded for a separate
decision: on a multicore CPU, for 8 HE₁ₘ modes at ≤ 0.95 P_cr, `modal_integral = :fixed`
with `nr = 128` is at least as accurate as the adaptive default on every metric checked
and 2.3× faster at `-t 8` (1.24× at `-t 1`); the default `nr = 64` is 4.5× faster but
2–3× less accurate in the weak modes.

## FFTW planning mode (`planmode.jl`, `runs_planmode.csv`; reviewer round 2)

Luna's default planning mode is `:patient`. On this machine PATIENT plans are **slower**
than ESTIMATE and MEASURE plans for batched and multi-dimensional transforms. Clean
microbenchmark (plan made on a copy of the input, executed on the original; one FFTW
thread; min / max over 3 processes):

| transform | estimate | measure | patient |
|---|---:|---:|---:|
| complex 1024 × 32 × 32, region (1,2,3) | 13.2 ms | 12.9 ms | 19.1–19.5 ms |
| real 4096 × 256, region 1 | 1.56 ms | 1.56 ms | 2.10–2.19 ms |
| real 16384 × 32, region 1 | 0.91 ms | 0.92 ms | 0.98–1.04 ms |
| real 1-D 131072 | 0.43 ms | 0.33 ms | 0.34 ms |
| complex 1-D 16384 | 0.106 ms | 0.118 ms | 0.118 ms |

Whole runs, best of 5 (`runs_planmode.csv`):

| case | J / FFTW | estimate | measure | patient |
|---|---|---:|---:|---:|
| 3-D 64 × 64 | 4 / 1 | 3.02 | 3.06 | 3.48 |
| 3-D 64 × 64 | 4 / 4 | 2.06 | 2.06 | 2.61 |
| 3-D 64 × 64 | 8 / 1 | 2.70 | 2.67 | 3.13 |
| 3-D 64 × 64 | 8 / 8 | 1.55 | 1.54 | 2.84 |
| 2-D χ⁽²⁾ | 4 / 1 | 0.42 | 0.44 | 0.44 |
| 2-D χ⁽²⁾ | 8 / 1 | 0.45 | 0.45 | 0.46 |

The effect is reproducible across processes (so a deterministic plan choice, not noise);
why FFTW's PATIENT search picks worse plans here (timing disturbed by the efficiency cores
or by Julia's task scheduler during planning are candidates) is not established.
`:measure` is within 5 % of the best mode in every whole run and up to 1.8× faster than
`:patient`; for the transforms alone it ranges from 23 % faster (real 131072) to 11 %
slower (complex 16384) than `:estimate`, and it plans 10–80× faster than PATIENT. It also qualifies the thread numbers: the 3-D 64 × 64 gain
from threads (5.4 → 1.7 s) is under `:estimate`/`:measure`; under the default `:patient`
the same run takes 2.6–2.9 s. The `:patient` FFT ratio table above is normalised to
PATIENT's own single-thread plans.

## Conclusions: what matters for CPU speed

In order of effect, from the measurements above (this machine):
1. **The transverse rule of multimode runs**: at equal or better accuracy on every metric
   (fixed `nr = 128` against the adaptive default), 2.3× with 8 threads and 1.24× with one,
   for 8 HE₁ₘ modes at 0.95 P_cr; the fixed rule also gains 3.3× from threads against the
   adaptive rule's 1.7×. Not a thread setting, and measured for one case only.
2. **BLAS threads against Julia threads**: multimode wants the cores the Julia pool leaves
   (`max(1, 8 − J)`; 8 BLAS threads cost up to 2.8×); radial wants 4 at `J = 4, 8` (8 costs
   up to 55 %, 1 up to 3×); at `J = 1` OpenBLAS's default (8) is right (radial 1024: 7.1 s
   against 20.2 s).
3. **FFTW threads**: only for large plans without GEMMs (≥ 2¹⁸ elements); otherwise 1.
   The current default (`4J` for every plan) costs up to 2× on small grids and up to 3.5×
   on radial runs together with BLAS 8, and `:patient` planning time grows with the FFTW
   count (25–64× that of one thread at 40 threads, i.e. `4J` for `J = 10`).
4. **Julia threads**: large gains for fixed-rule multimode (3.3×) and large free-space
   grids (3-D 64 × 64: 5.4 → 1.7 s, 3.2×, under `:estimate`/`:measure`; 2.6–2.9 s threaded
   under the default `:patient`); moderate for radial (1.2–1.5×, the GEMM already
   uses BLAS threads at `-t 1`) and long mode-averaged windows (8 ps: 7.9 → 5.1 s); none
   for small mode-averaged, GNLSE or small free-space runs, which the thresholds leave
   single-threaded. `-t 4` to `-t 8` (at most the performance cores) is the useful range
   here: `-t 4` was as fast as `-t 8` or faster in half the cases, `-t 10` never faster.
5. **FFTW planning mode**: Luna's default `:patient` gives plans up to 1.8× slower in whole
   3-D runs (1.45× for the transform alone) than `:measure`, which was within 5 % of the
   best mode in every whole run here. Not a thread setting; recorded for a separate decision
   (`:measure` as the default; the FFTW threshold follows it, see Decisions).
6. **What does not matter**: FFTW threads above `J`; any thread setting for small grids
   beyond not oversubscribing; host thread settings on a Metal run (≤ 30 %, through
   setup).

## Proposed rule ("auto", `fftw_threads = 0` and a new `blas_threads = 0`)

Revised after review round 1 (BLAS grid `runs_blasgrid.csv`, `:patient` runs
`runs_patient.csv`, 6-repetition noise runs `runs_noise.csv`). Per transform, at
`Luna.setup`, with `J = Threads.nthreads()`, `N` the elements of each FFT plan and `B₀`
OpenBLAS's own default count (8 here, whatever `J`):

| state | FFTW threads per plan | BLAS threads |
|---|---|---|
| device (Metal/CUDA) | 1 | unchanged |
| CPU, no GEMMs (mode-averaged, GNLSE, Cartesian free space) | `N ≥ 2¹⁸` ? `J` : 1 | unchanged |
| CPU, radial (large GEMMs) | 1 | `J = 1` ? unchanged : `B₀ ÷ 2` |
| CPU, multimode (small GEMMs next to the threaded plasma loop) | 1 | `J = 1` ? unchanged : `max(1, B₀ − J)` |

Why these forms:
- **Multimode BLAS `max(1, B₀ − J)`**: the BLAS pool gets the cores the Julia pool leaves.
  `modal_fixed` at `J=4`: BLAS 1/2/4/8 = 2.32/2.20/2.15/4.42 s (best 4); at `J=8`:
  1.74/2.05/2.16/4.14 s (best 1). At `J=1` BLAS makes no difference (≤ 1 %).
- **Radial BLAS `B₀ ÷ 2`**: radial 256 at `J=4`/`J=8`: BLAS 1/2/4/8 = 1.72/1.24/0.98/1.33
  and 1.66/1.23/0.97/1.50 s (best 4 at both); radial 1024: 17.8/10.2/6.44/5.86 and
  17.8/10.2/6.45/6.11 s (4 is within 10 % of the best). `B₀ − J` would give 1 at `J=8`
  (1.66 s, 71 % slower). At `J=1` the default 8 is right (7.1 vs 20.2 s).
- **FFTW threshold 2¹⁸ (raised from 2¹⁷)**: under Luna's default `:patient` planning the
  128k-element case (8 ps window) is not reliably faster with FFTW threads (`J=4`: F1 5.97,
  F4 6.16, F8 5.79, F16 5.10 s; `J=8`: F1 6.34, F8 10.31, F16 6.01 s; best of 3 with
  FFTW's in-memory wisdom reused, so mostly steady state), whereas the 512k case gains
  10–25 % (`J=4`: 3.52 → 2.64 s; `J=8`: 3.15 → 2.85 s) as it does under `:estimate`. The
  cost: under `:estimate` the 128k case loses the 9–20 % it gained with FFTW threads.
  The threshold is bracketed by one case at 128k (mixed) and one at 512k (gain), with no
  case in between, so 2¹⁸ is a conservative choice, not a measured crossover. The worst
  `:patient` result at 128k (`J=8`, F8 10.3 s) is a single erratic plan (F16 6.0 s, and at
  `J=4` F16 was the best), not a trend.
- The `B₀` forms are stated for this machine, where OpenBLAS's default (8) happens to equal
  the performance-core count. That default is not the physical core count in general (it
  depends on platform, build and Julia version), so an implementation should derive the
  core count explicitly or document the dependence; whether "cores left by the Julia pool"
  and "half the default" transfer to machines with SMT or uniform cores is untested.

Against the best configuration measured (all CSVs, best over repetitions), the revised
rule is within 10 % in every case at `J = 4, 8` (radial 1024 at `J=4` at the edge, 1.0995)
except the 128k mode-averaged case (1.20 at `J=4`, 1.09 at `J=8` under `:estimate`; 1.17
and 1.05 under `:patient`), which is the threshold trade-off above. Repeatability across
sessions: ≤ 5 % for the recommended (non-oversubscribed) configurations, up to ≈30 % for
oversubscribed ones (radial `J=4`/F16/B8: 3.33 vs 2.43 s) — oversubscription also makes
run times unstable.

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

1. The revised rule and its thresholds (FFT threshold 2¹⁸ elements; BLAS by geometry).
   Decided: no FFTW threads at `J = 1`.
2. The FFTW planning mode, which the FFT threshold depends on. Under `:measure`
   (`runs_measure.csv`, best of 5) the 128k case gains 18 % at `J=4` (5.92 → 4.87 s) and
   9 % at `J=8` (6.38 → 5.79 s) from FFTW threads, the 64k case is flat at `J=4` and 12 %
   slower at `J=8`, the 16k case 23–48 % slower: so **2¹⁷ with `:measure` as the default,
   2¹⁸ if `:patient` stays**. `:measure` is also the faster mode for multi-dimensional
   grids (previous section).
3. Setting BLAS threads from Luna changes global state that other code in the session sees.
   Options: set it at `setup` and log once (proposed), or only recommend it in the docs.
4. `Luna.tune_threads()` (fit the two thresholds per machine, stored in the scratch cache):
   implement now, or leave for later.
5. The docs could recommend `-t 4` to `-t 8` (not `-t auto`, which includes efficiency cores)
   for laptops; Luna cannot choose `J` itself.
