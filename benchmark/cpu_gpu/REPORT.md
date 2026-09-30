# CPU versus GPU: whole propagations (first round)

Apple M1 Pro (8 performance + 2 efficiency cores, 32 GB), Metal, Julia 1.13, branch
`gpu/34-cpu-vs-gpu` (Luna's thread defaults from `gpu/33-threads`). Raw rows:
`results/HW-WYQH90RQMT.csv`; the full table (every setting, accuracy columns) is
`results/HW-WYQH90RQMT_summary.md`, produced by `summarise.jl`. Run 2026-09-30 on mains
power under `caffeinate`, one setting at a time; 1 h 53 min in total.

## Cases

Every case is a full, physically meaningful propagation (`cases.jl`), timed without
compilation and without setup; the step count is the integrator's own (all attempted
steps). Strengths are set against the critical power (`Tools.Pcr`, 11.3 GW for Ar at
1 bar and 800 nm).

| case | geometry | physics | grid |
| --- | --- | --- | --- |
| `modeavg` | capillary, mode-averaged, 50 cm | field, Kerr + PPT plasma, 150 µJ | 1 ps window |
| `modeavg_long` | same, 20 cm | same | 8 ps window (≈2¹⁶ points) |
| `modal_adaptive` | capillary, 8 HE₁ₘ modes, 30 cm | field, Kerr + plasma, 0.95 P_cr | adaptive transverse rule |
| `modal_fixed64/128` | same | same | fixed rule, nr = 64 / 128 |
| `radial256/512/1024` | radial, focused w0 = 100 µm, 10 cm | field, Kerr + PPT plasma, 0.8 P_cr | N radial points |
| `free3d64/128/256` | 3-D Cartesian, focused w0 = 200 µm, 30 cm | envelope Kerr, 0.8 P_cr | N × N |

## Settings

| setting | what |
| --- | --- |
| `cpu64_t1_b1` | **baseline**: one Julia thread, one BLAS thread, Float64 |
| `cpu64_t1` | one Julia thread, OpenBLAS at its own default (8) — what `-t 1` gives by default |
| `cpu64_t4`, `cpu64_t8`, `cpu64_t10` | Luna's default threading at 4 / 8 / 10 Julia threads, Float64 |
| `cpu32_t8` | 8 threads, Float32 |
| `cpu64_t8_s100` | 8 threads, statistics 100 times per case instead of every step |
| `metal32` | Metal, Float32, default statistics (every step) |
| `metal32_s100` | Metal, statistics 100 times per case |

`cpu64_t8` and `metal32` were run three times (two for the large cases); the fastest is
reported, and the slowest was within 1–7 % of it except Metal `radial1024` (16 %) and the
adaptive rule on the CPU (12 %).

## Propagation time

| case | serial | CPU `-t 8` | CPU32 `-t 8` | Metal | Metal s100 | Metal s100 ÷ CPU `-t 8` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `modeavg` | 2.1 s | 2.0 s | 1.8 s | 5.4 s | 4.8 s | **0.42×** |
| `modeavg_long` | 6.6 s | 5.1 s | 4.9 s | 3.2 s | 2.4 s | 2.1× |
| `modal_adaptive` | 133 s | 73 s | — | — | — | — |
| `modal_fixed64` | 49 s | 15.2 s | 12.5 s | 9.4 s | 6.7 s | 2.3× |
| `modal_fixed128` | 92 s | 24.9 s | 21.6 s | 11.1 s | 8.7 s | 2.9× |
| `radial256` | 38 s | 14.4 s | 9.9 s | 3.4 s | 2.9 s | 5.0× |
| `radial512` | 130 s | 42.1 s | 27.1 s | 6.0 s | 5.4 s | 7.8× |
| `radial1024` | 430 s | 128 s | 75 s | 13.4 s | 11.9 s | **10.7×** |
| `free3d64` | 22.7 s | 7.6 s | 6.3 s | 2.4 s | 2.4 s | 3.1× |
| `free3d128` | 100 s | 29.3 s | 25.8 s | 8.2 s | 8.1 s | 3.6× |
| `free3d256` | not run | 118 s | 100 s | 33.9 s | 33.0 s | 3.6× |

Summed over the nine cases every setting ran: serial 871 s, CPU `-t 8` 269 s, CPU32
185 s, Metal 62.5 s, Metal s100 53.2 s.

## Accuracy

Final field against the serial Float64 run (`summarise.jl` columns; L2 = relative norm
of the difference, dB = worst spectral difference above −40 dB, ΔU_m = worst relative
mode energy):

- Every threaded Float64 run agrees with the serial one to ≤ 6e-11 in the field (FFTW and
  BLAS thread counts change rounding only); `modal_fixed64` is the exception at 1.9e-5,
  because the step sequence differs by one step (757 vs 758).
- Metal agrees to 1.5e-5 – 3.8e-4 in L2 and ≤ 0.05 dB (the long-window mode-averaged case;
  ≤ 0.02 dB elsewhere); the multimode mode energies to ≤ 1e-3. **CPU Float32 is the same
  size** (1.6e-5 – 3.6e-4, ≤ 0.02 dB, ΔU_m ≤ 1.5e-3), so the difference is the precision,
  not the device.
- Step counts: Float32 runs take the same number of steps to within 1–3.

## Findings

1. **A GPU is worth it for everything except short mode-averaged runs.** Metal is 2–11×
   faster than Luna's threaded CPU default, and the gap grows with the radial grid
   (5× → 11× from 256 to 1024 points). 3-D saturates at ≈3.6×. The 1 ps mode-averaged run
   (≈2¹³ points) is 2.4× *slower* on Metal; at 8 ps it is 2× faster, so the crossover for
   a single column lies between.
2. **Part of the GPU's gain is just Float32.** CPU Float32 is 1.1–1.7× faster than
   Float64 (most for radial). Against CPU Float32, Metal s100 is still 1.9–2.5× (multimode),
   3.4–6.3× (radial) and 2.6–3.2× (3-D).
3. **Per-step statistics cost the GPU up to 30 %.** The default statistics of multimode,
   radial and long mode-averaged runs include members with no device form, so the field
   is copied to the host every step (Luna warns once). Statistics 100 times per run make
   Metal 29 % faster for `modal_fixed64`, 22 % for `modal_fixed128`, 25 % for
   `modeavg_long` and 10–16 % for radial; nothing for 3-D. On the CPU the same change is
   worth ≤ 10 % (13 % for the adaptive rule).
4. **The fixed multimode rule on Metal is 8.4× faster than the CPU default (adaptive).**
   `modal_adaptive` 73 s vs `modal_fixed128` on Metal s100 8.7 s; `benchmark/threads`
   found nr = 128 at least as accurate as the adaptive rule for this case. On the CPU
   alone, fixed nr = 128 is already 2.9× faster than adaptive.
5. **CPU threading gains up to 3.7× over serial at `-t 8`**: most for multimode fixed
   (3.2–3.7×), radial (2.7–3.4×) and 3-D (3.0–3.4×); 1.8× for the adaptive rule; 1.3× for
   `modeavg_long` and nothing for the short mode-averaged run.
6. **`-t 1` is not serial by default.** With one Julia thread Luna leaves OpenBLAS at its
   own default (8 threads here), so the radial and multimode GEMMs still run in parallel:
   `radial1024` takes 144 s at `-t 1` against 430 s fully serial. The serial baseline
   needs `BLAS=1` (`Luna.set_blas_threads(1)`).
7. **The efficiency cores do not help.** `-t 10` is never faster than `-t 8` beyond noise
   and is slower for the adaptive rule (84 s vs 73 s) and 3-D (34 s vs 29 s). `-t 4` is
   within 5 % of `-t 8` for radial (GEMM-bound; BLAS gets 4 threads either way) and 1.2×
   faster than `-t 8` for `modeavg_long` (4.3 s vs 5.1 s; not investigated).

## Open for the next round

- Set of cases: add or drop anything? Candidates: radial 2048 (Metal only, the CPU would
  take ≈ 8 min at `-t 8`), a 2-D Cartesian χ⁽²⁾ case, envelope radial, Raman.
- Whether the benchmark should use `STATS_N=100` as the GPU default setting, and whether
  Luna should pick a statistics period automatically on a device (finding 3).
- CUDA: `run.jl` takes `DEVICE=cuda` (Float64 by default, `PRECISION=32` for single); a
  CUDA machine also runs its own CPU settings, since every comparison is within one
  machine. Suggested order on that machine: `BLAS=1 -t 1`, `-t <P-cores>`, `PRECISION=32
  -t <P-cores>`, `DEVICE=cuda`, `DEVICE=cuda PRECISION=32`, and the same two with
  `STATS_N=100`.
