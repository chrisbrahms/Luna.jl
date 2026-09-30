# gpu/34-cpu-vs-gpu: whole-propagation CPU versus GPU benchmark

Base: `gpu/33-threads`. Benchmark scripts and results only; no change to `src/`.

## What

`benchmark/cpu_gpu/`:
- `cases.jl`: full propagations set against the critical power — mode-averaged (1 ps and
  8 ps windows), 8-mode capillary at 0.95 P_cr with the adaptive and fixed (nr 64, 128)
  transverse rules, focused radial Kerr + plasma at 0.8 P_cr (256/512/1024 points) and
  focused 3-D envelope Kerr at 0.8 P_cr (64²/128²/256²).
- `run.jl`: one process per setting (`DEVICE` cpu/metal/cuda, `PRECISION`, `-t`, `BLAS`,
  `STATS_N`); times setup and propagation separately after a short compile run, takes the
  step count from the integrator, saves the final field (`fields/`, gitignored).
- `summarise.jl`: per-case table with best-of-repeats time, spread, speed-up over the
  serial baseline (`cpu64_t1_b1`) and over `cpu64_t8`, and accuracy against the baseline
  (field L2, spectral dB, mode energies).
- `REPORT.md`: first-round results on an M1 Pro with Metal; `results/` raw CSV and the
  full summary table.

## Results (M1 Pro, Metal)

Metal (Float32, statistics 100 times per run) against Luna's threaded CPU default (Float64,
`-t 8`): 0.42× for a 1 ps mode-averaged run, 2.1× at 8 ps, 2.3–2.9× for fixed-rule
multimode, 5–11× for radial (growing with N), ≈3.6× for 3-D; the fixed rule on Metal is
8.4× faster than the CPU's adaptive default for the hard multimode case. Accuracy against
serial Float64 is the same as CPU Float32's (≤ 4e-4 L2, ≤ 0.05 dB). Details, and the
findings on statistics cost, Float32 on the CPU, `-t 1` not being serial by default and
efficiency cores, in `benchmark/cpu_gpu/REPORT.md`.

## Not done

CUDA (the scripts take `DEVICE=cuda`; to be run on an NVIDIA machine).

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01JLyXeRJFXy3CczpvjZHJWW
