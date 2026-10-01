# Fixed versus adaptive transverse integral: accuracy and serial cost

Apple M1 Pro, serial CPU (`-t 1`, one BLAS and one FFTW thread), branch `gpu/35-modal-rules`.
Overnight run 2026-09-30 22:00 → 2026-10-01 06:14 (`queue.sh`, `overnight.txt`, logs in
`logs/`). Tables: `results/summary.md` (every run against each case's fixed reference),
`results/summary_vortex_adaptiveref.md` (the vortex against the uncapped adaptive
reference); raw rows `results/runs.csv` (accuracy runs, four concurrent processes) and
`results/timed.csv` (serial timings, one process alone).

**Question:** should Luna switch the multimode transverse integral from the adaptive rule
(`TransModal`, rtol 1e-3, at most `mfcn = 512` points per right-hand side) to the fixed rule
(`modal_integral=:fixed`, defaults nr = 64, nθ = 16) everywhere, including serial CPU runs?

**Answer: not as a plain switch of the default.** The fixed rule is as accurate and faster
for radially symmetric mode sets and for the θ-structured and Cartesian cases, but in the
vortex case — sharp radial structure at a self-compression point — it converges slowly in
r and needs about 4× the adaptive rule's time for the same accuracy. Separately, the
adaptive default's 512-point cap silently ruins its accuracy in every 2-D case where it
binds; that is the most important defect found.

## Cases

| case | modes | integral | what happens |
| --- | --- | --- | --- |
| `he1m_strong` | 8 HE₁ₘ, 125 µm, Ar 1 bar, 800 nm, 30 fs, 0.95 P_cr, 30 cm | radial | HE₁₂ reaches 34 %; adaptive refines 31–255 points |
| `he1m_weak` | 4 HE₁ₘ, 0.04 P_cr | radial | nothing spatial; adaptive stays at 31 |
| `vortex` | HE₂₁ (ϕ = 0, π/4) in quadrature + HE₂₂, HE₂₃ pairs; 1030 nm, 12 fs, 250 µJ, 0.4 bar Ar, 60 cm | full polar | self-compression to 1.65 fs; adaptive 187 → 527 (cap) near compression |
| `mixed` | HE₁₁ + 10 % HE₂₁ seed; HE₁₂, HE₂₁ ×2, TE₀₁, TM₀₁, HE₃₁ ×2; 30 cm | full polar | θ-dependent intensity; HE₂₁ 10 %, TE₀₁ 6 %, HE₁₂ 4 %; adaptive pinned at 527 |
| `rect` | Ag rectangular guide 300 × 150 µm, (n, m) ∈ {1,3,5} × {1,3}; 0.9 P_cr, 10 cm | Cartesian | (1,3) reaches 8 %; adaptive pinned at 527 |

For any combination of HE₂ₘ modes alone, |E|² is independent of θ (Luna's vector Marcatili
fields), so the θ integral of the vortex case is exact for any nθ; it tests r only, and its
fixed rules use nθ = 5 (F64x5, F64x8 and F64x16 agree exactly). `mixed` is the case which
tests θ.

## References and floors

- **he1m_strong**: F512 at prtol 1e-8; the uncapped adaptive rule at rtol 1e-5 agrees to
  8e-6 (global) / 4e-4 dB. Converged.
- **he1m_weak**: every rule agrees to 1e-14 at 1e-8; at 1e-6 everything sits on the same
  propagation floor (3e-6).
- **vortex**: F256x5 at 1e-7 (F512x5 was not reached). The fixed sequence converges slowly
  (below), so the reference is checked against the independent uncapped adaptive A1e-4 at
  1e-7: they agree to 6e-3 dB, which is the size of the step-control floor (4e-3 – 1e-2 dB
  between prtol 1e-6, 1e-7 and 1e-8). Differences below ≈ 1e-2 dB are not resolved here.
- **mixed**: F128x32 at 1e-7, validated only by the fixed sequence (the uncapped adaptive
  references did not finish: one ran 2.7 h before the deadline). r is converged at nr = 64
  (F64x16 and F128x16 agree); θ converges 1.7e-2 → 4.9e-3 → 1.2e-3 dB for nθ = 16, 24, 32.
  Floor: 5e-4 dB (F64x16 1e-7 vs 1e-8).
- **rect**: F128x64 at 1e-7; the uncapped adaptive rule (7 000 – 14 000 points) agrees to
  3.7e-3 dB, F64x64 to 4.5e-3 dB. Floor: ≈ 2e-3 dB at 1e-6.

## Results (prtol 1e-6, the default; serial times from the timed lane)

dB = worst spectral difference above −40 dB; permode = worst per-mode field difference over
modes above a 1e-3 energy share; both maxima over 30 saved z.

| case | rule | serial s | dB | permode |
| --- | --- | ---: | ---: | ---: |
| he1m_strong | **adaptive default** | 134 | 8.4e-3 | 1.1e-3 |
| | F64 (default) | 49 | 1.8e-2 | 1.2e-3 |
| | F128 | 91 | 3.8e-3 | 5.2e-4 |
| | F256 | 179 | 1.6e-3 | 1.8e-4 |
| he1m_weak | **adaptive default** | 8.0 | floor | floor |
| | F16 / F32 / F64 | 2.9 / 5.2 / 9.8 | floor | floor |
| vortex (vs adaptive ref) | **adaptive default** | 321 | 7.9e-3 | 4.7e-4 |
| | F64x16 (default) | 743 | 5.9e-2 | 1.6e-3 |
| | F64x5 | 198 | 5.9e-2 | 1.6e-3 |
| | F128x5 | ≈ 400 ¹ | 2.5e-2 | 7.5e-4 |
| | F256x5 | ≈ 1400 ¹ | 6.1e-3 | 8.3e-4 |
| mixed | **adaptive default** (capped) | 398 | 6.8e-2 | 5.9e-3 |
| | F32x12 | 179 | 2.6e-2 | 1.4e-3 |
| | F64x16 (default) | 442 | 1.6e-2 | 8.7e-4 |
| | F64x32 | ≈ 800 ¹ | 1.2e-3 ² | 6.6e-5 ² |
| rect | **adaptive default** (capped) | 107 | 0.46 | 7.2e-2 |
| | adaptive, cap lifted | 1074 | 3.5e-3 | 4.8e-4 |
| | F32x16 | 68 | 0.15 | 4.1e-2 |
| | F64x16 (default) | 121 | 0.15 | 3.4e-2 |
| | F64x32 | 223 | 1.3e-2 | 2.2e-3 |
| | F64x64 | ≈ 450 ¹ | 4.5e-3 ² | 1.3e-3 ² |

¹ From the concurrent accuracy run (≈ 5–10 % slower than serial) or scaled from the nearest
timed rule; not timed alone. ² At prtol 1e-7 (not run at 1e-6).

## Findings

1. **The adaptive default's `mfcn = 512` cap binds in every 2-D case and then fails
   silently.** `rect` stays at 527 points for the whole run and is 0.46 dB / 7 % per mode
   off; `mixed` is 0.07 dB off; the vortex reaches the cap only around the compression
   point. With the cap lifted, `rect` is accurate (3.5e-3 dB) but takes 10× longer, and
   `mixed` did not finish in 2.7 h. In 1-D (`he1m_strong`) the cap never binds (≤ 255).
2. **1-D (radially symmetric) sets: fixed wins.** F128 is more accurate than the adaptive
   default on every metric and 1.5× faster serially; the current default nr = 64 is 2.7×
   faster but about 2× less accurate in the spectrum. In the weak case every rule sits on
   the propagation floor and F32 is 1.5× faster than adaptive (F64 1.2× slower).
3. **θ-structured polar sets (`mixed`): fixed wins.** F32x12 is more accurate than the
   capped adaptive default and 2.2× faster; θ needs well beyond the 4h+1 exactness bound
   once plasma is present (nθ = 32 for 1e-3 dB with h = 2), while r is converged at 64.
4. **Cartesian (`rect`): fixed wins, but the default nθ = 16 is too few** for modes with
   m = 3: the y direction limits (F64x16 and F32x16 both 0.15 dB; F64x32 1.3e-2 dB; F64x64
   4.5e-3 dB). F64x32 is 4.8× faster than the uncapped adaptive rule for 3.7× its dB error;
   F64x64 matches it at ≈ 2.4× less time.
5. **Sharp radial structure (`vortex`): adaptive wins.** The fixed rule converges roughly
   linearly in nr (0.46 → 0.059 → 0.025 → 0.006 dB for nr = 32 … 256), consistent with a
   non-smooth plasma integrand near the compression point; the adaptive rule refines
   locally and reaches the step-control floor with a median of 289 points in 321 s. A fixed
   rule of the same accuracy needs nr ≈ 256 (≈ 1400 s, ≈ 4.4× slower); the Luna default
   F64x16 is 2.3× slower and 7.5× less accurate in the spectrum.
6. **The fixed default nθ = 16 is wasted on single-order sets** (vortex: F64x5 = F64x16 to
   every digit, 3.7× faster) **and too small for Cartesian sets with m ≥ 3**.

## Recommendations (no `src/` change on this branch)

1. Make the cap visible: warn when the adaptive rule hits `mfcn` (once per run, with the
   point count), and expose `mfcn` in `prop_capillary`; consider a larger default for
   `full=true`.
2. Do not switch the default to `:fixed` globally. If the fixed rule becomes the default
   for some class, the evidence supports radially symmetric sets (with nr = 128) and the
   Cartesian sets (with nθ ≥ 32 for m = 3), not polar sets that may develop sharp radial
   structure.
3. Choose the fixed rule's nθ from the mode set: `4h+1` (exact) for single-order polar sets,
   ≥ 2–3× that when plasma is present and orders mix, and in proportion to the highest m for
   Cartesian modes.
4. Worth trying: a hybrid polar rule, adaptive in r and fixed trapezoid in θ, which would
   combine finding 5 with findings 3 and 6.

## Limits

- One machine, serial only; the GPU comparison of the fixed rule is in `benchmark/cpu_gpu`.
- Vortex differences below ≈ 1e-2 dB are within the step-control floor and the reference
  agreement; F512x5 was not reached.
- The mixed reference is validated only by its own fixed sequence.
- Accuracy runs ran four at a time; only the timed lane's times are serial.
