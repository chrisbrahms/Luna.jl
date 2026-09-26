# Retire the Capillary `neff_grid`/`neff_β_grid` overloads

Base: `gpu/int-E2` (ced6b299). Plan §4.5, layer 3.

## Motivation

`Capillary` overloaded `LinearOps.neff_β_grid` (mode-averaged) and `LinearOps.neff_grid`
(multimode, via the `FixedCoreCollection` union) for fixed-radius `MarcatiliMode`s, caching
the waveguide part of the effective index so that a z-dependent (gradient) operator did not
re-evaluate the cladding index at every call. That mattered when the operator closure was
evaluated at every RK stage. Since `gpu/27` the default `linop_integral=:auto` either detects
a z-independent operator or tabulates Φ(z) once at setup, so the closure runs only at the
table nodes. The plan retires the overloads if tabulation is the default or the measured
cost is negligible; this branch measures it.

## Measurements

`benchmark/neff_overloads.jl` (run with `julia --project=. benchmark/neff_overloads.jl` on the parent commit; it deletes the overload methods itself), this machine (M1 Pro), 1 thread, He 0 → 1 bar gradient, 125 µm
core, 1 m, `RealGrid(800e-9, (150e-9, 4e-6), 400e-15)` (2049 frequency samples). Same
process, overload methods deleted between the two halves.

| | with overload | generic `Modes.neff` |
|---|---|---|
| mode-averaged `linop!`, per call | 6.4e-5 s | 3.1e-4 s (4.9×) |
| 4-mode `linop!`, per call | 7.8e-4 s | 7.8e-4 s (1.0×) |
| `prop_capillary`, mode-averaged, `:auto` | 0.47 s | 0.55 s |
| `prop_capillary`, mode-averaged, `:quadrature` | 1.23 s | 4.15 s |
| `prop_capillary`, 4 modes, `:auto` | 6.26 s | 6.16 s |

The operator values and every output `Eω` are bitwise identical with and without the
overloads (max difference 0). The multimode overload gave no speed-up at all, although the
per-call benchmark passes a `Tuple` and does reach it; its `AbstractArray` branch could never
match a `Vector` of modes anyway (the `where` inside the element type makes it invariant).
The per-call timings use the 2049-sample grid above; the `prop_capillary` runs use their own
4097-sample grid, where the mode-averaged setup cost of the generic path is about +0.08 s per
run (one timing each, of the order of the run-to-run noise), and it grows with `Nω` but not
with the fibre length.

The slowdown that does matter is on the `:quadrature` path, where the closure runs at every
stage (3.4× here, up to the 4.9× per-call ratio). That path is taken when requested
explicitly, or by `:auto` when the table is capped: the node budget is
`clamp(256 MiB ÷ (96 B · Nω), 2, 1024)` (`LinearOps.node_budget`), i.e. 682 nodes at
`Nω = 4097`, and this gradient needed 87 nodes, so a mode-averaged run falls back at
`Nω ≳ 32k` (picosecond windows with a UV limit), or at any `Nω` when the refinement needs
more nodes than the cap (tight `linop_tol`, sharp profiles, long fibres). The fallback is
logged at info level. On those runs a fixed-radius Marcatili gradient is now 3–5× slower per
operator evaluation than before; tapers, non-Marcatili modes and multimode runs are
unchanged. The removal is still taken: the overload served one mode type in one geometry,
duplicated the index formula (without the `real(n) < 1e-3` clamp of the generic one), and
the default path is now insensitive to it.

## Changes

- `src/Capillary.jl`: removed the `neff_β_grid` and `neff_grid` methods, `FixedCoreCollection`,
  the helpers only they used (`neff_wg`, the three-argument `neff(m, εco, nwg)`), and the
  `LinearOps` import. The comment on `zconstant` no longer mentions `FixedCoreCollection`.
- `test/test_linops.jl`: the "equivalence for fast z-dependent linops" testset compared the
  overload with the generic path through `Modes.delegated`; it now checks that a delegated
  mode gives the same operator (renamed, comments updated, the type-inequality assert dropped),
  and the `RealGrid` block now also checks `βfun!` against an independent `Modes.β` per
  frequency.
- `benchmark/neff_overloads.jl`: the measurement script. It is meaningful only on the parent
  commit (on this branch there is nothing to delete and both halves run the generic path).
- The generic `LinearOps.neff_grid`/`neff_β_grid` stay as the extension point their docstrings
  describe.

## Tests

- `test_linops.jl` 364/364, `test_capillary.jl` 181/181, `test_gradient.jl` 7/7,
  `test_tapers.jl` 2/2, `test_multimode.jl` 15/15.
- Regression gate against `gpu/int-E2` ced6b299: 506/506, every case and mode exactly 0 (fixed and adaptive, `Eω` and statistics). Gate run 178 s, baseline generation 253 s.

## Known gaps

None. The GPU path never used the overloads.
