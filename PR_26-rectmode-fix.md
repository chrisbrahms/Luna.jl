# The Cartesian in-domain test of the adaptive modal integral

Base: `gpu/int-E` at `5ec21f77`. Branch: `gpu/26-rectmode-fix`. GPU_PLAN.md §6 Group E2,
§9 decision 7.

## Motivation

`NonlinearRHS._points!` decides, for each transverse point the adaptive cubature driver
hands over, whether that point is inside the waveguide. For a Cartesian domain it read

```julia
inside &= !(x2 <= ll[2] || x1 >= ul[2])
```

The upper limit is tested against the wrong coordinate: `x1` rather than `x2`. The branch
is reached only for a Cartesian `Modes.dimlimits`, which in Luna means
`RectModes.RectMode`, whose limits are `(:cartesian, (-a, -b), (a, b))`. So for a guide
with `a > b`, every cubature point with `b <= x1 < a` was treated as outside the guide and
contributed nothing, and a strip of the domain was silently dropped from the transverse
integral of the nonlinear polarisation.

For `a <= b` the wrong test never fires (`x1 < a <= b`), and the missing upper test on
`x2` is harmless because the h-adaptive rule samples the open interior only. Nothing in
the test suite ran the adaptive rule on an `a > b` Cartesian domain: `test_rect_modes.jl`
uses 50 × 12.5 µm but only for mode properties, and `test_device.jl`'s 50 × 20 µm guide
only reaches `TransModalFixed`. The example
`examples/low_level_interface/rectangular/rectangular_modal.jl` does -- 50 × 10 µm, 18
modes, `full=true` -- so that example has been running with 40 % of the guide's area
excluded from the transverse integral. Measured: its in-domain point count at z = 0 goes
from 365 to 527, the nonlinear polarisation at z = 0 changes by 3.7e-1, and after 1.5 cm
the field changes by 1.1e-1 normalised over the whole field and by O(1) in the weak
antisymmetric modes. The example is not edited; it now gives the right answer.

It is a pre-existing bug, not something the GPU work introduced: the same condition is in
`evanescent` at `src/NonlinearRHS.jl:295`. `gpu/22-modal` preserved it deliberately, with
a `NOTE`, and its review reported it as known gap 6; §9 decision 7 is to fix it here, in
its own branch, with a regression case.

The fixed quadrature rule (`Modes.transverse_quadrature`, used by `TransModalFixed`) was
already correct: its Cartesian nodes are interior Gauss–Legendre nodes mapped affinely
onto the rectangle, and the mode matrix zeroes points by `Modes._outside`, which tests
`x2` against `ul[2]`. So before this branch the two drivers disagreed by 8.3 % on a
100 × 40 µm guide and by nothing at all on a square one, and `modal_integral=:fixed`
silently changed the answer of an `a > b` run.

## What changed

### `src/NonlinearRHS.jl`

`_points!`: `inside &= !(x2 <= ll[2] || x2 >= ul[2])`. The test is now the one
`Modes._outside` makes, so the adaptive driver, the mode matrix and the fixed quadrature
rule all agree on the rectangle. The `NOTE` that preserved the old condition is replaced
by a note saying what it was and when it changed.

The polar branch is deliberately left as it was, and the comment on `_points!` now says
so: there `r == 0` counts as outside, which `Modes._outside` does not do, and the
Clenshaw–Curtis rule `Cubature.pcubature_v` uses does sample `r == 0`. Both give the same
integral, because the Jacobian factor `pre` is zero there, so making them agree would be
churn on the one path every existing multimode case takes.

No other file in `src/` changes. `TransModalFixed`'s node set is unchanged, so no Metal or
JLArray work was needed.

### `test/test_multimode.jl`

A new testset, "Rectangular transverse integral" (9 tests), which checks the transverse
integral of a known function against its analytic value with both drivers.

For a single mode the Kerr projection is analytic. The field at a transverse point is
`Eₘ ê/√N`, so the polarisation is `ρ ε₀ γ₃ (Eₘ ê/√N)³` and its projection back onto the
mode is `(∫ê⁴dA / N²) ρ ε₀ γ₃ Eₘ³`. `NonlinearRHS.Erω_to_Prω!` returns the same windowed,
transformed and normalised polarisation at a single transverse point, so with `norm!` the
identity the ratio of the transform's output to `Erω_to_Prω!` at the centre of the guide,
where `ê = 1`, is `∫ê⁴dA / √N` — everything else (the time window, the oversampled
transforms, the spectral window, the density, `γ₃`) cancels. For the fundamental
rectangular mode `∫ê⁴dA = 9ab/16` and `N = ½√(ε₀/μ₀)ab`.

The test runs 100×40 µm (`a > b`), 40×100 µm (`a < b`) and 60×60 µm (`a == b`), with the
adaptive driver at `rtol=1e-6` and with `modal_integral=:fixed` at the default 64×16
nodes.

### `test/regression/`

A new case, `rect_modal_field` (`cases.jl`, `setup_rect_modal`): two `RectMode`s in a
100 × 40 µm argon-filled guide at 5 bar with silver cladding, `m = 1` and `m = 3` of the
`x` index (the coordinate the guide is wide in, which is the one the in-domain test is
applied to) with the same `y` index and the same polarisation, Kerr only, 5 µJ of 10 fs
at 800 nm over 3 cm, `full=true`. It is the matrix's only Cartesian transverse domain and
the only case which reaches the `:cartesian` branch of `_points!`.

`tolerances.jl` gains its entry and the header gains the paragraph explaining it;
`README.md` records it as the 23rd case, whose baseline no pre-`gpu/26` baseline directory
contains.

## Measurements

Machine: M1 Pro, Julia 1.13.0, `-t 1`, `set_fftw_mode(:estimate)`,
`set_fftw_threads(1)`, `BLAS.set_num_threads(1)`, FFTW wisdom disabled.

### The physics correction

Measured against a detached worktree of `5ec21f77` (the branch base, pre-fix) running the
same case definitions, both under the regression harness.

| mode | class | pre-fix vs post-fix difference | worst quantity |
| --- | --- | --- | --- |
| `:fixed` | `Eω` | 2.095e-01 | `Eω [component 2/2, save 2/11]` |
| `:fixed` | `stats` | 9.154e-01 | `stats/transverse_integral_error_rel` |
| `:adaptive` | `Eω` | 2.095e-01 | `Eω [component 2/2, save 2/11]` |
| `:adaptive` | `stats` | `Inf` | step count (73 steps, was 74) |

The `:fixed` statistics that move, in full: `transverse_integral_error_rel` 9.15e-01,
`transverse_integral_error_abs` 9.05e-01, `transverse_points` 2.06e-01, `fwhm_t_max`
7.24e-02, `mode_reconstruction_error` 5.56e-02, `peakintensity` 4.15e-02, `fwhm_r`
3.08e-02, `fwhm_t_min` 2.73e-02, `peakpower` 8.10e-03, `energy` 4.58e-03,
`peakpower_allmodes` 2.43e-03, `ω0` 1.93e-03, `fwhm_t_min_allmodes` /
`fwhm_t_max_allmodes` 1.44e-03.

In the `:adaptive` mode the step-size controller takes one step fewer post-fix (73 against
74), which the gate reports as its hard `step count` failure and which suppresses the
statistics comparison for that mode. That is the correct report: a different number of
steps is a different propagation, and the statistics are recorded per step.

The correction is not small and it is not confined to the weak mode: the total energy
moves by 0.46 %, the peak on-axis intensity by 4.1 %, and the second mode's `Eω` by 21 %.
The excluded strip `b ≤ x < a` is 30 % of the guide's area and held 17.5 % of the points
the driver evaluated; weighted by the fundamental's `ê⁴` it is 8.3 % of the integral, and
three centimetres of propagation turn that into the 21 %.

These baselines were generated from `5ec21f77` with *this branch's* `cases.jl`, in the
scratch worktree `Luna-gpu/baselines/gpu26-5ec21f77`; the post-fix side is this branch.

### A square guide is bit-identical

The same two-mode propagation, run pre-fix and post-fix in the fixed-step mode, with the
guide made square (`a = b = 60 µm`) and left wide (`100 × 40 µm`):

| guide | `Eω` bit-identical pre/post | max normalised difference of `Eω` | in-domain cubature points per RHS, pre → post |
| --- | --- | --- | --- |
| 60 × 60 µm | yes | 0.0 | 527 → 527 |
| 100 × 40 µm | no | 2.10e-01 (mode 2, save 2) | 435–447 → 527 |

The point count is `stats/transverse_points`, which is `TransModal.ncalls`, the number of
points the driver evaluated that the in-domain test accepted. Post-fix it is 527 at every
step; pre-fix it ranges over 435–447 across the steps (435 at the first), because the
driver subdivides differently when part of the integrand is zeroed. Roughly 90 of 527
points — the strip `b ≤ x < a` — were being thrown away.

### The analytic transverse integral

Largest relative deviation of `nl / Erω_to_Prω!(t, (0,0))` from `∫ê⁴dA/√N`, over the
frequency samples where the polarisation is above 1e-3 of its peak, post-fix:

| guide | adaptive (`rtol=1e-6`) | fixed (64 × 16 Gauss) |
| --- | --- | --- |
| 100 × 40 µm | 8.7e-10 | 1.2e-13 |
| 40 × 100 µm | 8.7e-10 | 1.1e-13 |
| 60 × 60 µm | 8.7e-10 | 4.4e-13 |

Run in the `5ec21f77` worktree the same three rows come out 8.2588e-02, 8.7e-10 and
8.7e-10 for the adaptive driver, and bit-identical for the fixed one. The first is the
analytic value of the bug: the `x` integral becomes `∫_{-a}^{b}cos⁴(πx/2a)dx`, which is
0.91741 of `3a/4` at `a = 2.5b`.

### The gate

```
LUNA_REGRESSION_SKIP=rect_modal_field LUNA_REGRESSION_BASE=b641025b \
    julia --project=$PWD -t 1 test/test_regression.jl
```

```
Regression gate
  baseline commit: b641025beea7b91646ecc1a85a44623e0e5db7ff
  cases:           22 of 23  (skipped: rect_modal_field)
...
Test Summary: | Pass  Total     Time
regression    |  466    466  2m36.5s
largest difference over all cases and modes: 0.000e+00
```

Adding `rect_modal_field` costs the gate about 14 s per run, +9 % on its wall time.

Every one of the 22 existing cases is exactly `0.000e+00` in both modes and both classes.
All of them are polar or mode-averaged, so none reaches the changed branch. The new case
is skipped because no pre-`gpu/26` baseline directory contains it; its own comparison is
the physics-correction table above.

### One-ulp sensitivity of the new case

```
case                     mode               Eω       stats  largest quantities
rect_modal_field         fixed       1.785e-15   2.586e-14  stats/transverse_integral_error_rel 2.59e-14, ...
rect_modal_field         adaptive    2.087e-14   7.252e-11  stats/ω0 7.25e-11, stats/fwhm_t_min 4.17e-11, ...
```

The adaptive step count (73) is reproducible under a one-ulp perturbation of the input,
so the case is admissible under the matrix's step-count rule. The tolerances are 100× the
above with the 1e-12 floor: `:fixed` 1.0e-12 / 2.6e-12, `:adaptive` 2.1e-12 / 7.3e-09.

### Run time of the new case

`:fixed` 3.2 s (20 steps), `:adaptive` 11.2 s (73 steps).

## Tests

Each of these is run in a fresh process with

```
julia --project=$PWD -t 1 -e 'using Luna, Test, LinearAlgebra
Luna.set_fftw_mode(:estimate); Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1); Luna.set_fftw_wisdom(false)
include("test/<file>")'
```

| file | result |
| --- | --- |
| `test/test_multimode.jl` | 15 / 15 pass (6 before this branch, 9 added) |
| `test/test_rect_modes.jl` | 95 / 95 pass |
| `test/test_modes.jl` | 724 / 724 pass |
| `test/test_device.jl` | 1148 / 1148 pass |
| `test/test_regression.jl` | 466 / 466 pass, every difference exactly zero |

`test/test_device.jl` needs `JLArrays`, which is not in the worktree's `Project.toml`
(only `Pkg.test()` adds it). Without it the file's five last testsets error with
`UndefVarError: JLSpec`, because the JLArray section is skipped but those testsets are
not guarded -- 461 pass, 5 error. That is pre-existing on `5ec21f77`; `test_device.jl` is
not touched by this branch. The 1148 above is with `JLArrays` made visible through a
side environment (`JULIA_LOAD_PATH="@:<env with JLArrays>:@stdlib"`), which is the only
way I ran it that exercises the JLArray path at all.

No Metal run: the change is host domain logic in the adaptive cubature driver, which is
host-and-`Float64`-only by construction, and `TransModalFixed`'s node set -- the only
multimode path that reaches a device -- is untouched.

## Known gaps and open questions

- `_points!`'s radial branch (`full=false`) sets `pre = 1.0` and `x2 = 0.0` for a
  Cartesian domain, i.e. it integrates along `x` only and ignores `y` entirely. Nothing
  guards against constructing a `TransModal` that way — `TransModalFixed` errors on
  `full=false` with a Cartesian domain, `TransModal` does not. Out of scope here; it is a
  missing argument check, not a wrong answer for any supported configuration.
- The `a > b` correction is not bounded by a small number: the transverse integral of the
  fundamental's Kerr term was 8.3 % low at an aspect ratio of 2.5, and the error grows
  with the aspect ratio. What that does to a propagation depends on how nonlinear it is --
  21 % in the second mode of the new case. Any existing result from a rectangular
  multimode run with `a > b` is wrong by some amount of that order. Square
  guides, all polar geometries (capillaries, ARR fibres, step-index fibres) and every
  mode-averaged run are unaffected, exactly and bit for bit.
- The new regression case is Kerr only. A rectangular guide with a plasma response would
  exercise the same domain logic through a batched response, but would double the case's
  run time for no extra coverage of the thing this branch changes.

## Differences from GPU_PLAN.md

None. §6 Group E2 asks for the fix, a regression case with an `a > b` `RectMode` guide,
the recorded pre/post delta and a note in the PR; §9 decision 7 is the decision this
implements.
