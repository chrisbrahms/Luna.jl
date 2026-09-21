# `gpu/22-modal`: the batched column evaluator and the fixed transverse quadrature

Base: `gpu/int-D` at `90826dc4`. Five commits.

## Motivation

Luna's multimode transform evaluated the transverse integral of the nonlinear polarisation
one point at a time. For every point of the adaptive cubature rule it did a matrix product
to synthesise the field there, an inverse Fourier transform of that one column, the
nonlinear responses on that column, a forward transform and a matrix product back. Two
consequences:

- the number of transform columns grew with the number of transverse points, and every one
  of them was a separate call;
- nothing about it could run on a device. It is host scalar code driven by `Cubature`,
  which returns `Vector{Float64}`, and a batched response (the plasma) was handed one
  column at a time, so its column threading never engaged.

This branch replaces that loop with one **batched column evaluator** shared by two modal
transforms — the adaptive one, which stays the default, and a new fixed-quadrature one,
which is the multimode transform that runs on a GPU (GPU_PLAN.md §4.4, §6 Group E, §9
decision 2).

## What changed

### `src/Modes.jl`

- **`mode_matrix`/`mode_matrix!`** — the batched form of what `to_space!` does at one
  point: the normalised transverse fields ([`Exy`](@ref)) of a mode collection at a list of
  coordinates, as `(nmodes, npol, npts)`. The layout is chosen so that `reshape(M, nmodes,
  npol*npts)` is the synthesis matrix of a Luna `(nto, npol, npts)` block (polarisation
  fastest). Points outside the modes' `dimlimits` are zero, as `to_space!` makes them; a
  negative polar radius is treated as outside rather than raised on, because the one caller
  which can produce one discards those points anyway. The normalisation constant `N` is
  evaluated once per mode instead of once per point.
- **`TransverseQuadrature`** and `transverse_quadrature`/`quadrature_nodes`/
  `quadrature_weights` — a fixed rule on the transverse plane, stored on a reference domain
  so that it maps to a `z`-dependent core radius: Gauss–Legendre (or Gauss–Kronrod) in r or
  x, a periodic trapezoid in θ or Gauss–Legendre in y, with an embedded coarse rule where
  one exists. Ported from the `jtravs/modal-fixed` fork (`Modes.jl:535-755`) with review;
  see "From the fork" below.
- **`zconstant`** and **`azimuthal_order`** traits, `false`/`nothing` by default.
  `Capillary.MarcatiliMode` overrides both — `zconstant` for a numeric core radius (the
  same type parameter `FixedCoreCollection` selects on), `azimuthal_order` as `|n-1|` for
  HE, 1 for TE/TM — and `Antiresonant.ZeisbergerMode`/`VincettiMode` delegate them to the
  `MarcatiliMode` they wrap, as they already do for `field`, `N` and `dimlimits`.

### `src/NonlinearRHS.jl`

- **`AbstractTransModal`**, **`ModalBlock`**, **`ModalRound`**, **`synthesise_responses!`**
  — the evaluator. The modal spectrum goes to the oversampled time domain once per
  right-hand side over all `nmodes` columns; one matrix product synthesises the field at
  every point of the current set; the responses see the whole `(nto, npol, npts)` block
  through the existing `Et_to_Pt!`. Synthesis is linear and diagonal in time, so
  synthesising before or after the transform gives the same field, and doing it after costs
  `nmodes` transform columns instead of one per point.
- **`TransModal`** rewritten on the evaluator; the per-point `pointcalc!` loop is deleted.
  The cubature driver still runs on the host and hands over a round of points at a time; a
  round is evaluated as one block, split into chunks of at most `maxbatch`
  (`MODAL_MAXBATCH`, 16). Each distinct round width gets its own buffers, forward plan
  *and* copy of any response which owns buffers, because a `Batched` response sizes its
  buffers to the block it is given; the widths are fixed by the rule (`pcubature_v` gives
  3, 2, 4, 8, …, `hcubature_v` 17, 34, …, independent of the integrand's dimension), so the
  dictionary stops growing after the first right-hand side. The projection is a broadcast
  writing straight into the driver's buffer, reinterpreted as the complex modal array.
- **`TransModalFixed`** — the same evaluator on a fixed quadrature rule. Because the
  weights are known in advance, the points are summed *before* the transform back (one
  matrix product with the weights folded in), so the number of transform columns does not
  grow with the number of nodes. Parametric in its array type, with grid mirrors, unit
  scaling, a residency assertion and the `Adapt`-free construction the other device
  transforms use. `update_matrices!` rebuilds the mode matrices when `z` moves unless
  `Modes.zconstant` says the profiles do not depend on it. `integral_error!` /
  `has_error_estimate` give the embedded Gauss–Kronrod error estimate on demand, with one
  further matrix product against the precomputed difference of the two weight sets.
- **`NormModal`/`norm_modal`** — a struct with a mirrored, precombined vector instead of a
  closure over `grid.ω`, with the unit scaling folded in and a `check_norm` method. The
  `Float64` values are exactly what the per-call expression produced.

### `src/Luna.jl`

The two `setup` methods for a mode collection become one `setup_modal`, which takes
`modal_integral` (`:adaptive` or `:fixed`), `device`, `precision`, `maxbatch` and the
quadrature keywords, builds the normalisation for the chosen device and precision, and
refuses `:adaptive` for anything but a host `Float64` run with a message naming
`modal_integral=:fixed`. `runscaling` and `save_modeinfo_maybe` cover the new transform.

### `src/Interface.jl`

`prop_capillary` takes `modal_integral`, `nr`, `nθ` and `kronrod`, records them in the
output file, and ignores them for mode-averaged propagation (which has no transverse
integral), as it already ignored `radial_integral_rtol`. The device sentinel resolves to
`Luna.settings["device"]` for a multimode call as well, but only with
`modal_integral=:fixed` and only when every response has a device kernel; an explicit
`device`/`precision` request goes through the same response check the mode-averaged path
uses. `_cpu_only!` is now only `prop_gnlse`'s.

`Stats.mode_reconstruction_error` is written against the adaptive transform — it
re-evaluates the transform at one transverse point and records the cubature's own error
estimate — so a `:fixed` run collects the other default statistics and not that one, and
asking for it explicitly through `stats_kwargs` is an error which says why.
`NonlinearRHS.integral_error!` is the fixed rule's own estimate and becomes a statistic in
`gpu/25`, after `Stats` is refactored by `gpu/24`.

### The three host-hardcodings §4.4 names

- `reset!`'s `::Array{ComplexF64,2}` annotation is gone; it takes any array.
- The host `Vector{Float64}` `Cubature` returns is still there, and it is now the
  *stated* reason the adaptive transform is host- and `Float64`-only: its constructor
  refuses any other `spec`/`scaling` and names the fixed rule.
- The concrete buffer fields of `TransModal` are therefore kept, not made parametric for
  a path the transform cannot take, and the docstring says why. The transform's own
  time-domain buffers are parametric in the element type, as they were before; the
  complex ones stay `Array{ComplexF64}` because the driver's result is reinterpreted into
  them. `TransModalFixed` is the parametric transform, and it is the one which runs on a
  device.

## Deviations from GPU_PLAN.md and from the brief

1. **Stages 1–3 landed in one commit** (`a5c2b373`) rather than three. The fixed transform
   shares the `ModalBlock`, the evaluator, the rewritten `Luna.setup` modal method and the
   docstring cross-references with the adaptive one; splitting them would have meant
   writing the modal `setup` twice and leaving dangling `@ref`s in between. The `Interface`
   wiring, the tests, the Metal tests and the benchmarks/docs are separate commits.
2. **Where the synthesis happens.** §4.4 describes "one `mul!` to synthesise" without
   saying in which domain. It is done in the *time* domain, on the modal field, which is
   the same field (the operation is linear and diagonal in time) and which removes the
   per-point inverse transform. The fork does the same.
3. **`maxbatch`.** §4.4 does not mention a cap on the round width. Without one, a run at a
   tight `radial_integral_rtol` would allocate a block, a plan and a copy of every
   buffer-owning response for every width up to `mfcn`. The cap is a keyword on
   `Luna.setup` and `TransModal`.
4. **`Modes.scale_invariant` is not ported.** The fork uses it to follow a taper by
   rescaling the mode matrices instead of re-evaluating them. Without it, a tapered
   `:fixed` run re-evaluates the mode fields on the host whenever `z` moves — which is what
   the adaptive rule does at *every point* of every round, so it is not a regression. See
   "Known gaps".
5. **The `x1 >= ul[2]` condition in the Cartesian branch of the point rule** looks like a
   typo for `x2 >= ul[2]` and is preserved exactly, with a comment saying so. Changing it
   would change the answer of every Cartesian-domain multimode run. Not fixed, per the
   brief.

## From the fork: taken, and what was changed

Taken with review from `scratchpad/modal-fixed`: `TransverseQuadrature`,
`transverse_quadrature`, `quadrature_nodes`, `quadrature_weights`, `mode_matrix`,
`zconstant`, `azimuthal_order`, the `TransModalFixed` evaluation order, `integral_error!`
and `has_error_estimate`.

Changed:

- the transform is built through `alloc`/`todevice`/`gridvectors` and asserts residency,
  instead of `Luna.device_zeros` and `Adapt.adapt`, and it carries the unit scaling, so it
  works in `Float32` (the fork is `Float64`/CUDA only);
- the responses go through `Nonlinear.rescale_responses` and `Et_to_Pt!` — the `gpu/12`
  protocol — instead of the fork's `batched_responses`/`apply_responses!`, so there is one
  response implementation, and the vector case is the fused two-broadcast path rather than
  a per-node column loop;
- `Wc` is stored as `Wd = Wc - Wp`, so the error estimate is one matrix product rather than
  a five-argument `mul!` (not every backend's `mul!` takes `α`/`β`);
- the adaptive transform is rebuilt on the same evaluator; the fork left it as it was, so
  the fork has two implementations of the multimode physics and this branch has one;
- `scale_invariant` and the `zconstant=nothing` keyword's fork defaults are not taken
  (see above); `zconstant` is still overridable per transform.

## Regression gate

`LUNA_REGRESSION_BASE=fa556e6f` (source-identical to the base `90826dc4`), M1 Pro, Julia
1.13.0, `-t 1`, `:estimate`, no wisdom, one FFTW and one BLAS thread:

**460 pass, 0 fail.** Every case is exactly `0.000e+00` except the two which go through a
modal transform:

| case | mode | Eω diff | Eω tol | stats diff | stats tol |
| --- | --- | ---: | ---: | ---: | ---: |
| `modeavg_field_vector` | fixed | 1.452e-15 | 1.0e-12 | 4.676e-13 | 8.0e-11 |
| `modeavg_field_vector` | adaptive | 3.086e-11 | 6.3e-09 | 5.159e-06 | 1.1e-03 |
| `multimode_field_plasma` | fixed | 6.913e-14 | 6.7e-12 | 5.417e-13 | 5.8e-11 |
| `multimode_field_plasma` | adaptive | 3.146e-07 | 7.3e-05 | 3.043e-04 | 7.3e-02 |

(`modeavg_field_vector` is elliptically polarised, so it is a two-mode, two-component
`TransModal` run despite the name.) Both worst quantities are
`stats/transverse_integral_error_*`, which are `HCubature`'s own error *estimates* for the
transverse integral: the adaptive rule subdivides differently when the integrand changes in
its last bits, so the estimate moves by far more than the integral does — which is why
`gpu/00-harness` gave those two cases their own tolerances. The accepted step counts are
unchanged in both cases and both modes.

Why they moved: the projection, the synthesis and the transforms are now batched, so the
order of accumulation is different. Where it cost nothing the order was kept — the
per-point projection broadcast sums the polarisation components in the same order the
matrix product did and applies the Jacobian factor to the sum, as the per-point code did,
and `norm_modal`'s precomputed vector holds exactly the values the per-call expression
produced.

`LUNA_REGRESSION_BASE=fdf8dbe3` (`evanescent`): **460 pass, 0 fail**, largest difference
2.319e-04. The extra non-zero rows there (`modeavg_field_plasma`, `modeavg_field_adk`, the
two radial cases, the two χ⁽²⁾ cases) are the drifts earlier branches of the project
recorded, unchanged by this one.

## Tests

| file | result | how |
| --- | --- | --- |
| `test/test_regression.jl` | 460 pass, 0 fail (both bases) | `LUNA_REGRESSION_BASE=fa556e6f julia --project=. -t 1 test/test_regression.jl` |
| `test/test_device.jl` | 715 pass, 0 fail, 39 testsets | `julia --project=<jlenv> -t 1 -e 'using Luna; include("test/test_device.jl")'` |
| `test/test_metal.jl` | 448 pass, 0 fail, 18 testsets | `julia --project=<metalenv> -t 1 -e 'using Luna, Metal; include("test/test_metal.jl")'` |
| `test_multimode`, `test_modes`, `test_polarisation_field`, `test_vectorplasma`, `test_stats`, `test_noise`, `test_interface`, `test_capillary`, `test_antiresonant` | 1336 pass, 0 fail, 47 testsets | `using Luna; Luna.set_fftw_mode(:estimate); Luna.set_fftw_threads(1); Luna.set_fftw_wisdom(false); include(...)` |

New coverage:

- **`test_modes.jl`** (137 new): the polar and Cartesian rules integrate their domain
  exactly with the fine *and* the embedded coarse weights; a Gauss rule of n nodes is exact
  for the polynomials it should be; the periodic trapezoid in θ is exact below `nθ` and
  wrong at `nθ`, which is where the `nθ ≥ 4h+1` rule comes from; the mode matrix is
  bit-for-bit the numbers `to_space!` leaves in the `ToSpace`, is zero outside `dimlimits`,
  and the shape checks and traits behave.
- **`test_device.jl`, host**: the fixed rule against the adaptive one (below); a tapered
  mode collection rebuilding its matrices and agreeing with the adaptive rule at the new
  position; the embedded error estimate, including the `NaN` when there is no embedded
  rule; the refusal of the adaptive integral for a `Float32` or scaled run.
- **`test_device.jl`, JLArray**: the evaluator against the host at 1e-10 for a field and an
  envelope grid, with and without plasma, with the residency of every matrix, buffer,
  mirror and response array asserted; four-mode Kerr and Kerr+plasma propagations end to
  end; the error estimate; the refusal of `:adaptive` on a device array.
- **`test_metal.jl`**: below.
- **`test_interface.jl`**: `:fixed` runs multimode in `Float32`, agrees with `:adaptive` to
  1e-6 on the strongest mode after a full propagation, collects the right statistics, and
  the keyword validation.

### The fixed rule against the adaptive one

One right-hand side, four HE₁ₘ modes of a 75 µm capillary in argon, `nr=64`:

| case | relative difference |
| --- | ---: |
| field grid, Kerr | 3.0e-16 |
| field grid, Kerr + plasma | 2.1e-14 (3.3e-12 with a tabulated PPT rate) |
| envelope grid, Kerr | 4.2e-16 |
| field grid, 2 modes, `:xy`, full 2-D, `nr=64`, `nθ=16` | 2.6e-08 |

For the smooth HE₁ₘ fields the Gauss rule is far more accurate than the adaptive rule at
its default 1e-3 tolerance, so what these numbers measure is the *adaptive* rule's error.
The full 2-D row is larger because `hcubature` in two dimensions stops much earlier. The
plasma rows are looser because the ionisation rate is a spline of `log(rate)` evaluated at
the local field, so the integrand it contributes is piecewise cubic in the field rather
than smooth.

Through the simple interface, after a 1 cm propagation with plasma (adaptive steps, so the
step sequences differ too): 2.6e-09 on the strongest mode, 3.7e-04 on the fourth, which
carries 2e-11 of the energy.

## Metal

M1 Pro, `Metal.jl` v1.11, Float32. Every number is normalised to the strongest mode over
the whole `(nω, nmodes)` array at the end of a 10-step fixed-step propagation, against
explicit host references built in the same call:

| case | Metal vs CPU F32 | Metal vs CPU F64 | CPU F32 vs F64 |
| --- | ---: | ---: | ---: |
| Kerr, 4 modes, `nr=32` | 3.8e-08 | 5.4e-07 | 5.4e-07 |
| Kerr + plasma, 4 modes, `nr=32` | 1.9e-08 | 4.4e-07 | 4.4e-07 |
| Kerr, 2 modes `:xy`, full 2-D, `nr=16`, `nθ=8` | 3.8e-08 | 5.4e-07 | 5.4e-07 |
| Kerr, 4 modes, `boundary=:rate` | 4.3e-06 | 3.4e-06 | 2.1e-06 |
| Kerr + plasma, 4 modes, 1600 fs | 5.7e-08 | 5.5e-07 | 5.5e-07 |

The device is never what limits the agreement: Metal against the same arithmetic on the CPU
is one to two orders tighter than `Float32` against `Float64`.

`test_metal.jl` also runs the whole thing through `prop_capillary_args` + a fixed-step
`Luna.run`, which brings in `boundary=:rate` and the default statistics, and checks that
`modal_integral=:adaptive` on Metal is refused with a message naming the fix.

## Benchmarks

M1 Pro, Julia 1.13.0, one FFTW thread, one BLAS thread, `:estimate`, no wisdom.

### The rewrite against the base, adaptive path

One right-hand side, argon at 0.1 bar, `Luna.setup(...; full=false, rtol=1e-3)` — the same
script run against a source export of the base commit and against this branch:

| case | state | base, `-t 1` | branch, `-t 1` | base, `-t 8` | branch, `-t 8` |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kerr, 4 modes, 400 fs | 4100 | 1.120 ms | **0.938 ms** | 1.121 ms | **0.938 ms** |
| Kerr, 4 modes, 1600 fs | 16388 | 4.868 ms | **4.107 ms** | 4.933 ms | **4.113 ms** |
| Kerr, 2 modes `:xy`, 400 fs | 4098 | 0.709 ms | **0.651 ms** | 0.714 ms | **0.652 ms** |
| Kerr, 2 modes `:xy`, 1600 fs | 16386 | 3.189 ms | **2.886 ms** | 3.193 ms | **2.894 ms** |
| Kerr+plasma, 4 modes, 400 fs | 4100 | 2.993 ms | 2.967 ms | 2.984 ms | **1.957 ms** |
| Kerr+plasma, 4 modes, 1600 fs | 16388 | 11.633 ms | 11.830 ms | 11.867 ms | **6.693 ms** |
| Kerr+plasma, 2 modes `:xy`, 400 fs | 4098 | 1.993 ms | 2.205 ms | 1.997 ms | **1.583 ms** |
| Kerr+plasma, 2 modes `:xy`, 1600 fs | 16386 | 8.158 ms | 8.956 ms | 8.279 ms | **5.581 ms** |

Read three things off it.

- **Kerr is 8–17 % faster at any thread count**, from the transform columns the batching
  removes.
- **Kerr + plasma is 1–10 % *slower* at one thread.** The plasma response now works on a
  block of up to `maxbatch` columns instead of one, and at 16388 samples × 2 polarisation
  components × 8 columns its five buffers no longer fit in cache. This is the one
  measured regression in the branch.
- **Kerr + plasma is 21–44 % faster at eight threads**, because the response's column
  threading (`Nonlinear.PLASMA_THREAD_MINLEN`, GPU_PLAN.md §4.9) now has more than one
  column to share out. The base is identical at one and eight threads, as it must be.

The regression matrix's two modal cases, through `benchmark/run.jl` (one thread, so the
worst case above): `multimode_field_plasma` rhs 45.4 → 43.9 ms, prop 14.30 → 14.20 s;
`modeavg_field_vector` rhs 36.4 → 35.4 ms, prop 13.66 → 13.65 s.

### `maxbatch`

Batching the points of a round, with the modal transform already hoisted out of the loop,
is worth about 9 %: one right-hand side, Kerr, 4 modes, 1600 fs — `nb=1` 4.53 ms, `nb=4`
4.12 ms, `nb=16` 4.12 ms, `nb=32` 4.11 ms. 16 is the default; the rule asks for 3, 2, 4 and
8 points per round at the default tolerance, so nothing is split there.

### The two rules, and the device (`benchmark/modal.jl`)

Four HE₁ₘ modes, argon at 0.1 bar, one right-hand side:

| | Kerr, 4100 | Kerr, 16388 | Kerr+plasma, 4100 | Kerr+plasma, 16388 |
| --- | ---: | ---: | ---: | ---: |
| CPU `Float64` `:adaptive` | 938 µs | 4.12 ms | 2.96 ms | 11.9 ms |
| CPU `Float64` `:fixed` `nr=32` | 188 µs | 844 µs | 2.13 ms | 8.45 ms |
| CPU `Float64` `:fixed` `nr=64` | 306 µs | 1.60 ms | 4.34 ms | 16.7 ms |
| CPU `Float32` `:fixed` `nr=64` | 165 µs | 714 µs | 3.82 ms | 15.4 ms |
| Metal `Float32` `:fixed` `nr=32` | 588 µs | 821 µs | 831 µs | 1.06 ms |
| Metal `Float32` `:fixed` `nr=64` | 410 µs | 733 µs | 956 µs | 1.60 ms |

A whole 10-step propagation, Kerr+plasma at 16388 samples: 756 ms (`Float64` adaptive), 975
ms (`Float32` fixed `nr=64`), **106 ms** (Metal, same rule).

- The fixed rule is **not** automatically cheaper on the CPU. The adaptive rule needs about
  17 transverse points for a smooth HE₁ₘ set, so a 64-node rule does roughly four times the
  work; the Kerr rows hide that (the transform count does not grow with `nr`, only the two
  matrix products and the responses do) and the plasma rows do not. `:fixed` is for the
  device, and for a fixed per-step cost; it is not a CPU optimisation, which is why the
  default is `:adaptive` (GPU_PLAN.md §9 decision 2).
- The GPU's advantage grows with the work per node: 1.0× for Kerr at 16388 samples, **9.6×
  for Kerr + plasma** against the same arithmetic on the CPU and **7.4× against the
  `Float64` adaptive default**. The small-grid Kerr rows are launch-bound, which is why
  `nr=32` can be slower than `nr=64` there.

## Known gaps and open questions

1. **A tapered mode collection with `:fixed` re-evaluates the mode fields on the host**
   whenever `z` moves — `nmodes × npol × npts` calls to `field`, plus one `Modes.N` per
   mode per `z`, and `N` is `@memoize`d, so the dictionary grows a entry per stage for a
   mode type without a closed-form `N`. The adaptive rule has exactly the same problem
   today (it evaluates `Exy` at every point of every round), so this is not a regression,
   but `Modes.scale_invariant` from the fork, or `gpu/23`'s tabulation, would fix it. A
   fixed-radius capillary — the overwhelmingly common case — is unaffected: `zconstant` is
   `true` and nothing is rebuilt.
2. **`nr` is not checked against the mode set.** The rule is not adaptive and will not say
   it is under-resolved. `nθ` *is* checked, against `Modes.azimuthal_order`, with a warning
   at setup. A residual-based check on `nr` (comparing the fine and coarse rules at setup,
   which `kronrod=true` already makes possible) would close this; it is not in this branch.
3. **`mode_reconstruction_error`, `transverse_points` and `transverse_integral_error_*` are
   not collected for a `:fixed` run.** `Stats.jl` is `gpu/24`'s and the Kronrod statistic is
   `gpu/25`'s, per the brief. `NonlinearRHS.integral_error!` exists and is tested; nothing
   calls it per step.
4. **Memory.** A round width holds its own block *and* its own copy of every
   buffer-owning response, so the adaptive transform's footprint is now roughly
   `Σ widths ≈ 2 × maxbatch` columns instead of one. At `maxbatch=16` and a 16384-sample
   oversampled grid with two polarisation components that is ~40 MB where it used to be
   ~1 MB. `maxbatch` is the dial; I did not try to share one buffer set across widths
   because `Nonlinear.rescale` sizes a batched response's buffers to the block it is given
   and a view would not satisfy its exact-size check.
5. **Planning inside the propagation.** A round width the driver has not asked for before
   is planned on the first right-hand side that needs it, not at setup, because the widths
   are not known in advance. `Utils.saveFFTwisdom()` is called when a width is added, so it
   is paid once per shape per machine, but with `:patient` planning the first step of a run
   is slower than the rest.
6. **`x1 >= ul[2]`** in the Cartesian branch of the point rule (deviation 5 above) looks
   like a bug. It is preserved, not fixed.
7. `prop_capillary` has no `maxbatch` keyword; it is on `Luna.setup` and on the transform.
   Nothing in the interface needs it, but a user with a very tight
   `radial_integral_rtol` cannot lower it from there.

## Conflicts to expect at `gpu/int-E`

| file | this branch | likely against |
| --- | --- | --- |
| `src/NonlinearRHS.jl` | the whole modal section replaced (base 404–623), `AbstractTransModal`/`ModalBlock`/`ModalRound` inserted before it | `gpu/20-radial-device` rewrites `TransRadial` and `FreeSpaceNorm`, which are *after* it; adjacent-hunk only |
| `src/Luna.jl` | the two modal `setup` methods replaced by `setup_modal`; `runscaling`, `save_modeinfo_maybe`, one comment in `run` | `gpu/23` and `gpu/24` both edit `Luna.run`; `gpu/20` edits the radial `setup` methods |
| `src/Interface.jl` | `prop_capillary_args` keywords and the device sentinel, the multimode `setup` methods, `_statskwargs` | `gpu/24` will touch the `Stats.default` call next to `_statskwargs` |
| `src/Stats.jl` | **not touched** | — |
| `test/test_device.jl`, `test/test_metal.jl` | new testsets appended before the last two | every Group E branch appends there |
| `docs/src/gpu.md`, `docs/src/developer/device_model.md` | the admonition, "What runs where", "Performance", and one new section each | every Group E branch edits the same three sentences |
