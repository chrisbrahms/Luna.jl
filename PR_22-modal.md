# `gpu/22-modal`: the batched column evaluator and the fixed transverse quadrature

Base: `gpu/int-D` at `90826dc4`. Ten commits.

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

`prop_capillary` takes `modal_integral`, `modal_nr`, `modal_nθ` and `modal_kronrod`,
records them in the output file, and ignores them for mode-averaged propagation (which has no transverse
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
| `test_rect_modes`, `test_polarisation`, `test_polarisation_env`, `test_tapers`, `test_gradient`, `test_mixtures`, `test_processing`, `test_output`, `test_linops` | 2563 pass, 0 fail, 28 testsets | as above |

The eight modal examples in `examples/low_level_interface` run (the propagation part of
each, up to the first plotting call): `basic_modal`, `basic_modal_env`,
`full_modal/basic_modal_full`, `full_modal/basic_modal_full_bothpolarisations`,
`rectangular/rectangular_modal`, `tapers/taper_modal`. Two do not, and neither is this
branch's doing — both fail identically on the base commit:
`polarisation/modal_nonvector_plasma.jl` uses `linop` two lines before it is defined, and
`polarisation/elliptical_env.jl` calls `Luna.setup` with a `normfun` positional argument
the signature has not had for a long time.

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
4.12 ms, `nb=16` 4.12 ms, `nb=32` 4.11 ms. 16 is the default; for this case the driver
asks for rounds of 3, 2, 4, 8 **and 16** points, so nothing is split at the default and
the dictionary of blocks is `[1, 2, 3, 4, 8, 16]`.

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
6. **`x1 >= ul[2]` in the Cartesian branch of the point rule** (deviation 5 above) is a
   pre-existing bug, deliberately not fixed here. Review round 1 worked out what it
   costs. The condition is reached only for a Cartesian domain with `full=true`, i.e.
   `RectModes.RectMode`, whose `dimlimits` is `(:cartesian, (-a, -b), (a, b))`:

   - for `a > b` — a wide, shallow guide — every cubature point with `b <= x1 < a` is
     treated as outside and contributes zero, so a strip of the guide is silently dropped
     from the transverse integral and the nonlinearity is underestimated. That is a real
     physics bug, and it predates this branch;
   - for `a <= b` the wrong test never fires (`x1 < a <= b` always), and the missing
     `x2 >= ul[2]` test is harmless because the h-adaptive rule samples the open interior
     only. Square guides — what `test_rect_modes.jl` and the examples use — are
     unaffected, which is why it has never been caught.

   The correct condition is `inside &= !(x2 <= ll[2] || x2 >= ul[2])`. Fixing it changes
   the answer of every `a > b` rectangular multimode run, so it wants its own commit with
   a regression case and an entry in the gate's tolerance table, not this one.
7. `prop_capillary` has no `maxbatch` keyword; it is on `Luna.setup` and on the transform.
   Nothing in the interface needs it, but a user with a very tight
   `radial_integral_rtol` cannot lower it from there.

## Conflicts to expect at `gpu/int-E`

| file | this branch | likely against |
| --- | --- | --- |
| `src/NonlinearRHS.jl` | the whole modal section replaced (base 404–623), `AbstractTransModal`/`ModalBlock`/`ModalRound` inserted before it | `gpu/20-radial-device` rewrites `TransRadial` and `FreeSpaceNorm`, which are *after* it; adjacent-hunk only |
| `src/Luna.jl` | the two modal `setup` methods replaced by `setup_modal`; `runscaling`, `save_modeinfo_maybe`, and (since review round 1) the stale comment in `run` at 689–695 | `gpu/23` and `gpu/24` both edit `Luna.run`; `gpu/20` edits the radial `setup` methods |
| `src/Interface.jl` | `prop_capillary_args` keywords and the device sentinel, the multimode `setup` methods, `_statskwargs` | `gpu/24` will touch the `Stats.default` call next to `_statskwargs` |
| `src/Stats.jl` | **not touched** | — |
| `test/test_device.jl`, `test/test_metal.jl` | new testsets appended before the last two | every Group E branch appends there |
| `docs/src/gpu.md`, `docs/src/developer/device_model.md` | the admonition, "What runs where", "Performance", and one new section each | every Group E branch edits the same three sentences |

## Changes after review round 1

`reviews/gpu-22-modal-1.md` — "approve with minor fixes", every load-bearing claim
reproduced independently. All ten findings are addressed here; one commit (`eacee07d`)
plus this file.

**Finding 1 (fix before `int-E`) — `has_error_estimate` for a Cartesian rule.** The θ
clause fired for any `full=true` rule, but a Cartesian domain's second coordinate is
Gauss–Legendre in y and has no embedded rule, so `Wd` was identically zero and
`integral_error!` reported a quadrature error of exactly zero. The clause is now
polar-only, and the error-estimate testset has a Cartesian case (two `RectMode`s) with
`kronrod=false` (no estimate, `Wd == 0`, `NaN`) and with `kronrod=true` (a finite,
non-zero estimate).

**Finding 3 (fix before `int-E`) — the wisdom write inside the right-hand side.**
`modalround!` no longer calls `Utils.saveFFTwisdom()`: that took the shared FFTW pid lock,
deleted the cache file and re-exported all wisdom from inside an RK45 step, which is the
lock several processes of a `Scans.runscan` share. The round widths are predictable —
`Cubature.pcubature_v` doubles a Clenshaw–Curtis rule and `hcubature_v` works in multiples
of 17, neither depending on the integrand, its dimension or the tolerance — so
`NonlinearRHS.modal_round_widths` computes them and the constructor builds them all,
inside `Luna.setup` where the wisdom is loaded and saved. At `maxbatch = 16` that is
`[1, 2, 3, 4, 8, 16]` for the radial integral and `[1, 2, 4, 6, 16]` for the full 2-D one.
A width which was not predicted is still built on demand, without writing wisdom. A new
testset checks the predicted list for both drivers and that a whole right-hand side adds
none of them. The memory is the same as before — the reviewer measured 33 MiB of round
buffers as the steady state, which is what is now allocated up front — except for a run
which converges in two rounds, which now allocates blocks it does not use.

**Finding 9 — keyword names.** `prop_capillary`'s `nr`, `nθ` and `kronrod` are now
`modal_nr`, `modal_nθ` and `modal_kronrod`: they are top-level keywords of a flat
interface, meaningful only for a mode collection with `modal_integral=:fixed`, and the
radial and free-space device branches will want the generic names. `Luna.setup` keeps
`nr`/`nθ`/`kronrod`, where the context leaves no room for confusion. `saveargs`, the
tests and both documentation pages follow.

**Finding 2 — the low-level `MethodError`.** `Stats.mode_reconstruction_error` is typed on
`TransModal`, so `Stats.default(..., modes, ...)` on a fixed transform raises a
`MethodError` unless `mode_error=false` is passed; only `prop_capillary` is protected.
`Stats.jl` is `gpu/24`'s, so the fix stays out of this branch, but the `TransModalFixed`
docstring, the `Luna.setup` modal docstring and the user page now say so.

**Findings 4, 5 — stale comments and two PR inaccuracies.** Fixed: `Luna.run`'s list of
which transforms can produce a device `Eω` (the PR's conflict table said this branch
already changed it — it did not, and now it does); three comments in `test_metal.jl`,
including a block that had been left 157 lines from its testset; and
`Nonlinear.PLASMA_THREAD_MINLEN`'s "a block with a single column — every mode-averaged and
modal transform", which this branch's own 8-thread numbers depend on being false. The
`maxbatch` paragraph above now says the driver asks for 16-point rounds too.

**Findings 6, 7, 8 — nits.** `Interface._statskwargs` accepts a `NamedTuple` as well as a
`Dict`. There is a comment on why the mode matrices are held in the field's element type
(complex with a zero imaginary part on an envelope grid: `mul!` needs matching element
types to reach a BLAS `gemm` on the host and the accelerated path on a device). Two of the
three coverage gaps are closed: the modified shot-noise model through `TransModalFixed` on
JLArray against the host (the noise buffer, `Emt_nl` and the `1/Eref` division), and the
multimode `:auto` sentinel on Metal hardware (`modes=4, modal_integral=:fixed` with no
`device` keyword resolves to the GPU, while the `:adaptive` default stays on the host and
`device=:cpu` opts out). The third — a multimode *envelope* propagation on a device — is
still one right-hand side on JLArray only.

**Finding 10** — the `x1 >= ul[2]` typo — is deliberately not fixed; the reviewer's
analysis of what it costs is in known gap 6 above.

### Tests after the fixes

| what | result |
| --- | --- |
| gate, `LUNA_REGRESSION_BASE=fa556e6f` | **460 pass, 0 fail**, the same two modal rows to every digit (1.452e-15 / 3.086e-11 and 6.913e-14 / 3.146e-07) |
| `test_device.jl` | **738 pass, 0 fail, 41 testsets** (was 715/39) |
| `test_metal.jl` | **452 pass, 0 fail, 18 testsets** (was 448; the four new assertions are the `:auto` sentinel ones in "Luna.set_device(:cpu) opts out") |
| `test_modes.jl`, `test_multimode.jl`, `test_interface.jl` | 724 + 6 + 349 pass, 0 fail |
