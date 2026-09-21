# The integral of a z-dependent linear operator

Branch `gpu/27-linop-integral`, base `gpu/int-E` (`5ec21f77`). GPU_PLAN.md §4.5, §4.6,
§6 Group E2, §9 decision 6, §11.

The stepper propagates the linear part of a step by `exp(Φ(t2) − Φ(t1))` with
`Φ(z) = ∫linop dz'`, which is exact. For a constant operator Luna has always formed that
exactly. For a z-dependent one it used `exp(linop(t2)·(t2 − t1))`, a one-point rule, which
is first order in the step size and — the part that matters — whose error is common to
both of the embedded Runge–Kutta solutions, so it cancels out of the error estimate and no
value of `rtol` responds to it. `gpu/23` measured it and made the exact propagator
available behind `tabulate_linop=true`. This branch removes the one-point rule, makes `Φ`
the interface the stepper is given, and adds a second and a third way to supply it.

**Physics change, intended (GPU_PLAN.md §9 decision 6).** Every taper and pressure
gradient result changes. Uniform fibres and every other geometry are exactly unchanged.
The regression gate is **443 pass, 23 fail**, the 23 being exactly `gradient_field_kerr`
and `taper_field_kerr`; every other case is `0.000e+00` in both classes and both modes.

## What changed

### `src/LinearOps.jl` — the interface and the two implementations

`AbstractIntegratedLinop` is a z-dependent operator supplied as its integral. A subtype
implements one of two styles, declared by `PhaseStyle`:

- `AbsolutePhase()` (the default): `phase!(out, op, z)` fills `out` with `Φ(z)`, less
  whatever straight line `secant(op)` reports (`nothing` by default). The propagator keeps
  two buffers, reads `Φ(t1)` once per step — the six stages share it — and `Φ(t2)` once per
  distinct `t2`, and adds `L̄·(t2 − t1)` back in the same broadcast as the exponential.
  This is gpu/23's `TabulatedLinop` arithmetic unchanged.
- `IncrementalPhase()`: `phasediff!(out, op, z1, z2)` fills `out` with `Φ(z2) − Φ(z1)`
  directly, and the propagator caches it on the `(t1, t2)` pair, which catches the same
  repeat (`fbar!` propagates forwards and back at one pair).

Every subtype also implements `derivative!(out, op, z)`, the operator itself. The
propagation never needs it; diagnostics do (a zero-dispersion wavelength, a linear
propagation applied to an input field), and a caller who replaced the closure with an
integrated operator has to be able to get the operator back. All three run per stage and,
on a device, `out` is on the device, so they are broadcasts and every scalar goes through
`Luna.scalar`.

`TabulatedLinop` is now `<: AbstractIntegratedLinop` with `secant` and `derivative!` added;
`derivative!` differentiates the same cubic Hermite `phase!` reads the integral back with,
so the operator it reports is consistent with the `Φ` the propagator uses and is exact at a
node. Nothing else about the table changed.

`QuadratureLinop(linop!, proto; tol, z0, order, maxevals)` is the `IncrementalPhase`
implementation: `QuadGK.quadgk!` over `[z1, z2]` in `Float64` on the host with an absolute
tolerance on the max norm, uploaded through a staging buffer in the state's element type.
One 15-point Gauss–Kronrod rule is the usual cost of a step. `quadgk` returns its best
estimate silently when it runs out of evaluations, so the error estimate is checked and
warned about; `maxevals` (10⁴) bounds a pathological integral — a gradient filled from
vacuum has a `√z` cusp, which Gauss–Kronrod resolves by bisecting towards it and the
table's own bisection handles cheaply. The integral is taken over the step and never from
a fixed origin: the accumulated `Φ` of a metre of fibre is hundreds of radians while the
difference over a step is a fraction of one, which is the same argument as the table's
secant subtraction.

`OffsetLinop(op, δ)` is an integrated operator with a constant added to `L`. A constant in
the operator is a straight line in its integral, so it goes in exactly and, for an
`AbsolutePhase` operator, for free — it is the secant.

### `src/RK45.jl`

`make_prop!(linop!, y0)` — the one-point rule and its host-buffer upload path — is deleted.
The method remains as an `ArgumentError` naming the two wrappers and the `Luna.run`
keyword, because a low-level caller of `RK45.solve_precon` would otherwise get a
`MethodError` on a closure. The constant `AbstractArray` method is untouched and keeps its
own docstring. `isdevice` is no longer imported there.

`RK45.make_prop!(::AbstractIntegratedLinop, y0)` is defined in `LinearOps`, as gpu/23's
was, because `RK45.jl` is included before `LinearOps.jl` and cannot name the type.

### `src/Luna.jl`

`run(...; linop_integral=:tabulated, linop_tol, tabulate_linop=nothing)`. A callable is
converted: `:tabulated` builds the operator table and the transform's `β`/`Aeff` tables as
before; `:quadrature` builds a `QuadratureLinop` and tabulates nothing, so `β` and `Aeff`
go back through the `HostMirror` per stage. An array and an `AbstractIntegratedLinop` are
used as they are, and neither gets value tables. `_linop_integral` validates the symbol and
carries the deprecation: `tabulate_linop` warns once and never overrides an explicit
`linop_integral`; both of its values mean `:tabulated`, because the propagator `false`
selected no longer exists. `linoptype` reports `quadrature` and `integrated` alongside
`constant`, `tabulated` and `variable`.

The tabulation block is now gated on the operator being z-dependent. Under gpu/23 a
constant operator still got two-node `β`/`Aeff` tables when the flag was on; with
tabulation as the default that would move every uniform case at rounding level for nothing.

### `src/Interface.jl`

`linop_integral` on `prop_capillary_args`, forwarded by `boundary_kwargs` and recorded by
`saveargs` in place of `tabulate_linop`, which stays as an accepted (deprecated) keyword.
`Aeff` is pre-tabulated for the statistics only when the operator is z-dependent
*and* `linop_integral === :tabulated`.

### `src/Boundaries.jl`

`addloss` and `addloss_k` gain methods for an `AbstractIntegratedLinop`, returning an
`OffsetLinop`: this is how the spectral and k-space absorbers reach an operator the caller
supplied already integrated, which they could not before (`addloss` would have called it
like a closure). `clampdecay` is not linear in the operator, cannot be pushed through an
integral, and raises with a message naming the fix. `Luna.run` integrates *after*
`Boundaries.setup`, so an operator Luna built itself is always clamped first; only a
caller-supplied integrated operator in free space can reach that error.

### `src/LinearOps.jl` — one behaviour change outside the interface

A `TabulatedScalar` read outside its span now calls its source callable instead of holding
its end value. `prop_capillary` tabulates `Aeff` over `[0, flength]` for the statistics,
and the statistics of the last accepted step are recorded a fraction of a step past the end
of the fibre. gpu/23 documented the held value there as a scale factor on a diagnostic;
with tabulation opt-in it was, but on the default path it is worth **2.824e-2** on
`stats/peakintensity` of the gate's taper case (adaptive) against the untabulated answer,
where the rest of that case's statistics move by 6e-4; with the source called instead it is
1.175e-3, in line with the rest. So the table now returns the true value there. Nothing inside a propagation can reach that branch —
`Luna.run` rebuilds the table over everything the stepper can ask about — so it costs one
host call per out-of-range diagnostic and nothing per stage. The *operator* table still
holds and warns, which is deliberate: the propagator adds the secant term whatever the
readback returns.

## What it does to the answer

0.1 m Ar capillary, Kerr only, `boundary=:none`, fixed steps, largest relative difference
of `Eω` per save from a 10240-step run of the old one-point rule:

| | 20 steps | 80 | 320 | 1280 |
| --- | --- | --- | --- | --- |
| gradient 0 → 1 bar, one-point (before) | 7.512e-2 | 1.956e-2 | 4.885e-3 | 1.115e-3 |
| gradient 0 → 1 bar, `:tabulated` | 1.581e-4 | 1.604e-4 | 1.608e-4 | — |
| gradient 0 → 1 bar, `:quadrature` | 1.580e-4 | 1.604e-4 | 1.607e-4 | — |
| taper 75 → 50 µm, one-point (before) | 4.096e-1 | 1.017e-1 | 2.478e-2 | 5.592e-3 |
| taper 75 → 50 µm, `:tabulated` | 7.990e-4 | 7.990e-4 | 7.988e-4 | — |
| taper 75 → 50 µm, `:quadrature` | 7.988e-4 | 7.988e-4 | 7.988e-4 | — |

The one-point rows are first order (ratios 3.84, 4.00, 4.38 and 4.03, 4.10, 4.43 for a 4×
step reduction); the other two do not move with the step count at all, and what is left of
them is mostly the reference's own error. The two new paths agree with each other far
below either:

| tabulated vs quadrature, same step count | 20 | 80 | 320 |
| --- | --- | --- | --- |
| gradient | 6.92e-8 | 6.27e-8 | 1.02e-7 |
| taper | 2.52e-7 | 2.01e-7 | 2.09e-8 |

and the difference from the removed rule at the same step count is:

| new vs one-point, same step count | 20 | 80 | 320 |
| --- | --- | --- | --- |
| gradient 0 → 1 bar | 7.528e-2 | 1.972e-2 | 5.046e-3 |
| taper 75 → 50 µm | 4.103e-1 | 1.025e-1 | 2.558e-2 |

(`scratchpad/scratch-27/conv.jl`, which reproduces gpu/23's method with the removed
propagator written as an `AbstractIntegratedLinop` — three lines, and it is the same
arithmetic bit for bit, including `exp(L·(t1 − t2))` for the backward direction.)

## Regression gate

`LUNA_REGRESSION_BASE=b641025b`, 22 cases in two run modes, the environment
`test/regression/README.md` prescribes (`:estimate`, one FFTW thread, one BLAS thread, no
wisdom, `julia -t 1`).

```
LUNA_REGRESSION_BASE=b641025b julia --project=$PWD -t 1 test/test_regression.jl
```

**443 pass, 23 fail.** Every case other than the two below is `0.000e+00` in both the `Eω`
and the statistics class, in both run modes. **These two are not re-baselined here** — the
integrator re-baselines at `gpu/int-E2`; the rows are the correction.

| case | mode | `Eω` diff | `Eω` tol | stats diff | stats tol | worst |
| --- | --- | --- | --- | --- | --- | --- |
| `gradient_field_kerr` | fixed | 7.968e-03 | 1.0e-12 | 6.733e-06 | 1.0e-12 | `Eω [save 11/11]` |
| `gradient_field_kerr` | adaptive | 7.275e-03 | 2.6e-07 | 4.899e-04 | 1.6e-04 | `Eω [save 11/11]` |
| `taper_field_kerr` | fixed | 9.228e-02 | 1.0e-12 | 5.085e-04 | 1.0e-12 | `Eω [save 11/11]` |
| `taper_field_kerr` | adaptive | 9.019e-02 | 8.7e-07 | 1.175e-03 | 3.1e-05 | `Eω [save 11/11]` |

Per quantity:

```
gradient_field_kerr  fixed     Eω 7.968e-03  energy 2.286e-08  peakintensity 6.733e-06
                               peakpower 6.733e-06  fwhm_t_max 5.451e-06
                               fwhm_t_min 5.451e-06  ω0 1.688e-10
gradient_field_kerr  adaptive  Eω 7.275e-03  density 4.206e-04  pressure 4.203e-04
                               zdw 4.899e-04
taper_field_kerr     fixed     Eω 9.228e-02  energy 4.831e-04  peakintensity 5.085e-04
                               peakpower 4.971e-04  fwhm_t_max 2.160e-05
                               fwhm_t_min 2.160e-05  ω0 2.492e-06
taper_field_kerr     adaptive  Eω 9.019e-02  energy 6.059e-04  peakintensity 1.175e-03
                               peakpower 6.180e-04  zdw 3.763e-04
```

(The gradient's adaptive `density`/`pressure`/`zdw` differences are the statistics of the
last accepted step, which the two runs take to slightly different `z`; the gate excludes
`stats/z` itself from the adaptive comparison and reported no step-count failure, so the
step *count* is the same in both.)

`stats/z` and `stats/dz` are excluded from the adaptive comparison by the gate, and no
step-count failure was reported, so the two runs took the same number of steps.

## Tests

`julia -t 1`, `Luna.set_fftw_mode(:estimate)`, one FFTW thread, one BLAS thread, no wisdom,
M1 Pro, Julia 1.13.0.

| file | result | how to run |
| --- | --- | --- |
| `test/test_linops.jl` | **338 pass, 0 fail** (68 of them new) | `julia --project=$PWD -t 1 -e 'using Luna; include("test/test_linops.jl")'` |
| `test/test_device.jl` | **1183 pass, 0 fail** (35 of them new) | as above, from an environment with `JLArrays` |
| `test/test_metal.jl` | **726 pass, 0 fail** | from an environment with `Metal`, `using Luna, Metal` first |
| `test/test_interface.jl` | **362 pass, 0 fail** | as `test_linops.jl` |
| `test/test_rk45.jl`, `test_gradient.jl`, `test_tapers.jl`, `test_boundaries.jl` | **all pass** | as `test_linops.jl` |
| `test/test_freespace.jl`, `test/test_radialgrid.jl` | **236 pass, 0 fail** between them | as `test_linops.jl` |

New in `test_linops.jl` (`@testset "integrated linear operators"`): the three sources of
`Φ` — the table, the quadrature and a caller-written analytic `Φ` — against a closed form
and against each other, on a synthetic `-im*(k0 + k1√z) - k2` operator, which is not an
arbitrary choice (a capillary filled by `Capillary.gradient` from `p₀ = 0` has a density
∝ √z and a dilute gas's propagation constant is linear in density, so that is the shape a
pressure gradient's operator has); `derivative!` for all three; one propagator method for
all three, its backward direction as the inverse, and its caches; the table and the
quadrature on Luna's own gradient and taper; the quadrature's cost (15 evaluations for a
step, none for a zero-length one); `OffsetLinop`, `Boundaries.addloss` on an integrated
operator and `clampdecay` refusing one; `RK45.make_prop!` refusing a callable; and
`linop_integral`'s validation and deprecation.

New in `test_device.jl`: `:quadrature` on `JLArray` for a gradient and a taper, matching
the host to 1e-10 and the tabulated run of the same propagation to 1e-6, with `β` shown
still going through the `HostMirror`; and a caller-supplied `AbstractIntegratedLinop`
(`L(z) = L0 + L1·z`, `Φ` in closed form, built from a real capillary operator) checked
against the closed form on host and device, propagated through `Luna.run` on both, and
`simulation_type/linop == "integrated"` with nothing tabulated for it.

`test_device.jl` also carries `OnePointLinop`, the removed propagator written as an
`AbstractIntegratedLinop`, which is what the convergence testset compares against and what
the "no host work per stage" testset uses for the untabulated side. It is a reference, not
a recommendation, and it is a good demonstration that the interface can express an
arbitrary discretisation of the linear step.

`test_metal.jl` runs the gradient with `:quadrature` on Metal — the operator is integrated
on the host in `Float64` and uploaded as `ComplexF32` — both in `metalgradientcase` (where
it is also the per-stage-host-work side of the call-counting check) and through
`prop_capillary`.

### Tests whose expectations changed

Three, all for the same reason: they compared a z-dependent operator which happens to be
**constant in z** against the constant-operator path elementwise, and were bit-identical
because the one-point rule evaluated the closure and used the result as it stood. A
two-node table holds that constant to within the rounding of a Simpson sum divided by the
span, which is enough for the adaptive controller to accept a slightly different step
somewhere.

- `test_gradient.jl` (`field` and `envelope`): `all(Eω_grad .≈ Eω_const)` elementwise
  becomes the normalised maximum difference per save, `< 1e-12`. **Measured: 7.1e-15.**
  The elementwise form was comparing spectral tails thirty orders below the peak.
- `test_tapers.jl` (`const vs afun`): the same change, same threshold.
- `test_rk45.jl`: `all(abs2.(Aarrp) .== abs2.(Aarrpf))` becomes `isapprox(..., rtol=1e-4)`.
  `Linfunc` is a constant operator written as a closure, now wrapped in a `TabulatedLinop`;
  ~1e-12 rad of phase difference over two soliton periods of an N = 5 soliton at
  `rtol=1e-8` moves the answer by of order `rtol`. The same test also gained
  `@test_throws ArgumentError RK45.make_prop!(Linfunc, Aω)`.

No test's *tolerance* was loosened to accommodate the new propagator; these three changed
the *metric* from elementwise to normalised-per-save, which is what the regression gate
uses.

One more, for the `TabulatedScalar` change: `test_interface.jl`'s memoised-cache test
asserted `t20 == t5` (the cache holds the tables' nodes and nothing the propagation added).
It is now `t20 <= t5 + 1`, because the statistics recorded past the end of the fibre call
`Modes.Aeff` at one `z` per run and the two runs end at different ones. Measured: 33 and
34, against 30-odd growing to 60-odd without a table.

## Metal

`metalenv-27` (Luna developed in, plus Metal, Test, GPUArraysCore, Adapt), M1 Pro,
`julia -t 1`. Largest relative difference of `Eω` per save.

`prop_capillary(125e-6, 0.1, :He, (1.0, 0.0); saveN=5, …)`, adaptive:

| | Metal vs CPU F32 | Metal vs CPU F64 | CPU F32 vs F64 |
| --- | --- | --- | --- |
| `:tabulated` | 3.715e-06 | 4.028e-06 | 1.266e-06 |
| `:quadrature` | 3.857e-06 | 4.001e-06 | 1.362e-06 |

Low-level, 0.1 m, 20 fixed steps, `boundary=:none`:

| | Metal vs CPU F32 | Metal vs CPU F64 | CPU F32 vs F64 |
| --- | --- | --- | --- |
| gradient `p₀ = 0`, `:tabulated` | 8.181e-07 | 1.065e-06 | 1.285e-06 |
| gradient `p₀ = 0`, `:quadrature` | 1.048e-06 | 9.884e-07 | 7.081e-07 |
| taper, `:tabulated` | 5.152e-07 | 1.875e-06 | 1.666e-06 |
| taper, `:quadrature` | 9.150e-07 | 5.458e-07 | 5.788e-07 |

The tabulated rows reproduce gpu/23's (3.72e-6 / 4.03e-6 and 8.18e-7) exactly, as they
must: that arithmetic is unchanged. The quadrature path sits at the same level, which is
the point — it is the same integral computed on the host and uploaded, so what is left is
the `Float32` rounding of the state.

## Cost

`benchmark/tabulated.jl` now compares `:quadrature` against `:tabulated`, which unlike its
previous comparison is like for like (the two compute the same integral to 1e-7 of each
other). The qualitative picture from gpu/23 is unchanged and the quadrature side is
*slower* than the one-point rule it replaced in the "no table" column, by roughly the ratio
of 15 operator evaluations to one per stage; `:tabulated` is the default for that reason.

The tabulated propagator costs one extra field-sized broadcast per `(t1, t2)` pair relative
to nothing, exactly as in gpu/23 — the `AbsolutePhase` propagator is that code with the
table replaced by the interface, so there is no per-step cost to the generalisation.
`PhaseStyle` and `isnothing(secant(op))` are resolved when the propagator is built.

## Known gaps and risks

- **A free-space or multimode run with a z-dependent operator now builds a table by
  default.** The table is `2·nnodes` copies of an operator which, for those geometries, is
  the size of the whole state. `TabulatedLinop` reports its size and warns above 256 MB
  (`LinearOps.TABLE_WARN_BYTES`), and `linop_integral=:quadrature` holds no table at all —
  but a 3-D free-space script which used to allocate nothing extra can now allocate a lot,
  and the warning is the only thing standing in front of that. Worth a decision at
  integration: whether the default should depend on the size of the operator.
- A caller-supplied `AbstractIntegratedLinop` cannot be used for a free-space propagation
  with any `boundary` mode, because the evanescent clamp is not linear in the operator.
  It raises rather than propagating something wrong. The callable path is unaffected.
- `QuadratureLinop` allocates on the host on every call: `QuadGK.quadgk!` allocates six
  buffers the size of the operator plus its segment buffer. Reusing them means holding a
  `QuadGK.InplaceIntegrand`, which is internal to that package, and the path is the slow
  one by construction; not done.
- `linop_tol` is still one knob for the operator (radians, absolute) and the value tables
  (relative), and it is now also the quadrature's tolerance. Splitting it would be easy.
- The `β` and `Aeff` value tables are still linear interpolants (gpu/23's gap), and are
  still tied to `:tabulated`.
- `secant(op)` takes no prototype, so an operator which returns a host array has it
  converted once by `_seclike` when the propagator is built. That is right for a
  caller-written operator and for `OffsetLinop`; `TabulatedLinop` already returns one on
  the state.

## Deviations from GPU_PLAN.md and the brief

- **The interface is `phase!` + `secant` *or* `phasediff!`, chosen by a `PhaseStyle`
  trait**, rather than `phase!` + `secant` alone. The brief specifies both `phase!(out, op,
  z)` as the interface *and* that `QuadratureLinop` computes `Φ(t2) − Φ(t1)` per call over
  `[t1, t2]`, and those two cannot both hold: `phase!` needs an origin, and integrating
  from a fixed origin at every stage costs O(z) evaluations and forms a small difference
  out of two large numbers — the very thing the table's secant subtraction exists to avoid
  (measured: `Φ` over 0.1 m of gradient is 3.4e1 rad, `ΔΦ` over a step 1.7e-1). There is
  still one public `RK45.make_prop!` method for `AbstractIntegratedLinop`; the trait
  selects between two internal closures, so the `AbsolutePhase` one is gpu/23's code with
  the same buffers, caches and fused exponent, and is bit-identical to it.
  `QuadratureLinop` implements `phase!` as well (as `phasediff!(z0, z)`), so all three
  sources can be compared through the same function in tests.
- **`RK45.make_prop!(linop!, y0)` is not deleted outright but raises.** Deleting it gives a
  `MethodError` listing every `make_prop!` method to a low-level caller of
  `RK45.solve_precon`, which is a worse diagnostic than a sentence naming the two wrappers
  and the keyword. No one-point propagator is left in `src`.
- **`OffsetLinop` and the `Boundaries` methods are not in the brief.** Without them a
  caller-supplied integrated operator is unusable with the default `boundary=:rate`,
  because `addloss` would call it like a closure. A constant added to `L` is a straight
  line in `Φ`, so it is four lines and exact.
- **The `TabulatedScalar` out-of-range change is not in the brief** (see above). It
  reverses a documented gpu/23 decision because the default changed underneath it.
- **`Stats.zdw_linop` needed no change.** Its methods dispatch on `linop::AbstractArray`
  (constant → `ConstZDW`) against everything else (z-varying → `zdw(mode_s)`, which
  evaluates `Modes.zdw` at each z and never calls the operator). An
  `AbstractIntegratedLinop` is not an `AbstractArray`, so it takes the z-varying method,
  which is correct. `Stats.jl` is therefore untouched, which also keeps it disjoint from
  `gpu/25`. `Fields.PropagatedField`'s `propagator!` is a user-supplied function with no
  relation to `linop`. What both would need if they did use the operator —
  `derivative!(out, op, z)` — is part of the interface and is tested.
- **The staging of the commits.** The brief asks for a commit per stage (interface,
  quadrature, keyword, tests, docs). The first three touch the same three functions in the
  same three files and were written together; they are one commit, with tests and docs
  following. The stages were still gated in order: the regression gate was run on the
  interface change alone (the same two cases, the same `Eω` numbers) before the keyword
  went in, and again on the finished tree.

## Files

```
src/LinearOps.jl     +~370  the interface, QuadratureLinop, OffsetLinop, TabulatedLinop's
                            secant/derivative!, the TabulatedScalar fallback
src/RK45.jl           ~-30  the one-point rule out, the error and the docstrings in
src/Luna.jl           ~+40  linop_integral, _linop_integral, linoptype, the run docstring
src/Interface.jl      ~+20  the keyword, saveargs, the aefftol gate
src/Boundaries.jl     ~+25  addloss/addloss_k/clampdecay for an integrated operator
src/NonlinearRHS.jl    ~8   docstrings only
src/Device.jl          ~5   docstring only
test/                +~330  see Tests
docs/                +~150  docs/src/gpu.md user section, developer/device_model.md
benchmark/            ~20   tabulated.jl compares :quadrature against :tabulated
```
