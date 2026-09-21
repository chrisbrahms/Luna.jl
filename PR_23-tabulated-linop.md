# The tabulated linear operator

Branch `gpu/23-tabulated-linop`, base `gpu/int-D` (`90826dc4`). GPU_PLAN.md §2 (the PR 440
bullet), §4.5 layer 2, §4.6, §6 Group E, §11.

A tapered or pressure-graded capillary is the one mode-averaged configuration which still
runs host code inside every stage of every step: the linear operator is a scalar loop over
`Modes.neff` evaluated into a host buffer and copied to the device (`RK45.make_prop!`), the
propagation constant `β(z)` is evaluated and uploaded through a `Luna.HostMirror` at every
right-hand side (`NonlinearRHS.NormModeAvg`), and the effective area `Aeff(z)` is a
memoised cubature which grows a `Dict` entry per distinct `z` for a non-Marcatili mode.
This branch replaces all three with tables over `z`, built at setup and resident wherever
the state is, behind the opt-in keyword `tabulate_linop=true`.

It also changes the discretisation of the linear step, deliberately and for the better: the
propagator becomes the exact `exp(∫linop dz)` over the step instead of a one-point rule.
That is why it is opt-in; the numbers are in "What it does to the answer" below.

**The regression gate is 0 failures and exactly `0.000e+00` on every case against both
`fa556e6f` (source-identical to the base) and `fdf8dbe3` (`evanescent`).** Tabulation is
off by default and nothing on the default path changed.

## What changed

### `src/LinearOps.jl` — the tabulation

Three new types and the bisection which builds them.

`TabulatedLinop(linop!, proto, z0, z1; tol, maxdepth, maxnodes)` holds the integrated
operator `Φ(z) = ∫_{z0}^{z} linop dz'` on adaptively placed nodes, as
`(size(linop)..., nnodes)` arrays on `proto`'s array type and precision, together with the
operator itself at the nodes. `phase!(out, tab, z)` reads it back with a cubic Hermite
interpolant in one broadcast over four node slices with four scalar weights;
`integrated!(out, tab, z)` gives the integral itself and is used only for checking.

`RK45.make_prop!(tab::TabulatedLinop, y0)` (the method lives in `LinearOps`, because
`RK45.jl` is included before it) is the propagator `y *= exp(Φ(t2) − Φ(t1))`. `Φ(t1)` is
read once per step — the six stages share one `t1` — and `Φ(t2)` once per distinct `t2`,
which is the same argument the untabulated propagator's last-`t2` cache rests on. Two
field-sized device buffers, one broadcast each, and no host traffic.

`TabulatedVector` (for `β(z)`) and `TabulatedScalar` (for `Aeff(z)`) tabulate quantities
which are needed as values rather than as integrals. The mode interface gives no derivative
of either, so they are interpolated linearly rather than with a Hermite, with the
tolerance taken relative to the largest value in the table (`β` is ~1e7 in SI units and
`Aeff` ~1e-8; one absolute tolerance cannot serve both). `TabulatedVector` owns its output
buffer and skips the readback when `z` has not changed.

**Node placement.** An interval is accepted when the interpolant that will actually be read
back agrees, at the interval's midpoint, with a directly computed value to the tolerance,
and is bisected otherwise. For the integral the comparison is Hermite-at-the-midpoint
against Simpson over the left half, which is PR 440's criterion
(`scratchpad/pr440/src/LinearOps.jl:494-680`) generalised from `∫imag(linop) dz` to the
full complex operator. Luna's z-dependent quantities are not smooth — `Capillary.gradient`
goes as `√(p₀² + z/L(p₁² − p₀²))`, which has a `1/√z` cusp in its derivative at the
entrance when `p₀ = 0`, and a multi-section fill has a derivative discontinuity at every
junction — and bisection costs a number of intervals proportional to the depth a kink is
resolved to rather than to a resolution imposed everywhere.

**The secant subtraction** is the one thing that is not in PR 440. `TabulatedLinop` stores

```
Φ̃(z) = Φ(z) − L̄·(z − z0),   L̄ = Φ(z1)/(z1 − z0)
```

and the propagator adds `L̄·(t2 − t1)` back in the same broadcast, formed from the step
length exactly as the constant propagator forms it. In exact arithmetic the result is
unchanged; what changes is the size of the stored numbers, which is what a `Float32` table
rounds. `make_linop` already subtracts the frame, so `max|Φ|` is 34 rad over 0.1 m of
gradient and 359 rad over 1 m — tens to hundreds of radians, not thousands — and `max|Φ̃|`
is 11.3 times smaller over the same span. Rounding `ΔΦ` to `Float32` both ways (review
round 1's measurement, `t1 = L/2`, one step of `L/20`):

| span | max‖Φ̃‖ | max‖Φ‖ | ‖ΔΦ‖ | error with the secant | without |
| --- | --- | --- | --- | --- | --- |
| 0.1 m | 3.26 | 35.9 | 1.71 | 2.1e-7 | 1.3e-6 |
| 1 m | 32.6 | 359 | 17.1 | 2.0e-6 | 1.2e-5 |
| 10 m | 326 | 3593 | 171 | 2.7e-5 | 1.0e-4 |

So it buys a factor of 4 to 11 in the rounding error of `Φ(t2) − Φ(t1)`, not the difference
between working and not working: the un-subtracted version would also have been usable at
these spans. It is kept because it costs nothing — a cubic Hermite reproduces a linear
function exactly, so the node placement and the interpolation error are unchanged — and
because it is exactly zero for a z-independent operator, which then stores nothing at all
(`tab.scale < 1e-6`, asserted in `test_linops.jl`). An earlier version of this description
claimed the un-subtracted table would be unusable in `Float32`; that was not measured and
is not true.

### `src/NonlinearRHS.jl` — the transform's own z-dependent quantities

`NormModeAvg`'s `β` branch goes through `_βdev(n.β, n.βfun!, z)`, which either stages and
uploads through the `HostMirror` as before or calls a `TabulatedVector`. `tabulate(transform,
z0, z1, tol, proto)` returns the transform with `aeff` and `norm!` tabulated; the generic
method returns it unchanged, so only the mode-averaged transform (the only device-capable
one so far) is affected.

### `src/Luna.jl` — where the tables are built

`run` takes `tabulate_linop=false` and `linop_tol=1e-6` and builds both tables after
`Boundaries.setup`, so the operator includes the absorber and the evanescent clamp, and
over `[z0, zmax + max(max_dz, init_dz)]` with the absorber's `max_dz` — `RK45.solve` runs
`while tn <= tmax` and so overshoots `zmax` by up to one step, and that step's stages are
what the last saved plane is interpolated from (GPU_PLAN.md review finding 9; PR 440's
fixed 5 % margin is not enough for `boundary=:none`). `init_dz` is in there because the
first step is taken at `init_dz` before `steplims!` can clamp it, and `Boundaries.setup`
only reduces it to `max_dz` for `boundary=:rate`. A constant operator is left alone: it
is already exact in the propagator. `linoptype` reports `"tabulated"` in the output's
`simulation_type` group.

### `src/Interface.jl`

`tabulate_linop` and `linop_tol` are declared on `prop_capillary_args`, recorded by
`saveargs` and forwarded to `Luna.run` by `boundary_kwargs`. No positional or existing
keyword changed.

## What it does to the answer

The untabulated propagator is `exp(linop(t2)·(t2 − t1))`, a one-point rule whose error is
first order in the step size. The tabulated one is the exact interaction-picture propagator
of the linear part. Both converge to the same solution; they do not converge at the same
rate. Largest relative difference of `Eω` per save from a 10240-step run of the same kind,
0.1 m Ar capillary, Kerr only, `boundary=:none`, fixed steps:

| | 20 steps | 80 | 320 | 1280 |
| --- | --- | --- | --- | --- |
| gradient 0 → 1 bar, default | 7.51e-2 | 1.96e-2 | 4.89e-3 | 1.12e-3 |
| gradient 0 → 1 bar, tabulated | 1.58e-4 | 1.60e-4 | 1.61e-4 | 1.61e-4 |
| taper 75 → 50 µm, default | 4.10e-1 | 1.02e-1 | 2.48e-2 | 5.59e-3 |
| taper 75 → 50 µm, tabulated | 7.99e-4 | 7.99e-4 | 7.99e-4 | 7.99e-4 |

The tabulated runs are at the well resolved answer from the coarsest step count tried, and
the 10240-step untabulated reference still differs from the tabulated one by 1.61e-4
(gradient) and 7.99e-4 (taper) — the residual of its own first-order sequence. So the
"difference from the exact path" asked for in the brief is, at a usable step count, mostly
the *untabulated* path's error:

| tabulated vs untabulated, same step count | 20 steps | 80 | 320 |
| --- | --- | --- | --- |
| gradient 1 → 0 bar | 8.99e-2 | 2.16e-2 | 5.28e-3 |
| gradient 0 → 1 bar (`p₀ = 0` entrance) | 7.53e-2 | 1.97e-2 | 5.05e-3 |
| taper 75 → 50 µm | 4.10e-1 | 1.03e-1 | 2.56e-2 |

`linop_tol` does not enter these numbers: at 20 steps, `1e-4`, `1e-6` and `1e-8` give
8.993e-2, 8.993e-2 and 8.993e-2 for the first row. The tolerance would have to be very
loose before the table's own error was visible next to the discretisation change.

The step-size controller does not close the gap on its own. The propagator's error is
common to both of the embedded Runge–Kutta solutions the error estimate is formed from, so
it cancels out of the estimate: an adaptive run at `rtol=1e-6` with `boundary=:rate`
differs by 6.86e-2 (gradient) and 3.35e-1 (taper) between the two paths.

**This is a change of discretisation, not of the equation.** A result computed with
`tabulate_linop=true` is not comparable element by element with a stored one computed
without it, which is why it is opt-in and why the gate keeps the default path at exactly
zero.

## Cost

Table sizes and build cost, 0.1 m, `linop_tol=1e-6`, 1025-point ω grid:

| case | linop nodes | linop! evaluations | β nodes | Aeff nodes |
| --- | --- | --- | --- | --- |
| gradient 1 → 0 bar | 55 | 217 | — | 2 |
| gradient 0 → 1 bar | 57 | 225 | 23 | 2 |
| taper 75 → 50 µm | 33 | 129 | — | 257 |

The cusp is resolved by bisection rather than by refining everywhere: the `p₀ = 0` entrance
puts 10 of its 57 nodes in the first 1 % of the fibre and its first interval is 1.3e-5 m,
against 3.3e-3 m for the smooth taper. `Aeff` of a linear taper needs 257 nodes because the
value tables are interpolated linearly and `Aeff ∝ a(z)²`; they are 257 `Float64`s.

### Benchmarks

`benchmark/tabulated.jl`, M1 Pro, Julia 1.13.0, 1 Julia thread, 1 FFTW thread, 1 BLAS
thread, `:estimate`, no wisdom. Ar, 0.1 m, 20 fixed steps, `boundary=:none`,
`linop_tol=1e-6`. `prop` is the whole `Luna.run`; the tabulated column therefore *includes*
the table build, which is listed separately.

| case | device | trange | state | prop | prop (tab) | speedup | table setup |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gradient | CPU Float64 | 400 fs | 1025 | 11.81 ms | 14.53 ms | 0.81× | 8.01 ms |
| gradient | CPU Float32 | 400 fs | 1025 | 10.84 ms | 13.52 ms | 0.80× | 7.95 ms |
| gradient | Metal Float32 | 400 fs | 1025 | 164.9 ms | 75.03 ms | **2.20×** | 9.83 ms |
| gradient | CPU Float64 | 1600 fs | 4097 | 46.41 ms | 57.99 ms | 0.80× | 31.49 ms |
| gradient | CPU Float32 | 1600 fs | 4097 | 42.42 ms | 53.05 ms | 0.80× | 32.48 ms |
| gradient | Metal Float32 | 1600 fs | 4097 | 199.96 ms | 108.53 ms | **1.84×** | 32.33 ms |
| taper | CPU Float64 | 400 fs | 1025 | 37.08 ms | 25.91 ms | 1.43× | 19.07 ms |
| taper | CPU Float32 | 400 fs | 1025 | 35.51 ms | 24.81 ms | 1.43× | 19.36 ms |
| taper | Metal Float32 | 400 fs | 1025 | 197.87 ms | 87.97 ms | **2.25×** | 20.47 ms |
| taper | CPU Float64 | 1600 fs | 4097 | 146.41 ms | 103.51 ms | 1.41× | 75.44 ms |
| taper | CPU Float32 | 1600 fs | 4097 | 142.27 ms | 98.47 ms | 1.44× | 75.83 ms |
| taper | Metal Float32 | 1600 fs | 4097 | 303.32 ms | 148.51 ms | **2.04×** | 77.86 ms |

Metal is 1.8–2.3× faster with tabulation *including* the table build; the stepping alone
goes from 164.9 ms to 65 ms for the 400 fs gradient, because the six host evaluations and
six uploads per step are gone. On the CPU the taper gains 1.4× outright and the gradient
loses 20 % at 20 steps: `Capillary.neff_β_grid` caches the waveguide term for a fixed-radius
mode, so a gradient's `linop!` is cheap (36 µs against 242 µs for the taper, which has to
redo the waveguide term at every `z`) and the table costs about as much to build as it
saves. The build is a fixed cost, so it is amortised by step count: measured separately
over 320 steps at 0.1 m, the gradient run goes from 0.18 s to 0.09 s and the taper from
0.57 s to 0.10 s.

A device run is slower than the CPU at these sizes either way — a mode-averaged single
column is launch-bound, which `benchmark/device.jl` already documented. What this branch
changes is that the gap is now the launch overhead alone, rather than the launch overhead
plus a host round trip per stage.

### On Metal

`prop_capillary(125e-6, 0.1, :He, (1.0, 0.0); ...)`, largest relative difference of `Eω`
per save:

| | Metal vs CPU Float32 | Metal vs CPU Float64 |
| --- | --- | --- |
| `tabulate_linop=false` | 9.94e-6 | 2.04e-4 |
| `tabulate_linop=true` | 3.72e-6 | 4.03e-6 |

The second row is the exit criterion. The Float64 column improves by fifty times because
the untabulated gradient's 2e-4 came from the adaptive controller taking slightly different
decisions in `Float32` — a pressure gradient's step sequence is sensitive to rounding — and
with an exact linear step there is much less for it to respond to. Low-level fixed-step
runs, 0.1 m, 20 steps, `boundary=:none`:

| | Metal vs CPU F32 | Metal vs CPU F64 | CPU F32 vs CPU F64 |
| --- | --- | --- | --- |
| gradient `p₀ = 0`, tabulated | 8.18e-7 | 1.07e-6 | 1.29e-6 |
| taper, tabulated | 5.37e-7 | 1.88e-6 | 1.67e-6 |
| gradient, untabulated (reference) | 7.09e-7 | | |

The tabulated device path agrees with the tabulated host path exactly as well as the
untabulated one does, and the `Float32` tables cost nothing measurable next to the
`Float32` state itself — which is what the secant subtraction is for.


## Regression gate

21 cases in two run modes, the environment `test/regression/README.md` prescribes
(`:estimate`, one FFTW thread, one BLAS thread, no wisdom, `julia -t 1`), run on the final
tree.

| baseline | result | largest difference |
| --- | --- | --- |
| `fa556e6f` (source-identical to the base `90826dc4`) | **460 pass, 0 fail** | **0.000e+00**, every case, both modes |
| `fdf8dbe3` (`evanescent`) | **460 pass, 0 fail** | 7.240e-05 (`multimode_field_plasma` adaptive statistics) |

The first row is the gate for this branch and it is exactly zero, as it has to be:
tabulation is opt-in and nothing on the default path changed. The second row is the
project's accumulated distance from `evanescent` and is identical, row for row, to what
`PR_int-D.md` recorded for the base — this branch adds nothing to it.

```
LUNA_REGRESSION_BASE=fa556e6f julia --project=$PWD -t 1 test/test_regression.jl
LUNA_REGRESSION_BASE=fdf8dbe3 julia --project=$PWD -t 1 test/test_regression.jl
```

The tabulated deltas are documented above ("What it does to the answer") and are
deliberately *not* gated: they are a different discretisation, not a rounding difference.

## Tests

All run with `Luna.set_fftw_mode(:estimate)`, `set_fftw_threads(1)`,
`BLAS.set_num_threads(1)` and `julia -t 1` on an M1 Pro under Julia 1.13.0.

| file | result | how to run |
| --- | --- | --- |
| `test/test_linops.jl` | **270 pass, 0 fail** (74 of them new) | `julia --project=$PWD -t 1 -e 'using Luna; include("test/test_linops.jl")'` |
| `test/test_device.jl` | **655 pass, 0 fail** (37 testsets) | as above, from an environment with `JLArrays` |
| `test/test_metal.jl` | **414 pass, 0 fail** (22 of them new) | from an environment with `Metal`, `using Luna, Metal` first |
| `test/test_interface.jl` | **351 pass, 0 fail** (11 of them new) | as `test_linops.jl` |
| `test/test_gradient.jl`, `test/test_tapers.jl` | **9 pass, 0 fail** between them | as `test_linops.jl` |

New in `test/test_linops.jl` (`@testset "tabulated linear operator"`, 62 tests): the table
against an operator whose integral is known in closed form and which has the same `√z`
awkwardness Luna's own have, at three tolerances; that a tighter tolerance gives more nodes
and a smaller error; the secant properties; the propagator against the closed form, its
backward direction as the inverse, and the `t1`/`t2` caches giving the same numbers as a
cold readback; the same tables built from Luna's own gradient and taper operators against a
reference table built to `1e-10`; where the nodes went (the `p₀ = 0` entrance resolves its
first interval to 1.3e-5 m and puts more than five nodes in the first 1 % of the fibre, the
taper does neither); a z-independent closure giving two nodes and reproducing the constant
propagator; and the `β`/`Aeff` tables against direct evaluation.

New in `test/test_device.jl`: the tables on `JLArray` for a gradient and a taper, with the
device run matching the host one to 1e-10 as the untabulated pair does; the "no host work
per stage" check; and the discretisation comparison above, as a convergence-rate assertion
rather than a fixed number.

New in `test/test_interface.jl`: that `prop_capillary(...; tabulate_linop=true)` hands
`Luna.setup` a table rather than the bare callable, so the statistics built from
`transform.aeff` hold one; and that the memoised `Modes.Aeff` cache stops growing with the
step count, measured on a delegated (non-Marcatili) tapered mode, whose `Aeff` really does
go through the memoised cubature. The cache holds 31 entries after a 6-step untabulated run
and 86 after a 17-step one; with tabulation it holds 33 — the two tables' nodes — after
both.

New in `test/test_metal.jl` (`@testset "tabulated operator on Metal"`, 22 tests): that
`Φ`, `dΦ`, the secant and the `β` table are `MtlArray{ComplexF32}`/`MtlArray{Float32}` —
which is the only place a stray `Float64` in them would show, since Metal's compiler
rejects one — that the readback and the propagator run as kernels, the Metal/CPU-Float32
agreement, and the same call-counting check. `prop_capillary` with `tabulate_linop=true` is
in `@testset "prop_capillary on Metal"`.

**The call-counting check** is how "no per-stage host work" is established. `CountCalls`
wraps `linop!`, `βfun!`, `aeff` and `densityfun` (which is called exactly once per
right-hand side, so it counts stages). Two runs of the same propagation with different step
counts but the same `max_dz` — so the tables span the same interval and are built from the
same evaluations — must give identical `linop`/`β`/`aeff` counts while the right-hand side
count moves. Measured on the gradient at 1 cm, between a 20-step run and one which took 38
accepted steps:

| host calls | `linop!` | `βfun!` | `Aeff` | right-hand sides |
| --- | --- | --- | --- | --- |
| `tabulate_linop=false` | 104 → 194 | 121 → 229 | 242 → 458 | 121 → 229 |
| `tabulate_linop=true` | 133 → 133 | 35 → 35 | 3 → 3 | 121 → 229 |

Without tabulation the counts scale with the stages; with it they do not move at all.


## Known gaps

- `Luna.run` tabulates into a transform of its own and does not modify the caller's. A
  low-level caller who builds statistics by hand from `transform.aeff` therefore keeps the
  untabulated callable in them; `prop_capillary` does not (review round 1, finding 3: it
  tabulates `Aeff` itself before `Stats.default` closes over it), and the recipe for a
  low-level caller is the same one line — pass a `LinearOps.TabulatedScalar` to
  `Luna.setup` as `aeff`.
- `prop_capillary`'s `Aeff` table spans `[0, flength]`, and the statistics of the last step
  are recorded a fraction of a step past the end of the fibre, where the table holds its
  end value. That is a scale factor on a diagnostic, and the propagation uses a table built
  over the full span, so nothing propagated is clamped.
- Tabulating a multimode or free-space operator is allowed and correct but expensive: the
  table is `2·nnodes` times the size of the operator, which for those geometries is the
  size of the whole state. The constructor now reports the size and warns above
  `LinearOps.TABLE_WARN_BYTES` (256 MB), naming `linop_tol` and the shape of the operator.
- The value tables (`β`, `Aeff`) are linear interpolants because no derivative is available
  from the mode interface. A quadratic through the bisection's own midpoint would cost
  nothing extra in evaluations and would cut the node count for a smooth quantity by an
  order of magnitude; not done here. It is why a linear taper's `Aeff` takes 257 nodes at
  `linop_tol=1e-6` — 257 `Float64`s, but also 257 cubatures at setup for a mode whose
  `Aeff` is not analytic.
- The relative tolerance of a value table is taken against its two endpoints, which is
  right for `β` and `Aeff` (monotonic in z over any fibre Luna describes) and would not be
  for a quantity small at both ends and large in between. Commented, not guarded.
- `linop_tol` is one knob for two different things (radians for the operator, relative for
  the values). Splitting it would be easy if anyone needs it.

## Changes after review round 1

The review (`scratchpad/reviews/gpu-23-tabulated-linop-1.md`, "approve with minor fixes")
reproduced the convergence table, the node and evaluation counts and the Metal figures
exactly, and found no defect in the numerics. Ten findings, all addressed in
`f75c5fcd`:

1. **A caller-supplied `norm!` made `tabulate_linop=true` throw.** `tabulate` read
   `t.norm!.aeff` before dispatch could help, so a normalisation that is not one of Luna's
   own — which `Luna.setup` documents and `NonlinearRHS.check_norm` accepts — aborted the
   run with a `FieldError`. `_normaeff` returns the field only for `NormModeAvg` and
   `NormModeAvgGNLSE`, the generic `tabulate` takes the `aeff` keyword the known ones take,
   and an unknown normalisation is now passed through unchanged. Tested in
   `test_device.jl` with a wrapper around the standard normalisation.
2. **The `Float32` justification for the secant subtraction was overstated.** Corrected in
   the docstring, `docs/src/developer/device_model.md` and the section above: `max|Φ|` is
   tens to hundreds of radians in the co-moving frame, and the subtraction buys a measured
   factor of 4 to 11 in the rounding of `ΔΦ`, not usability. The measured ratio of
   `max|Φ|` to `max|Φ̃|` over the same span is 11.3, not the 10.5 quoted before (which
   compared a table over `[0, 1.05L]` with `|Φ(L)|`).
3. **The `transform.aeff` gap is closed rather than documented.** `prop_capillary`
   tabulates `Aeff` over `[0, flength]` before `Stats.default` closes over it;
   `TabulatedScalar` keeps its source callable so `NonlinearRHS._aefftab` can rebuild a
   wider table for the propagation, which needs `Aeff` one step past the end of the fibre.
   The statistics keep the narrower table, whose end value is held past the fibre — a scale
   factor on a diagnostic. Tested by the memoised cache's size, above.
4. **No size guard.** The constructor reports the table's size in the `@info` line and
   warns above `LinearOps.TABLE_WARN_BYTES` (256 MB), naming `linop_tol` and the shape of
   the operator. `DEFAULT_MAXNODES` caps nodes, not bytes, and a 10-mode 8192-point
   operator at 57 nodes is about 260 MB.
5. **"Works for every mode type" was untested.** `test_linops.jl` now tabulates a 4-mode
   operator: `size(Φ) == (1025, 4, nnodes)`, the error within tolerance, agreement with a
   20000-panel midpoint quadrature, and a finite propagator. It is the only test of the
   `ndims` generality of `_stack`, `selectdim` and `phase!`.
6. **Going off the end of the table was silent.** `phase!` warns once: the propagator adds
   the secant term whatever the readback returns, so a step outside the table would
   propagate with the mean operator over the whole of it. The value tables hold their end
   value silently, which is what `prop_capillary`'s `Aeff` table does by design; the
   distinction is documented and tested.
7. `quiet` is documented.
8. A comment says why a value table's relative scale comes from its endpoints (`β` and
   `Aeff` are monotonic in z over any fibre Luna describes).
9. `prop_capillary`'s keyword documentation no longer implies a uniform fibre is unaffected:
   reading a constant off a two-node table is `(1-s)f + sf`, not bitwise `f`.
10. `docs/src/gpu.md` says what is left on the host — an interval lookup and four
    interpolation weights — instead of "nothing".

Every test file named above and the gate were rerun after the fixes: `test_linops.jl`
270 pass, `test_device.jl` 655 pass, `test_interface.jl` 351 pass, `test_metal.jl` 414 pass
(unchanged), `test_gradient.jl` + `test_tapers.jl` 9 pass, and the gate against `fa556e6f`
460 pass with `0.000e+00` on every case.

## Deviations from GPU_PLAN.md

- §4.5 says the tables are built "inside `Luna.run` after `Boundaries.setup`", which they
  are, but it does not mention the secant subtraction, which turned out to be necessary for
  `Float32` (see above).
- The plan describes the readback as "one broadcast (Hermite weights are four scalars per
  call)", which is what `phase!` is; the propagator then needs a second broadcast for the
  exponential, as the untabulated one does.
- `RK45.make_prop!`'s new method is defined in `LinearOps`, not in `RK45.jl`, because
  `RK45.jl` is included before `LinearOps.jl` and cannot name the type.
