# The device and precision model

How Luna runs the same propagation code on a host array, on a reduced-precision host array
and on a GPU. The user-facing page is [Running on a GPU](../gpu.md).

!!! note "Work in progress"
    Written for `gpu/10-device-model`, the first branch of the GPU work, and extended by
    each later branch. It is consolidated, with every branch's measurements, in
    `gpu/32-docs`.

## Where a run happens

[`Luna.DeviceSpec`](@ref)`{A, T}` is a singleton carrying the array type `A` (`Array`,
`MtlArray`, `CuArray`, `JLArray`) and the *real* element type `T`; the state is
`Complex{T}`. Because it is a singleton, everything derived from it — the buffer types of
a transform, which planner is used, whether a mirror is a copy or an alias — is resolved at
compile time.

`Luna.settings["device"]` holds what the user asked for: `:cpu` (also the meaning of an
absent key), `:auto`, `:metal`, `:cuda`, or a spec. [`Luna.device`](@ref) resolves it and
[`Luna.set_device`](@ref) sets it, validating immediately.

The building blocks in `src/Device.jl`:

| | |
| --- | --- |
| [`Luna.alloc`](@ref)`(spec, T, dims)` | a zero-filled buffer on the spec's array type |
| [`Luna.todevice`](@ref)`(spec, x)` | a host array moved and converted; the identity on `DeviceSpec(Array, Float64)` |
| [`Luna.tohost`](@ref)`(x)` | a host `Array` holding the same data |
| [`Luna.scalar`](@ref)`(x, s)` | a `Float64` scalar converted to `real(eltype(x))` |
| [`Luna.upload_like`](@ref)`(y, x)` | `x` in the array type and precision of `y` (the linear operator) |
| [`Luna.mask_like`](@ref)`(y, m)` | a boolean mask on `y`'s array type |
| [`Luna.assert_resident`](@ref)`(spec, arrays...)` | the structural residency check every transform makes |

## Kernel discipline

Per-step code consists only of broadcasts, `mul!`, planned FFTs applied with `mul!`,
`mapreduce`/`sum`/`maximum`, `accumulate!` and `copyto!`. No scalar indexing of
field-sized arrays. The rules, in the order they bite:

1. **One implementation per operation.** There is no CPU copy of any kernel. The backend
   trait `Utils.backend` exists, but only for FFT planning (FFTW flags and wisdom against
   the generic `AbstractFFTs` planners), for the host-copy fallbacks, and for residency
   checks. It never selects a kernel.
2. **Fused broadcasts are left-associated in the original order** wherever that costs
   nothing, so that the rounding-level differences from replacing a sequence of passes
   with one stay as small as they can be. In practice this has been free: the RK45 stage
   combines, the error estimate and the dense output are bit-identical to the loops they
   replace, and so is the whole default CPU path.
3. **No `Float64` anywhere in a kernel's arguments.** Metal's kernel compiler rejects any
   `double` which survives optimisation, and its kernel adaptor does not convert scalars.
   Scalars go through [`Luna.scalar`](@ref), vectors are mirrored once, structs are
   parametric in the real type. Rational literals inside a kernel body (`2/3`, `3/4`) are
   `Float64` and have to be converted or folded into a scalar factor.
4. **Residency is asserted structurally at construction.** Every transform checks that its
   buffers, mirrors and coefficient arrays are on the expected array type. `JLArrays` does
   not catch a mixed host/device broadcast, and on some backends it is a silent, enormous
   slowdown rather than an error.
5. **Inverse FFT plans are held explicitly** as `IFT` on every backend, and their `1/N` is
   folded into the scale factor of the oversampling copy
   ([`NonlinearRHS.to_time!`](@ref Luna.NonlinearRHS.to_time!)), one pass fewer. The
   folding is exact whenever `1/N` is a power of two, which covers every transform over
   the time axis alone (Luna's time grids are powers of two). A multi-axis free-space
   transform normalises by `1/(Nt·Nx·Ny)` and the transverse grids accept any `Nx`, `Ny`;
   for a non-power-of-two length the two routes differ at rounding level instead
   (measured: 5.6e-17 for length 24).
6. **Anything which cannot meet the above runs on the host through an explicit copy**, and
   says so once.
7. **Every branch which adds or changes a kernel runs the Metal hardware tests.** That is
   the only reliable detector of a stray `Float64`: `JLArray` and a `Float32` `Array` both
   promote silently.

## Grid mirrors

`Grid.RealGrid`/`Grid.EnvGrid` stay host `Float64` — they are metadata, and they are
written to output files. A host vector cannot be broadcast against a device array, so each
transform holds a [`Luna.GridVectors`](@ref) mirror of the vectors its kernels use (`ω`,
`ωwin`, `twin`, `towin`, and `sidx` as a `Bool` mask), adapted once at construction. On
`DeviceSpec(Array, Float64)` the mirror aliases the grid's own vectors: no copy, no extra
memory.

Quantities which host scalar code has to produce on every right-hand side — `β(z)`
for a taper or a pressure gradient — go through a
[`Luna.HostMirror`](@ref): a `Float64` host buffer, an optional staging buffer in the
device precision, and the device array. On the host in double precision the device array
*is* the host buffer and `upload!` does nothing. This is the fallback, and the path
`linop_integral=:quadrature` leaves `β` and `Aeff` on; `:tabulated` replaces it with
tables (below).

## The integrated linear operator

The stepper propagates the linear part of a step by `exp(Φ(t2) − Φ(t1))` with
`Φ(z) = ∫linop dz'`, which is exact, so a z-dependent operator reaches it as its integral
and never as the operator itself. [`LinearOps.AbstractIntegratedLinop`](@ref
Luna.LinearOps.AbstractIntegratedLinop) is that interface and
[`RK45.make_prop!`](@ref Luna.RK45.make_prop!) has one method for it. A constant operator
is an array and keeps its own exact method; a bare `linop!(out, z)` callable is refused,
and `Luna.run` converts one according to `linop_integral`.

`exp(linop(t2)·(t2 − t1))`, the one-point rule Luna used for a callable before gpu/27, is
removed. Its error is first order in the step and is common to both of the embedded
Runge–Kutta solutions, so it cancels out of the error estimate and the step controller
never responds to it; measured, it is 7.5e-2 relative on `Eω` at 20 steps over a 0.1 m
0 → 1 bar gradient and 4.1e-1 for a 75 → 50 µm taper. Gradient and taper results therefore
change, which is the one documented exception to the project's rounding-level regression
gate.

### The interface

A subtype implements one of two styles, declared by `LinearOps.PhaseStyle`:

- `AbsolutePhase()` (the default): `phase!(out, op, z)` fills `out` with `Φ(z)`, less
  whatever straight line `secant(op)` reports (`nothing` by default). The propagator keeps
  two buffers and forms the difference itself, reading `Φ(t1)` once per step — the six
  stages share it — and `Φ(t2)` once per distinct `t2`, and adds `L̄·(t2 − t1)` back in the
  same broadcast as the exponential.
- `IncrementalPhase()`: `phasediff!(out, op, z1, z2)` fills `out` with `Φ(z2) − Φ(z1)`
  directly, and the propagator caches it on the `(t1, t2)` pair, which catches the same
  repeat. This is for an operator which computes the integral rather than looking it up:
  integrating from a fixed origin at every stage would cost more and round worse, since
  the accumulated `Φ` of a metre of fibre is hundreds of radians and the difference over a
  step is a fraction of one.

Every subtype also implements `derivative!(out, op, z)`, the operator `L(z)` itself.
Nothing in the propagation needs it; the diagnostics do — a zero-dispersion wavelength, a
linear propagation applied to an input field — and a caller who replaced the closure with
an integrated operator has to be able to get the operator back.

`phase!`, `phasediff!` and `derivative!` run at every stage and, on a device run, `out` is
on the device, so they must be broadcasts or other kernels and every scalar must go through
`Luna.scalar`.

### The two implementations

[`LinearOps.QuadratureLinop`](@ref Luna.LinearOps.QuadratureLinop) is the
`IncrementalPhase` one: `QuadGK.quadgk!` over `[z1, z2]` in `Float64` on the host, with an
absolute tolerance on the max norm, uploaded through a staging buffer in the state's
element type. One 15-point Gauss–Kronrod rule is the usual cost of a step; `maxevals`
bounds a pathological one and the constructor's warning says when it was hit. It holds no
table and no state beyond its buffers, and it is the check on the table: the two agree to
1e-7 on a gradient and a taper at 20, 80 and 320 fixed steps.

[`LinearOps.OffsetLinop`](@ref Luna.LinearOps.OffsetLinop) wraps an integrated operator
with a constant added to `L`. That is a straight line added to `Φ`, so
`Boundaries.addloss` and `addloss_k` fold a spectral or k-space absorption rate into a
caller-supplied integrated operator exactly and, for an `AbsolutePhase` operator, for
free — it is the secant. `Boundaries.clampdecay` is not linear in the operator, cannot be
pushed through the integral, and raises.

## Tabulated z-dependent quantities

`linop_integral=:tabulated` — the default — replaces every host quantity the step would
otherwise evaluate with a table over `z`, built at setup and held on the state's array
type:

- [`LinearOps.TabulatedLinop`](@ref Luna.LinearOps.TabulatedLinop): the integrated operator
  `Φ(z) = ∫ linop dz'`, read back with a cubic Hermite interpolant. This is the
  `AbsolutePhase` implementation of the interface above; `derivative!` differentiates the
  same interpolant, so the operator it reports is consistent with the `Φ` the propagator
  uses and is exact at a node.
- [`LinearOps.TabulatedVector`](@ref Luna.LinearOps.TabulatedVector) for `β(z)` and
  [`LinearOps.TabulatedScalar`](@ref Luna.LinearOps.TabulatedScalar) for `Aeff(z)`, put in
  place by `NonlinearRHS.tabulate`. These are values rather than integrals and no
  derivative of either is available from the mode interface, so they are interpolated
  linearly.

The nodes are placed by bisection: an interval is accepted when the interpolant that will
be read back agrees, at the interval's midpoint, with a directly computed value to the
tolerance. That measures the error of the thing actually used, and it puts nodes at the
features Luna's operators have — the `1/√z` cusp of a gradient filled from vacuum, a
junction in a multi-section fill — without refining anywhere else. The design is PR 440's
`TabulatedUnitaryPhase` generalised to the full complex operator.

`TabulatedLinop` stores not `Φ` but its deviation from the secant through the ends of the
table, and the propagator adds `L̄·(t2 − t1)` back in the same broadcast. In exact
arithmetic that is the same number; in `Float32` it stores the operator's *variation* along
`z` rather than its accumulated phase. `make_linop` already subtracts the frame, so `max|Φ|`
is tens to hundreds of radians (34 rad over 0.1 m of gradient, 359 rad over 1 m) and the
deviation is 11.3 times smaller over the same span. Rounding `ΔΦ` to `Float32` both ways,
the subtraction buys a factor of 4 to 11 in its error — 2.1e-7 rad against 1.3e-6 at 0.1 m,
2.7e-5 against 1.0e-4 at 10 m. It is a reduction of the rounding error, not the difference
between working and not working, and it costs nothing: a cubic Hermite is exact on a linear
function, so the node placement and the interpolation error are unchanged, and a
z-independent operator ends up storing nothing at all.

The `β` and `Aeff` tables are tied to `:tabulated`: `:quadrature` and a caller-supplied
operator leave them evaluated per stage, as they were before tabulation existed. A constant
operator is not tabulated at all, and neither are its `β` and `Aeff` — reading a constant
off a two-node table is `(1-s)f + sf`, not `f`, and a uniform fibre is meant to be
bit-for-bit what it was.

`Luna.run` tabulates into a transform of its own and leaves the caller's object alone, so a
statistics function built from `transform.aeff` before the run would keep calling the
untabulated one. `prop_capillary` therefore tabulates `Aeff` itself, over `[0, flength]`,
before `Stats.default` closes over it; `NonlinearRHS._aefftab` rebuilds a wider table from
that one's `src` for the propagation, which needs `Aeff` up to one step past the end of the
fibre. A low-level caller who builds statistics by hand and wants the same has to pass a
`LinearOps.TabulatedScalar` to `Luna.setup` as `aeff`, which is all `prop_capillary` does.

The statistics of the last accepted step are recorded a fraction of a step past the end of
the fibre, i.e. outside the `[0, flength]` table `prop_capillary` built for them. A value
table read outside its span calls its source callable rather than holding its end value,
so that point is the same number the untabulated path gave; holding it is worth 2.8e-2 on
the peak intensity of the regression gate's taper case. Nothing inside a propagation can
reach that branch, since `Luna.run` rebuilds the table over everything the stepper can
ask about. The *operator* table holds instead, and warns: the propagator adds the secant
term whatever the readback returns, so a step outside it would propagate with the mean
operator over the whole table.

What this does not cover: a `linop!` which is genuinely discontinuous in `z` cannot be
tabulated to tolerance, and the bisection stops at its depth or node limit and warns.
Nor is the table free for a geometry whose operator is the size of the whole state — a
multimode, radial or free-space one — where it is `2·nnodes` copies of it; the constructor
reports its size and warns above `LinearOps.TABLE_WARN_BYTES` (256 MB), and
`linop_integral=:quadrature` is the way out.

## Unit scaling

Kernels never see SI-unit coefficients. The state is `e = E/E_ref` and the polarisation
buffer holds `p = P/(P_ref E_ref)`; see [`Luna.UnitScaling`](@ref).

`E_ref = P_ref = 1` for every `Float64` run, which makes the scaled and unscaled arithmetic
identical and is why the default CPU path did not move. For `Float32`, `E_ref` is the peak
of the time-domain input at `z0` rounded to a power of two — so that dividing by it is
exact in every linear operation — and `P_ref` is `ε₀`.

Each response's constants are combined with the right powers of `E_ref` and `P_ref` on the
host, in `Float64`, by [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale), and only the
result is converted to the device precision. The frequency-domain normalisation vectors are
precombined the same way, one vector per z-independent factor.

Dynamic-range audit, helium (the worst case in Luna's parameter range, since `γ₃(He)` is
the smallest):

| quantity | 1 bar | 0.3 bar |
| --- | ---: | ---: |
| `ρ ε₀ γ₃` (unscaled) | 1.1e-38 | 3.3e-39 |
| smallest normal `Float32` | 1.2e-38 | 1.2e-38 |
| scaled `γ₃` (`E_ref = 8192`, `P_ref = ε₀`) | — | 3.9e-34 |
| scaled `ρ ε₀ γ₃` as the kernel sees it | — | 2.5e-20 |

Unscaled, the 0.3 bar coefficient is subnormal and Metal flushes it to zero; scaled, it is
eighteen orders of magnitude above the subnormal threshold.

The unscaling happens in one place only, the output boundary.

## The response protocol

A nonlinear response is still a callable `resp!(out, E, ρ)` which accumulates into `out`.
What is new is a trait, [`Nonlinear.kind`](@ref Luna.Nonlinear.kind), which says *how*
[`NonlinearRHS.Et_to_Pt!`](@ref Luna.NonlinearRHS.Et_to_Pt!) evaluates it. `Et_to_Pt!`
walks the response list in order, takes the longest run of pointwise responses at a time
and evaluates it as one broadcast, and applies everything else one response at a time.
The first group written assigns into `Pt`; every later group folds `Pt` in as the leading
term of its own sum, so the sequence of additions each element sees is exactly the one the
per-response loop produced. Where a transform passes `idcs`, it must cover every column of
the block: the pointwise and batched paths act on the whole block and ignore it, and it is
the columnwise loop's range.

| kind | what the dispatcher does | what the response supplies | examples |
| --- | --- | --- | --- |
| [`Nonlinear.Pointwise`](@ref) | fuses it into one broadcast over the whole block with the other pointwise responses, no buffer | [`pointwise_kernel`](@ref Luna.Nonlinear.pointwise_kernel) (a `T -> T`), or [`pointwise_expr`](@ref Luna.Nonlinear.pointwise_expr) if it carries per-sample arrays | `KerrField`, `KerrEnv` (scalar field), `KerrEnvTHG` |
| [`Nonlinear.VectorPointwise`](@ref) | two broadcasts, one per polarisation component, each fused across the group | [`vector_kernel`](@ref Luna.Nonlinear.vector_kernel) (an `(ex, ey) -> SVector{2}`), or [`vector_expr`](@ref Luna.Nonlinear.vector_expr) | `KerrField`, `KerrEnv` (two-component field), `Chi2Field`, `Chi2Env` |
| [`Nonlinear.Batched`](@ref) | calls it once with the whole `(nt, npol, ncols)` block, through [`batched!`](@ref Luna.Nonlinear.batched!), which also hands it the unit scaling | [`batched!`](@ref Luna.Nonlinear.batched!) (or just the call operator, if its coefficients already carry the scaling), plus its own full-size buffers in the run's array type | `HostResponse`, `PlasmaCumtrapz`, `RamanPolarField`/`RamanPolarEnv`, `KerrFieldNoTHG` |
| [`Nonlinear.Columnwise`](@ref) | calls it once per column, on the host; refused on a device array | nothing — the default, so any callable works | anything user-written |

Three functions carry the units and the precision:

- [`Nonlinear.coefficients`](@ref Luna.Nonlinear.coefficients)`(r, ρ, scaling)` is the one
  place a response's physical constants, the density and the powers of `E_ref`/`P_ref`
  meet. It runs on the host, in `Float64`. [`Luna.polscale`](@ref)`(scaling, n)` gives the
  factor for a response of polynomial degree `n` in the field (`Eref^(n-1)/Pref`): cubic
  for Kerr, quadratic for χ⁽²⁾. A response with no single degree scales its intermediates
  itself.
- the kernel converts that scalar exactly once, with [`Luna.scalar`](@ref)`(E, c)`, so
  nothing reachable from a device kernel is a `Float64`.
- [`Nonlinear.rescale`](@ref Luna.Nonlinear.rescale) converts the *arrays* a response
  carries (a carrier phase, a tabulated coefficient) to the run's array type and
  precision, allocates whatever depends on the shape of the block, and is where a
  columnwise response is wrapped for a device run. A transform calls the **four**-argument
  form, `rescale(r, spec, scaling, Et)`, where `Et` is a prototype of the field block:
  `size(Et)` is `(nt, npol, ncols...)`, `eltype(Et)` says whether the field is real or
  complex and in what precision, and `similar(Et)` allocates a buffer of the run's array
  type. A [`Batched`](@ref Luna.Nonlinear.Batched) response implements that form, so that
  its buffers exist at construction and the transform's residency assertion can see them
  (rule 4 above). Everything else implements the three-argument form, which the
  four-argument fallback delegates to.

  **Which kinds may skip it.** A [`Pointwise`](@ref Luna.Nonlinear.Pointwise) or
  [`VectorPointwise`](@ref Luna.Nonlinear.VectorPointwise) response whose coefficients are
  all scalars needs no method at all: the fallback passes it through when
  [`resident_arrays`](@ref Luna.Nonlinear.resident_arrays) is empty, because
  `coefficients` is called per right-hand side with the run's scaling. A
  [`Batched`](@ref Luna.Nonlinear.Batched) response may **not**: it is handed the block,
  the density and the units and nothing else, so the fallback refuses it unless the run is
  unscaled and on the host. A response which carries an array and has no method is an
  error either way.

  The fallbacks also check, once per setup, that every non-`isbits` array field of a
  device-kind response is one `resident_arrays` names, and say which field is missing if
  not. A `StaticArrays` or `Rotations` matrix is `isbits`, travels inside the struct and
  is exempt.

**The response struct itself never enters a kernel.** Only the scalars the kernel captured
and the arrays `rescale` converted do. That is why `KerrField` keeps its physical
`Float64` `γ3` on a Metal run and the Metal tests check the *kernel's* element type rather
than the struct's.

It is also why a response which holds an FFT plan — [`Nonlinear.RamanPolarField`](@ref
Luna.Nonlinear.RamanPolarField), [`Nonlinear.KerrFieldNoTHG`](@ref
Luna.Nonlinear.KerrFieldNoTHG) — is moved by [`rescale`](@ref Luna.Nonlinear.rescale) and
**not** by an `Adapt.adapt_structure` rule. `rescale` constructs the response for the
target array type, which includes planning its transforms there; `Adapt` cannot replan,
so an `Adapt` rule over such a struct would produce plans that did not match their
buffers. A response whose only device-side state is an array (`KerrEnvTHG`'s carrier
phase) has an `Adapt` rule as well, because there the two agree.

A new scalar pointwise response is therefore three short methods (this is
`SquareResponse` in `test/test_device.jl`):

```julia
struct SquareResponse{T}
    c::T
end
Nonlinear.kind(::SquareResponse) = Nonlinear.Pointwise()
Nonlinear.coefficients(r::SquareResponse, ρ, scaling) = ρ*r.c*Luna.polscale(scaling, 2)
Nonlinear.pointwise_kernel(r::SquareResponse, E, ρ, scaling) =
    (fac = Luna.scalar(E, Nonlinear.coefficients(r, ρ, scaling)); e -> fac*e^2)
```

The kernel must capture nothing but `isbits` scalars: it is the body of a broadcast which
may be compiled for a GPU. A response which needs a per-sample array passes it as a second
broadcast argument by overriding `pointwise_expr` instead, as `KerrEnvTHG` does with its
carrier phase, and lists it in `resident_arrays` so the transform's residency assertion
covers it.

**Deviation from GPU_PLAN.md §4.3.** The plan describes the vector form as one broadcast
returning an `SVector{2}` written through a `reinterpret`ed `(2, nt, ncols)` view of the
output. Luna's buffers are `(nt, npol, ncols)` — the polarisation index is the *slow* axis
— so there is no contiguous leading axis of length 2 to reinterpret, and producing one
would mean a transpose and a buffer. The kernel contract is kept (the response returns an
`SVector{2}` of the two lab-frame components); the dispatcher materialises it into the two
output component views with one broadcast each, both fused across the whole pointwise
group.

For the two-component Kerr forms this is the same arithmetic, and the same number of
passes, as the pair of broadcasts they were already written as. It is not free in general:
because each output component is materialised by its own broadcast over the *same*
`vector_expr`, a genuinely coupled response evaluates its shared intermediates twice. For
the χ⁽²⁾ responses that is the lab-to-crystal rotation and the contracted field products,
which are evaluated once per output component. Dead-code elimination removes the unused
component of the returned `SVector`, and with it the row of the 3×6 contraction which
produced it, but not the work the two components share. That is the price of the layout;
a response for which it matters can override `vector_expr` to compute the shared part in
a form the compiler can hoist, or ask for a batched kind instead.

### The χ⁽²⁾ responses

[`Nonlinear.Chi2Field`](@ref) and [`Nonlinear.Chi2Env`](@ref) are the vector-pointwise
form in its intended shape, and the reason the kind exists. At one time sample the
response rotates the two lab-frame components into the crystal frame, forms the six
contracted second-order products, contracts them with the 3×6 tensor and rotates the
result back — a chain of small dense products which couples the components and nothing
else.

Everything it contracts with is an `isbits` static matrix: `SMatrix{3, 6, T}` for the
tensor, `Rotations.RotMatrix3{T}` for the two rotations, never an `MArray`. Static
matrices travel inside the closure a broadcast compiles, so the responses need no buffer
and carry no host array apart from `Chi2Env`'s carrier phase, which is a second broadcast
argument (as [`KerrEnvTHG`](@ref Luna.Nonlinear.KerrEnvTHG)'s is) and is listed in
`resident_arrays`. The four work vectors the old implementation allocated at construction
and wrote into per time sample are gone.

`coefficients` is `ε₀·polscale(scaling, 2)`: the responses are quadratic in the field, so
they carry one power of `E_ref`, against the Kerr responses' two. The density is ignored,
as it always was — a χ⁽²⁾ crystal is not a gas.

Both declare `kind` as `VectorPointwise()` *unconditionally* rather than only for
`Val(2)`. They have no one-component form, and a `Val{2}`-only declaration would leave
`kind(r)` at `Columnwise()` and with it `device_capable(r)` false. The scalar case is
therefore refused by the `Val(1)` check in `NonlinearRHS._scalarexpr` rather than
silently taken down the elementwise path.

`rescale` converts the crystal matrices to the run's real type and moves the carrier
phase. The kernel converts the matrices again, from whatever type the response holds to
`real(eltype(E))`. That is 27 host scalars per component broadcast, so 54 per right-hand
side — the two-broadcast layout doubles this as it doubles the shared arithmetic — and
the identity on a response `rescale` has already converted. It is there so that a kernel
built from a response which never went through `rescale` still carries no `Float64`
(rule 3 above), which is what a low-level caller assembling `Et_to_Pt!` by hand does.

**The guarantee is partial, and covers the matrices only.** An unrescaled `Chi2Field` is
entirely `isbits`, so it compiles and runs on a device. An unrescaled `Chi2Env` does not:
its carrier phase is still a host `Vector{ComplexF64}`, it enters the broadcast as a
host array, and Metal refuses to compile the kernel. That failure is loud and immediate,
not a silently wrong answer, but `Chi2Env` does need `rescale` — which is what
`resident_arrays` and the transform's residency assertion are for.

The χ⁽²⁾ transforms themselves — the Cartesian free-space ones — are host-only until
`gpu/21`, so the responses are exercised on a device block directly (`test_device.jl`,
`test_metal.jl`) rather than through a propagation.

## The host fallback

[`Nonlinear.HostResponse`](@ref Luna.Nonlinear.HostResponse) is what makes rule 6 of the
kernel discipline concrete for responses, and what GPU_PLAN.md §3 means by keeping Luna
hackable: a response written as a plain closure, in physical SI units, must not stop a
user from running on a GPU.

`rescale` wraps any [`Columnwise`](@ref Luna.Nonlinear.Columnwise) response in one as soon
as the run is not host `Float64` in physical units, and logs one `@info` line per wrapped
response at setup. The wrapper is `Batched`, so it receives the whole block, and at every
right-hand side it

1. copies the block to a host buffer in the run's element type, and multiplies by `E_ref`
   into a `Float64`/`ComplexF64` buffer — physical units, double precision, which is what
   the response was written for;
2. zeroes its own polarisation buffer and calls the response on it column by column,
   exactly as `Et_to_Pt!`'s `idcs` loop does;
3. multiplies by `1/(P_ref E_ref)` into the staging buffer and copies it back up, then
   adds it to the output.

Four buffers and two copies per right-hand side: `Eh` and `Ph` on the host in physical
`Float64`/`ComplexF64`, `stage` on the host in the run's element type, and `Pd` in the
run's array type. A scaled host run needs only the first two and leaves the others
`nothing`. The staging buffer exists because `copyto!` between a host array and a device
array does not convert the precision. All of them are allocated by the constructor, from
the block prototype `rescale` is given, so the per-call code is fully typed and `Pd` is
covered by the transform's residency assertion.

The simple interface deliberately does not use it: `prop_capillary` refuses an explicit
`device`/`precision` request whose responses are not all
[`device_capable`](@ref Luna.Nonlinear.device_capable), naming `device=:cpu`, rather than
silently producing a run slower than the CPU. `device_capable` means "has a kernel of its
own", i.e. `kind` is not `Columnwise`; it is about speed, not possibility.

## The plasma response

[`Nonlinear.PlasmaCumtrapz`](@ref Luna.Nonlinear.PlasmaCumtrapz) is the first response
with a device kernel which is not a single fused broadcast, and it is the one the rest of
this section's machinery exists for. It is [`Batched`](@ref Luna.Nonlinear.Batched): it
is handed the whole `(nt, npol, ncols...)` block and evaluates it as

1. one broadcast for the ionisation rate (below);
2. `fraction`, `J` and `P`, each a cumulative trapezoid integral along the time axis;
3. one `ifelse` broadcast for the ionisation-loss term;
4. one broadcast to accumulate into the output.

**The integrals are the reason for the trait.** A trapezoid integral written as a loop —
`out[i] = out[i-1] + δt(y[i-1] + y[i])/2`, which is what `Maths.cumtrapz!` does — has a
serial dependency along the time axis and cannot be a broadcast.
`Maths.cumtrapz_scan!` is the same integral as one
prefix scan and one broadcast,

```math
\mathrm{out}_i = δt\left(\sum_{k\le i} y_k - \frac{y_1 + y_i}{2}\right)
```

which is `accumulate!(+, out, y; dims=1)` followed by `@. out = δt*(out - (y + y1)/2)`.
Metal and CUDA both provide a native scan; `JLArrays` does not, so `test_device.jl`
supplies a host-backed one. The summation order is not the loop's, so the answer differs
at rounding level — measured end to end, over four propagations which ionise between 0.1
and 3 per cent of the gas, at 0.9e-15 to 4.3e-15 relative in `Eω`.

**The loss term** is `ifelse(abs(E) > 0, closs*rate*(1-fraction)/E, zero(E))`. `ifelse` is
a select, not a branch: both arms are evaluated, so the division by a zero field does
happen and produces an infinity, and the select discards it. That is the same value the
loop's `if abs(E[ii]) > 0` produced, and it is why the arm has to be an `ifelse` rather
than a short-circuit — a branch per element is what a device kernel cannot afford.

For a two-component field the term divides by `Em²`, and **the guard is on `Em²`, not on
`Em`**. In the wings of a pulse `Em` is many orders below its peak — `exp(-50)` is 1e-22
— so `Em²` is 1e-44, which is subnormal in `Float32` and which a device flushes to zero
while `Em` itself is still a normal number. The rate there is zero too, so the term
becomes `0/0`, and one NaN travels through the three scans which follow and destroys the
whole column. On the CPU the subnormal survives and `0/x` is `0`, so nothing but the
Metal hardware test sees this; it is the clearest example of why rule 7 of the kernel
discipline exists. In `Float64` the two guards differ only below 1e-162 V/m.

**Buffers.** The response owns `rate` and `fraction` (one polarisation component, so a
singleton second axis which broadcasts against both), `Em` for a two-component field, and
`J` and `P` at the size of the block. `P` holds the phase-modulation term until the last
scan and the polarisation after it; the two are never live at once, which is one
block-sized buffer fewer than the columnwise response held per column.

They are sized for the block, which for a radial or free-space transform is the whole
transverse grid and not the single column the constructor is given a prototype of, so
every transform calls
[`Nonlinear.rescale_responses`](@ref Luna.Nonlinear.rescale_responses) — the
four-argument [`rescale`](@ref Luna.Nonlinear.rescale) mapped over its response
collection, through the tuple-of-tuples a gas mixture is — at construction. On the
default CPU path every fallback returns the response unchanged, so this changes nothing
for anything else.

**Units.** The response has no single polynomial degree, so
[`Luna.polscale`](@ref) does not apply and
[`coefficients`](@ref Luna.Nonlinear.coefficients) returns four scalars instead of one:
`Eref` (the rate is evaluated at the physical field `Eref*e`), `e_ratio*Eref`,
`ionpot/Eref` and `ρ/(Pref*Eref)`. Each reduces to the physical constant at
`Eref == Pref == 1`.

**CPU threading.** A column of the block costs a rate evaluation and three prefix scans,
which is a large grain, so above `Nonlinear.PLASMA_THREAD_MINLEN` elements and more than
one column the host path shares the columns out with `Threads.@threads :dynamic`
(GPU_PLAN.md §4.9; `:dynamic` rather than `:static`, which throws when nested). Each task
takes one column of every buffer and runs the same `_plasma_block!` as the whole-block
call, so no column's arithmetic depends on the threading or on how many columns are
passed at once — `test_device.jl` asserts that as exact equality, not a tolerance. A
mode-averaged or modal transform has one column and is never threaded.

Measured on an M1 Pro, `batched!` on a 2048-sample block of columns in `Float64`:

| columns | 1 thread | 8 threads | Metal (`Float32`) |
| ---: | ---: | ---: | ---: |
| 1 | 36 µs | 36 µs | 587 µs |
| 16 | 538 µs | 95 µs | 653 µs |
| 128 | 4.51 ms | 580 µs | 620 µs |

The threading is close to linear from a handful of columns up. A single column — which is
every mode-averaged and modal run — is launch-bound on the GPU and is left alone by the
threading, as it should be.

### Dynamic range in `Float32`

The magnitude of every intermediate of the plasma response, in `Float32`, over 45
parameter combinations: He, Ne, Ar, Kr and Xe at 0.1, 1 and 10 bar, driven at 1e12, 1e14
and 1e16 W/cm² (800 nm, 10 fs, 1024 samples), with `E_ref` the peak field rounded to a
power of two and `P_ref = ε₀`. "Headroom" is the distance to the `Float32` limits
(`1.2e-38` normal, `3.4e38`); the last column is the smallest normal `Float32` divided by
the largest magnitude of the same array, i.e. the size of the error a flushed subnormal
could cause relative to that quantity's own peak.

| quantity | smallest non-zero | largest | headroom (low / high) | worst subnormal / own peak |
| --- | ---: | ---: | ---: | ---: |
| `e = E/E_ref` | 5.5e-23 | 1.3e+00 | 5e+15 / 3e+38 | 1.5e-38 |
| `E_ref` | 2.1e+09 | 2.8e+11 | 2e+47 / 1e+27 | 5.5e-48 |
| rate `W` [1/s] | 5.7e-27 | 1.9e+17 | 5e+11 / 2e+21 | 1.0e-43 |
| `cumsum(W)` | 5.7e-27 | 1.7e+19 | 5e+11 / 2e+19 | 2.4e-44 |
| `fraction` (integral) | 3.4e-43 | 2.0e+03 | 3e-05 / 2e+35 | 2.0e-28 |
| `fraction` (`1-exp`) | 6.0e-08 | 1.0e+00 | 5e+30 / 3e+38 | 1.8e-35 |
| `e_ratio*E_ref` | 6.1e+01 | 7.8e+03 | 5e+39 / 4e+34 | 1.9e-40 |
| phase term | 3.5e-23 | 7.7e+03 | 3e+15 / 4e+34 | 2.7e-38 |
| `cumsum(phase)` | 3.6e-06 | 2.9e+04 | 3e+32 / 1e+34 | 7.3e-39 |
| `J` (phase integral) | 4.9e-22 | 3.4e-12 | 4e+16 / 1e+50 | 6.4e-23 |
| loss term | 5.6e-45 | 9.1e-14 | 5e-07 / 4e+51 | 7.2e-16 |
| `J` (with loss) | 5.6e-45 | 3.4e-12 | 5e-07 / 1e+50 | 7.2e-16 |
| `cumsum(J)` | 5.6e-45 | 9.7e-11 | 5e-07 / 4e+48 | 5.9e-16 |
| `P` | 1.4e-45 | 1.1e-26 | 1e-07 / 3e+64 | 5.0e+00 |
| `ρ/(P_ref E_ref)` | 1.0e+24 | 1.4e+28 | 9e+61 / 2e+10 | 1.2e-62 |
| `out` | 1.4e-21 | 1.2e+00 | 1e+17 / 3e+38 | 6.2e-25 |

Nothing overflows and nothing is non-finite anywhere in the range: the largest
intermediate is `cumsum(W)` at 1.7e19, nineteen orders below `floatmax(Float32)`. The
scaling is what buys that — unscaled, the output coefficient `ρ/(P_ref E_ref)` would be
`ρ` itself and `P` would carry the whole of ε₀.

One intermediate is missing from the table because it only exists for a two-component
field: `Em²`, the square of the field magnitude the loss term divides by. It reaches
1e-44 in the wings of a pulse, is subnormal in `Float32`, and is the one place in the
response where flushing a subnormal to zero produces a NaN rather than a small error —
see the loss term above.

The bottom end is where it is worth looking. Four quantities reach the subnormal range,
and for three of them (`fraction`'s integral, the loss term, `J`) the subnormal values are
at least 1e15 times smaller than the peak of the array they are in, so flushing them to
zero — which Metal does — is below `Float32`'s own precision and cannot matter. `P` is the
exception: there are parameter combinations where the *whole* `P` array is subnormal, and
Metal produces zero for the plasma polarisation rather than a very small number.

| gas | pressure | intensity | peak `P` | peak plasma `out` / peak Kerr `out` |
| --- | ---: | ---: | ---: | ---: |
| He | 0.1, 1, 10 bar | 1e14 W/cm² | 2.3e-39 | 2.6e-07 |

That is the one case in the range where single precision loses the plasma term
altogether, and where it does, the term is 2.6e-7 of the Kerr term at the same
parameters — the size of a single `Float32` rounding of the Kerr term itself. The ratio
does not depend on pressure, since both terms are proportional to ρ. At 1e12 W/cm² helium
does not ionise at all (the rate is zero in both precisions) and at 1e16 W/cm² the plasma
term dominates and is nowhere near the subnormal range. A propagation which depends on a
plasma contribution that small needs `Float64`.

### Ionisation rates in a kernel

The rate is the part of the response which is not arithmetic on the field: it either
evaluates a formula from a dozen stored constants ([`Ionisation.IonRateADK`](@ref)) or
indexes a spline of `log(rate)` ([`Ionisation.IonRatePPTAccel`](@ref)). Both are captured
by the broadcast kernel rather than passed as broadcast arguments, because the spline is
indexed at a data-dependent position; `Adapt` rewrites the closure when the kernel is
compiled, which is what turns the rate's arrays into device pointers.

| | |
| --- | --- |
| [`Ionisation.device_capable`](@ref)`(ir)` | whether `ir` has a kernel: ADK, and a cached PPT rate on a uniformly spaced table |
| [`Ionisation.device_rate`](@ref)`(ir, spec)` | the same rate in `spec`'s precision and array type; `ir` itself on the default host path |
| [`Ionisation.ionrate!`](@ref)`(out, ir, E, Eref)` | the array-level evaluation, one broadcast on any array type |
| [`Ionisation.ratekernel`](@ref)`(ir, Eref)` | the kernel itself, `e -> W(Eref*e)` |

Three things had to change for this to compile for a GPU.

- **Every constant is parametric.** `IonRateADK` held nine `Float64` fields and
  `IonRatePPTAccel`'s `Emin`/`Emax` were `Float64`; a `Float64` struct field read inside
  a Metal kernel never compiles. So is the spline: `Maths.CSpline`'s arrays move with
  `Maths.todevice_spline`, and the index function it
  captured for a uniform axis, which was a closure over three `Float64`s, is now
  `Maths.UniformIndex`.
- **No error paths.** `Maths.spline_eval` is `CSpline`'s
  evaluation without the optional bounds check, whose `DomainError` message is built by
  string interpolation. Above the table `IonRatePPTAccel`'s kernel saturates at the
  table's last value; the host keeps today's error, raised once per call from a
  `maximum(abs, E)` check rather than once per element.
- **The field is reconstructed, not scaled away.** A response of polynomial degree `n`
  carries `Eref^(n-1)` in a coefficient; an ionisation rate has no degree, so the kernel
  computes `W(Eref*e)`. `Eref` is 1 for every `Float64` run, where `1*e` is exact.

What single precision costs the rate itself, over the field range where the rate is big
enough to matter (`W > 1e6` 1/s), measured against the same rate in `Float64`:

| gas | ADK | cached PPT table | field range [V/m] | stored `log(rate)` range |
| --- | ---: | ---: | ---: | ---: |
| He | 6.0e-06 | 6.4e-06 | 2.9e10 – 2.1e11 | −697 … +37 |
| Ne | 3.9e-06 | 6.6e-06 | 2.4e10 – 1.6e11 | −698 … +37 |
| Ar | 3.2e-06 | 6.3e-06 | 1.4e10 – 8.6e10 | −695 … +37 |
| Kr | 2.9e-06 | 6.1e-06 | 1.2e10 – 6.8e10 | −694 … +37 |
| Xe | 4.2e-06 | 6.8e-06 | 9.4e09 – 5.1e10 | −689 … +37 |

Both are a few times `Float32`'s own precision, from two places: the ADK exponent, whose
argument is order 10–100 so that a 1e-7 relative error in it becomes a 1e-5 one in the
rate, and the table's knots, which are ~1e11 apart by ~3e6 and so are resolved to ~5e-3
of an interval in `Float32` (the interpolation position, not the value).

Below that range both rates underflow to zero in `Float32` — the ADK rate for helium
first becomes non-zero at 8.0e9 V/m instead of the `Float64` threshold of 1.1e9 V/m,
where its value is 4e-304 1/s. Over a 100 fs pulse that is an ionisation fraction of
1e-290, so the two are the same answer.

**The stored table needs no offset.** GPU_PLAN.md §4.1 suggests holding the PPT table "in
`log` form with an offset so that `Float32` covers the range". It is already in `log`
form (`IonRatePPTAccel` splines `log(rate)`), and the range that puts in the table is
−708 to +37, which `Float32` holds with room to spare: the resolution is 6.1e-5 at −708
and 3.8e-6 at +37, i.e. a 6.1e-5 and 3.8e-6 relative error in the rate after `exp`. The
first of those is at a rate of 1e-308 1/s and is not a number anybody uses. An offset
would move the whole table towards zero and improve the unused end; it is not implemented,
because it would also change the `Float64` path, which is bit-identical to the one before
this branch.

A rate with no kernel — the direct [`Ionisation.IonRatePPT`](@ref), whose series
summation and `BigFloat` fallback are host code, a table which ended up on a
`Maths.FastFinder`, or a user's callable — is refused by
`device_rate` with a message naming the alternatives and `device=:cpu`.

## The Raman polarisation

[`Nonlinear.RamanPolarField`](@ref Luna.Nonlinear.RamanPolarField) and [`Nonlinear.RamanPolarEnv`](@ref Luna.Nonlinear.RamanPolarEnv) are
[`Batched`](@ref Luna.Nonlinear.Batched) for the same reason the plasma response is: the
convolution of the driving term with the Raman response function is a transform of the
whole column, not an operation on one sample. Per right-hand side each is

1. one broadcast for the driving term (`E²`, or `½|A|²` for an envelope and for
   `thg=false`, where `A` is the analytic signal);
2. one forward FFT along the time axis, over a **doubled** time grid — the driving term
   occupies the first half and the second is zero padding, which makes the
   multiplication in the frequency domain the full linear convolution rather than a
   circular one;
3. one broadcast for the product with the frequency-domain response function;
4. one inverse FFT;
5. one broadcast to multiply by the density and the field and accumulate into the output.

Both transforms are **batched over the block's columns**: one pair of FFTs per
right-hand side whatever the geometry, rather than the three per column the response used
to do. The plans are made by [`Utils.plan_ft`](@ref Luna.Utils.plan_ft) on the buffer
itself, so they are FFTW plans on the host and the backend's own on a device, and the
inverse plan is held unnormalised with its `1/N` folded into a scalar (below).

**The response function is host scalar code.** `r(h, ρ)` sums a few dozen damped
oscillators into a `Float64` host vector; that is not a kernel and does not need to be.
What changed is how often it runs: it used to be evaluated, and transformed, at *every*
right-hand side, and it is now keyed on the density, so a run at constant pressure
evaluates it once and a pressure gradient pays what it always did. Only the transformed
result crosses to the device, through a host staging buffer in the run's precision
(`copyto!` between a host and a device array does not convert).

### Splitting the coefficient (`_splitscale`)

The Raman constants are the smallest numbers in Luna. `K` in
[`Raman.RamanRespVibrational`](@ref Luna.Raman.RamanRespVibrational) is
`(4πε₀)²(dα/dQ)²/(4μΩ)`, around 1e-48 in SI units, and the frequency-domain response
function is around 1e-45. The scalar it is multiplied by — the time step, the unit
scaling and the power of two below — is around 1e-15. Their product is a perfectly
ordinary number, but **in `Float32` neither factor exists on its own**: 1e-45 is below
the smallest subnormal and a device flushes it to zero.

`_splitscale` divides the frequency-domain response function by a power of two chosen so
that the two factors land on either side of the square root of their product, i.e. it
splits the smallness evenly between them and gives each the widest margin against
underflow it can have. For the gases in the table below the exponent is between −54 and
−103.

Dividing by a power of two and multiplying by it again is exact, and an FFT of a
power-of-two-scaled input is the scaled FFT of that input, so **this changes no `Float64`
value**: the default CPU path is bit-for-bit what it was.

The `1/N` of the inverse transform goes on the *density*, at the end, rather than into
the frequency-domain scalar. Both are exact, but `N` is 2^18 or so for a typical grid and
applying it before the transform would cost the frequency-domain buffer five orders of
`Float32` headroom for nothing.

### Dynamic range in `Float32`

Every gas whose Raman response Luna can build — which is every gas `Interface.jl` turns
Raman on for by default, plus fused silica — at 0.1, 1 and 10 bar, with a 20 fs pulse at
800 nm of peak field 1e10 V/m (a few tens of µJ in a 75 µm capillary), `E_ref = 2^33` and
`P_ref = ε₀`. `UF` marks a quantity below the smallest normal `Float32` (1.2e-38).

`max |h(ω)|` is an unnormalised DFT sum over the doubled grid, so it depends on the number
of samples and on `δt`: the field rows are on the `grid.to` of
`Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)` (4096 samples, `δt` = 1.668e-16 s) and
the envelope rows on that of the matching `Grid.EnvGrid` (1024 samples, `δt` =
6.901e-16 s). The SiO₂ row uses the `scale` argument of
[`Raman.raman_response`](@ref Luna.Raman.raman_response), `0.18 ε₀ χ₃(SiO₂)` = 3.202e-34,
which is the size `prop_gnlse` supplies from `fr` and `n₂`.

| gas | kind | max \|h(t)\| | max \|h(ω)\| | split 2^m | max \|h(ω)\|/2^m | scalar | max \|product\| |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| H₂ | field | 1.3e-48 | 1.4e-45 `UF` | 2^-100 | 1.7e-15 | 1.1e-15 | 4.7e-29 |
| H₂ | envelope | 1.3e-48 | 3.5e-46 `UF` | 2^-102 | 1.8e-15 | 1.1e-15 | 1.1e-29 |
| N₂ | field | 6.9e-49 | 3.5e-46 `UF` | 2^-101 | 9.0e-16 | 5.5e-16 | 4.3e-29 |
| N₂ | envelope | 6.9e-49 | 8.5e-47 `UF` | 2^-103 | 8.7e-16 | 5.7e-16 | 1.0e-29 |
| D₂ | field | 1.0e-48 | 1.1e-45 `UF` | 2^-100 | 1.4e-15 | 1.1e-15 | 2.5e-29 |
| D₂ | envelope | 1.0e-48 | 2.6e-46 `UF` | 2^-102 | 1.3e-15 | 1.1e-15 | 6.8e-30 |
| CH₄ | field | 1.5e-48 | 2.4e-45 `UF` | 2^-99 | 1.5e-15 | 2.2e-15 | 1.9e-30 |
| CH₄ | envelope | 1.5e-48 | 5.9e-46 `UF` | 2^-101 | 1.5e-15 | 2.3e-15 | 4.6e-31 |
| SF₆ | field | 6.1e-49 | 9.9e-46 `UF` | 2^-100 | 1.3e-15 | 1.1e-15 | 5.5e-29 |
| SF₆ | envelope | 6.1e-49 | 2.5e-46 `UF` | 2^-102 | 1.3e-15 | 1.1e-15 | 1.4e-29 |
| N₂O | field | 4.0e-48 | 5.7e-45 `UF` | 2^-99 | 3.6e-15 | 2.2e-15 | 6.9e-28 |
| N₂O | envelope | 4.0e-48 | 1.4e-45 `UF` | 2^-101 | 3.5e-15 | 2.3e-15 | 1.7e-28 |
| SiO₂ | field | 1.4e-20 | 2.7e-18 | 2^-54 | 4.8e-02 | 7.7e-02 | 2.4e-01 |

O₂ is missing because its Raman parameters are incomplete (below). The pressure changes
the response function only through the dephasing time, so the rows at 0.1 and 10 bar are
within a few per cent of the ones shown and are left out; the density enters at the end,
on the output scalar. The transformed driving term is 87 (a field) or 21 (an envelope) in
these units, and the largest output is 1.0e-5.

Two things to read off. The unsplit response function is `UF` for **every** gas — without
the split a `Float32` Raman run gives exactly zero. And the worst case after the split,
CH₄ as an envelope, still has seven orders of headroom.

The last column scales as `E_ref²`, because the response is cubic in the field and the
driving term is `O(1)` by construction. At peak fields below about 1e8 V/m it approaches
the subnormal threshold — at which point the Raman polarisation, and every other
nonlinearity, is negligible anyway. `Luna.unitscaling` takes `E_ref` from the peak of the
input field, so the ratio in the table is what a run actually sees.

O₂ is missing because its Raman parameters are incomplete in
`PhysData.raman_parameters`: both the rotational and the vibrational linewidth are
`TODO`, so `Raman.raman_response(t, :O2)` raises a `FieldError` on any branch.

## The no-THG Kerr response

`Kerr_field_nothg(γ3, n)` builds a
[`Nonlinear.KerrFieldNoTHG`](@ref Luna.Nonlinear.KerrFieldNoTHG) rather than the closure it used to. Removing the
third-harmonic term means replacing `E³` with `|A|²E`, where `A` is the analytic signal
of the whole column, so this too is [`Batched`](@ref Luna.Nonlinear.Batched) rather than
pointwise. Its coefficient is the ordinary cubic one.

[`Nonlinear.AnalyticSignal`](@ref Luna.Nonlinear.AnalyticSignal) is the whole-block form of `Maths.plan_hilbert`: one
complex FFT along the time axis, one broadcast against a filter vector, one inverse FFT.
The host version keeps the mean, doubles the positive frequencies and zeroes the negative
ones with three slice assignments; the filter vector is the same three factors, which is
what makes it a kernel. The `1/N` of the inverse transform is folded into that vector, so
there is no separate normalisation pass. Folding is exact — Luna's time grids are powers
of two — and multiplying by an exact zero instead of assigning one differs only in the
sign of a zero, so the analytic signal is **bit-identical** to `Maths.plan_hilbert`'s and
so is the `Float64` propagation.

`RamanPolarField(t, r; thg=false)` uses the same transform for its driving term.

## The radial transform

[`NonlinearRHS.TransRadial`](@ref Luna.NonlinearRHS.TransRadial) is the first transform
with more than one transverse column, and the first one whose per-step work includes a
matrix multiply. Per right-hand side it is

1. one inverse FFT over the time axis, batched over the `(npol, nr)` columns;
2. one matrix multiply, k-space to real space;
3. the response protocol on the whole `(nto, npol, nr)` block;
4. one broadcast for the temporal apodisation;
5. one matrix multiply, real space back to k-space;
6. one forward FFT;
7. one broadcast for the frequency-domain normalisation.

**The Hankel step is one GEMM per direction**, on the block reshaped to `(nto·npol, nr)`
([`Grid.radial_matmul!`](@ref Luna.Grid.radial_matmul!)), where it used to be one `mul!`
per polarisation component on a `view`. On a device that is the difference between the
accelerated matrix multiply and a fallback: MPSGraph's matmul needs plain zero-offset
operands of equal element type, and a `view` is neither (it becomes an `MtlMatrixOperand`
and reaches only the native kernels, which for complex operands are scalar). The transform
therefore holds its own copies of [`Grid.RadialGrid`](@ref Luna.Grid.RadialGrid)'s `Tfwd`
and `Tbwd` **in the time-domain element type** -- `Float32`/`Float64` on a `RealGrid`,
complex on an `EnvGrid` -- rather than the grid's `Float64` ones.

For one polarisation component the reshaped operand is the same matrix, with the same
leading dimension, as the view was, so the CPU result is bit-identical; for two the
summation order may differ, and it is measured at 1e-14 relative in `test_device.jl`. The
regression gate's two radial cases are unchanged to the last bit
(`gpu/20-radial-device`).

`Grid.radial_matmul!` allows `out === A` and makes one copy when it is, which is how the
noise setup and the transverse collar use it; nothing per step does.

**The frequency-domain normalisation** is one fused broadcast over a precombined vector
`prefac = ωwin·(-iω)·Pref` and the normalisation array, written in the same association as
the expression it replaces (`pre ./ (2 .* norm)`), with the `2` converted by
[`Luna.scalar`](@ref). `Pref` is folded into `prefac` rather than into the normalisation,
which stays physical.

### The free-space normalisation

[`NonlinearRHS.FreeSpaceNorm`](@ref) is shared by the radial and both Cartesian free-space
transforms, and is device-capable for all three.

The isotropic fill is one broadcast over

- `n(ω)`, evaluated on the host by scalar code (a Sellmeier equation, or a user's
  function) into a `(Nω, Npol)` [`Luna.HostMirror`](@ref) and uploaded once per call. It is
  host code and does not need a kernel: it runs once per `z`, not once per element.
- mirrors of `ω` and `grid.sidx`, `(Nω,)`, and of `kperp2` and the k-space window,
  reshaped to `(1, 1, Nk...)` so that they broadcast against the `(Nω, Npol, Nk...)`
  output.

The kernel is `normfactor` with every constant converted to the element type: `c`, `μ₀`,
`κmax` and `ℓ` are captured scalars, and `βz` is written as
`βsq < 0 ? complex(0, -√-βsq) : complex(√βsq, 0)` rather than with an `im` literal, which
is a `Complex{Int}`. Out of band and at `ω = 0` it returns exactly 1, which is what the
loop it replaces assigned.

**The k-window mirror is rebuilt lazily.** `Boundaries.setup` calls
[`NonlinearRHS.reflength!`](@ref) *after* the transform exists, to set `ℓ`, `κmax` and the
k-space absorber profile; that invalidates the mirror, and the next `fillnorm!` refills it.
Nothing else the kernel broadcasts against changes after construction.

The **crystal-optics** variant, whose per-`(ω, kx)` root-finding for the internal angle is
host scalar code with no kernel, stays on the host and copies its result up through a
staging buffer. It is reached only from the Cartesian transforms, which are host-only until
`gpu/21`.

Because the normalisation is a positional argument of `Luna.setup` -- every low-level
radial script builds one before it knows what device the run will use --
[`NonlinearRHS.retarget`](@ref) rebuilds it for the run's spec, carrying over anything
`reflength!` has already set. `norm_radial`/`norm_free`/`norm_free2D` and the `const_`
variants also take a `spec` keyword directly. Anything which is not a `FreeSpaceNorm` (a
user's own `normfun(z)`) is accepted only for the default host `Float64` path.

### The transverse collar

[`Boundaries.RadialCollar`](@ref) holds adapted copies of the same two matrices in the
**complex spectral** type (it is applied to `Eω` directly, the collar being diagonal in ω),
the absorption rate and the radial integration weights in the state's real precision, and a
buffer the size of `Eω`. Per accepted step it is one `radial_matmul!` k→r, one `mapreduce`
over a lazy `Broadcasted` for the energy bookkeeping, one broadcast for the absorbing
multiply and one `radial_matmul!` r→k. The running totals are accumulated in `Float64` on
the host from each step's reduction.

### Measurements

Radial Kerr propagation, M1 Pro, Julia 1.13.0, `-t 1`, one FFTW thread, one BLAS thread,
`:estimate`, no wisdom; 20 fixed steps over 1 cm of argon at 1 bar on a 100 fs /
400--2000 nm grid, `boundary=:none` (`benchmark/radial.jl`).

| radial points | | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | --- | ---: | ---: | ---: |
| 64 | Hankel GEMM | 95.7 µs | 48.2 µs | 170.7 µs |
| | right-hand side | 365.9 µs | 208.8 µs | 494.0 µs |
| | propagation | 72.1 ms | 51.6 ms | 93.6 ms |
| 256 | Hankel GEMM | 1.395 ms | 703.6 µs | 287.1 µs |
| | right-hand side | 3.480 ms | 1.858 ms | 773.3 µs |
| | propagation | 556.0 ms | 341.8 ms | 114.6 ms |
| 1024 | Hankel GEMM | 22.23 ms | 11.12 ms | 592.6 µs |
| | right-hand side | 47.68 ms | 24.09 ms | 1.215 ms |
| | propagation | 6.51 s | 3.47 s | 210.6 ms |

The crossover is at 64--96 radial points against the `Float64` host and 96--128 against
the `Float32` one. At 1024 points Metal is 31 times the `Float64` host and 16 times the
`Float32` one, and the Hankel GEMM alone is 38 times faster -- which is the whole reason
this geometry is the one worth putting on a GPU.

## The Cartesian free-space transforms

[`NonlinearRHS.TransFree2D`](@ref Luna.NonlinearRHS.TransFree2D) (2-D, `(t, x)`) and
[`NonlinearRHS.TransFree`](@ref Luna.NonlinearRHS.TransFree) (3-D, `(t, x, y)`) are the
same transform in two dimensionalities, and their per-step body is written once
([`NonlinearRHS.freetransform!`](@ref Luna.NonlinearRHS.freetransform!)). They are
separate types because `Boundaries.spacegrid` and `Luna.setup` dispatch on them and
because their transform region differs. Per right-hand side:

1. one inverse FFT over the time axis **and** the transverse axes together, taking the
   state from `(ω, k⊥)` to `(t, r⊥)`;
2. the response protocol on the whole `(nto, npol, Nk...)` block;
3. one broadcast for the temporal apodisation;
4. one forward FFT the other way;
5. one broadcast for the frequency-domain normalisation.

There is no matrix multiply: the transverse transform *is* the FFT. The region is
`(1, 3)` in 2-D and `(1, 3, 4)` in 3-D -- axis 2 is polarisation and is not transformed --
and both are planned through [`Utils.plan_ft`](@ref Luna.Utils.plan_ft) on the run's array
type. Metal.jl's own tests cover `(1, 3)` and `(1, 4)` but not `(1, 3, 4)`; the region is
checked directly against FFTW in `test_metal.jl` ("multi-axis FFT plans on Metal"), real
and complex, forward and inverse, and it is correct.

**The frequency-domain normalisation** is the same fused `fsnorm!` broadcast over a
precombined `prefac = ωwin·(-iω)·Pref` as the radial transform uses.

**One buffer fewer: `Pωo === Eωo`.** The oversampled frequency-domain buffer does double
duty. `to_time!` writes the field into it, applies the inverse plan and never reads it
again; `to_freq!` then writes the nonlinear polarisation into the same array. The two
never appear as the input and the output of the same FFT call, which a device plan would
reject. Each Cartesian transform therefore holds **three** field-sized arrays -- `Eto`,
`Pto` and the shared frequency-domain buffer -- where it used to hold four, and
`test_device.jl` asserts the count ("the free-space transforms hold three field-sized
buffers").

### Memory

Free-space geometry is where device memory starts to matter, because the state has
`prod(Nk)` columns. Counted from the shapes for the 3-D example
(`examples/low_level_interface/freespace/full3D.jl`: field-resolved, 400--2000 nm,
0.2 ps, `nt = 512`/`nto = 1024`, `nω = 257`/`nωo = 513`, 128 x 128 transverse, one
polarisation) in `Float32`:

| item | Kerr | Kerr + plasma | Kerr + Raman |
| --- | ---: | ---: | ---: |
| transform `Eto`, `Pto` | 128 MB | 128 MB | 128 MB |
| transform `Eωo === Pωo` | 64 MB | 64 MB | 64 MB |
| normalisation `out` | 32 MB | 32 MB | 32 MB |
| stepper (`y`, `yn`, `yi`, `yerr`, 7 stages) | 353 MB | 353 MB | 353 MB |
| linear operator | 32 MB | 32 MB | 32 MB |
| absorber `Et` | 32 MB | 32 MB | 32 MB |
| response block buffers | -- | 256 MB | 384 MB |
| **total** | **0.63 GB** | **0.88 GB** | **1.00 GB** |

Three things follow.

- **The examples' 3-D plasma run fits a 16 GB device eighteen times over**, so the
  response block is not chunked. The total scales as `Nx·Ny`: 0.88 GB at 128 x 128 is
  3.5 GB at 256 x 256 and 14.0 GB at 512 x 512, which is where a 16 GB device runs out.
  Chunking the response block along the transverse axes would move that limit by less
  than it looks: the batched buffers are 29 % of the total and the stepper's eleven
  state-sized arrays, which cannot be chunked without changing `RK45`, are 39 %.
- **The aliasing is worth about 7 %** of a 3-D run (64 MB of 0.88 GB at 128 x 128), and
  more of a transform-dominated one.
- A `Float64` host run is exactly twice these numbers.

The table is what the propagation holds. `Luna.setup` allocates a little more, all of it
collectable once it returns: on the **host**, the `Float64` prototypes the input-field
plans are made against (`(nt, npol, Nk...)` and `(nt, 2, Nk...)`, 64 MB and 128 MB at this
grid) and the initial state before it is uploaded (64 MB), and on the **device** one
state-shaped time-domain block for the state's own plan (32 MB). The *oversampled* block
is not among them: `setup_free` hands the array it planned `FTo` against to the transform
as its `Eto` ([`NonlinearRHS.freebuffers`](@ref Luna.NonlinearRHS.freebuffers)) instead of
leaving it to the garbage collector.

### The transverse collar

[`Boundaries.CartesianCollar`](@ref Luna.Boundaries.CartesianCollar) needed no change to
run on a device. The Cartesian grids transform time and space together, so when
`RateAbsorber` applies the temporal collar the state is already in `(t, x[, y])` and the
transverse collar is applied in the same pass: one broadcast for
`exp(-α Δz/2)` over the rate mirrored with [`Luna.upload_like`](@ref), one `sum(abs2, ·)`
and one `mapreduce` over a lazy `Broadcasted` for the energy bookkeeping, and one
broadcast for the multiply. No transform, no index vectors, no host scalar code.

### Measurements

3-D free-space envelope Kerr, M1 Pro, Julia 1.13.0, `-t 1`, one FFTW thread, one BLAS
thread, `:estimate`, no wisdom; 10 fixed steps over 1 cm of argon at 1 bar on a 100 fs /
400--2000 nm envelope grid (`nω = 128`), `boundary=:none` (`benchmark/free.jl`).

| transverse grid | | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | --- | ---: | ---: | ---: |
| 32 x 32 | joint inverse FFT | 1.148 ms | 916.5 µs | 313.4 µs |
| | right-hand side | 3.680 ms | 2.747 ms | 465.6 µs |
| | one step | 46.15 ms | 39.77 ms | 3.297 ms |
| | propagation | 477 ms | 412 ms | 56.0 ms |
| 64 x 64 | joint inverse FFT | 5.246 ms | 3.887 ms | 424.3 µs |
| | right-hand side | 16.21 ms | 11.52 ms | 851.8 µs |
| | one step | 195.0 ms | 164.1 ms | 7.413 ms |
| | propagation | 2.03 s | 1.69 s | 109 ms |
| 128 x 128 | joint inverse FFT | 36.05 ms | 17.52 ms | 1.057 ms |
| | right-hand side | 93.23 ms | 50.47 ms | 2.688 ms |
| | one step | 1.004 s | 700.3 ms | 25.93 ms |
| | propagation | 10.56 s | 7.27 s | 352 ms |

There is no crossover to report: the smallest grid in the sweep already has 1024
transverse columns. Per step -- the figure to quote, since `BenchmarkTools` repeats it
many times and it reproduces to about 3 % -- Metal is **14 times** the `Float64` host at
32 x 32, 26 times at 64 x 64 and **39 times** at 128 x 128. End to end the propagation is
8.5, 18.7 and 30.0 times faster; it is lower because a propagation also does its setup,
its output and, on a device, the host copy of each saved field, none of which the GPU
helps with. The gap is the FFT: at 128 x 128 the joint inverse transform alone is 34 times
faster on the GPU.

The `fft`, `rhs` and `step` columns reproduce to a few per cent between runs. The
`propagation` column is a single `@elapsed` per sample and scatters more: across two full
sweeps and an independent single-size run the same rows came out within 10 % of each other
(Metal at 64 x 64: 108.8, 108.8 and 121.8 ms), so read it to two figures. `proptime`
discards a warm-up run and takes the minimum of five, which is what makes even that much
reproducible -- with two samples and no warm-up, review 1 of `gpu/21-free-device` measured
a factor of 2.2 on the same row.

## The output and statistics boundary

`Output.jl` stays device-unaware: `MemoryOutput`/`HDF5Output` know how to save an array
and a dictionary, nothing about where the array lives or what units it is in. Two things
make that possible.

**[`Output.willsave`](@ref)`(o, y, t, dt)`** answers whether calling `o` right now would
save at least one data point, without saving anything. It exists so that a wrapper can
decide, cheaply, whether it is worth doing something expensive (a device-to-host copy)
before the real call.

**[`Luna.ScaledOutput`](@ref)** is that wrapper. `Luna.run` constructs one around the
output handler whenever the state is on a device or the run is scaled (`E_ref != 1`,
i.e. every `Float32` run, host or device); on the default CPU `Float64` path the handler
is passed through unchanged. It holds two reusable host buffers in the state's element
type (so a `Float32` run saves `Float32`, not `Float64`):

- `ybuf`, the unscaled host copy of the per-step solution `y`, filled when the handler's
  statistics need it (`Output.nostats` does not, and neither does a statistics set which
  runs on the device — see "Statistics" below) or when an `HDF5Output`'s resume cache does
  (gated by `willsave`, since the cache is written only on a save step);
- `ibuf`, the unscaled host copy of a *saved* field, filled lazily inside the closure
  `Luna.run`'s `stepfun` already passes as `yfun` — so a step which does not save costs
  nothing, and an `HDF5Output`'s `while save` loop, which can call `yfun` several times
  in one step, gets a distinct buffer from `ybuf`'s (the two can be needed together with
  different contents: an `HDF5Output`'s cache write reads the step endpoint `y` after
  possibly writing several interpolated `yfun(ts)` saves earlier in the same call).

The stepper's own arrays are never modified: `ScaledOutput` only ever copies *into* its
own buffers before handing them to the wrapped output. Multiplying by `E_ref` happens on
the host, after the copy, and is skipped entirely when `E_ref == 1`.

`MemoryOutput`/`HDF5Output` allocate their solution array with `eltype(y)`, which by the
time it reaches them through `ScaledOutput` is the state's own precision (`ComplexF32` or
`ComplexF64`), not a hardcoded `ComplexF64` — a `Float32` run therefore saves `Float32`.
`HDF5Output`'s resume cache stores the same host, physical-unit array `ScaledOutput`
copies for it, so `check_cache`'s result is always in physical units; `Luna.run` rescales
it (dividing by the *new* run's `E_ref`, deterministic from the same input) before
uploading it back to the state's array type. `HDF5Output`'s `cachehash` includes
`eltype(y)`, so resuming a run in a different precision than the cache was written in is
refused rather than silently misinterpreted.

[`Output.PeriodicStats`](@ref) (`prop_capillary`'s/`prop_gnlse`'s `stats_period` keyword)
reduces how often the statistics run, by evaluating the wrapped function only every
`period`-th accepted step and returning `nothing` in between —
`MemoryOutput`/`HDF5Output` skip appending a `nothing` result rather than erroring on it.

## Statistics

`Stats.jl`'s default sets are callable structs, not closures, so that they can carry three
traits:

| trait | what it says |
| --- | --- |
| [`Stats.device_capable`](@ref Luna.Stats.device_capable) | the statistic can be evaluated on the state as the stepper holds it: a device array, in scaled units |
| [`Stats.needs_time`](@ref Luna.Stats.needs_time) | it reads the time-domain field `Et`, so the inverse transform has to be done |
| [`Stats.statlabel`](@ref Luna.Stats.statlabel) | what to call it in a diagnostic |

The fallbacks are `false`, `true` and the type's name, which is what a user-written closure
gets: it is handed a host copy in physical units and works exactly as it always has.

Each struct has two branches. The **host branch** is the code Luna has always run, kept
arithmetically identical to it. The **device branch** is reductions and broadcasts over the
state. Which one runs is **fixed at `Stats.prepare`** and carried in the statistic's
`ondevice` field; it is not read off the array the statistic is handed, and a statistic
called with an array that disagrees errors instead of computing something wrong. The
branch is chosen by array type, not by precision, so a `Float32` run on the host takes the
host branch and the `Float64` regression gate is exactly `0.000e+00` on every row.

One thing about a `Float32` *host* run did change: `plan_analytic` now allocates the
`RealGrid` analytic buffer with `similar(Eω)`, so it follows the state's precision, where
it used to be `ComplexF64` whatever the state was. Every statistic that reads `Et` in a
`Float32` host run therefore reduces a `Float32` analytic field; measured against the old
arithmetic that is ~1.2e-8 on `peakpower` and `peakintensity` and ~5.4e-8 on `fwhm_t`. It
is deliberate: it is what makes a Metal run and a CPU `Float32` run comparable, it is what
the `EnvGrid` method always did, and it is what lets the transform be planned for the
array type it will be applied to. The gate's 21 cases are all `Float64` and do not see it.

The device branches are:

| statistic | on a device |
| --- | --- |
| `ω0` | two fused reductions, instead of materialising `abs2.(Eω)` |
| `energy`, `energy_λ`, `energy_window` | the frequency integral as one weighted reduction; a window folds into the weights squared |
| `peakpower`, `peakintensity` | one `maximum`; the on-axis form takes the single mode's transverse field out of the reduction as a scalar |
| `fwhm_t` | `abs2` into a device buffer, one copy of *that* to the host, then the same root-finding on the same samples |
| `electrondensity` (mode-averaged) | one broadcast of the ionisation-rate kernel off the state, then two reductions |
| `density`, `pressure`, `core_radius`, `zdw`, `zdz!` | they read only `z` |
| `fwhm_r`, `mode_reconstruction_error` | **not device-capable**: they keep their algorithms on the host |

`Fields.energyfuncs(grid)[2]` integrates the spectral power density with
`NumericalIntegration`'s `SimpsonEven` (a `RealGrid`) or a plain `sum` (an `EnvGrid`).
`SimpsonEven` is the alternative extended Simpson rule — every sample with weight 1 except
the first four and last four — so as a weight vector it is one fused reduction.
`Stats.energy` takes the energy functional as an argument, so the weights are *checked*
against it on a probe spectrum at construction: a functional which is not the grid's own
gives a host-only statistic rather than a silently different quantity.

`electrondensity` only ever reads the end point of the cumulative ionisation integral, and
that end point is the trapezoid rule over the whole window, so the device branch is one
weighted reduction rather than a prefix scan followed by a scalar read of the last element
(which a device array refuses). The intensity conversion and the unit scaling are folded
into the field reference `Ionisation.ratekernel` multiplies each sample by, so the field
itself is never rescaled.

Every weighted reduction folds over a lazy `Broadcast.Broadcasted` — `RK45._zipreduce` for
a whole-array reduction, `Stats._zipreduce1` for one along the frequency axis — and never
over `w .* abs2.(Eω)`, which materialises a field-sized temporary before reducing. That is
the same reason the stepper norms use `_zipreduce` (see "Amendments" in GPU_PLAN.md). The
unweighted reductions use the `mapreduce` forms (`sum(abs2, Eω)`,
`maximum(abs2, Et; dims=1)`), which are already allocation-free.

On a device the state is **scaled** (`e = E/E_ref`), because nothing has unscaled it —
`ScaledOutput` unscales only into a host buffer, which is the copy these branches exist to
avoid. The device branches therefore apply `E_ref` themselves, always to the scalar result
of a reduction and in `Float64`, never to the field: in a `Float32` run the field is
deliberately of order one and the physical units, up to 1e20 apart, would not survive being
put back on it. `Stats.StatsContext` carries `E_ref`; `Stats.collect_stats` takes it as a
keyword and `Stats.default` reads it off the transform.

`Stats.plan_analytic` builds its buffers with `similar`/`Luna.alloc` and plans through
`Utils.plan_ft`/`Utils.plan_ift`, so the analytic signal is one inverse FFT on whatever the
state lives on, shared by every statistic which needs it. When no statistic reads `Et` the
transform is not *applied*, but it is still planned and its buffers still allocated, which
is what Luna has always done.

## Which array a statistics set is built for

`Stats.collect_stats` makes one decision, once, about the array the whole set will be
called with, and everything — buffers, mirrors, the inverse plan, each statistic's
`ondevice` flag — is built for it. `Luna.stats_device_capable` reports that decision, and
`ScaledOutput` asks it once at construction (`devstats`) to decide whether to copy the
state to the host. The two cannot disagree, which is the point: a set built for the device
and handed a host array would put a host field and a device buffer in the same broadcast,
which JLArrays tolerates silently and Metal refuses.

The device state is used when it is a device array, when every statistic in the set has a
device form, and when `stats_device` allows it. `stats_device` is `:auto` (default),
`:device` or `:host`; under `:auto` the device state is used only when it has at least
`Stats.STATS_DEVICE_MINLEN` elements. **Below that the copy is cheaper**, and measurably
so: on an M1 Pro through Metal every statistic which ends in a device-to-host transfer
costs ~400 µs regardless of the size of the state, and the default set makes six of them
plus one MPSGraph inverse FFT — 2.7–3.6 ms for 1025 to 16385 elements, against 0.34–1.31 ms
for one transfer plus the host branches. On the mode-averaged Kerr case that is a 1.99×
slower accepted step, so the size test restores the host path there, which is what
`gpu/int-D` did. Luna logs which path a device run took, once, at construction.

The rule is the size of the state and nothing else. An earlier version also took the device
path for any state with more than one column; the column sweep in `benchmark/stats.jl` does
not support that — on Metal the device path is still slower at 16 and 128 columns, because
`fwhm_t` copies the time-domain intensity to the host on either path and its per-column
root-finding is host work either way, so extra columns alone do not make it pay. The
threshold is where the transfer of the state reaches the fixed cost of the round trips,
which is a few million elements. Every device state Luna produces today is one
mode-averaged column, far below it. `gpu/int-E` should re-measure on real radial and
free-space device states, where the transfer is much larger relative to the round trips,
and may well lower the threshold.

The fix that would make the device path win on a single column is to stop making six round
trips: a two-phase protocol in which each statistic writes its scalar reductions into one
small device buffer and the collector transfers that buffer once per call, finishing the
arithmetic on the host. By the numbers above that would take the default set from ~2.8 ms
to ~0.6 ms. It is a change to the `(d, Eω, Et, z, dz)` contract every statistic — including
a user's — is written against, so it is not done here.

Two arrays, not one, on a save step: an `HDF5Output` with a resume cache writes the raw
per-step `y` into the file as well as passing it to its statistics function, and those need
different arrays when the statistics are on the device. `ScaledOutput` therefore passes the
state itself as `y` and the host copy separately, as `Output.jl`'s `cache_y` keyword, which
defaults to `y` for every caller that does not wrap the output. Without that split, every
save step of a device run with `filepath=` set would hand a host array to a device set.

When the host path is taken because a statistic has no device form at all, `ScaledOutput`'s
one-time warning names it (`userfuns[1]` for a user's own). When it is taken because of the
shape test, only the `@info` line at construction says so — there is nothing the user did
wrong to warn about.

## The extension and hook mechanism

`ext/LunaMetalExt.jl` and `ext/LunaCUDAExt.jl` are package extensions under
`[weakdeps]`/`[extensions]`. Each registers its spec and its four vendor operations
(`functional`, `synchronize`, `reclaim`, `memory_status`) from `__init__` through
[`Luna.register_device!`](@ref).

They are *function-valued hooks*, not methods: a method with the same signature as a stub
in Luna itself would overwrite it, which Julia rejects during precompilation and which
leaves the extension uncompilable. The hooks are called through `Base.invokelatest` — the
GPU package may have been loaded after the calling frame was compiled — and only at setup
or teardown, never per step.

Registering sets `settings["device"] = :auto` only if the key is absent, so an explicit
`set_device(:cpu)` is never overridden.

Luna's own precompile block runs during Luna's precompilation, never with a GPU package
loaded, so it always precompiles the CPU path.

## Testing

- **The regression gate** (`test/test_regression.jl`) is the contract for the default CPU
  path: a matrix of small propagations against a baseline generated from the branch's base
  commit in the same environment. See `test/regression/README.md`.
- **`test/test_device.jl`** runs the device code paths on `JLArrays` with
  `allowscalar(false)`, against the host. It is part of `Pkg.test()`. Blind spots:
  `JLArrays` interprets its kernels on the host, so a stray `Float64` or a mixed
  host/device broadcast passes. `test/test_boundaries.jl` has its own, smaller JLArray
  testset for `Boundaries.RateAbsorber`/`LegacyAbsorber` built directly (not through a
  full propagation), gated the same way.
- **`test/test_metal.jl`** is the hardware test, and the only thing which catches those.
  It is not part of the suite (Metal is never installed with Luna); it has its own CI job,
  which installs Metal into a separate environment.
- **`benchmark/device.jl`** times the same mode-averaged propagation on each device and
  precision, sweeping the time-grid size; **`benchmark/radial.jl`** does the same for the
  radial transform, sweeping the number of radial points; **`benchmark/free.jl`** does the
  same for the 3-D Cartesian transform, sweeping the transverse grid.
- **`benchmark/modal.jl`** times a four-mode propagation on each transverse integral,
  device and precision.

## The modal transforms

A multimode propagation evaluates the nonlinear polarisation at transverse points and
integrates it against each mode's transverse field. That integral used to be evaluated
one point at a time: for every point, a matrix product to synthesise the field there, an
inverse transform of that one column, the responses on it, a forward transform and a
matrix product back. Both modal transforms now share one **batched column evaluator**,
[`NonlinearRHS.synthesise_responses!`](@ref Luna.NonlinearRHS.synthesise_responses!):

1. the modal spectrum goes to the oversampled time domain **once per right-hand side**,
   over the `nmodes` columns at once. Synthesis is linear and diagonal in time, so it
   commutes with the transform: synthesising and then transforming gives the same field
   as transforming and then synthesising, and the second order needs `nmodes` transform
   columns instead of one per point;
2. one matrix product (`Et = Emt S`) synthesises the field at every point of the current
   set, with `S` the mode matrix
   ([`Modes.mode_matrix`](@ref Luna.Modes.mode_matrix)) reshaped so that its column order
   is the `(nto, npol, npts)` block's — polarisation fastest;
3. the responses see the whole block through
   [`NonlinearRHS.Et_to_Pt!`](@ref Luna.NonlinearRHS.Et_to_Pt!), exactly as a radial or
   free-space transform's columns do.

What happens next is the only thing the two transforms do differently, and it is what
decides which one can run on a device.

### `TransModal`: the adaptive rule

[`NonlinearRHS.TransModal`](@ref) is the default (`modal_integral=:adaptive`). Its driver
is `Cubature.pcubature_v`/`hcubature_v`, which chooses the transverse points itself and
needs the integrand **at each point separately**, so the polarisation has to be
transformed back per point and projected per point. The driver is host scalar code and
returns the integral and its error estimate as `Vector{Float64}`, which is why this
transform is host- and `Float64`-only and says so, naming the fixed rule, rather than
being made parametric for a path it could not take.

The driver hands over a *round* of points at a time — 3, 2, 4, 8, … for `pcubature_v`,
17, 34, … for `hcubature_v`, independent of the integrand's dimension — and a round is
evaluated as one block, or as a few blocks where the round is wider than `maxbatch`
([`NonlinearRHS.MODAL_MAXBATCH`](@ref Luna.NonlinearRHS.MODAL_MAXBATCH)). Each distinct
width has its own buffers, its own forward plan *and its own copy of any response which
owns buffers*, because a [`Batched`](@ref Luna.Nonlinear.Batched) response sizes its
buffers to the block it is given; that is what the cap is for. The set of widths is fixed
by the rule, so the dictionary holding them stops growing after the first right-hand
side.

The projection is a broadcast rather than a matrix product: `out[ω, m, i] = pre[i] Σₚ
Pω[ω, p, i] W[m, p, i]` writes straight into the driver's buffer, reinterpreted as the
complex modal array. There are one or two polarisation components, so the sum over `p` is
unrolled.

### `TransModalFixed`: the fixed rule

[`NonlinearRHS.TransModalFixed`](@ref) (`modal_integral=:fixed`) evaluates the same
integral on a fixed quadrature rule
([`Modes.TransverseQuadrature`](@ref Luna.Modes.TransverseQuadrature); Gauss–Legendre or
Gauss–Kronrod in r or x, a periodic trapezoid in θ or Gauss–Legendre in y). Because the
rule's weights are known in advance, the points are summed **before** the transform back:
one matrix product `Pmt = Pt Wp` with the weights folded into `Wp`, then the time window,
one batched transform over the `nmodes` columns, the spectral window and the
normalisation. The time and spectral windows are diagonal in time and in frequency, and
the projection is a sum over points at fixed time, so applying them after the projection
is the same operation in a different order — and it means the number of transform columns
does not grow with the number of nodes. Everything per step is a matrix product, a
batched transform or a broadcast, which is the kernel discipline, so this is the
multimode transform which runs on a device.

The mode matrices are rebuilt when `z` moves, unless
[`Modes.zconstant`](@ref Luna.Modes.zconstant) says the transverse profiles do not depend
on it. That trait is `false` by default and `true` for a `Capillary.MarcatiliMode` with a
numeric core radius (and for the `Antiresonant` modes wrapping one), which is the fixed-
radius case; a taper re-evaluates the mode fields on the host and uploads them, which is
what the adaptive rule does at every point anyway.

The rule carries an embedded coarse rule — the Gauss subset of a Kronrod rule in r
(`kronrod=true`), or every other node of the *polar* θ trapezoid (a Cartesian domain's
second coordinate is Gauss–Legendre, which has none) — and
[`NonlinearRHS.integral_error!`](@ref Luna.NonlinearRHS.integral_error!) turns it into
`P_coarse - P_fine` with one further matrix product against the precomputed difference of
the two weight sets. Nothing evaluates it per step; it becomes a statistic once `Stats`
has been refactored.

`Stats.mode_reconstruction_error` is the adaptive transform's: it re-evaluates the
transform at one transverse point (which is what
[`NonlinearRHS.Erω_to_Prω!`](@ref Luna.NonlinearRHS.Erω_to_Prω!) is for) and records the
cubature's own error estimate. `prop_capillary` turns it off for a `:fixed` run and
errors if it is asked for explicitly.

### What moved

Replacing a per-point loop with a matrix product and a batched transform changes the
order of the arithmetic, which is what the regression gate allows and measures. Against
`gpu/int-D`, the two cases which go through a modal transform moved by 1.5e-15
(`modeavg_field_vector`, two modes and two polarisation components) and 6.9e-14
(`multimode_field_plasma`, four modes) in `Eω` in the fixed-step mode; every other case
is exactly zero, and the adaptive step counts did not change.

The fixed rule is a different *discretisation*, so it agrees with the adaptive rule to
the accuracy of the quadrature and not to rounding: 3.0e-16 (Kerr), 2.1e-14 (Kerr and
plasma) and 4.2e-16 (envelope Kerr) on one right-hand side of a four-HE₁ₘ-mode capillary
with `nr=64`, which is the adaptive rule's error rather than the fixed rule's.
