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

Quantities which host scalar code still has to produce on every right-hand side — `β(z)`
for a taper or a pressure gradient, a user-supplied `linop!` — go through a
[`Luna.HostMirror`](@ref): a `Float64` host buffer, an optional staging buffer in the
device precision, and the device array. On the host in double precision the device array
*is* the host buffer and `upload!` does nothing. This is the interim arrangement until
`gpu/23` tabulates them.

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
| [`Nonlinear.Batched`](@ref) | calls it once with the whole `(nt, npol, ncols)` block, through [`batched!`](@ref Luna.Nonlinear.batched!), which also hands it the unit scaling | [`batched!`](@ref Luna.Nonlinear.batched!) (or just the call operator, if its coefficients already carry the scaling), plus its own full-size buffers in the run's array type | `HostResponse`; `PlasmaCumtrapz`, `RamanPolar*` after Group D |
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
`real(eltype(E))`. That is 27 host scalars once per right-hand side and the identity on a
response `rescale` has already converted; it is there so that a kernel built from a
response which never went through `rescale` still carries no `Float64` (rule 3 above),
which is what a low-level caller assembling `Et_to_Pt!` by hand does.

The χ⁽²⁾ transforms themselves — the free-space ones — are host-only until Group E, so the
responses are exercised on a device block directly (`test_device.jl`, `test_metal.jl`)
rather than through a propagation.

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
  statistics need it (`Output.nostats` does not) or when an `HDF5Output`'s resume cache
  does (gated by `willsave`, since the cache is written only on a save step);
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

Per-step statistics (`Stats.jl`) are unchanged host code: they do a full inverse FFT and
host reductions on whatever `y` `ScaledOutput` hands them, which on a device is a copy
made every accepted step. `ScaledOutput` warns once per propagation when this happens.
[`Output.PeriodicStats`](@ref) (`prop_capillary`'s/`prop_gnlse`'s `stats_period` keyword)
reduces how often that copy and the statistics themselves run, by evaluating the wrapped
function only every `period`-th accepted step and returning `nothing` in between —
`MemoryOutput`/`HDF5Output` skip appending a `nothing` result rather than erroring on it.
Device-capable statistics (skipping the copy) are `gpu/24`'s.

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
- **`benchmark/device.jl`** times the same propagation on each device and precision.
