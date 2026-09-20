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
