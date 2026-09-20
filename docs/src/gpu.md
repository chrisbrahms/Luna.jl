# Running on a GPU

Luna can run the heavy part of a propagation on a GPU. Neither Metal nor CUDA is a
dependency of Luna: they are weak dependencies, loaded through package extensions, so
`Pkg.add("Luna")` on a machine without either installs and runs the CPU version.

!!! warning "Work in progress: only mode-averaged Kerr propagation runs on a device"
    This page describes what the device model does as of `gpu/11-boundaries-output`.
    `prop_capillary` and `prop_gnlse` take `device` and `precision` keywords (below), and
    for mode-averaged propagation with Kerr responses (`modes` a single mode; no plasma,
    no Raman, no χ⁽²⁾) `prop_capillary` runs end to end on a device: the absorbing
    boundaries (`boundary=:rate`, the default) and the default per-step statistics both
    work now. Anything else -- multimode and radial propagation, `prop_gnlse`, and every
    response but Kerr -- is still host code; `Luna.setup`/`Luna.run` refuse a device or a
    reduced precision there rather than running it wrongly (or, for the simple interface,
    error with a message naming the actual limitation). Free space and multimode
    propagation, and the other nonlinear responses, follow in later branches; the page is
    completed in `gpu/32-docs`.

## Enabling it

```julia
using Luna
using Metal      # or: using CUDA
```

Loading the GPU package registers the backend and sets `Luna.settings["device"] = :auto`
**if the key is absent**, so a script which has not said anything about devices starts
using the GPU -- including a plain `prop_capillary(...)` call, for the cases it can run
on one (mode-averaged, Kerr only). To opt out:

```julia
Luna.set_device(:cpu)
```

An explicit setting is never overridden by loading a GPU package later, in either order.
[`Luna.device`](@ref) reports what the next run will use, and `Luna.setup` logs it once
per run:

```
[ Info: Propagating on MtlArray in Float32 precision.
```

`Luna.set_device` also accepts `:metal`, `:cuda` and a [`Luna.DeviceSpec`](@ref)
directly, which is how a non-default precision is requested:

```julia
Luna.set_device(Luna.DeviceSpec(CUDA.CuArray, Float32))
```

`Luna.setup`, `prop_capillary` and `prop_gnlse` all take `device` and `precision`
keywords which override the global setting for one call:

```julia
# Explicit, whatever Luna.settings["device"] says
out = prop_capillary(125e-6, 0.1, :He, 1.0; λ0=800e-9, τfwhm=10e-15, energy=100e-9,
                     plasma=false, raman=false, λlims=(200e-9, 3e-6), trange=200e-15,
                     device=:metal)

# Force the CPU for one call without touching the global setting
out = prop_capillary(...; device=:cpu)
```

`prop_gnlse` accepts the same two keywords for a uniform call signature, but is not
device- or reduced-precision-capable yet: it builds its own normalisation before the unit
scaling is known, so anything other than the default (the CPU, `Float64`) errors.
Multimode and radial propagation (`prop_capillary` with `modes` a collection) are the
same: the keywords are accepted and validated, and a request that cannot be honoured
errors naming the actual limitation, rather than running silently on the CPU or failing
with an unrelated `MethodError`.

### Scans

`Scans` workers load Luna on their own, so a scan needs the GPU package on every worker:

```julia
using Distributed
@everywhere using Luna, Metal
```

Without that, the workers see `:auto` with no registered device and run on the CPU. They
say so once:

```
[ Info: `:auto` requested but no GPU package is loaded in this process; running on the CPU.
```

## Precision

Metal refuses `Float64` arrays and its kernel compiler rejects `double` arithmetic, so a
Metal run is single precision. CUDA runs in double precision by default, so a CUDA result
is directly comparable with a CPU one.

Single precision is not just less accurate, it has less *range*, and SI units do not fit
in it: the combined Kerr coefficient `ρ ε₀ γ₃` for helium at 0.3 bar is 3.3e-39, below the
smallest normal `Float32` (1.2e-38), and GPUs flush subnormals to zero. Luna therefore
rescales a single-precision run: the state is `e = E/E_ref` with `E_ref` the peak of the
input field rounded to a power of two, the polarisation is measured in units of `ε₀E_ref`,
and every response's coefficients are combined with those references on the host in
`Float64` before being converted. See [`Luna.UnitScaling`](@ref). `E_ref` is exactly 1 for
every `Float64` run, so the default CPU path is arithmetically unchanged.

Measured on an Apple M1 Pro, for a mode-averaged Kerr propagation of a 20 fs, 1 µJ pulse
through 1 cm of helium-filled capillary, as the largest relative difference in `Eω` over
the saves:

| comparison | field-resolved | envelope |
| --- | ---: | ---: |
| Metal `Float32` vs CPU `Float32` | 2.7e-7 | 5.2e-8 |
| Metal `Float32` vs CPU `Float64` | 3.0e-7 | 7.3e-8 |
| CPU `Float32` vs CPU `Float64` | 3.5e-7 | 7.3e-8 |

At 0.3 bar -- the case whose unscaled coefficient does not exist in `Float32` at all --
Metal agrees with the `Float64` CPU path to 2.4e-7.

With the absorbing boundaries (`boundary=:rate`, the default) and the default
statistics, Metal vs CPU `Float32` for the same propagation is 3.9e-6 (field-resolved)
and 4.4e-6 (envelope) -- still three orders of magnitude inside the exit-criterion
tolerance, the extra digit coming from the collar broadcasts and the host round trip the
statistics need. `prop_capillary` itself, at the exit criterion's own parameters (100 nJ,
1 cm, Kerr only), agrees exactly: at that pulse energy the nonlinear phase is far below
`Float32`'s precision floor, so the right-hand side rounds to zero on both backends and
the propagation is the linear operator's `exp`, which has no backend-dependent
summation order.

## What runs where

Anything Luna has not yet made device-capable runs on the host. At the moment that means
plasma, Raman and χ⁽²⁾ responses, and the radial, free-space and multimode transforms:
`Luna.setup` refuses a device or a reduced precision for them, through the residency
checks every transform and response makes, rather than running them wrongly.

The absorbing boundaries (`boundary=:rate`, `:legacy` and `:none`) and the default
statistics *do* run with a device state, for the mode-averaged transform:

- `Boundaries.RateAbsorber`/`LegacyAbsorber` and the transverse collars are broadcasts
  and reductions over mirrored arrays (`Boundaries.jl`), like everything else per-step.
- Per-step statistics (`Stats.jl`) are still host code: the field is copied to the host
  every accepted step to compute them, and `Luna.run` warns once when this happens. Use
  `stats_period` to reduce how often they run (below), or `Output.nostats` to disable
  them; device statistics (computing them without the copy) are `gpu/24`'s.
- The output itself never sees a device array or a scaled one: `Luna.run` wraps it in
  `Luna.ScaledOutput`, which copies to the host and, for a `Float32` run, unscales, before
  handing it to `Output.MemoryOutput`/`HDF5Output`. A `Float32` run's saved field is
  therefore in physical units already, and is `ComplexF32` -- `eltype(y)` is what the
  output allocates with, not always `ComplexF64`.

A nonlinear response which is a plain closure (the way an ad hoc response is usually
written) still works on the default CPU path in double precision. On a device or in single
precision it is refused, because its coefficients cannot be rescaled and its body is
usually scalar code. The fallback which copies a column block to the host and runs it there
is `gpu/12`.

### `stats_period`

```julia
out = prop_capillary(...; stats_period=10)
```

Collects the default statistics every 10th accepted step instead of every step
(`Output.PeriodicStats`). The recorded statistics arrays are correspondingly shorter; the
saved field (`saveN`, `out["Eω"]`) is unaffected. Worth raising on a device, where
per-step statistics force a host copy every accepted step regardless of how often the
propagation actually saves the field.

## Performance

A GPU wins on many columns and large time grids; a mode-averaged single-column run is
launch-bound and will not speed up. `benchmark/device.jl` times the same propagation on
each device and precision, sweeping the grid size. Run it from an environment which has
Luna, BenchmarkTools and the GPU package.

## Running the hardware tests

```
julia --project=<env> -e '
  using Pkg
  Pkg.develop(path=".")
  Pkg.add(["Metal", "Test", "GPUArraysCore", "Adapt"])'
julia --project=<env> -t 1 -e 'using Luna, Metal; include("test/test_metal.jl")'
```

`<env>` needs Luna developed into it *and* the packages the test file imports directly,
which is the set the CI job installs. The file skips itself if Metal is not loaded or not
functional.

`test/test_device.jl`, which runs the same code paths on `JLArrays` without any GPU, is
part of `Pkg.test()`. On its own it needs
`Pkg.add(["JLArrays", "Test", "GPUArraysCore", "Adapt", "AbstractFFTs", "FFTW"])` in the
same way.
