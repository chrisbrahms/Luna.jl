# Running on a GPU

Luna can run the heavy part of a propagation on a GPU. Neither Metal nor CUDA is a
dependency of Luna: they are weak dependencies, loaded through package extensions, so
`Pkg.add("Luna")` on a machine without either installs and runs the CPU version.

!!! note "Work in progress"
    This page describes what the device model does as of `gpu/10-device-model`, the first
    branch of the GPU work. At this point only mode-averaged propagation with Kerr
    responses runs on a device, through the low-level interface, and with
    `boundary=:none`. The absorbing boundaries, the output wrapper and the simple
    interface (`prop_capillary`, `prop_gnlse`) follow in `gpu/11`; plasma, Raman, χ⁽²⁾,
    radial, free-space and multimode propagation follow after that. The page is completed
    in `gpu/32-docs`.

## Enabling it

```julia
using Luna
using Metal      # or: using CUDA
```

Loading the GPU package registers the backend and sets `Luna.settings["device"] = :auto`
**if the key is absent**, so a script which has not said anything about devices starts
using the GPU. To opt out:

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

`Luna.setup` takes `device` and `precision` keywords which override the global setting for
one run.

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

## What runs where

Anything Luna has not yet made device-capable runs on the host. At the moment that means
a device run has to use `boundary=:none` and an output which collects no per-step
statistics; `Luna.run` says so rather than falling back silently.

A nonlinear response which is a plain closure (the way an ad hoc response is usually
written) still works on the default CPU path in double precision. On a device or in single
precision it is refused, because its coefficients cannot be rescaled and its body is
usually scalar code. The fallback which copies a column block to the host and runs it there
is `gpu/12`.

## Performance

A GPU wins on many columns and large time grids; a mode-averaged single-column run is
launch-bound and will not speed up. `benchmark/device.jl` times the same propagation on
each device and precision, sweeping the grid size. Run it from an environment which has
Luna, BenchmarkTools and the GPU package.

## Running the hardware tests

```
julia --project=<env> -t 1 -e 'using Luna, Metal; include("test/test_metal.jl")'
```

where `<env>` is an environment with Luna developed into it and Metal added. The file
skips itself if Metal is not loaded or not functional. `test/test_device.jl`, which runs
the same code paths on `JLArrays` without any GPU, is part of `Pkg.test()`.
