# Running on a GPU

Luna can run the heavy part of a propagation on a GPU. Neither Metal nor CUDA is a
dependency of Luna: they are weak dependencies, loaded through package extensions, so
`Pkg.add("Luna")` on a machine without either installs and runs the CPU version.

!!! warning "Work in progress: everything but multimode propagation"
    This page describes what the device model does as of `gpu/21-free-device`.
    `prop_capillary` and `prop_gnlse` take `device` and `precision` keywords (below), and
    for mode-averaged propagation (`modes` a single mode) with the Kerr, plasma and
    Raman responses -- which is everything `prop_capillary` builds by default, in any gas
    -- it runs end to end on a device, including the absorbing boundaries
    (`boundary=:rate`, the default) and the default per-step statistics. All three
    free-space geometries -- radially symmetric, 2-D Cartesian and full 3-D, with the
    Kerr, plasma, Raman and χ⁽²⁾ responses -- do too, through the low-level interface,
    which is the only way to build one.
    Multimode propagation (`TransModal`) is still host code; `Luna.setup`/`Luna.run`
    refuse a device or a reduced precision for it rather than running it wrongly (or, for
    the simple interface, error with a message naming the actual limitation). It follows
    in `gpu/22-modal`; the page is completed in `gpu/32-docs`.

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
statistics, at the exit criterion's own parameters (100 nJ, 10 cm, Kerr only -- weakly
nonlinear: the spectrum barely broadens over that length) Metal agrees with a genuine CPU
`Float32` run to 4.4e-6 and with CPU `Float64` to 4.2e-6. With the nonlinearity clearly
visible (He at 5 bar, 300 µJ, ~×2.5 spectral broadening, fixed steps so the comparison is
of arithmetic rather than of step sequence) the agreement is 4.3e-6 -- comfortably inside
the 1e-4 the hardware tests assert throughout. (An earlier version of this page and of
`PR_11-boundaries-output.md` claimed the exit-criterion comparison was exactly `0.0` and
explained it as the nonlinear right-hand side rounding to zero below `Float32`'s precision
floor. That claim was wrong:
the comparison it was based on had no host reference at all -- both sides were resolving
to the GPU through `Luna.settings["device"] = :auto` -- so it was comparing Metal with
itself. Fixed in review round 1; see `PR_11-boundaries-output.md`'s "Changes after review
round 1".)

## What runs where

Anything Luna has not yet made device-capable runs on the host. At the moment that means
the multimode transform (`TransModal`) alone. The mode-averaged transform and all three
free-space ones -- `TransRadial`, `TransFree2D` and `TransFree` -- are device-capable, and
so is every nonlinear response Luna ships: the Kerr responses, the χ⁽²⁾ responses, the
plasma response and the Raman responses. `Luna.setup` refuses a device or a reduced
precision for a *transform* which cannot do it, through the residency checks each of them
makes, rather than running it wrongly. A *response* is not refused: it falls back to the
host copy described under "An ad hoc response on a device" below, which is correct and
slow.

`prop_capillary` and `prop_gnlse` never build a free-space run, so a free-space
propagation on a device is set up through the low-level interface:

```julia
using Luna, Metal
grid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 100e-15)
rg = Grid.RadialGrid(1e-3, 256)
nfunλ = PhysData.ref_index_fun(:Ar, 1.0)
nfun = (λ; z=0.0) -> nfunλ(λ)
linop = LinearOps.make_const_linop(grid, rg, nfun)
normfun = NonlinearRHS.const_norm_radial(grid, rg, nfun)
responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(:Ar)),)
inputs = Fields.GaussGaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6, w0=200e-6)
Eω, transform, FT = Luna.setup(grid, rg, z -> PhysData.density(:Ar, 1.0),
                               normfun, responses, inputs)  # device=:auto by default
output = Output.MemoryOutput(0, 1e-2, 11)
Luna.run(Eω, grid, linop, transform, FT, output; zmax=1e-2)
```

The `normfun` is built before the device is known, so `Luna.setup` moves it for you. On a
device or `Float32` run that means **replacing** it: `NonlinearRHS.retarget` builds a new
`FreeSpaceNorm` on the run's array type, the transform holds that one, and the object your
script still refers to is no longer part of the propagation — calling `reflength!` on it,
or reading its `out`, afterwards does nothing useful. Reach it through the transform
(`transform.normfun`), or build it on the device in the first place with the `spec` keyword
of `norm_radial` / `const_norm_radial`, in which case `retarget` returns it unchanged. On
the default host `Float64` path nothing is replaced.

The same applies to the two Cartesian free-space geometries: build the transverse grid
(`Grid.Free2DGrid` or `Grid.FreeGrid`), the normalisation (`const_norm_free2D` /
`const_norm_free`, or the crystal-optics pair from `PhysData.ref_index_fun_xy` for a
χ⁽²⁾ crystal) and the responses, and pass `device`/`precision` to `Luna.setup` in the same
way. The time and transverse axes are transformed together by one multi-axis FFT plan --
region `(1, 3)` in 2-D and `(1, 3, 4)` in 3-D -- on every backend.

This is the geometry where a GPU is worth using. A mode-averaged run has one transverse
column and is launch-bound; a free-space run has one per transverse grid point. Measured on an M1 Pro
(`benchmark/radial.jl`, 20 fixed steps over 1 cm of argon at 1 bar, a 100 fs / 400-2000 nm
grid, `boundary=:none`), wall time for the whole propagation:

| radial points | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | ---: | ---: | ---: |
| 64 | 72.1 ms | 51.6 ms | 93.6 ms |
| 96 | 125.9 ms | 86.4 ms | 83.5 ms |
| 128 | 192.0 ms | 126.8 ms | 102.8 ms |
| 256 | 556.0 ms | 341.8 ms | 114.6 ms |
| 1024 | 6.51 s | 3.47 s | 210.6 ms |

Metal passes the `Float64` host between 64 and 96 radial points and the `Float32` host
between 96 and 128; at 1024 it is 31 times the `Float64` host and 16 times the `Float32`
one. Below the crossover the CPU is faster and `device=:cpu` is the right answer.

A 3-D Cartesian run has `Nx*Ny` columns, so even the smallest useful grid is past the
crossover. Measured on the same machine (`benchmark/free.jl`, 10 fixed steps over 1 cm of
argon at 1 bar, a 100 fs / 400-2000 nm envelope grid, `boundary=:none`):

| transverse grid | columns | CPU `Float64` | CPU `Float32` | Metal `Float32` |
| ---: | ---: | ---: | ---: | ---: |
| 32 x 32 | 1024 | 478 ms | 412 ms | 61.1 ms |
| 64 x 64 | 4096 | 2.05 s | 1.69 s | 138 ms |
| 128 x 128 | 16384 | 10.96 s | 7.28 s | 521 ms |

Metal is 7.8 times the `Float64` host at 32 x 32 and 21 times at 128 x 128.

### Memory

A free-space state is `(nω, npol, Nk...)` and the transforms hold their buffers on the
oversampled time grid, so device memory, not speed, is what limits the grid. Each
transform holds **three** field-sized buffers -- the oversampled time-domain field and
polarisation, and one oversampled frequency-domain buffer used by both passes (`Pωo` and
`Eωo` are the same array). The stepper holds eleven state-sized ones, which is the largest
single item in any free-space run.

For the 3-D example (`examples/low_level_interface/freespace/full3D.jl`: a field-resolved
400-2000 nm, 0.2 ps grid on a 128 x 128 transverse grid, one polarisation) in `Float32`:

| | Kerr | Kerr + plasma | Kerr + Raman |
| --- | ---: | ---: | ---: |
| transform buffers | 192 MB | 192 MB | 192 MB |
| response buffers | 0 | 256 MB | 384 MB |
| stepper, operator, normalisation, absorber | 450 MB | 450 MB | 450 MB |
| **total** | **0.63 GB** | **0.88 GB** | **1.00 GB** |

so that run fits a 16 GB device with room to spare; the total scales as `Nx*Ny`, and
128 x 128 Kerr + plasma at 0.88 GB extrapolates to 3.5 GB at 256 x 256 and 14 GB at
512 x 512. `Float64` on the host is twice these numbers.

The absorbing boundaries (`boundary=:rate`, `:legacy` and `:none`) and the default
statistics *do* run with a device state:

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

### An ad hoc response on a device

A nonlinear response which is a plain closure `resp!(out, E, ρ)` (the way an ad hoc
response is usually written) works everywhere, including on a GPU and in single
precision. It is not run there: `Luna.setup` wraps it in a
[`Nonlinear.HostResponse`](@ref Luna.Nonlinear.HostResponse), which at every right-hand
side copies the whole field block to the host, converts it to physical SI units and
`Float64`, calls the response column by column exactly as the CPU path does, converts the
result back and adds it to the polarisation. One `@info` line at setup says which
response this applies to.

It is correct and slow. Two host copies and a host evaluation per right-hand side also
serialise the step on a GPU, so a propagation whose dominant nonlinearity goes through
this fallback will be slower than the same run on the CPU. It exists so that a response
you wrote yourself does not stop you using a device at all, not as a way to run one
quickly.

Because of that, the *simple* interface refuses instead of falling back: an explicit
`device` or `precision` request to `prop_capillary` with a response that has no device
kernel of its own -- which now means a response you wrote yourself, since every
response Luna ships has one -- is an error naming `device=:cpu`. A
call which does not mention `device` or `precision` is never affected — it stays on the
CPU as it always did. Use the low-level interface (`Luna.setup`/`Luna.run`) if you really
want the host fallback.

### Ionisation rates on a device

The plasma response evaluates its ionisation rate inside a kernel, which only two of
Luna's rates can do: [`Ionisation.IonRateADK`](@ref) (an analytic formula) and a cached
PPT rate, [`Ionisation.IonRatePPTAccel`](@ref)/`IonRatePPTCached`, whose table is
uniformly spaced — which every table Luna builds is. Those are what
`prop_capillary(...; plasma=true)` uses, so a default call needs no thought.

The direct [`Ionisation.IonRatePPT`](@ref), which sums a series and falls back to
`BigFloat`, cannot, and neither can a rate you wrote yourself. `Luna.setup` refuses those
on a device or in single precision with a message naming the alternatives; run on the CPU
with `device=:cpu`.

One behaviour differs on a device. The cached rate has a maximum field strength (twice
the barrier-suppression field), above which the CPU raises an error saying so. A device
kernel cannot raise one, so it returns the table's last value instead. The host path
keeps the error, which now comes from a check on the largest field in the block rather
than from every element.

**A weak plasma can vanish in single precision.** The plasma polarisation is built from
three cumulative integrals, and for a light gas at a moderate intensity the result is
small enough that the whole buffer is subnormal in `Float32`, which a GPU flushes to
zero. Helium at around 1e14 W/cm² is the case in Luna's usual range: there the plasma
term is about 2.6e-7 of the Kerr term — the size of a single `Float32` rounding of the
Kerr term itself — and `precision=Float32` or a Metal run drops it entirely. Use
`Float64` (`device=:cpu`, or CUDA, which runs in double precision) if a contribution that
small matters. At intensities where plasma actually shapes the pulse it is many orders
above the subnormal range and this does not arise. The developer guide has the full
dynamic-range audit.

### Raman on a device

A molecular gas runs on a device with no keywords of its own: the Raman polarisation
(`raman=true`, which `prop_capillary` sets by default for a gas that has one) is a
device-capable response like the Kerr and plasma ones. So is the no-THG Kerr response
`thg=false` selects for a field-resolved run.

Both do their work with FFTs along the time axis, batched over the transverse grid, so
the cost per step is a pair of transforms whatever the geometry rather than a pair per
column. The Raman response function itself is host scalar code and is evaluated only
when the density changes -- once, for a run at constant pressure.

**A pressure gradient is the bad case on a device.** With `pressure=(pin, pout)` the
density changes at every right-hand side, so every right-hand side pays a host
oscillator sum, a host FFT over the doubled time grid and a host-to-device copy. That
costs 0.32 ms at 16384 samples and 1.45 ms at 65536 against a Metal Raman right-hand side
of around 0.48 ms, so a differentially pumped capillary in a molecular gas runs two to
four times slower than the constant-pressure device run and gives back most of what the
GPU buys. Use `device=:cpu` for that configuration until the response function itself is
a kernel.

In single precision the Raman coefficients need care: the response function is around
1e-45 in SI units, far below what `Float32` can represent, and Luna splits the
coefficient so that no factor a kernel sees is subnormal. The developer guide has the
audit; there is nothing to set.

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
launch-bound and will not speed up. `benchmark/device.jl` times the mode-averaged
propagation on each device and precision, sweeping the grid size, and
`benchmark/radial.jl` does the same for the radial one, sweeping the number of radial
points, and `benchmark/free.jl` for the 3-D Cartesian one, sweeping the transverse grid
(the two tables under "What runs where" are their output). Run any of them from an
environment which has Luna, BenchmarkTools and the GPU package.

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
