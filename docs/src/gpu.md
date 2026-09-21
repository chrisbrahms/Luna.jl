# Running on a GPU

Luna can run the heavy part of a propagation on a GPU. Neither Metal nor CUDA is a
dependency of Luna: they are weak dependencies, loaded through package extensions, so
`Pkg.add("Luna")` on a machine without either installs and runs the CPU version.

!!! warning "Work in progress: mode-averaged and multimode propagation"
    This page describes what the device model does as of `gpu/22-modal`.
    `prop_capillary` and `prop_gnlse` take `device` and `precision` keywords (below), and
    for mode-averaged propagation (`modes` a single mode) with the Kerr, plasma and
    Raman responses -- which is everything `prop_capillary` builds by default, in any gas
    -- it runs end to end on a device, including the absorbing boundaries
    (`boundary=:rate`, the default) and the default per-step statistics. Multimode
    propagation does too, with `modal_integral=:fixed` (below).
    Anything else -- radial and free-space propagation, `prop_gnlse`, and the χ⁽²⁾
    responses -- is still host code; `Luna.setup`/`Luna.run` refuse a device or a
    reduced precision for a *transform* rather than running it wrongly, and fall back to
    the host for a *response* (or, for the simple interface, error with a message naming
    the actual limitation). Free space follows in a later branch; the page is completed
    in `gpu/32-docs`.

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
the radial and free-space transforms, and the *adaptive* multimode transform (see
"Multimode propagation" below). Every nonlinear response Luna ships
has a device kernel: the Kerr responses, the χ⁽²⁾ responses, the plasma response and the
Raman responses. The χ⁽²⁾ responses are only used by the free-space transforms, which are
not device-capable yet, so a χ⁽²⁾ propagation still runs on the host as a whole.
`Luna.setup` refuses a device or a reduced precision for the *transforms*, through the
residency checks each of them makes, rather than running them wrongly. A *response* is
not refused: it falls back to the host copy described under "An ad hoc response on a
device" below, which is correct and slow.

The absorbing boundaries (`boundary=:rate`, `:legacy` and `:none`) and the default
statistics *do* run with a device state, for the mode-averaged and the fixed-quadrature
multimode transforms:

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

### Multimode propagation

A multimode propagation evaluates the nonlinear polarisation at transverse points, then
integrates it against each mode's transverse field. `modal_integral` chooses how:

- `:adaptive` (the default) drives an adaptive cubature rule, which places transverse
  points itself until it reaches `radial_integral_rtol`. The driver
  (`Cubature.pcubature_v`/`hcubature_v`) is host scalar code and hands its results back
  as `Vector{Float64}`, so this cannot run on a device or in single precision; asking for
  one is an error naming `modal_integral=:fixed`.
- `:fixed` uses a fixed Gauss quadrature rule of `modal_nr` nodes along r (and, for the
  full 2-D integral, `modal_nθ` along θ; they are `nr` and `nθ` on `Luna.setup`). Every
  right-hand side costs the same, everything it does is a matrix product, a batched
  transform or a broadcast, and it runs wherever the mode-averaged transform does.

```julia
using Luna, Metal
prop_capillary(125e-6, 0.1, :Ar, 0.1;
               λ0=800e-9, energy=50e-6, τfwhm=20e-15,
               λlims=(200e-9, 3000e-9), trange=400e-15,
               modes=4, modal_integral=:fixed, modal_nr=64)
```

The two are different discretisations of the same integral, so they agree to the accuracy
of the quadrature rather than to rounding. For a set of HE₁ₘ modes, whose transverse
fields are smooth, the default 64-node rule is far more accurate than the adaptive rule
at its default 1e-3 tolerance: the two agree to 3e-16 (Kerr) and 2e-14 (Kerr and plasma)
on one right-hand side, which is the *adaptive* rule's error, not the fixed rule's.

Two things to check before using it:

- **`modal_nr` has to resolve the transverse structure of the highest mode.** The rule is
  not adaptive and will not tell you it is under-resolved. 64 is the default; a mode set
  reaching HE₁₈ or beyond wants more.
- **`modal_nθ` has to be at least `4h+1`** for modes of azimuthal order up to `h` (`h` is
  `|n-1|` for an HE\_{nm} mode), because the θ rule is a periodic trapezoid and the
  integrand of a cubic response projected back onto a mode reaches the harmonic `4h`.
  Luna warns at setup when it can tell that it is too small. An HE₁ₘ set has `h = 0`, so
  any number of θ nodes will do, and `full=false` (which `prop_capillary` picks for such a
  set) uses a single one.

`modal_integral=:fixed` does not collect the `mode_reconstruction_error`,
`transverse_points` and `transverse_integral_error_*` statistics: those describe the
adaptive rule's own behaviour. The fixed rule carries an embedded Gauss--Kronrod error
estimate instead (`modal_kronrod=true`, which rounds `modal_nr` up to an odd number),
which `NonlinearRHS.integral_error!` evaluates on demand; it becomes a statistic in a
later branch. At the low level, `Stats.default` needs `mode_error=false` for this
transform; `prop_capillary` does that for you.

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
launch-bound and will not speed up. `benchmark/device.jl` times the same propagation on
each device and precision, sweeping the grid size, and `benchmark/modal.jl` does the same
for a four-mode propagation on each transverse integral. Run them from an environment
which has Luna, BenchmarkTools and the GPU package.

Multimode is where a GPU has something to work on: the fixed rule evaluates the
responses at all `nr` transverse points at once. Measured on an M1 Pro, four HE₁ₘ modes
of a 75 µm capillary in argon, 10 fixed steps, one FFTW and one BLAS thread, one
right-hand side:

| | Kerr, 4100-sample state | Kerr, 16388 | Kerr+plasma, 4100 | Kerr+plasma, 16388 |
| --- | ---: | ---: | ---: | ---: |
| CPU `Float64` `:adaptive` (the default) | 938 µs | 4.12 ms | 2.96 ms | 11.9 ms |
| CPU `Float64` `:fixed`, `nr=64` | 306 µs | 1.60 ms | 4.34 ms | 16.7 ms |
| CPU `Float32` `:fixed`, `nr=64` | 165 µs | 714 µs | 3.82 ms | 15.4 ms |
| Metal `Float32` `:fixed`, `nr=64` | 410 µs | 733 µs | 956 µs | 1.60 ms |

Two things to read off. The fixed rule is not automatically cheaper on the CPU: the
adaptive rule needs only about 17 transverse points for a smooth HE₁ₘ set, so a 64-node
rule does roughly four times the work, which the Kerr rows hide (the transform count does
not grow with `nr`) and the plasma rows do not. And the GPU's advantage grows with the
work per node: 1.0× for Kerr at the small grid, 9.6× for Kerr and plasma at the large
one against the same arithmetic on the CPU, 7.4× against the `Float64` adaptive default.

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
