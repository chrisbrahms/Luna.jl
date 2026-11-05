# Running on a GPU

Luna can run the heavy part of a propagation on a GPU. Neither Metal nor CUDA is a
dependency of Luna: they are weak dependencies, loaded through package extensions, so
`Pkg.add("Luna")` on a machine without either installs and runs the CPU version.

!!! warning "Work in progress"
    This page describes what the device model does as of `gpu/int-E`.
    `prop_capillary` and `prop_gnlse` take `device` and `precision` keywords (below), and
    for mode-averaged propagation (`modes` a single mode) with the Kerr, plasma and
    Raman responses -- which is everything `prop_capillary` builds by default, in any gas
    -- it runs end to end on a device, including the absorbing boundaries
    (`boundary=:rate`, the default) and the default per-step statistics. All three
    free-space geometries -- radially symmetric, 2-D Cartesian and full 3-D, with the
    Kerr, plasma, Raman and χ⁽²⁾ responses -- do too, through the low-level interface,
    which is the only way to build one, and so does multimode propagation with
    `modal_integral=:fixed` (below).
    What is left on the host is `prop_gnlse`, the *adaptive* multimode transverse
    integral (`modal_integral=:adaptive`, the default, whose cubature driver is host
    scalar code), the per-step evaluation of the crystal-optics normalisation, and the
    `fwhm_r` statistic. `Luna.setup`/`Luna.run` refuse a device or a reduced precision
    for a *transform* which cannot do it rather than running it wrongly, and fall back
    to the host for a *response* (or, for the simple interface, error with a message
    naming the actual limitation). The page is completed in `gpu/32-docs`.

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
the *adaptive* multimode transform (`TransModal`, `modal_integral=:adaptive`, the default
-- see "Multimode propagation" below) and `prop_gnlse`. The mode-averaged transform, all
three free-space ones -- `TransRadial`, `TransFree2D` and `TransFree` -- and the
fixed-quadrature multimode transform (`TransModalFixed`, `modal_integral=:fixed`) are
device-capable, and so is every nonlinear response Luna ships: the Kerr responses, the
χ⁽²⁾ responses, the plasma response and the Raman responses. Two pieces of per-step work
are still host code inside otherwise device-capable runs: the crystal-optics
normalisation evaluates its refractive indices on the host and stages the result (the
propagation itself stays on the device), and the `fwhm_r` statistic has no device form.
`Luna.setup` refuses a device or a reduced precision for a *transform* which cannot do
it, through the residency checks each of them makes, rather than running it wrongly. A
*response* is not refused: it falls back to the host copy described under "An ad hoc
response on a device" below, which is correct and slow.

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

Free space is the geometry where a GPU is worth using. A mode-averaged run has one
transverse column and is launch-bound; a radial run has one per radial point, and a
Cartesian one has `Nx` or `Nx*Ny` of them.

The radial transform, measured on an M1 Pro (`benchmark/radial.jl`, 20 fixed steps over
1 cm of argon at 1 bar, a 100 fs / 400-2000 nm grid, `boundary=:none`), wall time for the
whole propagation:

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
| 32 x 32 | 1024 | 477 ms | 412 ms | 56.0 ms |
| 64 x 64 | 4096 | 2.03 s | 1.69 s | 109 ms |
| 128 x 128 | 16384 | 10.56 s | 7.27 s | 352 ms |

Metal is 8.5 times the `Float64` host at 32 x 32 and 30 times at 128 x 128 end to end;
measured per step, which reproduces to a few per cent where the propagation column
scatters by up to 10 %, it is 14 times and 39 times. Read the table to two figures.

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

Setting up costs a little more than the propagation holds, all of it collectable once
`Luna.setup` returns: the `Float64` prototypes its input-field plans are made against and
the initial state before it is uploaded (about 256 MB on the host at that grid), and one
state-shaped time-domain block on the device.

The absorbing boundaries (`boundary=:rate`, `:legacy` and `:none`) and the default
statistics *do* run with a device state, for every device-capable transform -- the
mode-averaged one, the three free-space ones and the fixed-quadrature multimode one:

- `Boundaries.RateAbsorber`/`LegacyAbsorber` and the transverse collars are broadcasts
  and reductions over mirrored arrays (`Boundaries.jl`), like everything else per-step.
- The default per-step statistics (`Stats.jl`) have a device form -- the energies, the
  peak power and intensity, the temporal FWHM, the electron density and the z-dependent
  quantities are reductions and broadcasts over the state where it is -- but **which path
  a run takes depends on the size of the state**, and Luna logs which one it chose:

  - a state with at least `Stats.STATS_DEVICE_MINLEN` elements computes its statistics on
    the device, with no copy;
  - a smaller one is copied to the host and the statistics are computed there, because the
    copy is cheaper. Each statistic that ends in a device-to-host transfer costs the same
    round trip whatever the size of the state -- about 400 µs on an M1 Pro through Metal --
    and the default set makes six of them, where the host path makes one transfer. On the
    mode-averaged Kerr case that is the difference between a 2.3 ms and a 4.6 ms accepted
    step.

  In practice the host path is taken for everything Luna produces today. The threshold
  was re-measured on `gpu/int-E` on real radial and 3-D Metal states (M1 Pro, Metal 1.11,
  one call of a set of `ω0`, `energy`, `peakpower`, `fwhm_t` and `density`, against one
  RK45 step of the same propagation):

  | state | elements | host + copy | device | one step |
  | --- | ---: | ---: | ---: | ---: |
  | radial, 256 radial points | 33024 | 2.03 ms | 333 ms | 3.5 ms |
  | radial, 1024 radial points | 132096 | 7.65 ms | 1.44 s | 10.5 ms |
  | 3-D, 64 x 64 | 524288 | 661 ms | 615 ms | 9.3 ms |
  | 3-D, 128 x 128 | 2097152 | 5.23 s | 4.88 s | 26.1 ms |

  so the device path does not pay at any of them and the threshold stays at 2^22. Two
  things to read off it. On a radial state the device path is two orders of magnitude
  worse, and all of it is `fwhm_t` (303 ms of the 333): its device branch reduces the
  *scaled* field, so the columns far off axis underflow to exactly zero in `Float32` and
  the host root-finding which follows is far slower on them than on the small but non-zero
  numbers the host path gives it. And on a many-column state both paths are dominated by
  that same per-column root-finding on the host -- 0.6 s at 4096 columns and 5 s at 16384,
  against a 9 ms and a 26 ms step -- so per-step statistics are not usable on a large
  free-space grid on either path. Raise `stats_period`, or use `Output.nostats`. To
  override the choice, `stats_kwargs=Dict(:stats_device => :device)` (or `:host`) reaches
  `Stats.collect_stats` through `prop_capillary`.

  `fwhm_r` and the modal reconstruction error have no device form at all and keep their
  algorithms on the host, so a multimode set is computed on a host copy whichever
  transverse integral it came from. A statistics function *you* write is host code too:
  the whole set is then
  computed on a host copy of the field on every step the statistics fire, and `Luna.run`
  warns once, naming it (`userfuns[1]`), when that happens. Use `stats_period` to reduce
  how often they run (below), or `Output.nostats` to switch them off.
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

### Tapers and pressure gradients (`linop_integral`)

```julia
out = prop_capillary(125e-6, 0.1, :He, (1.0, 0.0); ...)                             # :auto
out = prop_capillary(125e-6, 0.1, :He, (1.0, 0.0); ..., linop_integral=:quadrature) # or this
```

A capillary whose radius or pressure changes along `z` has a linear operator `L(z)` which
changes with it. The stepper propagates the linear part of a step by

```
exp(Φ(t2) − Φ(t1)),   Φ(z) = ∫ L dz',
```

which is exact, so what it needs is not the operator but its integral. `linop_integral`
says how that integral is obtained, and there are two built-in answers plus the option of
writing your own.

A `linop!` which ignores `z` is recognised first and used as the constant operator it is,
whichever setting you pick — a `Capillary.gradient` with the same pressure at both ends, a
taper function which returns a constant, `LinearOps.make_linop` given a refractive index
which does not depend on `z`. Such a run is unaffected by anything on this page.

**What changed.** Until this release Luna propagated a z-dependent operator with
`exp(L(t2)·(t2 − t1))`, a one-point rule. It is first order in the step size, and — this is
the part which matters in practice — its error is common to both of the embedded
Runge–Kutta solutions the step-size controller forms its error estimate from, so it cancels
out of the estimate and no value of `rtol` responds to it. **Every taper and pressure
gradient result therefore changes**, by the amount that rule was wrong:

| 0.1 m Ar capillary, fixed steps, difference from a 10240-step reference | 20 steps | 80 | 320 | 1280 |
| --- | --- | --- | --- | --- |
| pressure gradient, one-point rule (before) | 7.5e-2 | 2.0e-2 | 4.9e-3 | 1.1e-3 |
| pressure gradient, `:tabulated` or `:quadrature` (now) | 1.6e-4 | 1.6e-4 | 1.6e-4 | — |
| taper 75 → 50 µm, one-point rule (before) | 4.1e-1 | 1.0e-1 | 2.5e-2 | 5.6e-3 |
| taper, `:tabulated` or `:quadrature` (now) | 8.0e-4 | 8.0e-4 | 8.0e-4 | — |

(largest relative difference of `Eω` per save; gradient 0 → 1 bar, Kerr only,
`boundary=:none`. The reference is itself a 10240-step run of the old rule, so the 1.6e-4
and 8.0e-4 left in the second and fourth rows are mostly the reference's own remaining
error.) The new propagator is at the well resolved answer from the coarsest step count
tried; the old one needed about a thousand steps to get there. A result computed now is
therefore *not* comparable element by element with one stored before, and the difference is
larger than any tolerance — the regression gate records 8.0e-3 for its gradient case and
9.2e-2 for its taper case. A uniform fibre has a constant operator, which was always exact
in the propagator, and is bit-for-bit unchanged.

**`:tabulated`** evaluates `Φ` at setup on `z` nodes placed by adaptive
bisection and stores the table where the state lives, together with the propagation
constant `β(z)` and effective area `Aeff(z)` the mode-averaged nonlinear normalisation
needs. What is left inside the propagation is an interval lookup and four interpolation
weights per readback — scalar host arithmetic, no mode evaluation and no copy between host
and device. This is what makes a tapered or pressure-graded run go entirely on a GPU: the
alternative is six host evaluations of `Modes.neff` and six uploads per step, which
serialise a device run against the host.

`linop_tol` (default `1e-6`) sets how accurate the table is: absolute, in radians, for
`Φ`, and relative for `β` and `Aeff`. The bisection checks the interpolant it will actually
read back against a directly computed value at each candidate interval's midpoint, so a
kink — the `1/√z` cusp in the density of a gradient filled from vacuum, a junction in a
multi-section fill — costs a handful of extra intervals where it is rather than a finer
grid everywhere. A 0.1 m gradient takes about 60 nodes at the default tolerance and a taper
about 30.

The table costs `2·length(Eω)·nnodes` numbers in the state's precision, and about twice
that again on the host while it is being built. For a mode-averaged run that is a few
megabytes. For a multimode, radial or free-space operator, whose `linop` is the size of the
whole state, it is several times `nnodes` times the state; Luna warns before it starts
building if that could exceed 256 MB, and `:auto` bounds it instead of warning.

**`:quadrature`** integrates the operator over each step instead, by adaptive
Gauss–Kronrod quadrature on the host, and uploads the result. It holds no table, needs no
setup pass and makes no assumption about how the operator behaves in `z`, but it costs
about fifteen host evaluations of the operator per stage against one table readback, and
`β` and `Aeff` are evaluated per stage as well. It computes the same integral as the table
does — the two agree to 1e-7 of each other on the cases above — so it is the thing to check
a tabulated result against, and the thing to use when the table would be too large.

**`:auto` is the default**, and it is `:tabulated` with a memory bound: the node count is
capped so that the table's peak cost stays inside `linop_budget` (256 MB, counted at the
peak, which is dominated by the host copies held while the table is built), and if that cap
binds before `linop_tol` is met the run falls back to `:quadrature` and logs a line saying
so. A mode-averaged capillary never comes close — the gate's gradient and taper cases
tabulate to about 1 MB — so in practice `:auto` is `:tabulated` for everything
`prop_capillary` builds. It exists for a multimode, radial or free-space operator, which is
the size of the whole state, where the table is `2·nnodes` copies of it; if you see that
message, raising `linop_tol` so that the table fits is usually a better answer than the
fallback, because the fallback's ninety host evaluations of a state-sized operator per step
cost more than the memory did.

**Your own `Φ`.** If you know the integral in closed form, define a subtype of
`LinearOps.AbstractIntegratedLinop` and pass it to `Luna.run` in place of the callable; it
is used as it is, and neither `β` nor `Aeff` is tabulated for it. The interface is three
functions:

```julia
import Luna: LinearOps, scalar

struct MyLinop{aT} <: LinearOps.AbstractIntegratedLinop
    L0::aT   # L(z) = L0 + L1*z, so Φ(z) = L0*z + L1*z^2/2
    L1::aT
end

function LinearOps.phase!(out, op::MyLinop, z)      # Φ(z)
    a = scalar(out, z); b = scalar(out, z^2/2)
    L0, L1 = op.L0, op.L1
    @. out = L0*a + L1*b
end

function LinearOps.derivative!(out, op::MyLinop, z) # L(z) itself, for the diagnostics
    a = scalar(out, z)
    L0, L1 = op.L0, op.L1
    @. out = L0 + L1*a
end
```

`phase!` and `derivative!` must be broadcasts (or other kernels) over `out`'s array type,
and every scalar must go through `Luna.scalar`, which converts it to the state's real
element type; that is all it takes for the same type to run on a GPU. An operator which
can only give *differences* of `Φ` — a quadrature — declares
`LinearOps.PhaseStyle(::MyLinop) = LinearOps.IncrementalPhase()` and implements
`phasediff!(out, op, z1, z2)` instead of `phase!`. See the developer guide.

Two limits on a caller-supplied operator. An absorbing boundary adds a constant to the
operator, which Luna folds into the integral exactly, so `boundary=:rate` works with one;
but the free-space evanescent clamp is not linear in the operator, cannot be pushed through
an integral, and raises rather than propagating something wrong — for free space, pass the
callable and let `linop_integral` do the integrating.

### `stats_period`

```julia
out = prop_capillary(...; stats_period=10)
```

Collects the default statistics every 10th accepted step instead of every step
(`Output.PeriodicStats`). The recorded statistics arrays are correspondingly shorter; the
saved field (`saveN`, `out["Eω"]`) is unaffected.

This is still the lever it always was on a GPU. On the mode-averaged geometry the
statistics are computed from a host copy of the field (above), which costs a transfer and
the host reductions on every step they fire on -- about 0.3 ms against a 2.0 ms step for a
1025-point grid, i.e. 13% of the step, and more with plasma. Raising `stats_period`
removes that in proportion. It is worth raising further when you have added a statistics
function of your own, and further still on a multi-column state, where the statistics run
on the device but the field is large.

## Performance

A GPU wins on many columns and large time grids; a mode-averaged single-column run is
launch-bound and will not speed up. `benchmark/device.jl` times the mode-averaged
propagation on each device and precision, sweeping the grid size; `benchmark/radial.jl`
does the same for the radial one, sweeping the number of radial points;
`benchmark/free.jl` for the 3-D Cartesian one, sweeping the transverse grid (the two
tables under "What runs where" are their output); and `benchmark/modal.jl` for a
four-mode propagation on each transverse integral. Run any of them from an environment
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
