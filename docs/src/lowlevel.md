# The low-level interface
The [simple interface](@ref "The simple interface") builds a complete simulation from keyword arguments, but it only covers the common cases. For anything it cannot express — custom fills, pressure gradients, unusual waveguides, or free-space geometries — you assemble the individual pieces yourself and call [`Luna.setup`](@ref) and [`Luna.run`](@ref) directly. This is the low-level interface.

A simulation is built from a handful of independent, composable pieces:

- a **grid** ([`Grid.RealGrid`](@ref) or [`Grid.EnvGrid`](@ref)) defining the time/frequency axes,
- a **transverse basis** — a tuple of waveguide modes, or a spatial transform for free space,
- a **density function** `z -> number density`,
- a tuple of nonlinear **responses**,
- an **input** field,
- a **linear operator** (`linop`) carrying the dispersion and loss, and
- an **output**.

`Luna.setup` combines the first five into the initial frequency-domain field `Eω` and a `transform` that computes the nonlinear polarisation; `Luna.run` then integrates the propagation. The pieces must be chosen consistently — switching geometry or making the waveguide non-uniform requires several coordinated changes, which is exactly what this page illustrates.

This page works through three example scripts of increasing generality. The plotting/post-processing at the end of each script is omitted here, as it works identically to the [simple interface](@ref "The simple interface").

## A modal simulation
We start from [`examples/low_level_interface/basic_modal.jl`](https://github.com/LupoLab/Luna.jl/blob/master/examples/low_level_interface/basic_modal.jl), which propagates a pulse through a gas-filled capillary in two modes.

First, the physical parameters — the capillary core radius `a`, the gas and its (constant) fill pressure, the waveguide length, and the input pulse duration, central wavelength and energy:

```julia
using Luna

a = 13e-6
gas = :Ar
pres = 5
flength = 15e-2

τfwhm = 30e-15
λ0 = 800e-9
energy = 1e-6
```

The **modes** are a *tuple* of [`Capillary.MarcatiliMode`](@ref)s — here the `HE₁₁` and `HE₁₂` modes of the capillary. This tuple *is* the transverse basis of the simulation, so including two modes is what makes this a multi-mode simulation. Each mode is built from the core radius, gas and pressure, plus its mode indices `n`/`m` and family `kind`:

```julia
modes = (
    Capillary.MarcatiliMode(a, gas, pres, n=1, m=1, kind=:HE, ϕ=0.0, loss=false),
    Capillary.MarcatiliMode(a, gas, pres, n=1, m=2, kind=:HE, ϕ=0.0, loss=false)
)
```

The **grid** is a [`Grid.RealGrid`](@ref), so this is a field-resolved (carrier-resolved) simulation. Its arguments are the propagation length, the reference wavelength, the wavelength limits `(λmin, λmax)` of the frequency window, and the time window `trange`. See [Choosing the grid limits and number of modes](@ref) for how to pick these:

```julia
grid = Grid.RealGrid(flength, λ0, (160e-9, 3000e-9), 1e-12)
```

[`Fields.energyfuncs`](@ref) returns a pair of closures that compute the pulse energy from a time-domain or frequency-domain field; they are used by the statistics and normalisation below:

```julia
energyfun, energyfunω = Fields.energyfuncs(grid)
```

The **density function** maps propagation distance `z` to number density. For a uniform fill it is a closure that ignores `z` and returns the constant density from [`PhysData.density`](@ref). (This is the line that changes for a pressure gradient — see the [next section](@ref "Adding a pressure gradient").)

```julia
densityfun = let dens0=PhysData.density(gas, pres)
    z -> dens0
end
```

The nonlinear physics is selected by the **responses** tuple. Here it contains the Kerr effect ([`Nonlinear.Kerr_field`](@ref), with the gas's `γ3` coefficient) and photoionisation-driven plasma ([`Nonlinear.PlasmaCumtrapz`](@ref)), the latter needing an ionisation potential and rate ([`Ionisation.IonRatePPTCached`](@ref)):

```julia
ionpot = PhysData.ionisation_potential(gas)
ionrate = Ionisation.IonRatePPTCached(gas, λ0)

responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
             Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot))
```

The **input** is a transform-limited Gaussian pulse, [`Fields.GaussField`](@ref):

```julia
inputs = Fields.GaussField(λ0=λ0, τfwhm=τfwhm, energy=energy)
```

Now [`Luna.setup`](@ref) assembles the pieces. It applies the input to build the initial frequency-domain field `Eω`, and returns it along with the `transform` (which computes the nonlinear polarisation each step) and the FFT plan `FT`. Passing the `modes` tuple and a polarisation axis (`:y`) is what dispatches `setup` to the *modal* transform; `full=false` selects the faster mode-projection treatment appropriate when the modes share an axis of symmetry:

```julia
Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs, modes, :y; full=false)
```

The **linear operator** carries the dispersion and loss of each mode. [`LinearOps.make_const_linop`](@ref) builds a *constant* operator — valid here because the waveguide is uniform along its length, so the dispersion does not depend on `z`. (This is the second line that changes for a gradient.)

```julia
linop = LinearOps.make_const_linop(grid, modes, λ0)
```

Finally, set up the diagnostics with [`Stats.default`](@ref), an in-memory [`Output.MemoryOutput`](@ref) collecting 201 snapshots along the fibre, and run the propagation with [`Luna.run`](@ref):

```julia
statsfun = Stats.default(grid, Eω, modes, linop, transform; gas=gas, windows=((150e-9, 300e-9),))
output = Output.MemoryOutput(0, grid.zmax, 201, statsfun)

Luna.run(Eω, grid, linop, transform, FT, output, status_period=5)
```

After the run, `output` can be indexed like a dictionary (`output["Eω"]`, `output["z"]`) and passed to the [Processing.jl](@ref) and [Plotting.jl](@ref) functions exactly as for the simple interface.

## Adding a pressure gradient
A common variation is a **pressure gradient** — a fill that varies along the fibre, for example high pressure at the entrance falling to vacuum at the exit. This is *not* a one-line change: because the gas density now depends on `z`, both the density *and* the mode dispersion become `z`-dependent, and the linear operator can no longer be constant. Adapting the script above (following [`examples/low_level_interface/gradients/gradient_modal.jl`](https://github.com/LupoLab/Luna.jl/blob/master/examples/low_level_interface/gradients/gradient_modal.jl)) takes **three coordinated changes**.

**1. Build the density and core-index profiles with `Capillary.gradient`.** Instead of the constant-density closure, `Capillary.gradient` returns *both* a `z`-dependent core refractive-index function `coren` *and* the density function `densityfun` in one call. Here the pressure runs from `pres` at `z=0` to `0` at `z=L` (the function builds the physically-correct ``\sqrt{\,}`` pressure profile between the two endpoints):

```julia
L = 15e-2

coren, densityfun = Capillary.gradient(gas, L, pres, 0)
```

This single line replaces the entire constant-density `let` block from the uniform script.

**2. Build the modes from `coren`, not from `(gas, pres)`.** Because the gas density — and hence the refractive index the mode sees — now varies with `z`, each mode must be constructed from the `coren` function rather than a fixed gas and pressure. The mode's effective index then varies along the fibre:

```julia
modes = (
    Capillary.MarcatiliMode(a, coren, n=1, m=1, kind=:HE, ϕ=0.0, loss=false),
    Capillary.MarcatiliMode(a, coren, n=1, m=2, kind=:HE, ϕ=0.0, loss=false)
)
```

**3. Use a variable linear operator.** Since the dispersion now changes with `z`, `LinearOps.make_linop` (which recomputes the operator at each step) replaces `make_const_linop`:

```julia
linop = LinearOps.make_linop(grid, modes, λ0)
```

Everything else — the grid, energy functions, responses, input, `Luna.setup` call, output and `Luna.run` — is unchanged from the uniform case.

!!! warning "Change all three together"
    These three edits are interdependent. Making only one of them — for example switching to `make_linop` while still building the modes from a fixed pressure — will *not* raise an error; it will simply run a physically wrong (effectively uniform) simulation. When introducing a gradient, always update the density/index profile, the modes, *and* the linear operator together.

## A free-space simulation
For propagation in free space (or in a large-core waveguide treated as free space) there are no discrete waveguide modes. Instead the transverse coordinate is handled by a spatial transform, and the shape of the [`Luna.setup`](@ref) call changes accordingly. The following follows [`examples/low_level_interface/freespace/radial.jl`](https://github.com/LupoLab/Luna.jl/blob/master/examples/low_level_interface/freespace/radial.jl), which propagates a focusing Gaussian beam with radial symmetry.

Beyond `using Luna`, this example needs a few extra imports — the `Hankel` transform module (part of Luna) for the transverse basis, and `FFTW`/`NumericalIntegration` for the post-processing:

```julia
using Luna
import Luna.PhysData: wlfreq
import FFTW
import Luna: Hankel
import NumericalIntegration: integrate, SimpsonEven
```

The **grid** is a [`Grid.RealGrid`](@ref) as before. The transverse basis is now a *quasi-discrete Hankel transform* `q = Hankel.QDHT(R, N)` — a radial grid of `N` points out to radius `R` — which **replaces the mode tuple**:

```julia
grid = Grid.RealGrid(L, 800e-9, (400e-9, 2000e-9), 0.2e-12)
q = Hankel.QDHT(R, N, dim=3)
```

Because the spatial basis is now `q`, the energy functions take it as a second argument:

```julia
energyfun, energyfun_ω = Fields.energyfuncs(grid, q)
```

The density function and responses are as before (a constant fill and, here, the Kerr effect only):

```julia
densityfun = let dens0=PhysData.density(gas, pres)
    z -> dens0
end

responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
```

The **linear operator** for free space is built over `q` and takes a *bulk refractive-index function* (`PhysData.ref_index_fun(gas, pres)`) rather than a mode set and reference wavelength:

```julia
linop = LinearOps.make_const_linop(grid, q, PhysData.ref_index_fun(gas, pres))
```

Free-space geometries also require the nonlinear **normalisation** to be supplied explicitly — in the modal case `Luna.setup` constructed it internally, but here `NonlinearRHS.const_norm_radial` must be built and passed in:

```julia
normfun = NonlinearRHS.const_norm_radial(grid, q, PhysData.ref_index_fun(gas, pres))
```

The **input** is a *spatiotemporal* field, [`Fields.GaussGaussField`](@ref) — Gaussian in both time and space. In addition to the temporal parameters it carries a beam waist `w0` and a launch position `propz` (here the beam starts converging 0.3 m before its focus):

```julia
inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ, energy=energy, w0=w0, propz=-0.3)
```

The **setup** call has a different signature: `q` and `normfun` appear, while the `modes` tuple, the polarisation symbol and `full` are gone. It is this combination of arguments that dispatches `setup` to the radial free-space transform:

```julia
Eω, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)
```

Output and running are as before:

```julia
output = Output.MemoryOutput(0, grid.zmax, 51)
Luna.run(Eω, grid, linop, transform, FT, output)
```

One difference to note downstream: the output field now carries a radial dimension, with `output["Eω"]` shaped `(Nω, Npol, Nr, Nz)`. The stored field is in the Hankel-transformed (spatial-frequency) domain; apply the inverse transform `q \ Eω` to recover the field as a function of radius.

Other free-space geometries follow the same pattern with a different transverse basis — a 2D planar grid, or a full 3D `Grid.FreeGrid` — see the [Modal decompositions](@ref) page for the underlying theory and the other scripts in [`examples/low_level_interface/freespace/`](https://github.com/LupoLab/Luna.jl/tree/master/examples/low_level_interface/freespace).
