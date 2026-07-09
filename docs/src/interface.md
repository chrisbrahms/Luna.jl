# The simple interface
Luna's simple interface provides two high-level entry points that build and run a complete simulation from keyword arguments, returning an output object (a [`Output.MemoryOutput`](@ref) by default) that can be indexed like a dictionary and passed to the [Processing.jl](@ref) and [Plotting.jl](@ref) functions:

- [`prop_capillary`](@ref Interface.prop_capillary) for pulse propagation in a gas-filled hollow capillary fibre or hollow-core photonic-crystal fibre, and
- [`prop_gnlse`](@ref Interface.prop_gnlse) for propagation governed by the generalised nonlinear Schrödinger equation (e.g. in solid-core fibres).

Both functions assemble the grid, the waveguide modes, the gas fill, the input pulse, the nonlinear response terms and the output for you, so that a simulation can be set up in a single call. If you need full control over these individual pieces — custom waveguides, gas mixtures, or free-space geometries — use the low-level interface instead (see the `examples/` directory).

## Choosing the grid limits and number of modes
The parameters that most strongly affect both the accuracy and the cost of a simulation are the wavelength limits `λlims`, the time window `trange`, and the `modes`.

**Wavelength limits (`λlims`).** `λlims` sets the extent of the frequency grid, and must span every wavelength you expect to be present or generated during the propagation — not only the pump, but also dispersive waves, harmonics, Raman lines and so on — with some margin, because the edges of the window are smoothly apodised. For field-resolved simulations the *short*-wavelength limit is the expensive one: it fixes the highest frequency the grid has to resolve and hence the time step, so extending `λlims` further into the ultraviolet rapidly increases the number of samples. Choose the smallest range that comfortably contains the dynamics of interest.

**Time window (`trange`).** `trange` is the total width of the time grid. It sets the spectral resolution (which is finer for larger `trange`) and, crucially, must be wide enough to contain the pulse together with all of its temporal spreading — from dispersion, nonlinear reshaping and, in multi-mode simulations, inter-modal walk-off — over the *entire* length of the waveguide. If pulse energy reaches the edge of the window it is attenuated smoothly, so if you see unexplained structure at the temporal edges and/or unexplained energy loss, increase `trange`. A reasonable starting point is ~20x the initial pulse duration, but as the required grid size depends on the nonlinear optical dynamics, **always** check that you have captured everything within the time window. The number of time samples is roughly `trange` divided by the time step, so both a wider `trange` and a bluer `λlims` make the simulation larger and slower (the grid is rounded up so that the sample count is a power of two).

**Number of modes (`modes`).** By default [`prop_capillary`](@ref Interface.prop_capillary) uses mode-averaged propagation in the fundamental `:HE11` mode, which is fast and reasonably accurate at low to moderate intensity. Setting `modes` to a number `N` (the first `N` `HE₁ₘ` modes) or to a tuple of mode signifiers switches to full multi-mode propagation, which captures inter-modal coupling, self-focusing and the spatial dependence of photoionisation — but is substantially more expensive, because the cost grows steeply with the number of modes (each nonlinear step evaluates overlap integrals between them). Multi-mode propagation is generally far more accurate, especially for strong nonlinear reshaping of the driving pulse and when photionisation plays a role, and can additionally capture other effects such as imperfect in-coupling (see below), polarisation effects, and inter-modal four-wave mixing. In most circumstances, mode-averaged simulations should be at least spot-checked against a multi-mode equivalent to make sure that these effects do not make a significant difference. A good starting point is to use 4 modes (`modes=4`); increase the number as required to achieve convergence.

Note that `modes=:HE11` (mode-averaged) and `modes=1` (single-mode projection) are **not** equivalent except when only the Kerr effect is included — see the [Modal decompositions](@ref) page.

## Specifying the input pulse
For the common case of a single pulse in the fundamental mode, the pulse is defined directly through keyword arguments — `λ0`, `τfwhm` (or `τw`), `energy` (or `power`), `pulseshape` (`:gauss` or `:sech`), the spectral phase `ϕ`, and `polarisation`. Internally these are assembled into a [`Pulses.GaussPulse`](@ref) or [`Pulses.SechPulse`](@ref) launched into the lowest-order mode.

Note that the `polarisation` keyword argument also has an effect on the mode set used in the propagation when using the default (`modes=:HE11`):
- `polarisation=:linear` will run a *mode-averaged* simulation with linear polarisation
- `polarisation=:circular`, `:x`, `:y`, or an ellipticity value ``-1 \le \varepsilon \le 1`` is equivalent to also setting `modes=1`, i.e. switching from mode-averaged to spatially resolved but single-mode propagation. This is because the mode-averaged UPPE is derived for scalar fields (see [Single-mode guided](@ref)).

For anything more elaborate, pass a single pulse object — an [`Pulses.AbstractPulse`](@ref) — or a `Vector` of them through the `pulses` keyword argument. A pulse object carries its own `mode` and `polarisation`, so a `Vector` of pulses can launch several pulses at once in different modes, polarisations, wavelengths or delays. The available types cover an arbitrary temporal shape ([`Pulses.CustomPulse`](@ref)), Gaussian and sech² pulses ([`Pulses.GaussPulse`](@ref), [`Pulses.SechPulse`](@ref)), a measured spectrum ([`Pulses.DataPulse`](@ref)), a field taken from a previous Luna run ([`Pulses.LunaPulse`](@ref)), and a focused free-space beam ([`Pulses.GaussBeamPulse`](@ref)). When `pulses` is given, all pulse-shaping keyword arguments **except `λ0`** are ignored.

## Coupling a free-space beam with `GaussBeamPulse`
[`Pulses.GaussBeamPulse`](@ref) models the coupling of an ideal free-space Gaussian beam of ``1/e^2`` intensity radius `waist` into the waveguide. Rather than assuming that all of the energy enters the fundamental mode, it computes the overlap integral of the Gaussian beam with each waveguide mode and distributes the input energy across the modes with the correct relative amplitudes. The temporal shape is taken from an inner `timepulse` (any other pulse object, e.g. a [`Pulses.GaussPulse`](@ref)). This is the realistic way to account for imperfect input coupling — and the resulting higher-order-mode content — when a beam is focused into a capillary. The `Nmodes` argument limits how many modes are coupled to; the default `:all` uses every mode in the simulation, but restricting it can avoid numerical inaccuracies from weakly-overlapping high-order modes.

## Polarisation and degenerate modes
Each spatial mode of the waveguide can appear in the simulation either once — as a single, linearly polarised (``y``-axis) field — or twice, as a degenerate pair of orthogonal ``x`` and ``y`` polarisations. Which happens is decided automatically from two things:

- **The input polarisation.** With the default `polarisation=:linear` a single (``y``) polarisation is used. Requesting `:x`, `:y`, `:circular`, or an ellipticity value ``-1 \le \varepsilon \le 1`` requires both components, so every mode is then present twice.
- **The mode set.** Mode-averaged propagation, and mode sets containing only ``\mathrm{HE}_{1m}`` modes, are treated as linearly polarised. But if the modes include any non-``\mathrm{HE}_{1m}`` mode — a `TE`, `TM` or `EH` mode, or an `HE` mode with azimuthal index ``n > 1`` — both polarisation components are always required, because such modes inherently mix ``x`` and ``y``.

Whenever the vector treatment is triggered (by either condition), the number of propagated fields doubles and each mode appears twice in the output, once per polarisation.

## Envelope versus field-resolved simulations
By default `prop_capillary` solves the full **field-resolved** (carrier-resolved) equation on a [`Grid.RealGrid`](@ref). Setting `envelope=true` instead propagates the complex **envelope** on a [`Grid.EnvGrid`](@ref).

The field-resolved treatment resolves the optical carrier oscillation itself, so it captures all frequencies on the grid automatically — including third-harmonic generation, ``\chi^{(2)}`` mixing, self-steepening and optical shock formation — as well as photoionisation and plasma. The price is that the time step must resolve the highest frequency present (set by the short-wavelength end of `λlims` *multiplied by 3 to take into account THG*), which makes the grid finer. The envelope treatment does not resolve the carrier, so its grid can be much coarser and the simulation correspondingly cheaper — but it comes with restrictions: by default it **ignores third-harmonic generation** (it can be re-enabled, at the cost of a larger grid, via the `thg` option), and **plasma is not implemented** for envelope fields. Accordingly, the `thg` and `plasma` options default to *on* for field-resolved simulations and *off* for envelope ones. Use the envelope mode for narrowband, GNLSE-style problems where harmonics and ionisation are negligible, and the field-resolved mode for broadband or few-cycle dynamics, harmonic generation, and strong-field effects.

## Function reference
```@autodocs
Modules = [Interface]
```
## Pulse types for input to `prop_capillary`
```@autodocs
Modules = [Pulses]
```