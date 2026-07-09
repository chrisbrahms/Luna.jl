# Nonlinear responses

This page describes the nonlinear polarisation terms ``P_\mathrm{nl}`` that appear on the right-hand side of the UPPE (see [The numerical model](@ref "The numerical model")). Regardless of the [modal decomposition](@ref "Modal decompositions") in use, these terms are always evaluated in the real-space–time domain: the field is transformed to time (and, for spatially resolved geometries, to real space), each response is evaluated pointwise, and the result is transformed back. A simulation can combine several of these effects at once — they are supplied as a tuple of response functions (assembled from the keyword arguments by `Interface.makeresponse` for the simple interface) and simply summed to give the total ``P_\mathrm{nl}``.

Each response is written to be independent of the propagation geometry, so the same functions are used for mode-averaged, multi-mode, radial and full 3D simulations.

## Kerr effect
The Kerr effect is the third-order (``\chi^{(3)}``) nonlinearity responsible for self- and cross-phase modulation, four-wave mixing and third-harmonic generation. In Luna the susceptibility is built at runtime from the number density and a single-molecule hyperpolarisability, ``\chi^{(3)} = \rho\,\gamma_3``, where ``\gamma_3`` is provided by [`PhysData.γ3_gas`](@ref) and ``\rho`` is the gas density at the current position (allowing the nonlinearity to follow a pressure gradient). The macroscopic susceptibility can equivalently be obtained from the nonlinear index ``n_2`` via [`PhysData.χ3`](@ref), ``\chi^{(3)} = \tfrac{4}{3} n_2 \varepsilon_0 c\, n_0^2``.

Which form of the response is used depends on whether the simulation is field-resolved or envelope, and whether it is scalar (a single polarisation state) or vector (full polarisation resolution).

For a **field-resolved, vector** simulation the isotropic ``(\mathbf{E}\cdot\mathbf{E})\mathbf{E}`` contraction is used (implemented in the vector branch of [`Nonlinear.Kerr_field`](@ref)):
```math
\begin{align*}
\mathbf{P}_\mathrm{Kerr}(t, x, y, z) &= \varepsilon_0 \chi^{(3)} \left(\mathbf{E}(t, x, y, z) \cdot \mathbf{E}(t, x, y, z)\right) \mathbf{E}(t, x, y, z) \\[1em]
&= \varepsilon_0 \chi^{(3)} \begin{pmatrix}
           \left(E_x^2 + E_y^2\right) E_x \\[0.5em]
           \left(E_x^2 + E_y^2\right) E_y
\end{pmatrix}
\end{align*}
```

For a **field-resolved, scalar** simulation this reduces to the full cube of the real field:
```math
P_\mathrm{Kerr}(t, x, y, z) = \varepsilon_0 \chi^{(3)} \left(E(t, x, y, z)\right) ^3
```
Because ``E`` is the real, carrier-resolved field, this expression automatically contains both the self-phase-modulation term and the third-harmonic-generation (THG) term. If THG is not wanted (for example to compare against envelope models, or to avoid having to resolve it on the frequency grid), [`Nonlinear.Kerr_field_nothg`](@ref) instead evaluates only the ``\tfrac{3}{4}\varepsilon_0\chi^{(3)}\vert\mathcal{H}(E)\vert^2 E`` self-phase-modulation part, using a Hilbert transform ``\mathcal{H}`` to construct the analytic signal.

For an **envelope, scalar** simulation the rapidly oscillating THG term is dropped and only the self-phase-modulation term is kept, with the familiar factor of ``\tfrac{3}{4}`` (implemented in [`Nonlinear.Kerr_env`](@ref)):
```math
P_\mathrm{Kerr}(t, x, y, z) = \frac{3}{4}\varepsilon_0 \chi^{(3)} \left\vert\mathcal{E}(t, x, y, z)\right\vert^2\mathcal{E}(t, x, y, z)
```
If the third harmonic *is* required in an envelope simulation, [`Nonlinear.Kerr_env_thg`](@ref) restores it as an explicitly up-shifted term ``\propto \mathrm{e}^{2i\omega_0 t}\mathcal{E}^2\mathcal{E}`` (see Eq. 4 of Genty et al., *Opt. Express* **15**, 5382 (2007)).

For an **envelope, vector** simulation the response additionally carries the cross-polarisation coupling terms, with the characteristic ``\tfrac{2}{3}`` cross-phase and ``\tfrac{1}{3}`` phase-conjugate (coherent coupling) coefficients that arise from separating the co- and counter-rotating circular components:
```math
\mathbf{P}_\mathrm{Kerr}(t, x, y, z) = \frac{3}{4}\varepsilon_0 \chi^{(3)} \begin{pmatrix}
           (\mathcal{E}_x^2 + \frac{2}{3}\mathcal{E}_y^2)\mathcal{E}_x + \frac{1}{3}\mathcal{E}_x^*\mathcal{E}_y^2 \\[0.5em]
           (\mathcal{E}_y^2 + \frac{2}{3}\mathcal{E}_x^2)\mathcal{E}_y + \frac{1}{3}\mathcal{E}_y^*\mathcal{E}_x^2
\end{pmatrix}
```

## Second-order nonlinearity
The second-order (``\chi^{(2)}``) nonlinearity gives rise to second-harmonic generation, sum- and difference-frequency generation and optical rectification. For a scalar field it is
```math
P_2(t, x, y, z) = \varepsilon_0 \chi^{(2)} \mathbf{E}^2(t, x, y, z)
```
and more generally, with the second-order susceptibility tensor,
```math
P_{i} = \varepsilon_0 \chi^{(2)}_{ijk} E_j E_k
```

!!! note "Not yet implemented"
    The ``\chi^{(2)}`` response is documented here for completeness but is **not currently part of Luna's response set** — there is no corresponding response function in `Nonlinear.jl`. This section describes the intended model only.


## Photoionisation & plasma
When the field is intense enough to ionise the medium, the freed electrons contribute to the polarisation both through the energy lost to ionisation and through the current they carry. Luna models this with [`Nonlinear.PlasmaCumtrapz`](@ref) (adapted from Geissler et al., *Phys. Rev. Lett.* **83**, 2930 (1999)):
```math
P_\mathrm{ion}\!\left(t,x, y,z\right) = I_p\int_{-\infty}^t \!\!\mathrm{d}t'\frac{\partial_{t'} \rho_\mathrm{e}(t', x, y,z)}{E\!\left(t',x, y,z\right)} + \frac{e^2}{m_\mathrm{e}}\int_{-\infty}^{t}\!\!\mathrm{d}t'\int_{-\infty}^{t'} \!\!\mathrm{d}t'' \rho_\mathrm{e}(t'',x, y,z)E\left(t'',x, y,z\right)
```
The first term is the **ionisation loss**: it removes energy proportional to the ionisation potential ``I_p`` each time an electron is freed. The second term is the **Drude plasma current**, describing the free electrons accelerating in the field. The electron density is obtained by integrating the ionisation rate ``w`` over the field history:
```math
\rho_\mathrm{e}(t, x, y,z) = \rho(z)\left(1 - \exp \left\{-\int_{-\infty}^t\!\!\mathrm{d}t'w\!\left(\left\vert E\!\left(t',x, y,z\right)\right\vert\right)\right\}\right)
```
where ``\rho(z)`` is the neutral density. Numerically, the two nested time integrals in the equation above are evaluated as successive cumulative-trapezoidal integrations of the plasma current, which is why the response is named `PlasmaCumtrapz`. For a vector field the ionisation rate is computed from the field magnitude ``\vert\mathbf{E}\vert`` and the resulting polarisation is applied component-wise.

The ionisation rate ``w(\vert E\vert)`` is supplied by one of the rate models in [Ionisation.jl](@ref), for example the ADK rate [`Ionisation.IonRateADK`](@ref) or the PPT rate [`Ionisation.IonRatePPT`](@ref) (and its cached, spline-interpolated form [`Ionisation.IonRatePPTCached`](@ref), which is much faster for repeated runs). The same ``1 - \exp(-\int w\,\mathrm{d}t')`` electron fraction is available directly via [`Ionisation.ionfrac`](@ref).

## Raman response
Raman scattering from the vibration and rotation of molecules produces a *delayed* (non-instantaneous) nonlinear response: the polarisation at time ``t`` depends on the field at earlier times through a convolution with a response function ``h``. Luna models this with the non-rigid-rotor / single-vibrational-transition treatment of [Wahlstrand et al., *Phys. Rev. A* **92**, 063828 (2015)] as applied to gas-filled hollow-core fibres by [Gao et al., *Laser & Photonics Reviews* **16**, 2100426 (2022)]. In general the Raman polarisation has the form
```math
P^\mathrm{R} = E(t)\,N\,(4\pi\epsilon_0)^2\sum_i \kappa_i\int_{0}^{\infty}h(\nu_i, T_2^i, \tau)\, E(t-\tau)^2\,\mathrm{d}\tau\,,
```
where ``N`` is the number density and the sum runs over the contributing Raman transitions ``i`` (a single vibrational transition and/or many rotational transitions, see below). Each transition contributes a **damped-oscillator response function**
```math
h(\nu, T_2, t) = \sin (2\pi \nu t)\exp(-t/T_2)\,,\qquad T_2=(\pi \Delta\nu)^{-1}\,,
```
with ``\nu`` the transition frequency, ``T_2`` the dephasing (coherence) time and ``\Delta\nu`` the full-width at half-maximum linewidth of the transition. This single-oscillator response is implemented by [`Raman.RamanRespSingleDampedOscillator`](@ref); the convolution with ``E(t-\tau)^2`` is carried out (via FFT, on a doubled time grid to avoid truncation of the long-lived response) when the response is applied to the field by [`Nonlinear.RamanPolarField`](@ref) (field-resolved) or [`Nonlinear.RamanPolarEnv`](@ref) (envelope). As with the Kerr effect, the field-resolved response can optionally include or exclude the third-harmonic component of the ``E^2`` driving term.

### Vibrational Raman
Because the vibrational frequency ``\nu_v`` is usually much larger than the rotational frequencies — large enough that only the ground vibrational state is thermally populated at room temperature — scattering from molecular stretching is captured by a *single* vibrational transition:
```math
P^\mathrm{v} = E(t)\,N\,(4\pi\epsilon_0)^2\,\kappa_v\int_{0}^{\infty}h(\nu_v, T_2^v, \tau)\, E(t-\tau)^2\,\mathrm{d}\tau\,,
```
```math
\kappa_v = -\left(\frac{\partial \alpha }{\partial Q}\right)^2\frac{1}{8\pi \mu \nu_v}\,,
```
where ``\partial\alpha/\partial Q`` is the isotropic polarisability derivative with respect to the molecular stretch coordinate ``Q`` and ``\mu`` is the reduced molecular mass. This is implemented by [`Raman.RamanRespVibrational`](@ref).

### Rotational Raman
The rotational levels are closely spaced, so a whole manifold of levels is thermally populated and contributes to the response. The energy ``e_J`` and fractional thermal population ``\rho_J`` of the level with rotational quantum number ``J`` are
```math
e_J = hc\left[BJ(J+1)-DJ^2(J+1)^2\right]\,,
```
```math
\rho_J = \frac{q_J(2J+1)\exp (-e_J/k_B T)}{\sum_K q_K(2K+1)\exp (-e_K/k_B T)}\,,
```
where ``B`` and ``D`` are the rotational and centrifugal constants, ``T`` is the temperature and ``q_J`` is a nuclear-spin statistical weight that alternates between odd and even ``J``. Selection rules restrict rotational Raman transitions to pairs of levels with ``\Delta J = 2``, giving line frequencies ``\nu_r^J = (e_{J+2}-e_J)/h``. The rotational polarisation is a sum over all allowed transitions,
```math
P^\mathrm{r} = E(t)\,N\,(4\pi\epsilon_0)^2\sum_{J}\kappa_r^J\int_{0}^{\infty}h(\nu_r^J, T_2^r, \tau)\,E(t-\tau)^2\,\mathrm{d}\tau\,,
```
```math
\kappa_r^J = \frac{4\pi\Delta\alpha^2}{15 h}\frac{(J+1)(J+2)}{2J+3}\left(\frac{\rho_{J+2}}{2J+5}-\frac{\rho_J}{2J+1}\right)\,,
```
where ``\Delta\alpha`` is the molecular polarisability anisotropy. The per-line linewidth is pressure-broadened, ``\Delta\nu_r = (\pi T_2^r)^{-1} = \rho\, b_r``, with ``\rho`` the gas density (in amagat) and ``b_r`` an experimentally determined constant. This is implemented by [`Raman.RamanRespRotationalNonRigid`](@ref).

### Assembling the response
For a molecular gas, [`Raman.molecular_raman_response`](@ref) builds the rotational and/or vibrational parts and sums them into a [`Raman.CombinedRamanResponse`](@ref). The top-level [`Raman.raman_response`](@ref) dispatches on the material: molecular gases use the vibrational/rotational model above, while solids such as fused silica use the multi-mode broadened model [`Raman.RamanRespIntermediateBroadening`](@ref) (Hollenbeck & Cantrell, *J. Opt. Soc. Am. B* **19**, 2886 (2002)). In the simple interface, the fractional Raman contribution ``f_r`` to the third-order response of a glass is set by the `fr` keyword (default ``0.18``).

### Implementation
- [`Raman.raman_response`](@ref)
- [`Raman.molecular_raman_response`](@ref)
- [`Raman.RamanRespVibrational`](@ref)
- [`Raman.RamanRespRotationalNonRigid`](@ref)
- [`Raman.RamanRespSingleDampedOscillator`](@ref)
- [`Raman.RamanRespIntermediateBroadening`](@ref)
- [`Raman.CombinedRamanResponse`](@ref)
- [`Nonlinear.RamanPolarField`](@ref)
- [`Nonlinear.RamanPolarEnv`](@ref)

