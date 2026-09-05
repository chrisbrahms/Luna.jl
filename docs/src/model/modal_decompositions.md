# Modal decompositions
!!! note
    All Fourier transforms in this section are written in the optics sign convention. See [A note on sign conventions](@ref) for details on what is used in the code.

## Multi-mode guided
For propagation in waveguides taking into account multiple modes and the coupling between them, Luna uses the model laid out in [Kolesik and Moloney, *Nonlinear optical pulse propagation simulation: From Maxwell’s to unidirectional equations*](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.70.036604) and [Tani et al., *Multimode ultrafast nonlinear optics in optical waveguides: numerical modeling and experiments in kagomé photonic-crystal fiber*](http://josab.osa.org/abstract.cfm?URI=josab-31-2-311). This is implemented in [`NonlinearRHS.TransModal`](@ref). The electric field ``\mathbf{E}(t, \mathbf{r_\perp}, z)`` is expressed as the inverse Fourier transform in time and the superposition of waveguide modes in space. This means that the transverse wave vector ``\mathbf{k}_\perp`` turns into a modal index ``j`` (this transform is implemented in [`Modes.ToSpace`](@ref) and [`Modes.to_space!`](@ref)):
```math
\mathbf{E}(t, \mathbf{r_\perp}, z) = \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{d} \omega \sum_j \hat{\mathbf{e}}_j(\mathbf{r_\perp}, z) \tilde{A}_j(\omega, z) \mathrm{e}^{-i \omega t}\,,
```
where ``\hat{\mathbf{e}}_j(\mathbf{r_\perp}, z)`` is the orthonormal transverse field distribution of the ``j^{\mathrm{th}}`` mode and ``\tilde{A}_j(\omega, z)`` is the frequency-domain amplitude in mode ``j``. The mode fields ``\hat{\mathbf{e}}_j(\mathbf{r_\perp}, z)`` are taken to be independent of frequency but can depend on the propagation coordinate ``z`` (e.g. in tapered waveguides). They are vector quantities in general, since the polarisation direction of most waveguide modes varies across the mode profile (e.g. for ``\mathrm{TE}``, ``\mathrm{TM}`` and higher-order ``\mathrm{HE}_{nm}`` modes with ``n > 1``); only for modes with a uniform linear polarisation (such as the ``\mathrm{HE}_{1m}`` modes of a capillary) can they be reduced to a single scalar component, which Luna does where possible (see [Polarisation and degenerate modes](@ref)). The modes are normalised such that ``\vert \tilde{A}_j(\omega, z) \vert^2`` gives the spectral energy density in mode ``j`` (when also taking into account the normalisation of the FFT), and equivalently ``\vert A_j(t, z)\vert^2`` gives the instantaneous power. The forward transform to reciprocal space is simply the overlap integral of the total field with each mode combined with the Fourier transform in time:
```math
\tilde{A}_j(\omega, z) = \int_S \mathrm{d}^2\mathbf{r_\perp} \int_{-\infty}^\infty \mathrm{d} t\,\, \hat{\mathbf{e}}_j^*(\mathbf{r_\perp}, z) \cdot \mathbf{E}(t, \mathbf{r_\perp}, z) \mathrm{e}^{i \omega t}\,,
```
where ``S`` is the cross-sectional area of the waveguide. This transform is implemented in [`NonlinearRHS.TransModal`](@ref) for use within simulations and in [`Modes.overlap`](@ref) for decomposition of existing sampled fields. In both cases, the mode overlap integral is solved explicitly with a p-adaptive or h-adaptive cubature method.

The linear operator for a mode ``\mathcal{L}_j(\omega, z)`` is given by (see [`LinearOps.make_const_linop`](@ref))
```math
\mathcal{L}_j(\omega, z) = i\left(\beta_j(\omega, z) - \frac{\omega}{v}\right) - \frac{1}{2}\alpha_j(\omega, z)\,,
```
where ``\beta_j(\omega, z)`` is real-valued and describes the phase evolution of the mode, ``v`` is a chosen frame velocity (this is the same for all modes) and ``\alpha(\omega, z)`` (also real) describes the attenuation of the waveguide (i.e. ``1/\alpha`` is the ``1/\mathrm{e}`` power/energy loss length). This can also be expressed in terms of the *effective index* of the mode:
```math
\mathcal{L}_j(\omega, z) = i \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\,,
```
where ``c`` is the speed of light in vacuum and ``n_\mathrm{eff}`` is complex, ``n_\mathrm{eff} = n + i k``, with ``n`` describing the effective refractive index and ``k`` describing the attenuation.

More generally, the phase subtracted from ``\beta_j`` to move into the reference frame of the simulation is ``\varphi(\omega) = \beta_1(\omega - \omega_0) + \beta_0``, where ``\beta_1 = 1/v`` is the inverse group velocity of a reference mode at the reference wavelength. For field-resolved simulations (`Grid.RealGrid`), ``\omega_0 = \beta_0 = 0`` and this reduces to the ``\omega/v`` above. For envelope simulations (`Grid.EnvGrid`), by default (`thg=false`) ``\omega_0`` is the carrier frequency `grid.ω0` and ``\beta_0`` is the propagation constant of the reference mode at ``\omega_0``, so that the envelope is phase-stationary at the carrier. With `thg=true`, ``\omega_0 = \beta_0 = 0`` instead, so that only the group delay ``\beta_1\omega`` (a pure time shift) is subtracted; this is required for any nonlinear response which mixes the carrier, see [Reference frames and carrier-mixing nonlinearities](@ref).

With the modal power normalisation for ``\hat{\mathbf{e}}_j(\mathbf{r_\perp}, z)``, the normalisation factor ``N_{\mathrm{nl}}`` comes out as simply ``N_{\mathrm{nl}}=4``. The propagation equation, coupling the modes through the nonlinear polarisation, is therefore
```math
\partial_z \tilde{A}_j(\omega, z) = i \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\tilde{A}_j(\omega, z) + i\frac{\omega}{4} \tilde{\mathbf{P}}_\mathrm{nl}\,,
```
where ``\tilde{\mathbf{P}}_\mathrm{nl}`` is given by
```math
\tilde{\mathbf{P}}_\mathrm{nl} =  \int_S \mathrm{d}^2\mathbf{r_\perp} \int_{-\infty}^\infty \mathrm{d} t\,\, \hat{\mathbf{e}}_j^*(\mathbf{r_\perp}, z) \cdot \mathbf{P}_\mathrm{nl}\left[\mathbf{E}(t, \mathbf{r_\perp}, z)\right] \mathrm{e}^{i \omega t}
```
and ``\mathbf{E}(t, \mathbf{r_\perp}, z)`` is obtained from the set of ``\tilde{A}_j(\omega, z)`` as above.

The transverse coordinate ``\mathbf{r_\perp}`` for circular waveguides (e.g. hollow capillaries, optical fibres, and anti-resonant fibres) is in polar coordinates, ``\mathbf{r_\perp} = (r, \theta)``. For other waveguides (e.g. rectangular), it is Cartesian, ``\mathbf{r_\perp} = (x, y)``.

!!! note
    While ``\mathbf{r_\perp}`` can be given in either coordinate system, the **components** of the modal fields ``\hat{\mathbf{e}}_j(\mathbf{r_\perp}, z)`` are **always** given in Cartesian coordinates, i.e. the basis vectors for the polarisation of the field are always ``\mathbf{x}`` and ``\mathbf{y}``.

### Implementation
The modules and functions that define and implement this decomposition for different modes are
- [Modes.jl](@ref)
- [Capillary.jl](@ref)
- [RectModes.jl](@ref)
- [Antiresonant.jl](@ref)
- [`NonlinearRHS.TransModal`](@ref)
- [`NonlinearRHS.norm_modal`](@ref)
- [`LinearOps.make_const_linop`](@ref)
- [`LinearOps.make_linop`](@ref)


## Single-mode guided
In some situations, inter-mode coupling in a waveguide is negligible, so including several waveguide modes in the simulation unnecessarily slows down the computation. Simulating propagation in a single mode is trivially achieved by including only that single mode in both the forward and inverse transforms as defined above for [multi-mode propagation](@ref "Multi-mode guided"). For example, setting `modes=1` when calling `prop_capillary` achieves this and leads to a significant speed-up. However, in this simple implementation, the overlap integral between the nonlinear polarisation and the waveguide mode still needs to be calculated explicitly. We can make this unnecessary by making an assumption about the nonlinear polarisation.

If the nonlinear polarisation is *only due to third-order effects* like the Kerr effect or Raman scattering, we can express it as
```math
P_\mathrm{nl}\left(t, \mathbf{r}_\perp, z \right) = C\, E(t, \mathbf{r}_\perp, z)^3\,,
```
where ``C`` is a constant which depends on the specific effect (e.g. for the Kerr effect, ``C`` becomes ``\varepsilon_0 \chi^{(3)}`` with ``\chi^{(3)}`` the third-order susceptibility of the nonlinear medium) and we have switched to *explicitly real-valued* and *scalar* fields to make the notation simpler; the same result can be obtained with vector fields and more algebra. Expanding the field in terms of its modal content as above, this turns into
```math
P_\mathrm{nl}\left(t, \mathbf{r}_\perp, z \right) = C\, \Big[\sum_j \hat{e}_j(\mathbf{r_\perp}, z) A_j(t, z)\Big]^3\,,
```
where we have simply carried out the time-domain inverse Fourier transform to obtain ``A_j(t, z)``. For a single mode, this simplifies greatly to
```math
P_\mathrm{nl}\left(t, \mathbf{r}_\perp, z \right) = C\, \hat{e}_0(\mathbf{r_\perp}, z)^3 A(t, z)^3\,.
```
Now we can explicitly calculate the overlap integral with the single mode we are considering:
```math
\begin{align*}
P_\mathrm{nl}(t, z) &=  CA(t, z)^3\times\int_S \mathrm{d}^2\mathbf{r_\perp} \, \hat{e}_0^*(\mathbf{r_\perp}, z) \hat{e}_0(\mathbf{r_\perp}, z)^3\\[1em]
&= CA(t, z)^3\int_S \mathrm{d}^2\mathbf{r_\perp} \, \hat{e}_0(\mathbf{r_\perp}, z)^4\,,
\end{align*}
```
where ``A`` is now the modal amplitude in the single mode. In the second step we have made use of the fact that we are considering real-valued fields and hence ``\hat{e}_0(\mathbf{r_\perp}, z)`` is also real.

The mode normalisation in Luna is chosen such that the absolute value squared of the modal field amplitudes ``A_j(t, z)`` is the instantaneous power. For this to be fulfilled, we need
```math
\frac{1}{2} c \varepsilon_0 \int_S \mathrm{d}^2\mathbf{r_\perp} \left\vert \hat{e}_j(\mathbf{r_\perp}, z) \right\vert^2 = 1\,.
```
This, in turn, means that the *effective area* of the mode,
```math
A_{\mathrm{eff}, j}(z) = \frac{\left(\int_S \mathrm{d}^2\mathbf{r_\perp} \,\left\vert \hat{e}_j(\mathbf{r}_\perp, z)\right\vert^2\right)^2}{\int_S \mathrm{d}^2\mathbf{r_\perp}  \,\left\vert \hat{e}_j(\mathbf{r}_\perp, z)\right\vert^4}\,,
```
for the single mode we are considering is
```math
A_\mathrm{eff} = \Big(\frac{1}{4} c^2 \varepsilon_0^2 \int_S \mathrm{d}^2\mathbf{r_\perp} \, \hat{e}_0(\mathbf{r_\perp}, z)^4 \Big)^{-1}\,.
```
Note that ``A_\mathrm{eff}`` is **independent of the normalisation**, because the overall power of ``\hat{e}_j`` and any constants inside it is the same in the numerator and denominator. From this we can see that we can replace the integral expression in the projection of ``P_\mathrm{nl}`` with
```math
\int_S \mathrm{d}^2\mathbf{r_\perp} \, \hat{e}_0(\mathbf{r_\perp}, z)^4 = \frac{4}{\varepsilon_0^2c^2 A_\mathrm{eff}}\,.
```
We can now write down the **single-mode UPPE:**
```math

\partial_z \tilde{A}(\omega, z) = i \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\tilde{A}(\omega, z)
 + i\frac{\omega}{4} C\frac{4}{\varepsilon_0^2c^2 A_\mathrm{eff}} \int_{-\infty}^\infty \mathrm{d} t\, A(t, z)^3 \mathrm{e}^{i \omega t}\,.
```
Crucially, the effective area depends on the mode shape ``\hat{e}_0(\mathbf{r_\perp}, z)``, but only needs to be calculated *once* (assuming the cross-section of the waveguide does not change along its length). When only third-order nonlinear effects are present, the single-mode UPPE is exactly equivalent to explicitly solving the projection integral, but *much* faster. However, it has two important drawbacks:

1. For obvious reasons, inter-modal coupling mediated by the nonlinearity is completely ignored.
2. The equation only works for third-order nonlinear effects, and hence photoionisation and plasma dynamics cannot be modelled in this way.

We can derive a different single-mode equation which can treat other nonlinear effects approximately. As written above, in the single-mode UPPE only the modal amplitude ``A`` appears. We now define a re-scaled **mode-averaged field** ``E_\mathrm{av}`` through
```math
A(t, z) = \sqrt{\frac{1}{2}\varepsilon_0 c A_\mathrm{eff}}E_\mathrm{av}(t, z)\,.
```
Note that, because ``\vert A(t, z)\vert^2`` has units of power, ``E_\mathrm{av}`` is in fact an electric field (with units of ``\mathrm{V/m}``). We can also define the *mode-averaged intensity* by
```math
I_\mathrm{av} = \frac{1}{2}\varepsilon_0 c \vert E_\mathrm{av}\vert^2 = \frac{\vert A(t, z)\vert^2}{A_\mathrm{eff}}\,.
```
Plugging in the definition of ``E_\mathrm{av}``, the UPPE reads
```math
\begin{align*}
\sqrt{\frac{1}{2}\varepsilon_0 c A_\mathrm{eff}}\partial_z \tilde{E}_\mathrm{av}(\omega, z) &= i \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\sqrt{\frac{1}{2}\varepsilon_0 c A_\mathrm{eff}}\tilde{E}_\mathrm{av}(\omega, z)\\[1em]
&\qquad + i\frac{\omega}{4} C\Big(\frac{1}{2}\varepsilon_0 c A_\mathrm{eff}\Big)^{\frac{3}{2}}\frac{4}{\varepsilon_0^2c^2 A_\mathrm{eff}} \int_{-\infty}^\infty \mathrm{d} t\, E_\mathrm{av}(t, z)^3 \mathrm{e}^{i \omega t}\,.
\end{align*}
```
Cancelling the various constants, we arrive at the **mode-averaged field UPPE**
```math
\partial_z \tilde{E}_\mathrm{av}(\omega, z) = i \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\tilde{E}_\mathrm{av}(\omega, z) + i\frac{\omega}{4} \frac{2}{\varepsilon_0 c} \int_{-\infty}^\infty \mathrm{d} t\, P_\mathrm{nl}\left[E_\mathrm{av}(t, z)\right] \mathrm{e}^{i \omega t}\,.
```
This now includes only a single inverse Fourier transform to obtain ``E_\mathrm{av}(t, z)`` followed by the calculation of ``P_\mathrm{nl}`` and then a forward transform. This equation is **still only valid for third-order responses** but we have now written it for an arbitrary polarisation ``P_\mathrm{nl}\left[E_\mathrm{av}(t, z)\right]``. Because ``E_\mathrm{av}`` is (a version of) the actual electric field, we can calculate arbitrary polarisation contributions, including the photoionisation and plasma term. However, this is still a significant approximation. For example, the mode-averaged intensity is approximately half of the on-axis intensity for the fundamental mode of a capillary fibre. Due to the exponential scaling of strong-field ionisation with intensity, this means that the peak ionisation fraction can be underestimated significantly. Only fully mode-resolved (and multi-mode) propagation can accurately model that situation.

The mode-averaged field UPPE as written above is very useful, but the scaling from ``A`` to ``E_\mathrm{av}`` changes the normalisation: ``\vert E_\mathrm{av}(t, z) \vert^2`` no longer gives the instantaneous power. To remain consistent with modal propagation simulations (e.g. for data analysis), Luna internally uses the same normalisation for both, which leads to a "hybrid" equation. The propagating quantity (and hence the simulation output) is ``A(z, t)`` and we switch to ``E_\mathrm{av}(z, t)`` to calculate the nonlinear polarisation. This leads to the appearance of an additional factor of in the equation
```math
\begin{align*}
\left(\frac{1}{2}\varepsilon_0 c A_\mathrm{eff}\right)^{-\frac{1}{2}}\partial_z \tilde{A}(\omega, z) &= i\left(\frac{1}{2}\varepsilon_0 c A_\mathrm{eff}\right)^{-\frac{1}{2}} \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\tilde{A}(\omega, z) + i\frac{\omega}{4} \frac{2}{\varepsilon_0 c} \int_{-\infty}^\infty \mathrm{d} t\, P_\mathrm{nl}\left[E_\mathrm{av}(t, z)\right] \mathrm{e}^{i \omega t}\\[1em]

\Rightarrow \partial_z\tilde{A}(\omega, z) &= i \left(\frac{\omega}{c} n_\mathrm{eff}(\omega, z) - \frac{\omega}{v}\right)\tilde{A}(\omega, z) + i\frac{\omega}{4}\sqrt{\frac{2A_\mathrm{eff}}{\varepsilon_0 c} }\int_{-\infty}^\infty \mathrm{d} t\, P_\mathrm{nl}\left[E_\mathrm{av}(t, z)\right] \mathrm{e}^{i \omega t}\,.
\end{align*}
```



## Free-space propagation
Luna can also propagate pulses in free space (e.g. in a gas cell or a bulk nonlinear crystal) in three geometries, which differ in the transverse basis used to represent the field: radial symmetry, one transverse dimension (``x``-``z``), and full ``(3+1)``-dimensional propagation. Much of the machinery is shared between them, so it is described once here before the individual geometries.

### Polarisation
Unlike the modal decompositions above, in which the polarisation state is carried by the modes, the free-space transforms always represent the field with explicit Cartesian components. In all three geometries, the propagating array has the polarisation as its **second** dimension, i.e. it has the shape ``(N_\omega, N_\mathrm{pol}, N_\perp\ldots)``, where ``N_\mathrm{pol}`` is either 1 or 2. With ``N_\mathrm{pol}=2`` the two components are ``(E_x, E_y)``, and nonlinear responses which couple the polarisation components (e.g. the vector Kerr effect or the [second-order nonlinearity](@ref "Second-order nonlinearity")) act on both at once. With ``N_\mathrm{pol}=1`` the field is a scalar with a fixed linear polarisation, which is taken to be the ``y`` component of the input field.

Which of the two is used is decided by the refractive-index function which is passed to the linear operator and the normalisation function: if it returns a single refractive index, the simulation is scalar; if it returns two, ``(n_x, n_y)``, both polarisation components are propagated, each with its own dispersion (see below). The spatio-temporal input fields (e.g. [`Fields.GaussGaussField`](@ref)) always contain both components, with the linear polarisation angle set by their `θ` keyword argument (measured from the ``y`` axis, so ``\theta = 0`` is ``y``-polarised), and `Luna.setup` selects the appropriate component(s).

### Linear operator
In free space there are no modes, so the propagation constant of each plane-wave component follows directly from the dispersion relation. For polarisation component ``p`` with refractive index ``n_p(\omega, z)`` and transverse wave vector ``\mathbf{k}_\perp`` (see the individual geometries below for what ``\mathbf{k}_\perp`` is), the longitudinal wave vector is
```math
k_{z,p}(\omega, \mathbf{k}_\perp, z) = \sqrt{\left(\frac{\omega}{c}n_p(\omega, z)\right)^2 - \vert\mathbf{k}_\perp\vert^2}\,,
```
and the linear operator is
```math
\mathcal{L}_p(\omega, \mathbf{k}_\perp, z) = i\left[k_{z,p}(\omega, \mathbf{k}_\perp, z) - \varphi(\omega)\right]\,,
```
where ``\varphi(\omega) = \beta_1(\omega - \omega_0) + \beta_0`` is the reference-frame phase, with exactly the same conventions as for the guided cases above: ``\beta_1`` is the inverse group velocity at the reference wavelength `grid.referenceλ` (calculated from the *last* refractive index returned by the index function, i.e. ``n_y`` if there are two), and ``\omega_0`` and ``\beta_0`` are zero for field-resolved simulations and for envelope simulations with `thg=true`, and equal to `grid.ω0` and ``\omega_0 n(\omega_0)/c`` for envelope simulations with `thg=false` (the default). This operator includes diffraction, since ``k_z`` depends on ``\mathbf{k}_\perp``. For evanescent components, i.e. those for which ``k_{z,p}^2 < 0``, the phase evolution is dropped and instead the component is attenuated, ``\mathcal{L}_p = -i\varphi(\omega) - \min\left(\sqrt{\vert k_{z,p}^2\vert}, 200\,\mathrm{m}^{-1}\right)``, where the cap avoids the numerical problems described in [The interaction picture, a.k.a. pre-conditioned Runge-Kutta](@ref).

The constant (``z``-invariant) operator is created by [`LinearOps.make_const_linop`](@ref)`(grid, spacegrid, nfun[, thg])` from a function `nfun(λ)` returning the refractive index (or indices) as a function of wavelength. For ``z``-dependent operators (e.g. pressure gradients), [`LinearOps.make_linop`](@ref)`(grid, spacegrid, nfun[, thg])` takes `nfun(ω; z)` instead and returns a closure which fills the operator at each ``z``. In both cases the same function also works for [birefringent crystals](@ref "Birefringent crystals") if it returns two indices. Recall from [A note on sign conventions](@ref) that in the code, the fields use the mathematics sign convention, so the operators actually created carry the opposite sign.

### Normalisation
In all free-space geometries, the normalisation factor ``N_\mathrm{nl}`` is not constant but depends on the propagation angle of each plane-wave component. For polarisation ``p`` at frequency ``\omega`` and transverse wave vector ``\mathbf{k}_\perp`` it is
```math
N_\mathrm{nl}(\omega, \mathbf{k}_\perp, z) = \frac{2 k_{z,p}(\omega, \mathbf{k}_\perp, z)}{\mu_0 \omega}\,,
```
which turns the UPPE into the familiar form
```math
\partial_z E_p(\omega, \mathbf{k}_\perp, z) = i k_{z,p} E_p(\omega, \mathbf{k}_\perp, z) + \frac{i \mu_0 \omega^2}{2 k_{z,p}}P_{\mathrm{nl}, p}(\omega, \mathbf{k}_\perp, z)\,.
```
For paraxial components (``\vert\mathbf{k}_\perp\vert \ll k_z``) this reduces to ``N_\mathrm{nl} = 2n_p\varepsilon_0 c``, which is the same impedance factor as in the [mode-averaged field UPPE](@ref "Single-mode guided"), while for oblique components it applies the corresponding angular-spectrum correction. For evanescent components the normalisation factor is set to unity to avoid division by zero; these components are strongly attenuated by the linear operator in any case. The normalisation functions (`NonlinearRHS.norm_radial` etc., see the individual geometries) return ``k_{z,p}/(\mu_0\omega)``; the factor of 2 is applied inside the transform.

### Inputs and energy
The spatio-temporal input fields in [Fields.jl](@ref) (e.g. [`Fields.GaussGaussField`](@ref)) work for all three geometries and are created directly in real space on the transverse grid, normalised to the requested energy using [`Fields.energyfuncs`](@ref)`(grid, spacegrid)`, and then transformed to reciprocal space. The keyword argument `propz` linearly propagates the input in vacuum before the simulation starts, so that e.g. a focusing beam can be launched by placing the waist at `propz` before the start of the simulation. [`Fields.energyfuncs`](@ref) also returns the functions to calculate the energy of a propagated field in the time or frequency domain, including the normalisation of the FFT and of the transverse transform, and should be used to calculate energies from simulation output in the correct units.

## Radially symmetric free-space
For free-space propagation of a radially symmetric beam, the transverse basis is the continuum of Bessel beams ``J_0(k_\perp r)``, so the generalised transverse spatial frequency ``\mathbf{k}_\perp`` becomes the (scalar) radial spatial frequency ``k_\perp``. The field is written as
```math
\mathbf{E}(t, r, z) = \frac{2\pi}{\left(2\pi\right)^3} \int_{-\infty}^{\infty} \int_{0}^{\infty}  \mathbf{E}(\omega, k_\perp, z) J_0(k_\perp r) \mathrm{e}^{-i\omega t} k_\perp\mathrm{d} k_\perp  \mathrm{d}\omega\,,
```
where the ``k_\perp\,\mathrm{d}k_\perp`` measure and the leading ``2\pi`` (from the azimuthal integration) are the conventions of the Hankel transform. In Luna this transverse transform is a zeroth-order quasi-discrete Hankel transform (`Hankel.QDHT` from the [Hankel.jl](https://github.com/LupoLab/Hankel.jl) package), and the radial UPPE is implemented in [`NonlinearRHS.TransRadial`](@ref). On each step the nonlinear polarisation is evaluated by transforming the field from reciprocal space to real-space–time (an inverse Fourier transform ``\omega\to t`` and an inverse Hankel transform ``k_\perp\to r``), calculating ``\mathbf{P}_\mathrm{nl}`` pointwise, and transforming back (``r\to k_\perp`` and ``t\to\omega``), followed by the ``i\omega/N_\mathrm{nl}`` prefactor of the UPPE with ``\vert\mathbf{k}_\perp\vert^2 = k_\perp^2`` in the [normalisation](@ref "Normalisation") and [linear operator](@ref "Linear operator").

The propagating array is three-dimensional, ``(N_\omega, N_\mathrm{pol}, N_{k_\perp})``, and the output of a simulation saved for ``N_z`` steps has the shape ``(N_\omega, N_\mathrm{pol}, N_{k_\perp}, N_z)``. Because the polarisation axis comes second, the `Hankel.QDHT` **must** be created with `dim=3`, e.g. `Hankel.QDHT(R, N; dim=3)` for a spatial window of radius `R` with `N` radial points. The field in real space is recovered by the inverse Hankel transform, `q \ Eωk`. Note that the Hankel transform does not describe the field on axis (``r = 0``) directly; `Hankel.onaxis` calculates it from the reciprocal-space field.

Since the radial symmetry does not allow for any physical process which distinguishes between ``x`` and ``y`` (e.g. spatial walk-off), a two-polarisation simulation in this geometry treats both components as separate radially symmetric fields with their own refractive index which are coupled only through the nonlinear response.

### Implementation
- [`NonlinearRHS.TransRadial`](@ref)
- [`NonlinearRHS.norm_radial`](@ref)
- [`NonlinearRHS.const_norm_radial`](@ref)
- [`LinearOps.make_const_linop`](@ref)
- [`LinearOps.make_linop`](@ref)
- [`Fields.energyfuncs`](@ref)

## Two-dimensional free-space
For propagation with only one transverse dimension ``x``, i.e. a field which is invariant along ``y`` (or, more practically, a field for which the ``y`` dependence is not relevant), the transverse basis is the one-dimensional plane-wave spectrum ``\mathrm{e}^{ik_x x}`` and ``\mathbf{k}_\perp = k_x``. The field is
```math
\mathbf{E}(t, x, z) = \frac{1}{\left(2\pi\right)^2}\int_{-\infty}^{\infty} \int_{-\infty}^{\infty}  \mathbf{E}(\omega, k_x, z) \mathrm{e}^{i\left(k_x x - \omega t\right)}\mathrm{d} k_x  \mathrm{d} \omega\,.
```
The ``x`` grid and its spatial frequencies ``k_x`` are held by [`Grid.Free2DGrid`](@ref), and the transform is implemented in [`NonlinearRHS.TransFree2D`](@ref). A single combined two-dimensional FFT handles both the transverse spatial axis and time at once: the field is transformed ``(\omega, k_x)\to(t, x)``, the nonlinear polarisation is evaluated pointwise, and the result is transformed back before applying the ``i\omega/N_\mathrm{nl}`` prefactor with ``\vert\mathbf{k}_\perp\vert^2 = k_x^2`` (see [`NonlinearRHS.norm_free2D`](@ref) and [`NonlinearRHS.const_norm_free2D`](@ref)). The propagating array is three-dimensional, ``(N_\omega, N_\mathrm{pol}, N_{k_x})``, the output has the shape ``(N_\omega, N_\mathrm{pol}, N_{k_x}, N_z)``, and the field in real space is recovered by an inverse FFT along the third dimension, `FFTW.ifft(Eωk, 3)`.

This geometry is much cheaper than full 3D propagation and, unlike radial symmetry, allows for physical processes which break the symmetry between ``x`` and ``y``. Its main use is the simulation of [birefringent crystals](@ref "Birefringent crystals") with spatial walk-off, which occurs in the plane containing the optic axis.

!!! note
    Since the field is invariant along ``y``, "energy" in this geometry (as calculated by [`Fields.energyfuncs`](@ref) or requested from an input field) is actually the energy *per unit length along* ``y`` in ``\mathrm{J/m}``. To represent a beam with a given pulse energy and a Gaussian profile with ``1/\mathrm{e}^2`` radius ``w_0`` along ``y``, the energy has to be divided by ``\int\mathrm{e}^{-2y^2/w_0^2}\mathrm{d}y = \sqrt{\pi/2}\,w_0`` (see `examples/low_level_interface/freespace/free2D_bbo.jl`).

### Implementation
- [`Grid.Free2DGrid`](@ref)
- [`NonlinearRHS.TransFree2D`](@ref)
- [`NonlinearRHS.norm_free2D`](@ref)
- [`NonlinearRHS.const_norm_free2D`](@ref)
- [`LinearOps.make_const_linop`](@ref)
- [`LinearOps.make_linop`](@ref)

## Three-dimensional free-space
For full ``(3+1)``-dimensional free-space propagation without any symmetry, the transverse basis is the two-dimensional plane-wave spectrum ``\mathrm{e}^{i(k_x x + k_y y)}``, so ``\mathbf{k}_\perp = (k_x, k_y)``. The field is
```math
\mathbf{E}(t, x, y, z) = \frac{1}{\left(2\pi\right)^3}\int_{-\infty}^{\infty} \int_{-\infty}^{\infty}  \int_{-\infty}^{\infty}  \mathbf{E}(\omega, k_x, k_y, z) \mathrm{e}^{i\left(k_x x + k_y y - \omega t\right)}\mathrm{d} k_x \mathrm{d} k_y \mathrm{d} \omega\,,
```
with the standard three-fold inverse-Fourier-transform measure ``1/(2\pi)^3``. The Cartesian ``x``/``y`` grid and its spatial frequencies ``k_x``/``k_y`` are held by [`Grid.FreeGrid`](@ref), and the transform is implemented in [`NonlinearRHS.TransFree`](@ref). Here a single combined three-dimensional FFT handles both transverse spatial axes *and* time at once: the field is transformed ``(\omega, k_x, k_y)\to(t, x, y)``, the nonlinear polarisation is evaluated pointwise, and the result is transformed back before applying the ``i\omega/N_\mathrm{nl}`` prefactor with ``\vert\mathbf{k}_\perp\vert^2 = k_x^2 + k_y^2`` (see [`NonlinearRHS.norm_free`](@ref) and [`NonlinearRHS.const_norm_free`](@ref)). The propagating array is four-dimensional, ``(N_\omega, N_\mathrm{pol}, N_{k_x}, N_{k_y})``, the output has the shape ``(N_\omega, N_\mathrm{pol}, N_{k_x}, N_{k_y}, N_z)``, and the field in real space is recovered by an inverse FFT along the third and fourth dimensions, `FFTW.ifft(Eωk, (3, 4))`.

### Implementation
- [`Grid.FreeGrid`](@ref)
- [`NonlinearRHS.TransFree`](@ref)
- [`NonlinearRHS.norm_free`](@ref)
- [`NonlinearRHS.const_norm_free`](@ref)
- [`LinearOps.make_const_linop`](@ref)
- [`LinearOps.make_linop`](@ref)

## Birefringent crystals
The two-polarisation free-space geometries can describe linear propagation in a uniaxial birefringent crystal, e.g. for the simulation of [second-order nonlinear processes](@ref "Second-order nonlinearity") like second-harmonic generation. The crystal is cut such that the optic axis lies in the ``x``-``z`` plane at an angle ``\theta`` to the propagation direction ``z``. The ``y`` polarisation is then the ordinary wave with refractive index ``n_o(\omega)``, while the ``x`` polarisation is the extraordinary wave, whose refractive index depends on the angle ``\theta'`` between its wave vector and the optic axis:
```math
n_e(\omega, \theta') = \left[\frac{\cos^2\theta'}{n_o(\omega)^2} + \frac{\sin^2\theta'}{n_e(\omega)^2}\right]^{-1/2}\,,
```
where ``n_e(\omega)`` on the right-hand side is the principal extraordinary index (the index for ``\theta' = \pi/2``). [`PhysData.ref_index_fun_uniax`](@ref) creates this function, and [`PhysData.ref_index_fun_xy`](@ref)`(material, θ)` creates the pair of functions `(nfunx, nfuny)` for the two polarisations for a given cut angle ``\theta``. The Sellmeier expansions for the principal axes of the supported crystals are available through [`PhysData.ref_index`](@ref) with the `axis` keyword argument (`:o` and `:e` for uniaxial crystals, `:x`, `:y` and `:z` for LBO).

In the simplest approximation, the angle is the same for all plane-wave components, ``\theta' = \theta``. This is the only option available in the [radially symmetric geometry](@ref "Radially symmetric free-space"): passing an index function which returns ``(n_e(\omega, \theta), n_o(\omega))`` (see `examples/low_level_interface/freespace/radial_xypol_bbo.jl`) results in two radially symmetric fields with different dispersion, but neglects spatial walk-off.

In the [two-dimensional](@ref "Two-dimensional free-space") and [three-dimensional](@ref "Three-dimensional free-space") geometries, the dedicated methods of [`LinearOps.make_const_linop`](@ref), [`NonlinearRHS.const_norm_free2D`](@ref) and [`NonlinearRHS.const_norm_free`](@ref) which take a tuple `(nfunx, nfuny)` go one step further. Here, `nfunx(λ, δθ)` must return ``n_e(\omega, \theta + \delta\theta)``, where ``\delta\theta`` is the angle between the internal propagation direction of a plane-wave component and ``z``. For each ``\omega`` and ``k_x``, this angle is found from the conservation of the transverse wave vector,
```math
k_x = \frac{\omega}{c}n_e(\omega, \theta + \delta\theta)\sin\delta\theta\,,
```
which is solved numerically by [`PhysData.crystal_internal_angle`](@ref). The longitudinal wave vector for the ``x`` polarisation is then calculated from ``n_e(\omega, \theta+\delta\theta)`` and ``k_x`` (and ``k_y``) as in the [linear operator](@ref "Linear operator") above, and the ``y`` polarisation uses ``n_o(\omega)`` throughout. Because the refractive index of the extraordinary wave now depends on the direction of its wave vector, its wave-vector surface is tilted, which makes the direction of energy flow (the group velocity) differ from the direction of the wave vector: **spatial walk-off** is included automatically. The angle ``\delta\theta`` only depends on ``k_x``, since the optic axis lies in the ``x``-``z`` plane; ``k_y`` enters only through ``\vert\mathbf{k}_\perp\vert^2``. Note that these crystal operators always subtract only the group delay ``\beta_1\omega`` (calculated from `nfuny`) as the reference phase, for both field-resolved and envelope simulations, since the second-order nonlinearity they are designed for requires this reference frame (see [Reference frames and carrier-mixing nonlinearities](@ref)).

### Implementation
- [`PhysData.ref_index_fun_uniax`](@ref)
- [`PhysData.ref_index_fun_xy`](@ref)
- [`PhysData.crystal_internal_angle`](@ref)
- [`LinearOps.make_const_linop`](@ref)
- [`NonlinearRHS.const_norm_free2D`](@ref)
- [`NonlinearRHS.const_norm_free`](@ref)
