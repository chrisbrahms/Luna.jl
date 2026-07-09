# Modal decompositions
!!! note
    All Fourier transforms in this section are written in the optics sign convention. See [A note on sign conventions](@ref) for details on what is used in the code.

## Multi-mode guided
For propagation in waveguides taking into account multiple modes and the coupling between them, Luna uses the model laid out in [Kolesik and Moloney, *Nonlinear optical pulse propagation simulation: From Maxwell’s to unidirectional equations*](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.70.036604) and [Tani et al., *Multimode ultrafast nonlinear optics in optical waveguides: numerical modeling and experiments in kagomé photonic-crystal fiber*](http://josab.osa.org/abstract.cfm?URI=josab-31-2-311). This is implemented in [`NonlinearRHS.TransModal`](@ref). The electric field ``\mathbf{E}(t, \mathbf{r_\perp}, z)`` is expressed as the inverse Fourier transform in time and the superposition of waveguide modes in space. This means that the transverse wave vector ``\mathbf{k}_\perp`` turns into a modal index ``j`` (this transform is implemented in [`Modes.ToSpace`](@ref) and [`Modes.to_space!`](@ref)):
```math
\mathbf{E}(t, \mathbf{r_\perp}, z) = \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{d} \omega \sum_j \hat{\mathbf{e}}_j(\mathbf{r_\perp}, z) \tilde{A}_j(\omega, z) \mathrm{e}^{-i \omega t}\,,
```
where ``\hat{\mathbf{e}}_j(\mathbf{r_\perp}, z)`` is the orthonormal transverse field distribution of the ``j^{\mathrm{th}}`` mode and ``\tilde{A}_j(\omega, z)`` is the frequency-domain amplitude in mode ``j``. The mode fields ``\hat{\mathbf{e}}_j(\mathbf{r_\perp}, z)`` are taken to be independent of frequency but can depend on the propagation coordinate ``z`` (e.g. in tapered waveguides). They can be vector quantities if polarisations other than purely lineary ``x``- or ``y``-polarisations need to be taken into account. The modes are normalised such that ``\vert \tilde{A}_j(\omega, z) \vert^2`` gives the spectral energy density in mode ``j`` (when also taking into account the normalisation of the FFT), and equivalently ``\vert A_j(t, z)\vert^2`` gives the instantaneous power. The forward transform to reciprocal space is simply the overlap integral of the total field with each mode combined with the Fourier transform in time:
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
In some situations, inter-mode coupling in a waveguide is negligible, so including several waveguide modes in the simulation unnecessarily slows down the computation. Simulating propagation in a single mode is trivially achieved by including only that single mode in both the forward and inverse transforms as defined above for [multi-mode propagation](#multi-mode-guided). For example, setting `modes=1` when calling `prop_capillary` achieves this and leads to a significant speed-up. However, in this simple implementation, the overlap integral between the nonlinear polarisation and the waveguide mode still needs to be calculated explicitly. We can make this unnecessary by making an assumption about the nonlinear polarisation.

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


## Radially symmetric free-space
For free-space propagation of a radially symmetric beam, the transverse basis is no longer a set of waveguide modes but the continuum of Bessel beams ``J_0(k_\perp r)``, so the generalised transverse spatial frequency ``\mathbf{k}_\perp`` becomes the (scalar) radial spatial frequency ``k_\perp``. The field is written as
```math
\mathbf{E}(t, r, z) = \sum_{s\in\{x,y\}}\hat{\mathbf{e}}_s\, \frac{2\pi}{\left(2\pi\right)^3} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}  E_s(\omega, k_\perp, z) J_0(k_\perp r) \mathrm{e}^{-i\omega t} k_\perp\mathrm{d} k_\perp  \mathrm{d}\omega\,,
```
where ``E_s`` is the amplitude of polarisation component ``s``, and the ``k_\perp\,\mathrm{d}k_\perp`` measure and the leading ``2\pi`` (from the azimuthal integration) are the conventions of the Hankel transform. In Luna this transverse transform is a zeroth-order quasi-discrete Hankel transform (`Hankel.QDHT`), and the whole radial UPPE is implemented in [`NonlinearRHS.TransRadial`](@ref). On each step the nonlinear polarisation is evaluated by transforming the field from reciprocal space to real-space–time (an inverse Fourier transform ``\omega\to t`` and an inverse Hankel transform ``k_\perp\to r``, the latter applied to each polarisation component in turn), calculating ``\mathbf{P}_\mathrm{nl}`` pointwise on the ``(E_x, E_y)`` field at each radial point, and transforming back (``r\to k_\perp`` and ``t\to\omega``), followed by the ``i\omega/N_\mathrm{nl}`` prefactor of the UPPE.

The normalisation factor here is not constant but depends on the propagation angle of each plane-wave component, and on its polarisation. For a component at frequency ``\omega`` and radial spatial frequency ``k_\perp`` in polarisation ``s`` it is (see [`NonlinearRHS.norm_radial`](@ref) and its ``z``-independent form [`NonlinearRHS.const_norm_radial`](@ref))
```math
N_\mathrm{nl}(\omega, k_\perp, z) = \frac{\beta_s(\omega, k_\perp, z)}{\mu_0 \omega}\,,\qquad \beta_s = \sqrt{\left(\frac{\omega}{c}n_s(\omega, z)\right)^2 - k_\perp^2}\,,
```
where ``\beta_s`` is the longitudinal wave vector for polarisation ``s``. On axis (``k_\perp = 0``) this reduces to the ``n_s\varepsilon_0 c`` impedance factor familiar from the guided cases, while for oblique components it applies the corresponding angular-spectrum correction. Evanescent components (``\beta_s^2 \le 0``) are guarded and do not propagate. The propagating array is three-dimensional in ``(\omega, N_\mathrm{pol}, k_\perp)``.

### Implementation
- [`NonlinearRHS.TransRadial`](@ref)
- [`NonlinearRHS.norm_radial`](@ref)
- [`NonlinearRHS.const_norm_radial`](@ref)

## Two-dimensional (planar) free-space
For planar propagation with a single transverse dimension ``x`` (a field that is uniform along ``y``, e.g. a sheet beam in the ``x``–``z`` plane), the transverse basis is the one-dimensional plane-wave spectrum ``\mathrm{e}^{i k_x x}``. The Cartesian ``x`` grid and its spatial frequencies ``k_x`` are held by [`Grid.Free2DGrid`](@ref), and the transform is implemented in [`NonlinearRHS.TransFree2D`](@ref):
```math
\mathbf{E}(t, x, z) = \sum_{s\in\{x,y\}}\hat{\mathbf{e}}_s\,\frac{1}{\left(2\pi\right)^2}\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}  E_s(\omega, k_x, z) \mathrm{e}^{i\left(k_x x - \omega t\right)}\mathrm{d} k_x \mathrm{d} \omega\,.
```
A combined FFT over ``(\omega, k_x)\to(t, x)`` (skipping the polarisation axis) transforms to real space, where ``\mathbf{P}_\mathrm{nl}`` is evaluated pointwise, before the inverse transform and the ``i\omega/N_\mathrm{nl}`` prefactor. The normalisation is the same ``\beta_s/(\mu_0\omega)`` as the radial case, with ``\beta_s = \sqrt{(\omega n_s/c)^2 - k_x^2}`` (implemented in `NonlinearRHS.norm_free2D` / `NonlinearRHS.const_norm_free2D`). The propagating array is three-dimensional in ``(\omega, N_\mathrm{pol}, k_x)``. Note that both ``x`` **and** ``y`` polarisation components exist even though there is only one transverse *spatial* dimension — this geometry is what makes birefringent-crystal ``\chi^{(2)}`` interactions cheap to simulate.

### Implementation
- [`Grid.Free2DGrid`](@ref)
- [`NonlinearRHS.TransFree2D`](@ref)

## Three-dimensional free-space
For full ``(3+1)``-dimensional free-space propagation without any symmetry, the transverse basis is the two-dimensional plane-wave spectrum ``\mathrm{e}^{i(k_x x + k_y y)}``, so ``\mathbf{k}_\perp = (k_x, k_y)``. The field is
```math
\mathbf{E}(t, x, y, z) = \sum_{s\in\{x,y\}}\hat{\mathbf{e}}_s\,\frac{1}{\left(2\pi\right)^3}\int_{-\infty}^{\infty} \int_{-\infty}^{\infty}  \int_{-\infty}^{\infty}  E_s(\omega, k_x, k_y, z) \mathrm{e}^{i\left(k_x x + k_y y - \omega t\right)}\mathrm{d} k_x \mathrm{d} k_y \mathrm{d} \omega\,,
```
with the standard three-fold inverse-Fourier-transform measure ``1/(2\pi)^3``. The Cartesian ``x``/``y`` grid and its spatial frequencies ``k_x``/``k_y`` are held by [`Grid.FreeGrid`](@ref), and the transform is implemented in [`NonlinearRHS.TransFree`](@ref). Here a single combined FFT handles both transverse spatial axes *and* time at once (over dimensions ``(\omega, k_y, k_x)``, again skipping the polarisation axis): the field is transformed ``(\omega, k_y, k_x)\to(t, y, x)``, the nonlinear polarisation is evaluated pointwise on the ``(E_x, E_y)`` field, and the result is transformed back before applying the ``i\omega/N_\mathrm{nl}`` prefactor.

The normalisation takes the same form as the radial case, now with the two-dimensional transverse spatial frequency ``k_\perp^2 = k_x^2 + k_y^2`` (see [`NonlinearRHS.norm_free`](@ref) and [`NonlinearRHS.const_norm_free`](@ref)):
```math
N_\mathrm{nl}(\omega, k_x, k_y, z) = \frac{\beta_s(\omega, k_x, k_y, z)}{\mu_0 \omega}\,,\qquad \beta_s = \sqrt{\left(\frac{\omega}{c}n_s(\omega, z)\right)^2 - k_x^2 - k_y^2}\,,
```
with the same evanescent-component guard. The propagating array is four-dimensional in ``(\omega, N_\mathrm{pol}, k_y, k_x)``.

### Implementation
- [`Grid.FreeGrid`](@ref)
- [`NonlinearRHS.TransFree`](@ref)
- [`NonlinearRHS.norm_free`](@ref)
- [`NonlinearRHS.const_norm_free`](@ref)

## Polarisation in free space
Free-space geometries can propagate a two-component transverse vector field ``\mathbf{E} = E_x\hat{\mathbf{x}} + E_y\hat{\mathbf{y}}``. Every free-space field, nonlinear-polarisation, linear-operator and normalisation array therefore carries a **polarisation axis** (of length ``N_\mathrm{pol}``) immediately after the frequency/time axis, giving the shape ``(N_\omega, N_\mathrm{pol}, N_\perp\ldots)``. ``N_\mathrm{pol}`` is either

- **1** — a scalar simulation along a single (``y``) polarisation, or
- **2** — a full ``(E_x, E_y)`` vector simulation.

Which one is used is determined by the refractive-index function passed to the linear operator: if it returns a single index the run is scalar, and if it returns a pair ``(n_x, n_y)`` the run is polarisation-resolved. Because the two components can be given different indices, birefringent media are supported (for uniaxial crystals ``n_x`` is evaluated at the internal angle, see below and [`LinearOps.make_const_linop`](@ref)). The input polarisation state can be set by a rotation angle ``\theta`` on the spatiotemporal input field (see [`Fields.GaussGaussField`](@ref)), which splits the pulse into ``E_x = E\sin\theta`` and ``E_y = E\cos\theta`` (so the default ``\theta = 0`` is ``y``-polarised). Vector nonlinear responses — the second-order response [`Nonlinear.Chi2Field`](@ref), and the vector Kerr and plasma responses — then couple the two components; this is what makes e.g. birefringent-crystal ``\chi^{(2)}`` interactions possible in free space. In the equations above, the sum over ``s \in \{x, y\}`` collapses to a single term for a scalar (``N_\mathrm{pol} = 1``) run.

## Birefringent crystals
One motivation for polarisation-resolved free-space propagation is nonlinear frequency conversion in birefringent crystals (for example ``\chi^{(2)}`` second-harmonic generation, implemented by [`Nonlinear.Chi2Field`](@ref)). These require a genuine two-component field, but they also complicate the *linear* operator, because in a uniaxial crystal the **extraordinary refractive index depends on the propagation angle** relative to the optic axis.

This matters in free space specifically because the transverse decomposition spreads the field over a range of plane-wave directions: a component with transverse wavevector ``k_x`` is a plane wave travelling at an angle ``\theta_i`` to the ``z`` axis, with ``k_x = (\omega/c)\sin\theta_i``. For the ordinary (``y``) polarisation the index is angle-independent, so a single ``n_y(\omega)`` suffices as in the isotropic case. For the extraordinary (``x``) polarisation, **each transverse wavevector sees a different index**, so ``n_x`` cannot be tabulated as a function of ``\omega`` alone.

Two helpers in [PhysData.jl](@ref) provide the angle-dependent index:

- [`PhysData.ref_index_fun_xy`](@ref)`(material, θ)` returns the pair ``(n_x, n_y)`` for a crystal cut at angle ``\theta`` to the optic axis. ``n_y(\lambda)`` is the ordinary index; ``n_x(\lambda, \delta\theta)`` is the extraordinary index at propagation angle ``\theta + \delta\theta``, where ``\delta\theta`` is the offset of a given plane-wave component from the nominal cut angle.
- [`PhysData.crystal_internal_angle`](@ref)`(n_x, ω, k_x)` recovers that offset ``\delta\theta`` from the transverse wavevector. Refraction at the crystal surface conserves the transverse wavevector, so the external ``k_x`` must equal the internal ``(\omega/c)\,n_x(\lambda, \delta\theta)\sin\delta\theta``. Because ``n_x`` itself depends on ``\delta\theta``, this is a transcendental equation and is solved numerically.

The birefringent linear operator is assembled by the tuple-valued methods of [`LinearOps.make_const_linop`](@ref) (and its ``z``-dependent counterpart), which take ``\mathrm{nfuns} = (n_x, n_y)``. For each frequency ``\omega`` and transverse wavevector, the operator fills the two polarisation columns separately:
```math
\beta_y = \sqrt{\left(\tfrac{\omega}{c} n_y(\omega)\right)^2 - k_\perp^2}\,,\qquad
\beta_x = \sqrt{\left(\tfrac{\omega}{c} n_x\!\big(\omega, \delta\theta(k_x)\big)\right)^2 - k_\perp^2}\,,
```
where ``\delta\theta(k_x)`` comes from [`PhysData.crystal_internal_angle`](@ref) for that ``k_x`` (the optic axis is taken to lie in the ``x``–``z`` plane, so only ``k_x`` tilts a component towards or away from it). The reference-frame velocity is taken from the ordinary index ``n_y``. Everything else — the transforms and normalisation — is exactly as for the isotropic free-space cases above, now simply carrying two polarisation columns with different longitudinal wavevectors.