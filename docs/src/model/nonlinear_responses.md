# Nonlinear responses

This page describes the nonlinear polarisation terms ``P_\mathrm{nl}`` that appear on the right-hand side of the UPPE (see [The numerical model](@ref "The numerical model")). Regardless of the [modal decomposition](@ref "Modal decompositions") in use, these terms are always evaluated in the real-space–time domain: the field is transformed to time (and, for spatially resolved geometries, to real space), each response is evaluated pointwise, and the result is transformed back. A simulation can combine several of these effects at once — they are supplied as a tuple of response functions (assembled from the keyword arguments by `Interface.makeresponse` for the simple interface) and simply summed to give the total ``P_\mathrm{nl}``.

Each response is written to be independent of the propagation geometry, so the same functions are used for mode-averaged, multi-mode, radial, two-dimensional and full 3D simulations. In every geometry, a response receives the time-domain field at one transverse position as an array of shape ``(N_t, N_\mathrm{pol})``, where ``N_\mathrm{pol}`` is 1 for scalar and 2 for two-polarisation simulations (``(E_x, E_y)``), and adds its contribution to the polarisation array of the same shape. The third argument is the number density of the medium (as returned by `densityfun`), which is how pressure gradients enter the nonlinearity. For solid media the number density is not a useful quantity, so the macroscopic susceptibility is passed directly to the response and a unit density is used (`densityfun = z -> 1`).

## Kerr effect
The Kerr effect is the third-order (``\chi^{(3)}``) nonlinearity responsible for self- and cross-phase modulation, four-wave mixing and third-harmonic generation. In Luna the susceptibility is built at runtime from the number density and a single-molecule hyperpolarisability, ``\chi^{(3)} = \rho\,\gamma_3``, where ``\gamma_3`` is provided by [`PhysData.γ3_gas`](@ref) and ``\rho`` is the gas density at the current position (allowing the nonlinearity to follow a pressure gradient). The macroscopic susceptibility can equivalently be obtained from the nonlinear index ``n_2`` via [`PhysData.χ3`](@ref), ``\chi^{(3)} = \tfrac{4}{3} n_2 \varepsilon_0 c\, n_0^2``. For solids (glasses and the nonlinear crystals listed in `PhysData.crystal`), `PhysData.χ3(material)` returns the macroscopic susceptibility directly, so the response is created as e.g. `Kerr_field(PhysData.χ3(:BBO))` and used with unit density.

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
If the third harmonic *is* required in an envelope simulation, [`Nonlinear.Kerr_env_thg`](@ref) restores it as an explicitly up-shifted term ``\propto \mathrm{e}^{2i\omega_0 t}\mathcal{E}^2\mathcal{E}`` (see Eq. 4 of Genty et al., *Opt. Express* **15**, 5382 (2007)). This requires a grid which can hold the third harmonic (`Grid.EnvGrid(...; thg=true)` with `λ_lims` extending below ``\lambda_0/3``) and, because the response mixes the carrier, a linear operator using the carrier-transparent reference frame, see [Reference frames and carrier-mixing nonlinearities](@ref). In the simple interface, `prop_capillary(...; envelope=true, thg=true)` takes care of both.

For an **envelope, vector** simulation the response additionally carries the cross-polarisation coupling terms, with the characteristic ``\tfrac{2}{3}`` cross-phase and ``\tfrac{1}{3}`` phase-conjugate (coherent coupling) coefficients that arise from separating the co- and counter-rotating circular components:
```math
\mathbf{P}_\mathrm{Kerr}(t, x, y, z) = \frac{3}{4}\varepsilon_0 \chi^{(3)} \begin{pmatrix}
           (\mathcal{E}_x^2 + \frac{2}{3}\mathcal{E}_y^2)\mathcal{E}_x + \frac{1}{3}\mathcal{E}_x^*\mathcal{E}_y^2 \\[0.5em]
           (\mathcal{E}_y^2 + \frac{2}{3}\mathcal{E}_x^2)\mathcal{E}_y + \frac{1}{3}\mathcal{E}_y^*\mathcal{E}_x^2
\end{pmatrix}
```

## Second-order nonlinearity
The second-order (``\chi^{(2)}``) nonlinearity gives rise to three-wave mixing: second-harmonic generation (SHG), sum- and difference-frequency generation (SFG/DFG) and optical rectification. It only exists in non-centrosymmetric media, i.e. in Luna it is used for propagation in nonlinear crystals (see [Birefringent crystals](@ref) for the linear propagation in this case). Luna includes this response in two variants: [`Nonlinear.Chi2Field`](@ref) for real (carrier-resolved) fields on a `Grid.RealGrid`, and [`Nonlinear.Chi2Env`](@ref) for complex envelope fields on a [`Grid.EnvGrid`](@ref). The real-field variant follows directly from the constitutive relation for the second-order polarisation and needs no further approximation. The envelope variant is derived from it, and this section documents that derivation in detail, including the choices which make it exact for all field content the envelope grid can represent, and the constraints it places on the reference frame used by the linear operator.

Both variants require a **two-polarisation** simulation, i.e. the field passed to them must have both ``x`` and ``y`` components (they throw an error otherwise), and both ignore the density argument, so they are used with `densityfun = z -> 1`. Since these responses are evaluated pointwise on the time-domain field, they act on the **total** field, so all three-wave mixing processes—SHG, SFG/DFG between all spectral components, optical rectification, and their cascades—are included automatically, limited only by the frequency window of the grid.

### The real-field vector response
The starting point is the second-order nonlinear polarisation in the frame of the crystal axes,

```math
P_i(t) = \varepsilon_0 \sum_{jk} \chi^{(2)}_{ijk} E_j(t) E_k(t)\,,
```

where ``i,j,k \in \{x,y,z\}`` run over the three crystal axes and ``\chi^{(2)}_{ijk}`` is the second-order susceptibility tensor. Since ``\chi^{(2)}_{ijk}`` is symmetric in ``j \leftrightarrow k`` (the two driving fields are physically indistinguishable), the pair ``jk`` can be contracted into a single index, giving the standard 3×6 contracted (Voigt) matrix with column order

```math
jk \in [xx,\; yy,\; zz,\; yz,\; xz,\; xy]\,,
```

in which the mixed products appear with a factor of 2. The polarisation is then a matrix-vector product

```math
\mathbf{P}(t) = \varepsilon_0\, \chi^{(2)}\, \mathbf{s}(t)\,,
\qquad
\mathbf{s} = \big[E_x^2,\; E_y^2,\; E_z^2,\; 2E_yE_z,\; 2E_xE_z,\; 2E_xE_y\big]\,,
```

where ``\chi^{(2)}`` is now the 3×6 matrix (twice the ``d`` matrix commonly tabulated for nonlinear crystals). The product vector ``\mathbf{s}`` is computed by [`Nonlinear.field_products!`](@ref). The 3×6 matrices for supported crystals (in SI units, m/V) are available from [`PhysData.χ2`](@ref).

The simulation stores the field in the *lab* frame, in which the pulse propagates along ``z`` and carries two transverse polarisation components ``(E_x, E_y)``. The crystal is rotated relative to this frame by the polar angle ``\theta`` (about ``y``) and the azimuthal angle ``\phi`` (about ``z``), so the field must be rotated into the crystal frame before applying ``\chi^{(2)}``, and the resulting polarisation rotated back:

```math
\mathbf{P}_\mathrm{lab}(t) = \varepsilon_0\, R_\mathrm{lab}\, \chi^{(2)}\, \mathbf{s}\big(R_\mathrm{cr} \mathbf{E}_\mathrm{lab}(t)\big)\,,
```

where ``R_\mathrm{cr}`` rotates lab-frame vectors into the crystal frame and ``R_\mathrm{lab} = R_\mathrm{cr}^{-1}`` (in the code, `toCrystal = RotMatrix(RotZY(-ϕ, -θ))` and `toLab = RotMatrix(RotYZ(θ, ϕ))`, and the product ``R_\mathrm{lab}\chi^{(2)}`` is precomputed). The angle ``\theta`` is the same cut angle which enters the linear propagation in a [birefringent crystal](@ref "Birefringent crystals"), since the optic axis of a uniaxial crystal is its ``z`` axis. Two approximations are made, both consequences of the unidirectional, transverse field description:

- the lab-frame longitudinal field is neglected on input, ``E_{z,\mathrm{lab}} = 0`` (the *crystal*-frame ``E_z`` is generally non-zero after rotation and does contribute to the products);
- the lab-frame longitudinal component of the polarisation, ``P_{z,\mathrm{lab}}``, is discarded.

Applied sample-by-sample in the time domain, this is exactly [`Nonlinear.Chi2Field`](@ref), constructed as `Chi2Field(θ, ϕ, χ2)`.

### Envelope conventions
!!! note
    Following the [A note on sign conventions](@ref) section, all equations in the remainder of this section are written in the convention of the *simulated* fields, i.e. the mathematics/FFT convention in which the analytic signal is proportional to ``\mathrm{e}^{+i\omega_0 t}``. This means every expression here maps one-to-one onto the source code (for example the carrier phase array `C = exp.(1im*ω0.*t)` in [`Nonlinear.Chi2Env`](@ref)). In the optics convention the signs of all exponents are flipped.

For envelope propagation, `Luna` stores the complex envelope ``\mathbf{A}(t)`` defined by

```math
\mathbf{E}(t) = \mathrm{Re}\left[\mathbf{A}(t)\,\mathrm{e}^{i\omega_0 t}\right]
             = \frac{1}{2}\left(\mathbf{A}(t)\,\mathrm{e}^{i\omega_0 t} + \mathbf{A}^*(t)\,\mathrm{e}^{-i\omega_0 t}\right)\,,
```

where ``\omega_0`` is the carrier frequency of the grid (`grid.ω0`). The quantity ``\mathbf{A}\,\mathrm{e}^{i\omega_0 t}`` is the *analytic signal* of the field: its spectrum is twice the positive-frequency half of the real field's spectrum, and zero at negative frequencies. Consequently ``|\mathbf{A}|`` equals the amplitude of the real field. On a [`Grid.EnvGrid`](@ref) the envelope's spectral content at offset frequency ``\nu`` corresponds to the absolute frequency ``\omega = \omega_0 + \nu``, and `grid.ω` stores these absolute frequencies. The spectral apodisation window `grid.ωwin`, which is applied to the field and to the nonlinear polarisation at every step, is zero for all absolute frequencies below half the lower edge of the requested wavelength window, so any content at negative absolute frequencies is removed immediately.

The task is therefore: given the real-field response above, find the envelope ``\tilde{\mathbf{P}}(t)`` of the nonlinear polarisation, defined in the same way, ``\mathbf{P}(t) = \mathrm{Re}[\tilde{\mathbf{P}}(t)\,\mathrm{e}^{i\omega_0 t}]``, with ``\tilde{\mathbf{P}}`` equal to the demodulated analytic signal of ``\mathbf{P}``.

### Deriving the envelope response
Since the response is built entirely from the pair products ``E_j E_k``, it is enough to find the envelope of one such product. Substituting the envelope representation for both factors:

```math
E_j(t) E_k(t) = \frac{1}{4}\Big[
      \underbrace{A_j A_k\, \mathrm{e}^{2i\omega_0 t}}_{\text{SFG}}
    + \underbrace{A_j A_k^* + A_j^* A_k}_{\text{DFG}}
    + \underbrace{A_j^* A_k^*\, \mathrm{e}^{-2i\omega_0 t}}_{\text{negative frequencies}}
\Big]\,.
```

The three groups of terms live in different spectral regions (assuming, as always for an envelope description, that the field content lies at positive absolute frequencies):

1. **The SFG term** ``A_jA_k\mathrm{e}^{2i\omega_0t}`` contains only positive absolute frequencies (each factor sits at ``\omega_0 + \nu > 0``, and the sum of two positive frequencies is positive). Its contribution to the analytic signal of the product is therefore obtained simply by doubling: ``\tfrac{1}{2}A_jA_k\mathrm{e}^{2i\omega_0t}``. Demodulating by ``\mathrm{e}^{i\omega_0 t}`` gives the envelope contribution ``\tfrac{1}{2}A_jA_k\,\mathrm{e}^{+i\omega_0 t}``, which shifts the envelope spectrum *up* by ``\omega_0``: content at offsets ``\nu_1,\nu_2`` lands at absolute frequency ``2\omega_0 + \nu_1 + \nu_2``. This term carries sum-frequency mixing, and in particular SHG.
2. **The DFG term** ``A_jA_k^* + A_j^*A_k`` is real, so its spectrum is symmetric about ``\omega = 0``: every difference-frequency pair appears once at ``+(\nu_1-\nu_2)`` and once, conjugated, at ``-(\nu_1-\nu_2)``. The analytic signal keeps twice the positive-frequency half. Rather than explicitly filtering, the implementation keeps the *whole* real product with a factor ``\tfrac{1}{2}\times 2 = 1`` relative to the ``\tfrac14`` above, demodulated *down* by the carrier: ``\tfrac{1}{2}(A_jA_k^* + A_j^*A_k)\,\mathrm{e}^{-i\omega_0 t}``. Content at offsets ``\nu_1, \nu_2`` now lands at absolute frequency ``\nu_1 - \nu_2`` *and* at ``\nu_2 - \nu_1``: for each pair, one partner lies at a positive absolute frequency—on the grid, with exactly the correct analytic-signal coefficient because it appears only once in the sum—and the other lies at a negative absolute frequency, where it is removed by the spectral apodisation (`grid.ωwin`) applied at every propagation step. The same applies to the optical-rectification content (``j`` and ``k`` content at the *same* frequency), which lands at absolute frequencies near zero.
3. **The negative-frequency term** ``A_j^*A_k^*\mathrm{e}^{-2i\omega_0t}`` contains only negative absolute frequencies and contributes nothing to the analytic signal; it is dropped entirely.

The envelope of each pair product is therefore

```math
(E_j E_k)_\mathrm{env} = \frac{1}{2} A_j A_k\, \mathrm{e}^{+i\omega_0 t}
                       + \frac{1}{2}\left(A_j A_k^* + A_j^* A_k\right) \mathrm{e}^{-i\omega_0 t}\,.
```

#### Checking the coefficients
The factors of ``\tfrac12`` are easy to get wrong, so it is worth verifying them against the real-field result for a concrete case: a fundamental and its second harmonic co-propagating on one grid, ``A = A_1 + A_2\,\mathrm{e}^{i\omega_0 t}`` (the second harmonic sits at offset ``+\omega_0``, i.e. absolute frequency ``2\omega_0``). Directly from the real field ``E = \mathrm{Re}[A_1\mathrm{e}^{i\omega_0t}] + \mathrm{Re}[A_2\mathrm{e}^{2i\omega_0t}]``:

- the ``2\omega_0`` content of ``E^2`` is ``\tfrac14 A_1^2\mathrm{e}^{2i\omega_0t} + \mathrm{c.c.}``, whose analytic signal is ``\tfrac12A_1^2\mathrm{e}^{2i\omega_0t}``—the **SHG** drive;
- the ``\omega_0`` content of ``E^2`` is ``\tfrac12 A_2A_1^*\mathrm{e}^{i\omega_0t} + \mathrm{c.c.}``, whose analytic signal is ``A_2A_1^*\mathrm{e}^{i\omega_0t}``—the **back-conversion** drive.

The envelope expression reproduces both exactly: the SFG term gives ``\tfrac12A_1^2`` at offset ``+\omega_0`` (absolute ``2\omega_0``), and the DFG term contains ``A_2A_1^*\mathrm{e}^{i\omega_0t}`` once (from ``A_jA_k^*`` with ``j`` on the second harmonic), which after the ``\mathrm{e}^{-i\omega_0t}`` demodulation lands at offset ``0`` (absolute ``\omega_0``) with unit coefficient—while its conjugate partner ``A_1A_2^*`` lands at absolute ``-\omega_0``, off the grid. Back-conversion is *not* an optional extra: it is required for energy-conserving SHG dynamics (pump depletion and re-conversion), which is why the DFG term is always included in the response.

#### The contracted product vectors
Carrying the two carrier phases through the contracted notation gives two 6-vectors built from the *crystal-frame* envelope components. The SFG vector is complex,

```math
\mathbf{s}_\mathrm{SFG} = \big[A_x^2,\; A_y^2,\; A_z^2,\; 2A_yA_z,\; 2A_xA_z,\; 2A_xA_y\big]\,,
```

while the DFG vector is purely real (each entry is of the form ``A_jA_k^* + A_j^*A_k = 2\,\mathrm{Re}(A_jA_k^*)``),

```math
\mathbf{s}_\mathrm{DFG} = \big[|A_x|^2,\; |A_y|^2,\; |A_z|^2,\; 2\,\mathrm{Re}(A_yA_z^*),\; 2\,\mathrm{Re}(A_xA_z^*),\; 2\,\mathrm{Re}(A_xA_y^*)\big]\,,
```

and the full envelope polarisation, including the rotations, is

```math
\tilde{\mathbf{P}}_\mathrm{lab}(t) = \varepsilon_0\, R_\mathrm{lab}\, \chi^{(2)} \left[
      \frac{1}{2}\mathrm{e}^{+i\omega_0 t}\,\mathbf{s}_\mathrm{SFG}(t)
    + \mathrm{e}^{-i\omega_0 t}\,\mathbf{s}_\mathrm{DFG}(t)
\right],
```

evaluated with ``\mathbf{A}_\mathrm{cr} = R_\mathrm{cr}\mathbf{A}_\mathrm{lab}`` and, as for the real-field response, ``A_{z,\mathrm{lab}} = 0`` on input and ``\tilde{P}_{z,\mathrm{lab}}`` discarded. This is exactly what [`Nonlinear.env_products!`](@ref) and the [`Nonlinear.Chi2Env`](@ref) callable compute, with the two phase factors combined into a single precomputed carrier array (``\mathrm{e}^{-i\omega_0t}`` is obtained as the conjugate of ``\mathrm{e}^{+i\omega_0t}``). Since the response is evaluated on the oversampled time grid, `Chi2Env` is constructed as `Chi2Env(θ, ϕ, χ2, grid.ω0, grid.to)`. Note that since the rotations and the ``\chi^{(2)}`` matrix are real, they commute with taking the envelope, so the same precomputed ``R_\mathrm{lab}\chi^{(2)}`` is used as in the real-field response.

Because the derivation above is exact for all field content representable on the grid—the only discarded terms are those the grid cannot hold anyway—the envelope response agrees with the real-field response to numerical precision within the grid's frequency window. This is verified directly in `test/test_chi2.jl`, which compares the two responses on the same two-colour, two-polarisation field, and in `test/test_freespace.jl`, where type I SHG in BBO simulated with `RealGrid`/`Chi2Field` and `EnvGrid`/`Chi2Env` through the full propagation pipeline agrees in the converted second-harmonic energy.

### Grid requirements
The SFG term shifts spectral content *up* by ``\omega_0``, so for the envelope response an [`Grid.EnvGrid`](@ref) must fulfil **both** of the following:

- The wavelength limits `λ_lims` must extend below ``\lambda_0/2``, so that the second harmonic lies within the frequency window. Otherwise it is simply removed by the apodisation window `grid.ωwin` at every step. (The same is of course true for a `Grid.RealGrid` and the real-field response.)
- The grid must be created with `thg=true`. This samples the oversampled time grid `grid.to`, on which the response is evaluated, three times more densely than the coarse grid, so that the products of the field with itself and their up-shift by ``\omega_0`` do not alias back into the frequency window. Without oversampling, mixing products which lie above the top of the coarse frequency grid (e.g. the sum frequency of the second harmonic and the fundamental) alias to low absolute frequencies within the window and corrupt the simulation.

In addition, the linear operator must use the reference frame described in the next section.

## Reference frames and carrier-mixing nonlinearities
Envelope propagation in `Luna` is performed in a moving and, optionally, rotating frame: the linear operator subtracts a reference phase ``\varphi(\omega)`` from the propagation constant (see [Modal decompositions](@ref)), so the stored field relates to the physical one by ``E'(\omega, z) = E(\omega, z)\,\mathrm{e}^{i\varphi(\omega)z}`` (in the optics sign convention). This choice interacts with nonlinearities which *mix the carrier*, such as THG and all χ⁽²⁾ processes, in an important way.

Consider a nonlinear product of ``n`` field factors (conjugated factors counting negatively) at frequencies ``\omega_1, \ldots, \omega_n``, generating polarisation at ``\omega_\Sigma = \sum \pm\,\omega_i``. Computed from the *stored* fields, the product carries the frame phase ``\sum \pm\varphi(\omega_i)``, whereas the frame source at ``\omega_\Sigma`` should carry ``\varphi(\omega_\Sigma)``. For an affine reference phase ``\varphi(\omega) = \beta_1\omega + c_0`` the discrepancy is

```math
\sum \pm\,\varphi(\omega_i) - \varphi(\omega_\Sigma) = (n_+ - n_- - 1)\, c_0\,,
```

independent of the frequencies involved, where ``n_\pm`` count the unconjugated/conjugated factors. The group-delay part ``\beta_1\omega`` always cancels—it is a pure time shift—but a constant offset ``c_0`` survives whenever ``n_+ - n_- \neq 1`` and acts as a *spurious phase mismatch* accumulating at ``(n_+ - n_- - 1)c_0`` per unit length:

| process | factors | spurious mismatch |
|:--|:--|:--|
| SPM/XPM (Kerr, ``\lvert A\rvert^2A``) | ``n_+=2,\ n_-=1`` | ``0`` |
| SHG/SFG (``A^2``) | ``n_+=2,\ n_-=0`` | ``c_0`` |
| DFG/back-conversion (``AA^*``) | ``n_+=1,\ n_-=1`` | ``-c_0`` |
| THG (``A^3``) | ``n_+=3,\ n_-=0`` | ``2c_0`` |

The standard co-rotating envelope frame subtracts ``\varphi(\omega) = \beta_1(\omega - \omega_0) + \beta_0`` with ``\beta_0 = \beta(\omega_0)``, so that the envelope is phase-stationary at the carrier. Its constant offset is ``c_0 = \beta_0 - \beta_1\omega_0 = (n_\phi - n_g)\,\omega_0/c``, the difference of phase and group index at the carrier. This is perfectly fine for the Kerr effect (as well as Raman scattering and the plasma response, none of which mix the carrier), but *not* for χ⁽²⁾ or THG: in BBO at 800 nm, ``c_0 \approx -1.9\times10^5\,\mathrm{m}^{-1}``, comparable to or larger than real phase-mismatch scales, which visibly suppresses (or otherwise falsifies) SHG. Setting ``\beta_0 = 0`` instead is even worse: then ``c_0 = -\beta_1\omega_0``, of order ``10^7\,\mathrm{m}^{-1}``, which both destroys harmonic generation entirely and forces the adaptive integrator to take extremely small steps to resolve the spurious oscillation of the source term.

The only frame transparent to *all* instantaneous nonlinearities is ``c_0 = 0``, i.e. a subtracted phase strictly linear in the absolute frequency:

```math
\varphi(\omega) = \beta_1\,\omega\,,
```

a pure group-delay (time-shift) frame. This still co-moves with the pulse exactly as the co-rotating frame does (the two differ only by the constant ``c_0``, a global phase rotation ``\mathrm{e}^{ic_0z}`` of the envelope which affects no intensity or spectrum). It is also the frame which is always used for field-resolved simulations on a `Grid.RealGrid`, where there is no carrier to subtract.

In practice, this frame is selected as follows:

- All linear operators in [LinearOps.jl](@ref) ([`LinearOps.make_const_linop`](@ref) and [`LinearOps.make_linop`](@ref) for mode-averaged, multi-mode and free-space propagation) accept a `thg` argument. On a `Grid.RealGrid` it must be `true` (the default). On a `Grid.EnvGrid`, `thg=false` (the default) selects the co-rotating frame, and `thg=true` selects the carrier-transparent frame ``\varphi(\omega) = \beta_1\omega``. The internal function `LinearOps.getω0` documents the convention.
- The [birefringent-crystal](@ref "Birefringent crystals") operators (the methods of [`LinearOps.make_const_linop`](@ref) taking a tuple of index functions) use the carrier-transparent frame unconditionally for both grid types, since they are designed for χ⁽²⁾ propagation.
- In the simple interface, `prop_capillary(...; envelope=true, thg=true)` passes the `thg` flag through to the linear operator automatically.

Whenever [`Nonlinear.Chi2Env`](@ref) or [`Nonlinear.Kerr_env_thg`](@ref) is used, the linear operator **must** use this frame. That the combination of `thg=true` grid, response and frame reproduces the field-resolved result is checked in the test suite ("THG: real vs envelope" in `test/test_freespace.jl` and `test/test_interface.jl`, and "BBO SHG: real vs envelope" in `test/test_freespace.jl`).

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
