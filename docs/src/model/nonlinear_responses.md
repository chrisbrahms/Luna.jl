# Nonlinear responses

## Kerr effect
Real fields vector
```math
\begin{align*}
\mathbf{P}_\mathrm{Kerr}(t, x, y, z) &= \varepsilon_0 \chi^{(3)} \left(\mathbf{E}(t, x, y, z) \cdot \mathbf{E}(t, x, y, z)\right) \mathbf{E}(t, x, y, z) \\[1em]
&= \varepsilon_0 \chi^{(3)} \begin{pmatrix}
           \left(E_x^2 + E_y^2\right) E_x \\[0.5em]
           \left(E_x^2 + E_y^2\right) E_y
\end{pmatrix}
\end{align*}
```

Real fields scalar
```math
P_\mathrm{Kerr}(t, x, y, z) = \varepsilon_0 \chi^{(3)} \left(E(t, x, y, z)\right) ^3
```

Envelope fields scalar
```math
P_\mathrm{Kerr}(t, x, y, z) = \frac{3}{4}\varepsilon_0 \chi^{(3)} \left\vert\mathcal{E}(t, x, y, z)\right\vert^2\mathcal{E}(t, x, y, z)
```

Envelope fields vector
```math
\mathbf{P}_\mathrm{Kerr}(t, x, y, z) = \frac{3}{4}\varepsilon_0 \chi^{(3)} \begin{pmatrix}
           (\mathcal{E}_x^2 + \frac{2}{3}\mathcal{E}_y^2)\mathcal{E}_x + \frac{1}{3}\mathcal{E}_x^*\mathcal{E}_y^2 \\[0.5em]
           (\mathcal{E}_y^2 + \frac{2}{3}\mathcal{E}_x^2)\mathcal{E}_y + \frac{1}{3}\mathcal{E}_y^*\mathcal{E}_x^2
\end{pmatrix}

```
## Photoionisation & plasma
```math
P_\mathrm{ion}\!\left(t,x, y,z\right) = I_p\int_{-\infty}^t \!\!\mathrm{d}t'\frac{\partial_{t'} \rho_\mathrm{e}(t', x, y,z)}{E\!\left(t',x, y,z\right)} + \frac{e^2}{m_\mathrm{e}}\int_{-\infty}^{t}\!\!\mathrm{d}t'\int_{-\infty}^{t'} \!\!\mathrm{d}t'' \rho_\mathrm{e}(t'',x, y,z)E\left(t'',x, y,z\right)
```

```math
\rho_\mathrm{e}(t, x, y,z) = \rho(z)\left(1 - \exp \left\{-\int_{-\infty}^t\!\!\mathrm{d}t'w\!\left(\left\vert E\!\left(t',x, y,z\right)\right\vert\right)\right\}\right)
```
## Raman response