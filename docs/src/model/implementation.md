# Implementation details

## The interaction picture, a.k.a. pre-conditioned Runge-Kutta
The basic form of the UPPE in Luna,
```math
\partial_z E(\omega, \mathbf{k}_\perp, z) = \mathcal{L}(\omega, \mathbf{k}_\perp, z)E(\omega, \mathbf{k}_\perp, z) + \frac{i\omega}{N_{\mathrm{nl}}} P_{\mathrm{nl}}(\omega, \mathbf{k}_\perp, z)\,,
```
is a general partial differential equation. However, by taking into account the discretisation in time and space, this turns into a system of coupled first-order ordinary differential equations of the general form
```math
\frac{\partial \overrightarrow{y}}{\partial t} = f(\overrightarrow{y}, t)\,,
```
where $\overrightarrow{y}$ is the (vector of) function(s) to be solved for, and $t$ is the independent variable. Note the difference in notation: in Luna, ``\overrightarrow{y}`` is the field ``E(\omega, \mathbf{k}_\perp, z)`` (evaluated at discrete values of ``\omega`` and ``\mathbf{k}_\perp``) and the "time" ``t`` is actually the distance ``z``. This most general form can be directly solved numerically with various methods. However, that would not make direct use of the fact that *without the nonlinear part*, the UPPE has an analytical solution:
```math
E(\omega, \mathbf{k}_\perp, z_2) = \mathrm{e}^{\int_{z_1}^{z_2}\mathcal{L}(\omega, \mathbf{k}_\perp, z')\mathrm{d}z}E(\omega, \mathbf{k}_\perp, z_1)\,.
```
We therefore re-write the general form (and drop the vector notation for readability) as
```math
\frac{\partial y}{\partial t} = \mathcal{L}(t)y(t) + f(y, t)\,,
```
We now define the "preconditioned" or "interaction-picture" function ``\bar{y}``:
```math
\bar{y}(t) = \exp\left[-\int_0^t\mathcal{L}(t')\mathrm{d} t'\right]y(t) = \mathrm{e}^{-L(t)} y(t)\,,
```
where we have also defined the linear propagator $L(t)$ for brevity:
```math
    L(t) \equiv \int_0^t\mathcal{L}(t')\mathrm{d} t'\,.
```
The left-hand side of the differential equation can thus be expressed as
```math
\begin{align*}
    \frac{\partial y}{\partial t} &= \mathrm{e}^{L(t)}\left[\frac{\partial \bar{y}}{\partial t} + \frac{\partial L}{\partial t}\bar{y}(t)\right]\\
    &= \mathrm{e}^{L(t)}\left[\frac{\partial \bar{y}}{\partial t} + \mathcal{L}(t)\bar{y}(t)\right]\,.
\end{align*}
```
Also substituting for $y(t)$ using $\bar{y}(t)$ on the right-hand side, the whole equation then becomes
```math
    \mathrm{e}^{L(t)}\left[\frac{\partial \bar{y}}{\partial t} + \mathcal{L}(t)\bar{y}(t)\right] = \mathrm{e}^{L(t)}\mathcal{L}(t)\bar{y}(t) + f(y, t)\,.
```
The two terms proportional to $\mathcal{L}(t)$ cancel, leaving us with
```math
    \frac{\partial \bar{y}}{\partial t} = \mathrm{e}^{-L(t)}f(y, t)\,.
```
By further defining the pre-conditioned function $\bar{f}(\bar{y}, t)$ as
```math
    \bar{f}(\bar{y}, t) = \mathrm{e}^{-L(t)}f\left(\mathrm{e}^{L(t)}\bar{y}, t\right)\,,
```
we recover a generic first-order ODE which now includes the pre-conditioner:
```math
    \frac{\partial \bar{y}}{\partial t} = \bar{f}(\bar{y}, t)\,.
```