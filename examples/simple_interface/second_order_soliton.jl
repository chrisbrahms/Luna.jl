using Luna
import PyPlot: plt
import FFTW

γ = 0.1
β2 = -1e-26
N = 2
τ0 = 100e-15
τfwhm = (2*log(1 + sqrt(2)))*τ0
fr = 0.18
P0 = N^2*abs(β2)/((1-fr)*γ*τ0^2)
Ld = τ0^2/abs(β2)
flength = π*Ld/2

λ0 = 1030e-9
λlims = [450e-9, 8000e-9]
trange = 8e-12

Eω, grid, linop, transform, FT, soliton = Interface.prop_gnlse_args(γ, flength, [0.0, 0.0, β2]; λ0, τfwhm, power=P0, pulseshape=:sech, λlims, trange,
                    raman=false, shock=false, fr, shotnoise=false)

Luna.run(Eω, grid, linop, transform, FT, soliton)

##
z = soliton["z"]
f, If = Processing.getIω(soliton, :f)

t, Et = Processing.getEt(soliton)

Pt = abs2.(Et)

f .-= PhysData.wlfreq(λ0)/2π

##
fig = plt.figure()
fig.set_size_inches(5.5, 2.3)
gs = fig.add_gridspec(1, 2, wspace=0.75, bottom=0.2, top=0.95, left=0.13, right=0.87)
gss = gs[1].subgridspec(1, 2, width_ratios=[1, 0.05], wspace=0.05)
ax = fig.add_subplot(gss[1, 1])
ax.pcolormesh(1e15t, z/Ld, Pt'/maximum(Pt); rasterized=true)
ax.set_xlim(-220, 220)
ax.set_xlabel("\$\\tau\$ (fs)", labelpad=0)
ax.set_ylabel("\$z/L_D\$", labelpad=0)
cax = fig.add_subplot(gss[1, 2])
fig.colorbar(ax.collections[1], cax=cax, orientation="vertical")
cax.set_ylabel("\$\\left\\vert A(z, \\tau)\\right\\vert^2\$", labelpad=0)

gss = gs[2].subgridspec(1, 2, width_ratios=[1, 0.05], wspace=0.05)
ax = fig.add_subplot(gss[1, 1])
img = ax.pcolormesh(1e-12f, z/Ld, If'/maximum(If); rasterized=true, norm=plt.matplotlib.colors.LogNorm())
img.set_clim(1e-4, 1)
ax.set_xlim(-25, 25)
ax.set_xlabel("\$\\Omega/2\\pi\$ (THz)", labelpad=0)
ax.set_ylabel("\$z/L_D\$", labelpad=0)
cax = fig.add_subplot(gss[1, 2])
fig.colorbar(ax.collections[1], cax=cax, orientation="vertical")
cax.set_ylabel("\$\\left\\vert \\tilde{A}(z, \\Omega)\\right\\vert^2\$", labelpad=0)
