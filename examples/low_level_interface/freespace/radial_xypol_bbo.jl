using Luna
import FFTW
import Luna: Hankel
import PyPlot: plt
import NumericalIntegration: integrate

λ0 = 1030e-9
τfwhm = 250e-15
w0 = 50e-6
energy = 1e-6

thickness = 100e-6
material = :BBO

R = 4*w0
N = 2^6

grid = Grid.RealGrid(thickness, λ0, (300e-9, 4e-6), 500e-15)
q = Hankel.QDHT(R, N, dim=3)

θ = deg2rad(23.3717)
ϕ = deg2rad(30)


responses = (Nonlinear.Chi2Field(θ, ϕ, PhysData.χ2(material)),
             Nonlinear.Kerr_field(PhysData.χ3(material)))

nfun = PhysData.ref_index_fun_uniax(material)
function nfunreal(λ; z=0)
    # n_e, n_o
    real(nfun(λ, θ)), real(nfun(λ, 0))
end
linop = LinearOps.make_const_linop(grid, q, nfunreal)

normfun = NonlinearRHS.const_norm_radial(grid, q, nfunreal)
densityfun = z -> 1

inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0)
# inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0)
Eω, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)

##
output = Output.MemoryOutput(0, grid.zmax, 101)
Luna.run(Eω, grid, linop, transform, FT, output; init_dz=1e-6)

##
z = output["z"]
Eωk = output["Eω"] # (ω, pol, k, z)

ωprefac = 2π*PhysData.c*PhysData.ε_0/2 * 2π/(grid.ω[end]^2)

Eωr = q \ Eωk # (ω, pol, r, z)
Etr = FFTW.irfft(Eωr, 2*(length(grid.ω)-1), 1) # (t, pol, r, z)
Etr = Maths.hilbert(Etr)
Iωr = abs2.(Eωr) # (ω, pol, r, z)
Itr = 0.5*PhysData.c*PhysData.ε_0*abs2.(Etr) # (t, pol, r, z)

Irxy = dropdims(sum(Iωr; dims=1); dims=1) # (pol, r, z)
Iωxy = dropdims(Hankel.integrateR(Iωr, q; dim=3); dims=3)*ωprefac # (ω, pol, z)
Ir = dropdims(sum(Iωr; dims=(1, 2)); dims=(1, 2)) # (r, z)

Itxy = dropdims(Hankel.integrateR(Itr, q; dim=3); dims=3) # (t, pol, z)

Eω0 = dropdims(Hankel.onaxis(Eωk, q); dims=q.dim) # (ω, pol, z)
Et0 = FFTW.irfft(Eω0, 2*(length(grid.ω)-1), 1) # (t, pol, z)
Et0 = Maths.hilbert(Et0)
It0 = 0.5*PhysData.c*PhysData.ε_0*abs2.(Et0)

et, eω = Fields.energyfuncs(grid, q)

energy_x = eω(Eωk[:, 1, :, end])
energy_y = eω(Eωk[:, 2, :, end])

ω = grid.ω

##
I1 = 0.94*energy/τfwhm / (π*w0^2)
d31 = 0.16e-12
d22 = -2.3e-12
deff = d31*sin(θ) - d22*cos(θ)*sin(3ϕ)

ω3 = PhysData.wlfreq(λ0/2)
n1 = real(nfun(λ0, 0))
n2 = n1
n3 = real(nfun(λ0/2, θ))

χ2eff = 2deff

I3 = 2*χ2eff^2*ω3^2/(n1*n2*n3*PhysData.ε_0 * PhysData.c^3) * I1^2 * z.^2

##
plt.figure()
plt.plot(z*1e6, I3*1e-4; label="Calculated")
plt.plot(z*1e6, dropdims(maximum(It0[:, 1, :]; dims=1); dims=1)*1e-4; label="Simulated")
# plt.plot(z*1e6, dropdims(maximum(Itr[:, 1, 1, :]; dims=1); dims=1)*1e-4; label="Simulated")

##
plt.figure()
plt.pcolormesh(z*1e3, q.r*1e6, Ir)
plt.xlabel("Distance (mm)")
plt.ylabel("radius (μm)")

##
plt.figure()
plt.subplot(1, 2, 1)
plt.pcolormesh(z*1e3, q.r*1e6, Irxy[1, :, :])
plt.subplot(1, 2, 2)
plt.pcolormesh(z*1e3, q.r*1e6, Irxy[2, :, :])
plt.xlabel("Distance (mm)")
plt.ylabel("radius (μm)")
plt.suptitle("Frequency domain")
##
plt.figure()
plt.subplot(1, 2, 1)
plt.pcolormesh(z*1e3, grid.t*1e15, Itxy[:, 1, :])
plt.title("X Pol")
plt.xlabel("Distance (mm)")
plt.ylabel("Time (fs)")
plt.subplot(1, 2, 2)
plt.pcolormesh(z*1e3, grid.t*1e15, Itxy[:, 2, :])
plt.title("Y Pol")
plt.xlabel("Distance (mm)")
plt.suptitle("Time domain")

##
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(grid.t*1e15, Itxy[:, 1, 1]; label="input X")
plt.plot(grid.t*1e15, Itxy[:, 2, 1]; label="input Y")
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(grid.t*1e15, Itxy[:, 1, end]; label="output X")
plt.plot(grid.t*1e15, Itxy[:, 2, end]; label="output Y")
plt.xlabel("Time (fs)")
plt.legend()

##
mm = 2π*maximum(Iωxy)*1e12
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogy(ω*1e-12/2π, 2π*1e12*Iωxy[:, 1, 1]; label="input X")
plt.semilogy(ω*1e-12/2π, 2π*1e12*Iωxy[:, 2, 1]; label="input Y")
plt.ylim(1e-5mm, 5mm)
plt.xlim(200, 800)
plt.legend()
plt.subplot(2, 1, 2)
plt.semilogy(ω*1e-12/2π, 2π*1e12*Iωxy[:, 1, end]; label="output X")
plt.semilogy(ω*1e-12/2π, 2π*1e12*Iωxy[:, 2, end]; label="output Y")
plt.semilogy(ω*1e-12/2π, 2π*1e12*Iωxy[:, 2, 1]; c="C1", linestyle="--", label="input Y")
plt.ylim(1e-5mm, 5mm)
plt.xlim(200, 800)
plt.xlabel("Frequency (THz)")
plt.ylabel("SED (J/Hz)")
plt.legend()