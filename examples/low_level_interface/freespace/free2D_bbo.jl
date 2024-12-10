using Luna
import FFTW
import Luna: Hankel
import PyPlot: plt
import NumericalIntegration: integrate

λ0 = 1030e-9
τfwhm = 250e-15
w0 = 50e-6
energy = 1e-6 / 8

thickness = 1000e-6
material = :BBO

R = 4*w0
N = 2^6

grid = Grid.RealGrid(thickness, λ0, (300e-9, 4e-6), 500e-15)
xgrid = Grid.Free2DGrid(R, N)

θ = deg2rad(23.3717)
ϕ = deg2rad(30)


responses = (Nonlinear.Chi2Field(θ, ϕ, PhysData.χ2(material)),
             Nonlinear.Kerr_field(PhysData.χ3(material)))

nfun = PhysData.ref_index_fun_uniax(material)
function nfunreal(λ; z=0)
    # n_e, n_o
    real(nfun(λ, θ)), real(nfun(λ, 0))
end
nfunx(λ, δθ; z=0) = real(nfun(λ, θ+δθ))
nfuny(λ; z=0) = real(nfun(λ, 0))
linop = LinearOps.make_const_linop(grid, xgrid, nfunx, nfuny)

normfun = NonlinearRHS.const_norm_free2D(grid, xgrid, nfunx, nfuny)
densityfun = z -> 1
##
inputs = Fields.GaussGaussField(;λ0, τfwhm, energy=energy/(sqrt(π/2)*w0), w0)
# inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0)
Eω, transform, FT = Luna.setup(grid, xgrid, densityfun, normfun, responses, inputs)

##
output = Output.MemoryOutput(0, grid.zmax, 101)
Luna.run(Eω, grid, linop, transform, FT, output; init_dz=1e-6)

##
z = output["z"]
Eωk = output["Eω"] # (ω, pol, k, z)
x = xgrid.x

ωprefac = 2π*PhysData.c*PhysData.ε_0/2 * 2π/(grid.ω[end]^2)

Eωr = FFTW.ifft(Eωk, 3) # (ω, pol, x, z)
Etr = FFTW.irfft(Eωr, 2*(length(grid.ω)-1), 1) # (t, pol, x, z)
Etr = Maths.hilbert(Etr)
Iωr = abs2.(Eωr) # (ω, pol, x, z)
Itr = 0.5*PhysData.c*PhysData.ε_0*abs2.(Etr) # (t, pol, x, z)

Irxy = dropdims(sum(Iωr; dims=1); dims=1) # (pol, x, z)
# Iωxy = dropdims(integrate(x, Iωr; dims=3); dims=3)*ωprefac # (ω, pol, z)
Ir = dropdims(sum(Iωr; dims=(1, 2)); dims=(1, 2)) # (x, z)

# Itxy = dropdims(integrate(x, Itr; dims=3); dims=3) # (t, pol, z)

Eω0 = Eωr[:, :, length(x)÷2+1, :]
Et0 = FFTW.irfft(Eω0, 2*(length(grid.ω)-1), 1) # (t, pol, z)
Et0 = Maths.hilbert(Et0)
It0 = 0.5*PhysData.c*PhysData.ε_0*abs2.(Et0)

et, eω = Fields.energyfuncs(grid, xgrid)

energy_x = eω(Eωk[:, 1, :, 1])*sqrt(π/2)*w0
energy_y = eω(Eωk[:, 2, :, 1])*sqrt(π/2)*w0

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
# plt.plot(z*1e6, I3*1e-4; label="Calculated")
plt.plot(z*1e6, dropdims(maximum(It0[:, 1, :]; dims=1); dims=1)*1e-4; label="Simulated")
# plt.plot(z*1e6, dropdims(maximum(Itr[:, 1, 1, :]; dims=1); dims=1)*1e-4; label="Simulated")

##
plt.figure()
plt.pcolormesh(z*1e3, x*1e6, Ir)
plt.xlabel("Distance (mm)")
plt.ylabel("X (μm)")

##
plt.figure()
plt.subplot(1, 2, 1)
plt.pcolormesh(z*1e3, x*1e6, Irxy[1, :, :])
plt.subplot(1, 2, 2)
plt.pcolormesh(z*1e3, x*1e6, Irxy[2, :, :])
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

##
plt.figure()
# plt.plot(grid.ω, imag(linop[:, 1, 32]))
plt.plot(grid.ω, imag(linop[:, 2, 32]))
plt.axvline(PhysData.wlfreq(λ0))
plt.axvline(PhysData.wlfreq(λ0/2))
# plt.ylim(-1e6, 1e6)