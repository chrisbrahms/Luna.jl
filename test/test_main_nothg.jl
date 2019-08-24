import Luna
import Luna: Grid, Maths, Capillary, PhysData, Nonlinear, Ionisation, Modes, Output, Stats
import Logging
import NumericalIntegration: integrate, SimpsonEven
Logging.disable_logging(Logging.BelowMinLevel)

import DSP.Unwrap: unwrap

import PyPlot:pygui, plt

a = 13e-6
gas = :Ar
pres = 5

τ = 30e-15
λ0 = 800e-9

grid = Grid.RealGrid(15e-2, 800e-9, (160e-9, 3000e-9), 1e-12)

m = Capillary.MarcatilliMode(a, gas, pres)

energyfun = Modes.energy_mode_avg(m)

function gausspulse(t)
    It = Maths.gauss(t, fwhm=τ)
    ω0 = 2π*PhysData.c/λ0
    Et = @. sqrt(It)*cos(ω0*t)
end

β1const = Capillary.dispersion(m, 1; λ=λ0)
βconst = zero(grid.ω)
βconst[2:end] = Capillary.β(m, grid.ω[2:end])
βconst[1] = 1
βfun(ω, m, n, z) = βconst
frame_vel(z) = 1/β1const
αfun(ω, m, n, z) = log(10)/10 * 2

dens0 = PhysData.density(gas, pres)
densityfun(z) = dens0

normfun = Modes.norm_mode_average(grid.ω, βfun)

transform = Modes.trans_mode_avg(grid)

ionpot = PhysData.ionisation_potential(gas)
ionrate = Ionisation.ionrate_fun!_ADK(ionpot)

responses = (Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas),length(grid.to)),
             Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot))

in1 = (func=gausspulse, energy=1e-6, m=1, n=1)
inputs = (in1, )

x = Array{Float64}(undef, length(grid.t))
FT = FFTW.plan_rfft(x, 1, flags=FFTW.MEASURE)

statsfun = Stats.collect_stats((Stats.ω0(grid), ))
output = Output.MemoryOutput(0, grid.zmax, 201, (length(grid.ω),), statsfun)

linop = Luna.make_linop(grid, βfun, αfun, frame_vel)
Luna.run(grid, linop, normfun, energyfun, densityfun, inputs,
         responses, transform, FT, output)

ω = grid.ω
t = grid.t

zout = output.data["z"]
Eout = output.data["Eω"]

Etout = FFTW.irfft(Eout, length(grid.t), 1)

Ilog = log10.(Maths.normbymax(abs2.(Eout)))

idcs = @. (t < 30e-15) & (t >-30e-15)
to, Eto = Maths.oversample(t[idcs], Etout[idcs, :], factor=16)
It = abs2.(Maths.hilbert(Eto))
Itlog = log10.(Maths.normbymax(It))
zpeak = argmax(dropdims(maximum(It, dims=1), dims=1))

Et = Maths.hilbert(Etout)
energy = zeros(length(zout))
for ii = 1:size(Etout, 2)
    energy[ii] = energyfun(t, Etout[:, ii], 1, 1)
end

pygui(true)
plt.figure()
plt.pcolormesh(ω./2π.*1e-15, zout, transpose(Ilog))
plt.clim(-6, 0)
plt.colorbar()

plt.figure()
plt.pcolormesh(to*1e15, zout, transpose(It))
plt.colorbar()
plt.xlim(-30, 30)

plt.figure()
plt.plot(zout.*1e2, energy.*1e6)
plt.xlabel("Distance [cm]")
plt.ylabel("Energy [μJ]")

plt.figure()
plt.plot(to*1e15, Eto[:, 121])
plt.xlim(-20, 20)