# 1. mode-averaged field, Kerr + PPT plasma (default prop_capillary physics for Ar)
modeavg() = prop_capillary(125e-6, 0.5, :Ar, 1.0; λ0=800e-9, τfwhm=30e-15, energy=150e-6,
                           λlims=(150e-9, 4e-6), trange=1e-12, shotnoise=false,
                           status_period=1e9)
# 2. four-mode fixed-rule multimode, Kerr + plasma
modal() = prop_capillary(125e-6, 0.05, :Ar, 1.0; λ0=800e-9, τfwhm=30e-15, energy=150e-6,
                         λlims=(150e-9, 4e-6), trange=0.5e-12, shotnoise=false, modes=4,
                         modal_integral=:fixed, status_period=1e9)
# 3. radial free space, Kerr + plasma
function radial()
    gas = :Ar; pres = 1.2; λ0 = 800e-9; L = 0.1
    grid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    q = Grid.RadialGrid(4e-3, 256)
    densityfun = let d = PhysData.density(gas, pres); z -> d; end
    ionrate = Ionisation.IonRatePPTCached(gas, λ0)
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
                 Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, PhysData.ionisation_potential(gas)))
    linop = LinearOps.make_const_linop(grid, q, PhysData.ref_index_fun(gas, pres))
    normfun = NonlinearRHS.const_norm_radial(grid, q, PhysData.ref_index_fun(gas, pres))
    inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=20e-15, energy=20e-6, w0=200e-6, propz=-0.02)
    Eω, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)
    statsfun = Stats.default(grid, Eω, linop, transform; gas=gas)
    output = Output.MemoryOutput(0, L, 11, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=L, status_period=1e9)
    output
end
# 4. mode-averaged, long window (nt ≈ 2^16)
modeavg_long() = prop_capillary(125e-6, 0.2, :Ar, 1.0; λ0=800e-9, τfwhm=30e-15, energy=150e-6,
                                λlims=(150e-9, 4e-6), trange=8e-12, shotnoise=false,
                                status_period=1e9)
# 5. radial, 1024 points
function radial_big()
    gas = :Ar; pres = 1.2; λ0 = 800e-9; L = 0.05
    grid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    q = Grid.RadialGrid(4e-3, 1024)
    densityfun = let d = PhysData.density(gas, pres); z -> d; end
    ionrate = Ionisation.IonRatePPTCached(gas, λ0)
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
                 Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, PhysData.ionisation_potential(gas)))
    linop = LinearOps.make_const_linop(grid, q, PhysData.ref_index_fun(gas, pres))
    normfun = NonlinearRHS.const_norm_radial(grid, q, PhysData.ref_index_fun(gas, pres))
    inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=20e-15, energy=20e-6, w0=200e-6, propz=-0.02)
    Eω, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)
    statsfun = Stats.default(grid, Eω, linop, transform; gas=gas)
    output = Output.MemoryOutput(0, L, 11, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=L, status_period=1e9)
    output
end
