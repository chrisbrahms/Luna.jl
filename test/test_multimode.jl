import Test: @test, @testset, @test_throws

@testset "Radial" begin
    # mode average and radial integral for single mode and only Kerr should be identical
    using Luna
    import LinearAlgebra: norm
    a = 13e-6
    gas = :Ar
    pres = 5
    τ = 30e-15
    λ0 = 800e-9
    grid = Grid.RealGrid(800e-9, (160e-9, 3000e-9), 1e-12)
    m = Capillary.MarcatiliMode(a, gas, pres, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    energyfun, energyfunω = Fields.energyfuncs(grid)

    dens0 = PhysData.density(gas, pres)
    densityfun(z) = dens0
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun!, frame_vel, αfun = LinearOps.make_const_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                               Stats.ω0(grid),
                               Stats.energy(grid, energyfunω))
    output = Output.MemoryOutput(0, 5e-2, 201, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output, status_period=5, zmax=5e-2)

    modes = (
         Capillary.MarcatiliMode(a, gas, pres, n=1, m=1, kind=:HE, ϕ=0.0, loss=false),
    )
    energyfun, energyfunω = Fields.energyfuncs(grid)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs,
                                modes, :y; full=false)
    outputr = Output.MemoryOutput(0, 5e-2, 201, statsfun)
    linop = LinearOps.make_const_linop(grid, modes, λ0)
    Luna.run(Eω, grid, linop, transform, FT, outputr, status_period=10, zmax=5e-2)

    Iω = abs2.(output.data["Eω"])
    Iωr = abs2.(dropdims(outputr.data["Eω"], dims=2))
    @test norm(Iω - Iωr)/norm(Iω) < 0.003
end

@testset "Full" begin
    # mode average and full integral for single mode and only Kerr should be identical
    using Luna
    import LinearAlgebra: norm
    a = 13e-6
    gas = :Ar
    pres = 5
    τ = 30e-15
    λ0 = 800e-9
    grid = Grid.RealGrid(800e-9, (160e-9, 3000e-9), 1e-12)
    m = Capillary.MarcatiliMode(a, gas, pres, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    energyfun, energyfunω = Fields.energyfuncs(grid)

    dens0 = PhysData.density(gas, pres)
    densityfun = let dens0=PhysData.density(gas, pres)
        z -> dens0
    end
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun!, frame_vel, αfun = LinearOps.make_const_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                               Stats.ω0(grid),
                               Stats.energy(grid, energyfunω))
    output = Output.MemoryOutput(0, 5e-2, 201, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output, status_period=5, zmax=5e-2)

    modes = (
         Capillary.MarcatiliMode(a, gas, pres, n=1, m=1, kind=:HE, ϕ=0.0, loss=false),
    )
    energyfun, energyfunω = Fields.energyfuncs(grid)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs,
                                modes, :y; full=true)
    outputf = Output.MemoryOutput(0, 5e-2, 201, statsfun)
    linop = LinearOps.make_const_linop(grid, modes, λ0)
    Luna.run(Eω, grid, linop, transform, FT, outputf, status_period=10, zmax=5e-2)

    Iω = abs2.(output.data["Eω"])
    Iωf = abs2.(dropdims(outputf.data["Eω"], dims=2))
    @test norm(Iω - Iωf)/norm(Iω) < 0.003
end

@testset "Full, LP11 vs TM01" begin
    # TM01 and LP11 (which is TM01+HE21) have the same dispersion and a ratio of effective
    # areas of 2/3 -- mode averaged propagation in TM01 with 2/3x Aeff should be identical
    # to full modal propagation in LP11, as long as only Kerr is considered.
    using Luna
    import LinearAlgebra: norm
    a = 13e-6
    gas = :Ar
    pres = 5
    τ = 30e-15
    λ0 = 800e-9
    energy = 1e-6
    grid = Grid.RealGrid(800e-9, (160e-9, 3000e-9), 0.4e-12)
    m = Capillary.MarcatiliMode(a, gas, pres; kind=:TM, n=0, m=1, loss=false)
    aeff(z) = 2/3*Modes.Aeff(m, z=z) # Aeff of LP11 is 2/3x smaller than that of TM01
    energyfun, energyfunω = Fields.energyfuncs(grid)

    dens0 = PhysData.density(gas, pres)
    densityfun = let dens0=PhysData.density(gas, pres)
        z -> dens0
    end
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun!, frame_vel, αfun = LinearOps.make_const_linop(grid, m, λ0)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=energy)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                               Stats.ω0(grid),
                               Stats.energy(grid, energyfunω))
    output = Output.MemoryOutput(0, 5e-2, 201, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output, status_period=5, zmax=5e-2)

    # vertical LP11 double-lobe
    modes = (
        Capillary.MarcatiliMode(a, gas, pres, n=0, m=1, kind=:TM, ϕ=0.0, loss=false),
        Capillary.MarcatiliMode(a, gas, pres, n=2, m=1, kind=:HE, ϕ=-π/4, loss=false),
    )

    field = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=energy/2)
    inputs = ((mode=1, fields=(field,)), (mode=2, fields=(field,)))
    # inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=energy)
    Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs,
                                   modes, :xy; full=true)
    outputf = Output.MemoryOutput(0, 5e-2, 201, statsfun)
    linop = LinearOps.make_const_linop(grid, modes, λ0)
    Luna.run(Eω, grid, linop, transform, FT, outputf, status_period=10, zmax=5e-2)

    Iω = abs2.(output["Eω"])
    Iωf = dropdims(sum(abs2.(outputf["Eω"]); dims=2), dims=2)
    @test norm(Iω - Iωf)/norm(Iω) < 0.0006
end

@testset "FieldInputs" begin
    using Luna
    import LinearAlgebra: norm
    a = 13e-6
    gas = :Ar
    pres = 5
    τ = 30e-15
    λ0 = 800e-9
    grid = Grid.RealGrid(800e-9, (160e-9, 3000e-9), 1e-12)
    modes = (
         Capillary.MarcatiliMode(a, gas, pres, n=1, m=1, kind=:HE, ϕ=0.0, loss=false),
         Capillary.MarcatiliMode(a, gas, pres, n=1, m=2, kind=:HE, ϕ=0.0, loss=false)
    )
    aeff(z) = Modes.Aeff(m, z=z)
    energyfun, energyfunω = Fields.energyfuncs(grid)
    dens0 = PhysData.density(gas, pres)
    densityfun(z) = dens0
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω_single, transform, FT = Luna.setup(grid, densityfun, responses, inputs,
                                modes, :y; full=false)
    inputs = (Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6),)
    Eω_tuple, transform, FT = Luna.setup(grid, densityfun, responses, inputs,
                                modes, :y; full=false)
    inputs = ((mode=1, fields=(Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6),)),)
    Eω_tuple_of_namedtuples, transform, FT = Luna.setup(grid, densityfun, responses, inputs,
                                modes, :y; full=false)
    @test Eω_single ≈ Eω_tuple
    @test Eω_single ≈ Eω_tuple_of_namedtuples
end

@testset "Nonlinear coupling" begin
    using Luna
    import LinearAlgebra: norm
    a = 125e-6
    gas = :Ar
    pres = 0.167
    flength = 3

    τfwhm = 10e-15
    λ0 = 800e-9
    energy = 150e-6

    # HE11 and LP11 (consisting of HE21 + TM01)
    modes = (
        Capillary.MarcatiliMode(a, gas, pres, n=0, m=1, kind=:TM, ϕ=0.0),
        Capillary.MarcatiliMode(a, gas, pres, n=2, m=1, kind=:HE, ϕ=float(-π/4)),
        Capillary.MarcatiliMode(a, gas, pres, n=1, m=1, kind=:HE, ϕ=float(0.0)),
    )

    grid = Grid.RealGrid(λ0, (160e-9, 3000e-9), 0.4e-12)

    densityfun = let dens0=PhysData.density(gas, pres)
        z -> dens0
    end
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)

    field = Fields.GaussField(;λ0, τfwhm, energy=energy/2)
    inputs = ((mode=1, fields=(field,)), (mode=2, fields=(field,)))

    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, modes, :xy; full=true)
    nl = similar(Eω)

    transform(nl, Eω, 0.0);
    # nonlinear polarisation in LP11
    Inl12 = dropdims(sum(abs2.(nl[:, 1:2]); dims=2); dims=2)
    # nonlinear polarisation in HE11
    Inl3 = abs2.(nl[:, 3])
    # LP11 should not couple nonlinearly to HE11
    @test norm(Inl3)/norm(Inl12) < 1e-32
end
@testset "Rectangular transverse integral" begin
    using Luna
    import Luna: NonlinearRHS, RectModes
    import Luna.PhysData: ε_0, μ_0

    #= For a single mode the transverse integral of the Kerr polarisation is analytic.
       The field at a transverse point is Eₘ ê/√N, with ê the mode profile and N its
       normalisation, so the Kerr polarisation is ρ ε₀ γ₃ (Eₘ ê/√N)³ and its projection
       back onto the mode is

           ∫ dA (ê/√N) ρ ε₀ γ₃ (Eₘ ê/√N)³ = (∫ê⁴ dA / N²) ρ ε₀ γ₃ Eₘ³ .

       `NonlinearRHS.Erω_to_Prω!` returns the same windowed, transformed and normalised
       polarisation at a single transverse point, so with `norm!` the identity the ratio
       of the transform's output to `Erω_to_Prω!` at the centre of the guide, where ê = 1,
       is ∫ê⁴ dA / N² · N^(3/2) = ∫ê⁴ dA / √N. Everything else -- the time window, the
       oversampled transforms, the spectral window, the density and γ₃ -- cancels.

       For the fundamental mode of a rectangular guide of half-widths a and b,
       ∫ê⁴ dA = 9ab/16 and N = ½√(ε₀/μ₀)ab.

       The `a > b` guide is the one this measures. The Cartesian in-domain test of the
       adaptive driver (`NonlinearRHS._points!`) used to read `x1 >= ul[2]` rather than
       `x2 >= ul[2]`, which treats every point with b ≤ x < a as outside the guide: for
       a = 2.5b that drops 21 % of ∫ê⁴ dA and the ratio comes out 8.3 % low. =#
    "∫ê⁴dA/√N for the fundamental mode of a rectangular guide of half-widths `a`, `b`."
    kerrfactor(a, b) = (9*a*b/16)/sqrt(0.5*sqrt(ε_0/μ_0)*a*b)

    λ0 = 800e-9
    gas = :Ar
    grid = Grid.RealGrid(λ0, (400e-9, 2e-6), 100e-15)
    densityfun = let dens=PhysData.density(gas, 1.0)
        z -> dens
    end
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=10e-15, energy=1e-6)

    "The nonlinear polarisation of one `RectMode` at z = 0, and the transform that made it."
    function rectnl(a, b; kwargs...)
        modes = (RectModes.RectMode(a, b, gas, 1.0, :Ag),)
        Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs, modes, :x;
                                       full=true, norm! = identity, kwargs...)
        nl = similar(Eω)
        transform(nl, Eω, 0.0)
        nl, transform
    end

    for (a, b) in ((100e-6, 40e-6), (40e-6, 100e-6), (60e-6, 60e-6))
        nl, transform = rectnl(a, b; rtol=1e-6, mfcn=20_000)
        # the same quantity at the centre of the guide, where the mode profile is 1
        Prω = copy(NonlinearRHS.Erω_to_Prω!(transform, (0.0, 0.0)))
        # only where the polarisation is not window taper or numerical dust
        idcs = findall(x -> abs(x) > 1e-3*maximum(abs, Prω), Prω[:, 1])
        @test length(idcs) > 10
        @test all(isapprox.(nl[idcs, 1] ./ Prω[idcs, 1], kerrfactor(a, b), rtol=1e-7))

        # the fixed quadrature rule integrates the same domain
        nlf, _ = rectnl(a, b; modal_integral=:fixed)
        @test all(isapprox.(nlf[idcs, 1] ./ Prω[idcs, 1], kerrfactor(a, b), rtol=1e-10))
    end
end
