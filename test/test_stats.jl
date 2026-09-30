using Luna
import Test: @test, @testset
import FFTW
import HCubature: hquadrature
import SpecialFunctions: besselj
import FunctionZeros: besselj_zero

@testset "On-axis intensity statistics" begin
    # Manually normalise the field distribution of HE11 to find scaling factor between
    # power and intensity
    a = 13e-6
    unm = besselj_zero(0, 1)
    E(r) = besselj(0, unm*r/a)
    norm, err = hquadrature(0, a) do r
        2π*r * abs2(E(r))
    end

    energy = 1e-6
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

    statsfun = Stats.default(grid, Eω, m, linop, transform; gas=gas, onaxis=true)
    output = Output.MemoryOutput(0, 15e-2, 201, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output, status_period=5, zmax=15e-2)

    @test all(output["stats"]["peakintensity"] .≈ output["stats"]["peakpower"]/norm)
end

#= The fixed transverse quadrature rule's embedded error estimate, as
   `Stats.transverse_integral_error` records it per accepted step. What it should record
   is `NonlinearRHS.integral_error!` of the same state, which is computed here directly
   from the transform. =#
@testset "the transverse integral error statistic" begin
    λ0 = 800e-9
    grid = Grid.RealGrid(λ0, (200e-9, 3000e-9), 400e-15)
    modes = Tuple(Capillary.MarcatiliMode(75e-6, :Ar, 0.1; n=1, m=mi, loss=false)
                  for mi in 1:2)
    ρ = PhysData.density(:Ar, 0.1)
    dens = z -> ρ
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:Ar)),)
    inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=50e-6)
    linop = LinearOps.make_const_linop(grid, modes, λ0)

    @testset "kronrod=$kronrod" for (kronrod, nr) in ((true, 33), (false, 32))
        Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, modes, :y;
                                       modal_integral=:fixed, nr, kronrod)
        @test transform isa NonlinearRHS.TransModalFixed
        @test NonlinearRHS.has_error_estimate(transform) == kronrod
        sf = Stats.collect_stats(grid, Eω, Stats.transverse_integral_error(transform))
        d = sf(Eω, 0.0, 1e-4)
        @test d["transverse_points"] == nr

        # the same quantity, computed here from the transform
        nl = similar(Eω)
        transform(nl, Eω, 0.0)
        err = NonlinearRHS.integral_error!(transform)
        if kronrod
            rms = sqrt(sum(abs2, err)/length(err))
            @test d["transverse_integral_error_abs"] ≈ rms
            @test d["transverse_integral_error_rel"] ≈
                  rms/sqrt(sum(abs2, nl)/length(nl))
            @test 0 < d["transverse_integral_error_rel"] < 1e-3
        else
            @test all(isnan, err)
            @test isnan(d["transverse_integral_error_abs"])
            @test isnan(d["transverse_integral_error_rel"])
        end
        #= The right-hand side is evaluated only when there is an estimate to record:
           without one the statistic must not cost an evaluation per step. =#
        stat = Stats.transverse_integral_error(transform)
        fill!(stat.nl, NaN)
        dstat = Dict{String, Any}()
        stat(dstat, Eω, nothing, 0.0, 1e-4)
        @test all(isnan, stat.nl) == !kronrod
        @test dstat["transverse_points"] == nr
        # no mode reconstruction error: that is the adaptive rule's diagnostic
        @test !haskey(d, "mode_reconstruction_error")
        # ... and it is what `Stats.default` puts in the set for this transform
        sfd = Stats.default(grid, Eω, modes, linop, transform)
        dd = sfd(Eω, 0.0, 1e-4)
        @test haskey(dd, "transverse_integral_error_rel")
        @test !haskey(dd, "mode_reconstruction_error")
        @test Stats.host_statistics(sfd) ⊇ ["TransverseIntegralError"]
    end

    #= The adaptive transform keeps the statistic it always had, which records the
       reconstruction error as well as the cubature's own estimate. =#
    Eω, atr, FT = Luna.setup(grid, dens, resp, inputs, modes, :y;
                             modal_integral=:adaptive)
    @test atr isa NonlinearRHS.TransModal
    da = Stats.default(grid, Eω, modes, linop, atr)(Eω, 0.0, 1e-4)
    @test haskey(da, "mode_reconstruction_error")
    @test haskey(da, "transverse_integral_error_rel")
end

#= `Stats.default` for a radial or free-space state, against the same quantities computed
   here by hand from the state. The state is in transverse reciprocal space, so "by hand"
   means the inverse transverse transform (`Grid.to_rspace`, `FFTW.ifft`) and then the
   definition of each statistic. =#

"The pieces of a radial or free-space propagation, without running it."
function freesetup(; geom=:radial, GT=Grid.RealGrid, N=16, R=400e-6, gas=:Ar, pres=1.0,
                   λ0=800e-9, w0=100e-6, energy=1e-9, plasma=false)
    grid = GT === Grid.RealGrid ?
        Grid.RealGrid(λ0, (200e-9, 3000e-9), 100e-15) :
        Grid.EnvGrid(λ0, (200e-9, 3000e-9), 100e-15)
    sg = geom === :radial ? Grid.RadialGrid(R, N) :
         geom === :free2d ? Grid.Free2DGrid(R, N) : Grid.FreeGrid(R, N)
    ρ = PhysData.density(gas, pres)
    dens = z -> ρ
    resp = GT === Grid.RealGrid ?
        Any[Nonlinear.Kerr_field(PhysData.γ3_gas(gas))] :
        Any[Nonlinear.Kerr_env(PhysData.γ3_gas(gas))]
    if plasma
        push!(resp, Nonlinear.PlasmaCumtrapz(grid.to, grid.to,
                                             Ionisation.IonRateADK(gas),
                                             PhysData.ionisation_potential(gas)))
    end
    nfun = PhysData.ref_index_fun(gas, pres)
    linop = LinearOps.make_const_linop(grid, sg, nfun)
    normfun = geom === :radial ? NonlinearRHS.const_norm_radial(grid, sg, nfun) :
              geom === :free2d ? NonlinearRHS.const_norm_free2D(grid, sg, nfun) :
                                 NonlinearRHS.const_norm_free(grid, sg, nfun)
    inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy, w0)
    Eω, transform, FT = Luna.setup(grid, sg, dens, normfun, Tuple(resp), inputs)
    (; grid, sg, Eω, transform, FT, linop)
end

"The state in transverse real space, which every hand computation below starts from."
rspace(sg::Grid.RadialGrid, Eω) = Grid.to_rspace(sg, Eω; dim=3)
rspace(::Grid.Free2DGrid, Eω) = FFTW.ifft(Eω, 3)
rspace(::Grid.FreeGrid, Eω) = FFTW.ifft(Eω, (3, 4))

"The on-axis spectral field, from the state in real space."
function onaxisfield(sg::Grid.RadialGrid, Eω, Er)
    dropdims(Grid.onaxis(sg, Eω; dim=3); dims=2)
end
onaxisfield(sg::Grid.Free2DGrid, Eω, Er) = Er[:, 1, argmin(abs.(sg.x))]
onaxisfield(sg::Grid.FreeGrid, Eω, Er) =
    Er[:, 1, argmin(abs.(sg.x)), argmin(abs.(sg.y))]

"The analytic time-domain field of an on-axis spectrum, as `Stats.plan_analytic` makes it."
function analytic(grid::Grid.RealGrid, Eω0)
    FFTW.ifft(vcat(2 .* Eω0[1:end-1], zeros(ComplexF64, length(grid.ω) - 1)))
end
analytic(grid::Grid.EnvGrid, Eω0) = FFTW.ifft(Eω0)

@testset "$geom statistics, $GT" for geom in (:radial, :free2d, :free3d),
                                    GT in (Grid.RealGrid, Grid.EnvGrid)
    s = freesetup(; geom, GT)
    sf = Stats.default(s.grid, s.Eω, s.linop, s.transform;
                       gas=:Ar, windows=((600e-9, 1000e-9),))
    d = sf(s.Eω, 0.0, 1e-4)

    _, energyω = Fields.energyfuncs(s.grid, s.sg)
    @test d["energy"] ≈ [energyω(dropdims(selectdim(s.Eω, 2, 1:1); dims=2))]
    @test d["density"] ≈ PhysData.density(:Ar, 1.0)
    @test d["pressure"] ≈ 1.0
    @test d["z"] == 0.0 && d["dz"] == 1e-4

    Er = rspace(s.sg, s.Eω)
    Eω0 = onaxisfield(s.sg, s.Eω, Er)
    Et0 = analytic(s.grid, Eω0)
    @test d["ω0"] ≈ Maths.moment(s.grid.ω, abs2.(Eω0))
    @test d["peakintensity"] ≈ PhysData.c*PhysData.ε_0/2 * maximum(abs2, Et0)
    @test d["fwhm_t_max"] ≈ Maths.fwhm(s.grid.t, abs2.(Et0); method=:linear, minmax=:max)
    @test d["fwhm_t_min"] ≈ Maths.fwhm(s.grid.t, abs2.(Et0); method=:linear, minmax=:min)

    # the windowed energy is a part of the total
    @test 0 < d["energy_600.00nm_1000.00nm"][1] <= d["energy"][1]

    # the beam size and the collar fraction, from the transverse fluence profile
    P = dropdims(sum(abs2, Er; dims=(1, 2)); dims=(1, 2))
    mask = Boundaries.rprofile(s.sg, Boundaries.DEFAULT_RCOLLAR) .< 1
    if geom === :radial
        Psym = vcat(reverse(P), sum(abs2, Eω0), P)
        @test d["fwhm_r"] ≈ Maths.fwhm(Grid.rsymmetric(s.sg), Psym;
                                       method=:linear, minmax=:max)
        @test d["collar_energy_fraction"] ≈
              sum((s.sg.wr .* P)[mask])/sum(s.sg.wr .* P)
    elseif geom === :free2d
        @test d["fwhm_x"] ≈ Maths.fwhm(s.sg.x, P; method=:linear, minmax=:max)
        @test d["collar_energy_fraction"] ≈ sum(P[mask])/sum(P)
    else
        ix, iy = argmin(abs.(s.sg.x)), argmin(abs.(s.sg.y))
        @test d["fwhm_x"] ≈ Maths.fwhm(s.sg.x, P[:, iy]; method=:linear, minmax=:max)
        @test d["fwhm_y"] ≈ Maths.fwhm(s.sg.y, P[ix, :]; method=:linear, minmax=:max)
        @test d["collar_energy_fraction"] ≈ sum(P[mask])/sum(P)
    end
    # the beam is well inside the grid, so almost nothing is in the collar
    @test 0 <= d["collar_energy_fraction"] < 1e-6

    #= The beam profile is the only member with no device form, so it is what makes the
       set host-only; without it every statistic in the set has one. =#
    @test Stats.host_statistics(sf) == ["BeamProfile"]
    sfd = Stats.default(s.grid, s.Eω, s.linop, s.transform; beam_profile=false)
    @test isempty(Stats.host_statistics(sfd))
    # ... and turning it off changes nothing about the rest
    dd = sfd(s.Eω, 0.0, 1e-4)
    for k in keys(dd)
        @test (k, all(isapprox.(dd[k], d[k]))) == (k, true)
    end
end

#= The on-axis electron density, and the two-component field which drives the rate
   through its quadrature sum. =#
@testset "on-axis electron density" begin
    s = freesetup(; geom=:radial, N=32, R=200e-6, w0=40e-6, energy=20e-6, plasma=true)
    sf = Stats.default(s.grid, s.Eω, s.linop, s.transform)
    d = sf(s.Eω, 0.0, 1e-4)
    Eω0 = dropdims(Grid.onaxis(s.sg, s.Eω; dim=3); dims=2)
    Et0 = analytic(s.grid, Eω0)
    rate = Ionisation.IonRateADK(:Ar)
    frac = similar(s.grid.t)
    rate(frac, real(Et0))
    @test d["peak_ionisation_rate"] ≈ maximum(frac)
    Maths.cumtrapz!(frac, s.grid.t[2] - s.grid.t[1])
    @. frac = 1 - exp(-frac)
    @test d["electrondensity"] ≈ frac[end]*PhysData.density(:Ar, 1.0)
    @test d["electrondensity"]/PhysData.density(:Ar, 1.0) > 1e-8
end

#= End to end: a radial propagation with the default statistics recorded at every
   accepted step. The energy is conserved by a Kerr-only run, the beam diffracts, and
   nothing is NaN. =#
@testset "a radial propagation with the default statistics" begin
    s = freesetup(; geom=:radial, N=32, R=1e-3, w0=200e-6, energy=1e-9)
    sf = Stats.default(s.grid, s.Eω, s.linop, s.transform; gas=:Ar)
    out = Output.MemoryOutput(0, 0.02, 5, sf)
    Luna.run(s.Eω, s.grid, s.linop, s.transform, s.FT, out; zmax=0.02,
             status_period=100)
    st = out["stats"]
    @test length(st["z"]) > 3
    @test all(isfinite, st["energy"])
    @test all(isfinite, st["fwhm_r"])
    @test all(isfinite, st["ω0"])
    e = st["energy"][1, :]
    @test maximum(abs, e .- e[1])/e[1] < 1e-3
    @test st["fwhm_r"][end] > st["fwhm_r"][1] # the beam diffracts
end
