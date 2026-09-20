using Luna
import Luna.Capillary: besselj, get_unm
import Luna.Modes: hquadrature
import Test: @test, @testset, @test_throws
import Logging

logger = Logging.SimpleLogger(stdout, Logging.Warn)
old_logger = Logging.global_logger(logger)

@testset "Polarisation" begin
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, τfwhm=10e-15, energy=1e-12, trange=400e-15,
              λlims=(200e-9, 4e-6), shotnoise=false)
    lin = prop_capillary(args...; polarisation=:linear, kwargs...)
    @testset "linear, single mode" begin
        o1 = prop_capillary(args...; polarisation=0.0, kwargs...)
        o2 = prop_capillary(args...; polarisation=:linear, modes=1, kwargs...)
        @test o1["Eω"][:, 1, :] ≈ o2["Eω"][:, 1, :]
        @test size(o1["Eω"], 2) == 2*size(o2["Eω"], 2)
    end
    @testset "linear, multimode" begin
        o1 = prop_capillary(args...; polarisation=0.0, modes=4, kwargs...)
        o2 = prop_capillary(args...; polarisation=:linear, modes=4, kwargs...)
        @test o1["Eω"][:, 1:2:end, :] ≈ o2["Eω"][:, :, :]
    end
    @testset "x/y" begin
        o1 = prop_capillary(args...; polarisation=:x, modes=4, kwargs...)
        o2 = prop_capillary(args...; polarisation=:y, modes=4, kwargs...)
        @test o1["Eω"][:, 1:2:end, :] ≈ o2["Eω"][:, 2:2:end, :]
        @test all(iszero, o1["Eω"][:, 2:2:end, 1])
        @test isapprox(o1["stats"]["energy"][1, 1], kwargs.energy; rtol=1e-4)
        @test isapprox(o2["stats"]["energy"][2, 1], kwargs.energy; rtol=1e-4)
        @test all(iszero, o2["Eω"][:, 1:2:end, 1])
    end
    @testset "circular, $modes" for modes in (:HE11, :HE12, 1, 2)
        o1 = prop_capillary(args...; modes, polarisation=:circular, kwargs...)
        o2 = prop_capillary(args...; modes, polarisation=1.0, kwargs...)
        @test o1["Eω"] == o2["Eω"]
        o3 = prop_capillary(args...; modes, polarisation=-1.0, kwargs...)
        @test o2["Eω"][:, 1:2:end, :] ≈ o3["Eω"][:, 1:2:end, :]
        @test o2["Eω"][:, 2:2:end, :] ≈ -o3["Eω"][:, 2:2:end, :]
    end
    @testset "elliptical, $modes" for modes in (:HE11, :HE12, 1, 2)
        ε = 0.5
        o1 = prop_capillary(args...; polarisation=ε, kwargs...)
        @test ε^2*sum(abs2.(o1["Eω"][:, 1:2:end, :])) ≈ sum(abs2.(o1["Eω"][:, 2:2:end, :]))
    end
end

##
@testset "Peak power vs energy" begin
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, τfwhm=10e-15, shotnoise=false,
              trange=500e-15, λlims=(200e-9, 4e-6),
              saveN=51, plasma=false)
    shape_fac = ((:gauss, kwargs.τfwhm*sqrt(pi/log(16))),
                 (:sech, 2*kwargs.τfwhm/(2*log(1 + sqrt(2)))))
    pp = 1e8
    @testset "$pol polarisation" for pol in (:linear, :circular)
        @testset "modes: $modes" for modes in (:HE11, :HE12, 1, 2)
            @testset "$(sf[1])" for sf in shape_fac
                s, f = sf
                op = prop_capillary(args...; pulseshape=s, power=pp, modes, polarisation=pol, kwargs...)
                oe = prop_capillary(args...; pulseshape=s, energy=f*pp, modes, polarisation=pol, kwargs...)
                @test Processing.energy(op) ≈ Processing.energy(oe)
                if s == :gauss
                    # stretch by factor of sqrt(2) with GDD
                    # peak power drops by 1/sqrt(2) but energy is the same
                    φ2 = Tools.τfw_to_τ0(kwargs.τfwhm, :gauss)^2
                    op = prop_capillary(args...; pulseshape=s, power=pp/sqrt(2),
                                                modes, polarisation=pol, ϕ=[0, 0, φ2],
                                                kwargs...)
                    oe = prop_capillary(args...; pulseshape=s, energy=f*pp,
                                        modes, polarisation=pol, ϕ=[0, 0, φ2],
                                        kwargs...)
                    @test Processing.energy(op) ≈ Processing.energy(oe)
                end
            end
        end
    end
end

##
prop!(Eω, grid) = nothing # do-nothing propagator
@testset "Input into higher-order modes" begin
    @testset "propagator $prop" for prop in (nothing, prop!)
        args = (100e-6, 0.1, :He, 1)
        kwargs = (λ0=800e-9, shotnoise=false, trange=250e-15, λlims=(200e-9, 4e-6),
                saveN=51, plasma=false, propagator=prop)
        pkwargs = (τfwhm=10e-15, energy=1e-12, λ0=800e-9)
        @testset "input into $m, mode average" for m in (:HE11, :HE12, :TE01, :TE02, :TM01)
            ip = Pulses.GaussPulse(;mode=m, pkwargs...)
            o = prop_capillary(args...; pulses=ip, modes=m, kwargs...)
            @test Processing.energy(o)[1] ≈ pkwargs.energy
        end
        @testset "input into $m, modal" for (midx, m) in enumerate((:HE11, :HE12, :HE13, :HE14))
            ip = Pulses.GaussPulse(;mode=m, pkwargs...)
            o = prop_capillary(args...; pulses=ip, modes=4, kwargs...)
            @test Processing.energy(o)[midx, 1] ≈ pkwargs.energy
        end
        @testset "input into $m, modal circular" for (midx, m) in enumerate((:HE11, :HE12, :HE13, :HE14))
            ip = Pulses.GaussPulse(;mode=m, polarisation=:circular, pkwargs...)
            o = prop_capillary(args...; pulses=ip, modes=4, kwargs...)
            @test Processing.energy(o)[2midx-1, 1] ≈ pkwargs.energy/2
            @test Processing.energy(o)[2midx, 1] ≈ pkwargs.energy/2
        end
        modes = (:HE21, :HE22, :HE23, :HE24)
        @testset "input into $m, modal" for (midx, m) in enumerate(modes)
            ip = Pulses.GaussPulse(;mode=m, pkwargs...)
            o = prop_capillary(args...; pulses=ip, modes, kwargs...)
            @test Processing.energy(o)[midx, 1] ≈ pkwargs.energy
        end
        @testset "input into $m, modal circular" for (midx, m) in enumerate(modes)
            ip = Pulses.GaussPulse(;mode=m, polarisation=:circular, pkwargs...)
            o = prop_capillary(args...; pulses=ip, modes, kwargs...)
            @test Processing.energy(o)[2midx-1, 1] ≈ pkwargs.energy/2
            @test Processing.energy(o)[2midx, 1] ≈ pkwargs.energy/2
        end
        modes = (:TE01, :TE02, :TE03, :TE04, :TM01, :TM02, :TM03, :TM04)
        @testset "input into $m, modal" for (midx, m) in enumerate((modes))
            ip = Pulses.GaussPulse(;mode=m, pkwargs...)
            o = prop_capillary(args...; pulses=ip, modes=modes, kwargs...)
            @test Processing.energy(o)[midx, 1] ≈ pkwargs.energy
        end
    end
end

##
@testset "multiple inputs" begin
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, shotnoise=false, trange=250e-15, λlims=(200e-9, 4e-6),
              saveN=51, plasma=false)
    p1 = (λ0=800e-9, energy=1e-12, τfwhm=10e-15)
    p2 = (λ0=400e-9, energy=2e-12, τfwhm=30e-15)
    modes = (:HE11, :HE12, :HE13, :HE14)
    @testset "first pulse into $m1" for (idx1, m1) in enumerate(modes)
        ip1 = Pulses.GaussPulse(;mode=m1, p1...)
        @testset "second pulse into $m2" for (idx2, m2) in enumerate(modes)
            ip2 = Pulses.GaussPulse(;mode=m2, p2...)
            o = prop_capillary(args...; pulses=[ip1, ip2], modes=4, kwargs...)
            if idx1 == idx2
                @test Processing.energy(o)[idx1, 1] ≈ p1.energy + p2.energy
            else
                @test Processing.energy(o)[idx1, 1] ≈ p1.energy
                @test Processing.energy(o)[idx2, 1] ≈ p2.energy
                @test isapprox(Processing.fwhm_t(o)[idx1, 1], p1.τfwhm, rtol=1e-3)
                @test isapprox(Processing.fwhm_t(o)[idx2, 2], p2.τfwhm, rtol=1e-3)
            end
        end
    end
end

##
@testset "propagators" begin
    # passing ϕ keyword argument and an equivalent propagator function should yield
    # the same result.
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, shotnoise=false, trange=250e-15, λlims=(200e-9, 4e-6),
              saveN=51, plasma=false)
    p = (λ0=800e-9, energy=1e-12, τfwhm=10e-15)
    ϕ = [0, 0, 10e-30, 100e-45]
    function prop!(Eω, grid)
        Fields.prop_taylor!(Eω, grid, ϕ, 800e-9)
    end
    pp = Pulses.GaussPulse(;p..., propagator=prop!)
    pt = Pulses.GaussPulse(;p..., ϕ)
    op = prop_capillary(args...; pulses=pp, kwargs...)
    ot = prop_capillary(args...; pulses=pt, kwargs...)
    @test isapprox(Processing.fwhm_t(op)[1], Processing.fwhm_t(ot)[1], rtol=1e-3)
    @test Processing.energy(op)[1] ≈ p.energy
    @test Processing.energy(ot)[1] ≈ p.energy

    op = prop_capillary(args...; p..., propagator=prop!, kwargs...)
    ot = prop_capillary(args...; p..., ϕ, kwargs...)
    @test isapprox(Processing.fwhm_t(op)[1], Processing.fwhm_t(ot)[1], rtol=1e-3)
    @test Processing.energy(op)[1] ≈ p.energy
    @test Processing.energy(ot)[1] ≈ p.energy
end

##
@testset "GaussBeamPulse" begin
    function overlap(n, m, kind, w0)
        unm = get_unm(n, m, kind)

        mode(ρ) = ρ >= 1 ? 0.0 : besselj(0, unm*ρ)
        beam(ρ) = exp(-ρ^2/w0^2)

        sqrt_numerator, _ = hquadrature(0, 1) do ρ
            mode(ρ)*beam(ρ)*ρ
        end

        den1, _ = hquadrature(0, 1) do ρ
            abs2(mode(ρ))*ρ
        end

        den2, _ = hquadrature(0, 10) do ρ
            abs2(beam(ρ))*ρ
        end

        sqrt_numerator/sqrt((den1*den2))
    end
    Nmodes = 16
    ovlp = overlap.(1, 1:Nmodes, :HE, 0.64)
    gauss_overlaps = abs2.(ovlp)
    phases = angle.(ovlp)

    a = 100e-6
    args = (a, 0.1, :He, 1)
    kwargs = (λ0=800e-9, shotnoise=false, trange=250e-15,
              λlims=(200e-9, 4e-6), saveN=51, plasma=false, loss=false)
    p = (λ0=800e-9, energy=1e-12, τfwhm=10e-15)
    gpl = Pulses.GaussPulse(;p...)
    gpc = Pulses.GaussPulse(;polarisation=:circular, p...)
    pulse = Pulses.GaussBeamPulse(0.64*a, gpl)
    Eω, grid, linop, transform, FT, o = Interface.prop_capillary_args(args...; pulses=pulse, modes=Nmodes, kwargs...)
    Luna.run(Eω, grid, linop, transform, FT, o; zmax=args[2])
    @testset for m in 1:Nmodes
        @test Processing.energy(o)[m, 1] ≈ p.energy * gauss_overlaps[m]
    end

    # testing internals of GaussBeamPulse separately
    # do we get the same overlap integrals?
    modes = transform.ts.ms
    k = 2π/kwargs[:λ0]
    gauss = Fields.normalised_gauss_beam(k, pulse.waist)
    ovlps = [Modes.overlap(mi, gauss) for mi in modes]
    @test all(ovlps .≈ ovlp)

    # circular polarisation
    pulse = Pulses.GaussBeamPulse(0.64*a, gpc)
    o = prop_capillary(args...; pulses=pulse, modes=Nmodes, kwargs...)
    @testset for m in 1:Nmodes
        @test Processing.energy(o)[2m, 1] ≈ p.energy * gauss_overlaps[m]/2
        @test Processing.energy(o)[2m-1, 1] ≈ p.energy * gauss_overlaps[m]/2
    end

    # two GaussBeamPulses
    pulse1 = Pulses.GaussBeamPulse(0.64*a, gpl)
    gpl2 = Pulses.GaussPulse(;ϕ=[0, 100e-15], p...)
    pulse2 = Pulses.GaussBeamPulse(0.64*a, gpl2)
    o = prop_capillary(args...; pulses=[pulse1, pulse2], modes=Nmodes, kwargs...)
    @testset for m in 1:Nmodes
        @test Processing.energy(o)[m, 1] ≈ 2p.energy * gauss_overlaps[m]
    end
end

##
@testset "Defaults" begin
    @testset "Envelope propagation: $env" for env in [false, true]
        @testset "Defaults for $gas" for gas in PhysData.gas
            # Check that prop_capillary has working defaults for all gases
            gas == :Air && continue
            a = 100e-6
            flength = 0.1
            pressure = 1
            kwargs = (λ0=800e-9, energy=1e-12, τfwhm=10e-15, shotnoise=false, trange=250e-15,
                      λlims=(200e-9, 4e-6), saveN=51, envelope=env)
            prop_capillary(a, flength, gas, pressure; kwargs...)
            @test true
        end
    end
end

##
@testset "LunaPulse" begin
    # single-mode
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(200e-9, 4e-6), shotnoise=false, modes=:HE11)
    e1 = 10e-9
    o1 = prop_capillary(args...; energy=e1, kwargs...)
    eo1 = Processing.energy(o1)[end]

    # change wavelength limits to force re-gridding
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=:HE11)

    # Defining nothing
    p = Pulses.LunaPulse(o1)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # energy should be the same as out of stage 1
    ei2 = Processing.energy(o2)[1]
    @test ei2 ≈ eo1

    # Defining the overall energy
    e2 = 5e-9
    p = Pulses.LunaPulse(o1; energy=e2)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # energy should be e2 as defined
    ei2 = Processing.energy(o2)[1]
    @test ei2 ≈ e2

    # Defining the energy
    es = 0.5
    p = Pulses.LunaPulse(o1; scale_energy=es)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # total energy should be es times output energy
    ei2 = Processing.energy(o2)[1]
    @test ei2 ≈ eo1*es


    # multi-mode
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(200e-9, 4e-6), shotnoise=false, modes=4)
    e1 = 10e-9
    o1 = prop_capillary(args...; energy=e1, kwargs...)
    eo1 = Processing.energy(o1)[:, end]

    # change wavelength limits to force re-gridding
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=4)

    # Defining nothing
    p = Pulses.LunaPulse(o1;)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # total energy should be the same as out of stage 1
    ei2 = Processing.energy(o2)[:, 1]
    @test sum(ei2) ≈ sum(eo1)
    # relative energy in each mode should be the same
    @test ei2 ./ sum(ei2) ≈ eo1 ./ sum(eo1)

    # Defining the overall energy
    e2 = 5e-9
    p = Pulses.LunaPulse(o1; energy=e2)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # total energy should be e2 as defined
    ei2 = Processing.energy(o2)[:, 1]
    @test sum(ei2) ≈ e2
    # relative energy in each mode should be the same
    @test ei2 ./ sum(ei2) ≈ eo1 ./ sum(eo1)

    # Defining the overall energy scale
    es = 0.5
    p = Pulses.LunaPulse(o1; scale_energy=es)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # total energy should be es times output energy
    ei2 = Processing.energy(o2)[:, 1]
    @test ei2 ≈ eo1*es
    # relative energy in each mode should be the same
    @test ei2 ./ sum(ei2) ≈ eo1 ./ sum(eo1)

    # Defining mode dependent energy scale
    es = [1, 0.75, 0.5, 0.25]
    p = Pulses.LunaPulse(o1; scale_energy=es)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    # total energy should be weighted sum from before
    ei2 = Processing.energy(o2)[:, 1]
    @test ei2 ≈ eo1 .* es

    # make sure coupling to fewer modes throws an error
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=2)
    p = Pulses.LunaPulse(o1;)
    @test_throws ErrorException o2 = prop_capillary(args...; pulses=p, kwargs...)

    # make sure coupling to *different* modes throws an error
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=(:HE21, :HE22, :HE23, :HE24))
    p = Pulses.LunaPulse(o1;)
    @test_throws ErrorException o2 = prop_capillary(args...; pulses=p, kwargs...)

    # coupling to more modes should work fine
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=6)
    p = Pulses.LunaPulse(o1;)
    o2 = prop_capillary(args...; pulses=p, kwargs...)
    ei2 = Processing.energy(o2)[:, 1]
    @test sum(ei2) ≈ sum(eo1)
    @test all(ei2[5:end] .== 0)

    # check that adding a LunaPulse with additional inputs works
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=4)
    e2 = 1e-9
    p = Pulses.LunaPulse(o1;)
    p2 = Pulses.GaussPulse(;λ0=400e-9, τfwhm=30e-15, energy=e2)
    o2 = prop_capillary(args...; pulses=[p, p2], kwargs...)
    ei2 = Processing.energy(o2)[:, 1]
    @test sum(ei2) ≈ sum(eo1) + e2

    # check that adding a LunaPulse with additional multi-mode input works
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=8)
    e2 = 1e-9
    p = Pulses.LunaPulse(o1;)
    p2 = Pulses.GaussPulse(;λ0=400e-9, τfwhm=30e-15, energy=e2)
    gp2 = Pulses.GaussBeamPulse(0.64*args[1], p2)
    o2 = prop_capillary(args...; pulses=[p, gp2], kwargs...)
    ei2 = Processing.energy(o2)[:, 1]
    # need higher tolerance here since the gaussian beam overlap introduces a bit of error
    @test isapprox(sum(ei2), sum(eo1) + e2, rtol=1e-3)

    # two LunaPulses
    kwargs2 = (λ0=400e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=4)
    o12 = prop_capillary(args...; energy=e1, kwargs2...)
    eo12 = Processing.energy(o12)[:, end]
    kwargs = (λ0=800e-9, τfwhm=10e-15, trange=400e-15, λlims=(150e-9, 4e-6), shotnoise=false, modes=8)
    e2 = 1e-9
    p = Pulses.LunaPulse(o1)
    p2 = Pulses.LunaPulse(o12)
    o2 = prop_capillary(args...; pulses=[p, p2], kwargs...)
    ei2 = Processing.energy(o2)[:, 1]
    # need higher tolerance here since the gaussian beam overlap introduces a bit of error
    @test sum(ei2) ≈ sum(eo1) + sum(eo12)

end

##
@testset "Temperature" begin
    # test that changing temperature changes results
    args = (100e-6, 0.1, :He, 1)
    kwargs = (λ0=800e-9, τfwhm=10e-15, energy=1e-12, trange=400e-15,
              λlims=(200e-9, 4e-6), shotnoise=false)
    o = prop_capillary(args...; temperature=300, kwargs...)
    o2 = prop_capillary(args...; temperature=300, kwargs...)
    o3 = prop_capillary(args...; temperature=400, kwargs...)
    @test o["Eω"] == o2["Eω"]
    @test o["Eω"] ≠ o3["Eω"]

    # test with kerr/plasma off (only Raman depends on temperature here)
    args = (100e-6, 0.1, :H2, 1)
    kwargs = (λ0=800e-9, τfwhm=10e-15, energy=1e-12, trange=400e-15,
              λlims=(200e-9, 4e-6), shotnoise=false, kerr=false, plasma=false)
    o = prop_capillary(args...; temperature=300, kwargs...)
    o2 = prop_capillary(args...; temperature=300, kwargs...)
    o3 = prop_capillary(args...; temperature=400, kwargs...)
    @test o["Eω"] == o2["Eω"]
    @test o["Eω"] ≠ o3["Eω"]

    Eω, _, _, t300, _, _ = Interface.prop_capillary_args(args...; temperature=300, kwargs...)
    _, _, _, t400, _, _ = Interface.prop_capillary_args(args...; temperature=400, kwargs...)
    Raman300 = t300.resp[1]
    Raman400 = t400.resp[1]
    NonlinearRHS.to_time!(t300.Eto, Eω, t300.Eωo, t300.IFT)
    ρ = t300.densityfun(0) # note: using same density to compare only NL response
    Pto300 = zero(t300.Eto)
    Pto300_2 = zero(t300.Eto)
    Pto400 = zero(t300.Eto)
    Raman300(Pto300, t300.Eto, ρ)
    Raman300(Pto300_2, t300.Eto, ρ)
    Raman400(Pto400, t300.Eto, ρ)
    @test Pto300 ≠ Pto400
    @test Pto300 == Pto300_2
end

##
#=
THG in a gas-filled capillary (mode-averaged): real and envelope fields must agree.
This tests the full envelope+thg chain through prop_capillary: the thg=true grid, the
Kerr_env_thg response, and the carrier-transparent (group-delay-only) linear operator
reference frame, which prop_capillary passes through to LinearOps. With a frame that
subtracts a constant carrier phase instead, the THG term acquires a spurious phase
mismatch 2β1ω0 and the third harmonic (almost) vanishes.
=#
@testset "THG: real vs envelope" begin
    args = (125e-6, 1e-3, :Ar, 2)
    kwargs = (λ0=800e-9, energy=1e-6, τfwhm=20e-15, trange=100e-15,
              λlims=(220e-9, 2000e-9), shotnoise=false, plasma=false, raman=false,
              thg=true, saveN=3)
    or = prop_capillary(args...; envelope=false, kwargs...)
    oe = prop_capillary(args...; envelope=true, kwargs...)

    etot_r = Processing.energy(or)[end]
    etot_e = Processing.energy(oe)[end]
    @test etot_r ≈ etot_e rtol=1e-6

    thgband = (240e-9, 300e-9)
    ethg_r = Processing.energy(or; bandpass=thgband)[end]
    ethg_e = Processing.energy(oe; bandpass=thgband)[end]
    @test ethg_r > 1e5*eps(etot_r) # THG was actually generated
    @test ethg_r ≈ ethg_e rtol=0.05
end

##
@testset "device, precision and stats_period keywords" begin
    import Luna: DeviceSpec, Output
    args = (125e-6, 1e-2, :He, 1.0)
    kwargs = (λ0=800e-9, energy=100e-9, τfwhm=10e-15, trange=400e-15,
              λlims=(300e-9, 2000e-9), shotnoise=false, plasma=false, raman=false,
              saveN=5)

    # An untouched call is unaffected: `device`'s default (Luna.device_request()) is
    # :cpu unless a GPU package is loaded, exactly what prop_capillary always gave.
    oref = prop_capillary(args...; kwargs...)
    odefault = prop_capillary(args...; kwargs..., device=Luna.device_request())
    @test odefault["Eω"] == oref["Eω"]
    @test eltype(oref["Eω"]) === ComplexF64

    # precision=Float32 on the CPU: unscaled automatically, saved as ComplexF32
    o32 = prop_capillary(args...; kwargs..., precision=Float32)
    @test eltype(o32["Eω"]) === ComplexF32
    for idx in axes(oref["Eω"], 2)
        h = oref["Eω"][:, idx]
        d = ComplexF64.(o32["Eω"][:, idx])
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-5
    end

    # device=<a DeviceSpec>: the low-level array-type test (JLArray, Metal) is
    # test_device.jl's/test_metal.jl's; here only the plumbing, at Float64 so the
    # comparison is exact bar the FFT/reduction order the array type changes.
    ohost = prop_capillary(args...; kwargs..., device=DeviceSpec(Array, Float64))
    @test ohost["Eω"] == oref["Eω"]

    # stats_period shortens the recorded statistics but not the propagation itself
    o3 = prop_capillary(args...; kwargs..., stats_period=3)
    @test length(o3["stats"]["z"]) < length(oref["stats"]["z"])
    @test o3["z"] == oref["z"]

    # a pressure gradient still runs (the z-dependent operator branch)
    ograd = prop_capillary(125e-6, 1e-2, :He, (1.0, 0.0); kwargs...)
    @test size(ograd["Eω"]) == size(oref["Eω"])

    #= A response with no device kernel is refused for an explicit `device` *or*
       `precision` request, naming the fix. Review round 1 of gpu/12-response-traits,
       finding 9: the `precision`-only path was untested, and without it
       `precision=Float32` would run the response through `Nonlinear.HostResponse` at
       every step instead of erroring. Raman is the case since gpu/13-plasma gave the
       plasma response a kernel; it gets one of its own in gpu/14-raman. =#
    ramankw = (λ0=800e-9, energy=100e-9, τfwhm=10e-15, trange=400e-15,
               λlims=(300e-9, 2000e-9), shotnoise=false, plasma=false, raman=true,
               saveN=3)
    ramanargs = (args[1], args[2], :N2, 0.5)
    err = try
        prop_capillary(ramanargs...; ramankw..., precision=Float32)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("device=:cpu", err.msg)
    @test occursin("RamanPolarField", err.msg)
    # the message names the type, not its several hundred characters of parameters
    @test !occursin("FFTW", err.msg)
    @test_throws ErrorException prop_capillary(ramanargs...; ramankw...,
                                               device=DeviceSpec(Array, Float32))
    # ... and the same call with neither keyword is unaffected
    oraman = prop_capillary(ramanargs...; ramankw...)
    @test eltype(oraman["Eω"]) === ComplexF64

    #= Plasma is device-capable since gpu/13-plasma, so the default field-resolved
       response set of a non-Raman gas -- Kerr and plasma -- now follows an explicit
       `precision`/`device` request instead of being refused. =#
    plasmakw = (λ0=800e-9, energy=100e-9, τfwhm=10e-15, trange=400e-15,
                λlims=(300e-9, 2000e-9), shotnoise=false, plasma=true, raman=false,
                saveN=3)
    oplasma = prop_capillary(args...; plasmakw...)
    @test eltype(oplasma["Eω"]) === ComplexF64
    oplasma32 = prop_capillary(args...; plasmakw..., precision=Float32)
    @test eltype(oplasma32["Eω"]) === ComplexF32
    @test size(oplasma32["Eω"]) == size(oplasma["Eω"])
    @test maximum(abs, oplasma32["Eω"][:, end] .- oplasma["Eω"][:, end])/
          maximum(abs, oplasma["Eω"][:, end]) < 1e-4

    # multimode/radial propagation is not device-capable: refused, not silently ignored
    @test_throws ErrorException prop_capillary(args...; kwargs..., modes=4,
                                               device=DeviceSpec(Array, Float32))
    # ... but unaffected at the default device
    om = prop_capillary(args...; kwargs..., modes=4)
    @test size(om["Eω"], 2) == 4

    # prop_gnlse: same keywords, but only the CPU/Float64 default is honoured
    gargs = (0.1, 1e-3, [0.0, 0.0, -1e-26])
    gkwargs = (λ0=835e-9, τfwhm=100e-15, power=1e3, pulseshape=:sech,
               λlims=(450e-9, 2e-6), trange=1e-12, saveN=3, raman=false,
               shotnoise=false)
    ognlse = prop_gnlse(gargs...; gkwargs..., device=Luna.device_request())
    @test_throws ErrorException prop_gnlse(gargs...; gkwargs...,
                                           device=DeviceSpec(Array, Float32))
    ognlse3 = prop_gnlse(gargs...; gkwargs..., stats_period=2)
    @test length(ognlse3["stats"]["z"]) < length(ognlse["stats"]["z"])
end

##
Logging.global_logger(old_logger)
