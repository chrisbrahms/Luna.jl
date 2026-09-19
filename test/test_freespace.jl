#= here we test that everything runs without throwing errors for every combination of:
    1. Real and envelope grids
    2. Radial (QDHT), 2D cartesian and 3D cartesian spatial grids
    3. (For envelope) THG on/off
    4. Constant pressure and pressure gradient
We also check that the spatial linear operators work correctly by testing the focusing
of a Gaussian beam. (This cross-checks LinearOps vs Fields.prop!)
=#
using Luna
import Luna.PhysData: wlfreq
import Luna: Hankel, FFTW
Luna.set_fftw_mode(:estimate)
import LinearAlgebra: norm
import Test: @test, @testset

#= 1 mm rather than 0.6 mm: the transverse absorbing collars (r > 0.9R radially, |x| > R
   on the Cartesian grids) absorb the wings of the beam coherently, and a wing amplitude
   of 1e-2 at the collar shifts the focal profile by ~1e-3, which the Cartesian tests
   resolve. At 1 mm the input beam (w = 280 μm) has an amplitude below 1e-4 there. =#
R = 1.0e-3
Nr = 64
Nx = 64
Ny = 32
gas = :Ar
pressure = 1

λ0 = 800e-9
w0 = 200e-6
τfwhm = 20e-15
energy = 1e-12
# 0.15 m rather than 0.3 m so that the input beam stays clear of the collars, see R above
L = 0.15

rgrid = Grid.RealGrid(λ0, (400e-9, 2000e-9), 0.2e-12)
egrid = Grid.EnvGrid(λ0, (400e-9, 2000e-9), 0.2e-12)
q = Hankel.QDHT(R, Nr, dim=3)
xygrid = Grid.FreeGrid(R, Nx, R, Ny)
xgrid = Grid.Free2DGrid(R, Nx)

getshape(grid, q::Hankel.QDHT, pol) = (length(grid.ω), pol ? 2 : 1, q.N)
getshape(grid, sg::Grid.Free2DGrid, pol) = (length(grid.ω), pol ? 2 : 1, length(sg.x))
getshape(grid, sg::Grid.FreeGrid, pol) = (length(grid.ω), pol ? 2 : 1, length(sg.x), length(sg.y))

makekerr(grid::Grid.RealGrid, thg) = Nonlinear.Kerr_field(PhysData.γ3_gas(gas))
makekerr(grid::Grid.EnvGrid, thg) = thg ? Nonlinear.Kerr_env_thg(PhysData.γ3_gas(gas), grid.ω0, grid.to) : Nonlinear.Kerr_env(PhysData.γ3_gas(gas))

makeconstnorm(grid, q::Hankel.QDHT, nfunλ) = NonlinearRHS.const_norm_radial(grid, q, nfunλ)
makeconstnorm(grid, sg::Grid.Free2DGrid, nfunλ) = NonlinearRHS.const_norm_free2D(grid, sg, nfunλ)
makeconstnorm(grid, sg::Grid.FreeGrid, nfunλ) = NonlinearRHS.const_norm_free(grid, sg, nfunλ)

makenorm(grid, q::Hankel.QDHT, nfunω) = NonlinearRHS.norm_radial(grid, q, nfunω)
makenorm(grid, sg::Grid.Free2DGrid, nfunω) = NonlinearRHS.norm_free2D(grid, sg, nfunω)
makenorm(grid, sg::Grid.FreeGrid, nfunω) = NonlinearRHS.norm_free(grid, sg, nfunω)

function testfocus(q::Hankel.QDHT, Eω, w0)
    Eωfoc = Eω[:, :, :, end]
    Eωr = q \ Eωfoc
    Ir = dropdims(sum(abs2.(Eωr); dims=(1, 2)); dims=(1, 2))
    Ir_analytical = Maths.gauss.(q.r, w0/2)
    @test Ir/norm(Ir) ≈ Ir_analytical/norm(Ir_analytical) rtol=0.01
end

function testfocus(sg::Grid.Free2DGrid, Eω, w0)
    Eωfoc = Eω[:, :, :, end]
    Eωx = FFTW.ifft(Eωfoc, 3)
    Ix = dropdims(sum(abs2.(Eωx); dims=(1, 2)); dims=(1, 2))
    Ix_analytical = Maths.gauss.(sg.x, w0/2)
    @test Ix/norm(Ix) ≈ Ix_analytical/norm(Ix_analytical) rtol=0.001
end

function testfocus(sg::Grid.FreeGrid, Eω, w0)
    Eωfoc = Eω[:, :, :, :, end]
    Eωxy = FFTW.ifft(Eωfoc, (3, 4))
    Iy = dropdims(sum(abs2.(Eωxy); dims=(1, 2, 3)); dims=(1, 2, 3))
    Ix = dropdims(sum(abs2.(Eωxy); dims=(1, 2, 4)); dims=(1, 2, 4))
    Ix_analytical = Maths.gauss.(sg.x, w0/2)
    Iy_analytical = Maths.gauss.(sg.y, w0/2)
    @test Ix/norm(Ix) ≈ Ix_analytical/norm(Ix_analytical) rtol=0.001
    @test Iy/norm(Iy) ≈ Iy_analytical/norm(Iy_analytical) rtol=0.001
end

function runprop_const(grid, sg, thg, pol)
    nfunλ = PhysData.ref_index_fun(gas, pressure)
    if pol
        nfun = (λ; z=0.0) -> (nfunλ(λ), nfunλ(λ))
    else
        nfun = (λ; z=0.0) -> nfunλ(λ)
    end

    linop = LinearOps.make_const_linop(grid, sg, nfun, thg)
    dens0 = PhysData.density(gas, pressure)
    densityfun(z) = dens0

    responses = (makekerr(grid, thg),)
    normfun = makeconstnorm(grid, sg, nfun)

    inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0, propz=-L)

    Eω, transform, FT = Luna.setup(grid, sg, densityfun, normfun, responses, inputs)
    output = Output.MemoryOutput(0, L, 11)
    Luna.run(Eω, grid, linop, transform, FT, output; init_dz=0.1, zmax=L)
    output["Eω"]
end

function runprop_grad(grid, sg, thg, pol)
    nfunλ = PhysData.ref_index_fun(gas, pressure)
    if pol
        nfun = (λ; z=0.0) -> (nfunλ(λ), nfunλ(λ))
    else
        nfun = (λ; z=0.0) -> nfunλ(λ)
    end
    nfunω = (ω; z) -> nfun(wlfreq(ω); z)

    linop = LinearOps.make_linop(grid, sg, nfunω, thg)
    dens0 = PhysData.density(gas, pressure)
    densityfun(z) = dens0

    responses = (makekerr(grid, thg),)
    normfun = makenorm(grid, sg, nfunω)

    inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0, propz=-L)

    Eω, transform, FT = Luna.setup(grid, sg, densityfun, normfun, responses, inputs)
    output = Output.MemoryOutput(0, L, 11)
    Luna.run(Eω, grid, linop, transform, FT, output; init_dz=0.1, zmax=L)
    output["Eω"]
end

@testset "Constant pressure: $(typeof(grid)), $(typeof(sg)), pol = $pol, thg = $thg" for sg in (q, xgrid, xygrid),
                                                            pol in (false, true),
                                                            thg in (false, true),
                                                            grid in (rgrid, egrid)
    if grid isa Grid.RealGrid && ~thg
        continue
    end
    Eω = runprop_const(grid, sg, thg, pol)

    testfocus(sg, Eω, w0)
end
##
@testset "Gradient pressure: $(typeof(grid)), $(typeof(sg)), pol = $pol, thg = $thg" for sg in (q, xgrid, xygrid),
                                                            pol in (false, true),
                                                            thg in (false, true),
                                                            grid in (rgrid, egrid)
    if grid isa Grid.RealGrid && ~thg
        continue
    end
    Eω = runprop_grad(grid, sg, thg, pol)

    testfocus(sg, Eω, w0)
end
##
#=
Test the χ⁽²⁾ response in propagation: type I SHG in BBO (miniature version of
examples/low_level_interface/freespace/free2D_bbo.jl), comparing real and envelope fields.
The input (y polarisation, the ordinary axis) converts to second harmonic on the
x polarisation (extraordinary axis).
=#
@testset "BBO SHG: real vs envelope" begin
    λ0_bbo = 800e-9
    τfwhm_bbo = 30e-15
    w0_bbo = 20e-6
    energy_bbo = 10e-9
    thickness = 30e-6
    θ = deg2rad(29.2) # type I phase-matching angle
    ϕ = deg2rad(30) # ϕ for type I phase-matching

    bboxgrid = Grid.Free2DGrid(4w0_bbo, 2^5)
    nfuns = PhysData.ref_index_fun_xy(:BBO, θ)

    makechi2(grid::Grid.RealGrid) = Nonlinear.Chi2Field(θ, ϕ, PhysData.χ2(:BBO))
    makechi2(grid::Grid.EnvGrid) = Nonlinear.Chi2Env(θ, ϕ, PhysData.χ2(:BBO), grid.ω0, grid.to)

    function runshg(grid)
        linop = LinearOps.make_const_linop(grid, bboxgrid, nfuns)
        normfun = NonlinearRHS.const_norm_free2D(grid, bboxgrid, nfuns)
        densityfun = z -> 1 # unity density: we're considering a solid
        responses = (makechi2(grid),)
        inputs = Fields.GaussGaussField(;λ0=λ0_bbo, τfwhm=τfwhm_bbo,
                                         energy=energy_bbo/(sqrt(π/2)*w0_bbo), w0=w0_bbo)
        Eω, transform, FT = Luna.setup(grid, bboxgrid, densityfun, normfun, responses, inputs)
        output = Output.MemoryOutput(0, thickness, 5)
        Luna.run(Eω, grid, linop, transform, FT, output; init_dz=1e-6, zmax=thickness)
        eω = Fields.energyfuncs(grid, bboxgrid)[2]
        Eωk = output["Eω"] # (ω, pol, k, z)
        dropdims(mapslices(eω, Eωk; dims=(1, 3)); dims=(1, 3)) # energy: (pol, z)
    end

    rgrid_bbo = Grid.RealGrid(λ0_bbo, (250e-9, 2e-6), 120e-15)
    egrid_bbo = Grid.EnvGrid(λ0_bbo, (250e-9, 2e-6), 120e-15; thg=true)

    energy_r = runshg(rgrid_bbo)
    energy_e = runshg(egrid_bbo)

    for energy_out in (energy_r, energy_e)
        # input is on the y polarisation (index 2)
        @test energy_out[2, 1] ≈ sum(energy_out[:, 1]) rtol=1e-6
        # second harmonic grows on the x polarisation (index 1)
        @test energy_out[1, end] > 1e3*max(energy_out[1, 1], eps())
        @test energy_out[1, end] > 1e-4*energy_out[2, 1]
    end
    # real and envelope fields solve the same physics
    @test energy_r[1, end] ≈ energy_e[1, end] rtol=0.05
    @test energy_r[2, end] ≈ energy_e[2, end] rtol=0.05
end
##
#=
Test THG in propagation: real vs envelope fields in argon (radial symmetry). This checks
the carrier-transparent (`thg=true`) envelope reference frame in the generic linear
operators: a frame which subtracts a constant carrier phase adds a spurious phase mismatch
2β1ω0 to the THG term, suppressing the third harmonic and forcing tiny steps.
=#
@testset "THG: real vs envelope" begin
    λ0thg = 800e-9
    ω0 = wlfreq(λ0thg)
    pr_thg = 2
    Lthg = 1e-3

    rgridt = Grid.RealGrid(λ0thg, (220e-9, 2000e-9), 100e-15)
    egridt = Grid.EnvGrid(λ0thg, (220e-9, 2000e-9), 100e-15; thg=true)

    makethg(grid::Grid.RealGrid) = Nonlinear.Kerr_field(PhysData.γ3_gas(gas))
    makethg(grid::Grid.EnvGrid) = Nonlinear.Kerr_env_thg(PhysData.γ3_gas(gas), grid.ω0, grid.to)

    function runthg(grid)
        nfunλ = PhysData.ref_index_fun(gas, pr_thg)
        nfun = (λ; z=0.0) -> real(nfunλ(λ))
        linop = LinearOps.make_const_linop(grid, q, nfun, true) # thg=true for both grids
        dens0 = PhysData.density(gas, pr_thg)
        responses = (makethg(grid),)
        normfun = NonlinearRHS.const_norm_radial(grid, q, nfun)
        inputs = Fields.GaussGaussField(;λ0=λ0thg, τfwhm=20e-15, energy=1e-6, w0=100e-6)
        Eω, transform, FT = Luna.setup(grid, q, z -> dens0, normfun, responses, inputs)
        output = Output.MemoryOutput(0, Lthg, 3)
        Luna.run(Eω, grid, linop, transform, FT, output; init_dz=1e-5, zmax=Lthg)
        eω = Fields.energyfuncs(grid, q)[2]
        Eωk = output["Eω"][:, 1, :, end] # (ω, k)
        thgband = @. 2.9ω0 < grid.ω < 3.1ω0
        Em = copy(Eωk)
        Em[.!thgband, :] .= 0
        eω(Em), eω(Eωk)
    end

    ethg_r, etot_r = runthg(rgridt)
    ethg_e, etot_e = runthg(egridt)

    @test etot_r ≈ etot_e rtol=1e-6
    @test ethg_r > 1e5*eps(etot_r) # THG was actually generated
    @test ethg_r ≈ ethg_e rtol=0.05
end
##
# counts the accepted steps of a run, forwarding everything to a real output
mutable struct CountingOutput{oT}
    out::oT
    n::Int
end
CountingOutput(o) = CountingOutput(o, 0)
(c::CountingOutput)(Eω, z, dz, interp) = (c.n += 1; c.out(Eω, z, dz, interp))
(c::CountingOutput)(args...; kwargs...) = c.out(args...; kwargs...)
Base.getindex(c::CountingOutput, k) = c.out[k]

#= Evanescent channels: a transverse grid fine enough that k⊥,max exceeds k(ω) at the
   long-wavelength end of the band, and an input with a hard edge to populate those
   channels. The exact, uncapped decay rate reaches ~4e6 /m here; without the clamp and the
   matching source taper the stepper would collapse its step to ~1/κ or produce NaN. =#
@testset "evanescent channels" begin
    Re = 100e-6
    qe = Hankel.QDHT(Re, 128, dim=3)
    gride = Grid.RealGrid(800e-9, (400e-9, 4000e-9), 100e-15)
    ℓe = 2e-3/Boundaries.DEFAULT_N
    nfunλ = PhysData.ref_index_fun(gas, pressure)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    ωs = gride.ω[gride.sidx]
    k = [real(nfunλ(wlfreq(ω)))*ω/PhysData.c for ω in ωs]
    κ = sqrt.(max.(-(k.^2 .- (qe.k.^2)'), 0))
    @test maximum(κ) > 1e6 # the evanescent region exists and is stiff
    function runedge(aperture)
        linop = LinearOps.make_const_linop(gride, qe, nfun, true)
        dens0 = PhysData.density(gas, pressure)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
        normfun = NonlinearRHS.const_norm_radial(gride, qe, nfun)
        inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy=1e-9, w0=30e-6)
        Eω, transform, FT = Luna.setup(gride, qe, z -> dens0, normfun, responses, inputs)
        if aperture
            Eωr = qe \ Eω
            Eωr[:, :, qe.r .> 40e-6] .= 0
            Eω = qe * Eωr
        end
        output = CountingOutput(Output.MemoryOutput(0, 2e-3, 3))
        Luna.run(Eω, gride, linop, transform, FT, output; zmax=2e-3)
        output.n, output["Eω"][:, 1, :, end]
    end
    nsmooth, Es = runedge(false)
    nedge, Ee = runedge(true)
    @test !any(isnan, Ee)
    @test nedge <= 2nsmooth
    # channels the clamp touches decay by at least exp(-30) per reference length
    strong = κ .> Boundaries.MAX_αℓ/(2ℓe)
    @test count(strong) > 0
    Eb = Ee[gride.sidx, :]
    @test sum(abs2, Eb[strong])/sum(abs2, Eb) < 1e-12
end
##
#= Nothing reaches any edge: the absorbers must not change the answer. =#
@testset "rate and none agree away from the edges" begin
    function runfocus(boundary)
        nfunλ = PhysData.ref_index_fun(gas, pressure)
        nfun = (λ; z=0.0) -> nfunλ(λ)
        linop = LinearOps.make_const_linop(rgrid, q, nfun, true)
        dens0 = PhysData.density(gas, pressure)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
        normfun = NonlinearRHS.const_norm_radial(rgrid, q, nfun)
        inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0, propz=-L)
        Eω, transform, FT = Luna.setup(rgrid, q, z -> dens0, normfun, responses, inputs)
        output = Output.MemoryOutput(0, L, 3)
        Luna.run(Eω, rgrid, linop, transform, FT, output; init_dz=0.1, boundary, zmax=L)
        output["Eω"][:, 1, :, end]
    end
    Er = runfocus(:rate)
    En = runfocus(:none)
    @test norm(Er - En)/norm(En) < 1e-4
end
##
#= The transverse collar: a beam diverging into the aperture of the QDHT, whose E(R) = 0
   wall otherwise reflects it back into the beam. With the collar the inner part of the
   profile stays close to the analytic Gaussian; without it the reflection spoils it. =#
@testset "transverse collar" begin
    Ld = 0.1
    w0d = 30e-6
    gridd = Grid.RealGrid(λ0, (400e-9, 2000e-9), 0.2e-12)
    zR = π*w0d^2/λ0
    wL = w0d*sqrt(1 + (Ld/zR)^2)
    @test wL > 0.8R # the beam really does reach the wall
    function rundiverge(boundary)
        nfunλ = PhysData.ref_index_fun(gas, pressure)
        nfun = (λ; z=0.0) -> nfunλ(λ)
        linop = LinearOps.make_const_linop(gridd, q, nfun, true)
        dens0 = PhysData.density(gas, pressure)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
        normfun = NonlinearRHS.const_norm_radial(gridd, q, nfun)
        inputs = Fields.GaussGaussField(;λ0, τfwhm, energy, w0=w0d)
        Eω, transform, FT = Luna.setup(gridd, q, z -> dens0, normfun, responses, inputs)
        output = Output.MemoryOutput(0, Ld, 3)
        Luna.run(Eω, gridd, linop, transform, FT, output; init_dz=0.1, boundary, zmax=Ld)
        Eωr = q \ output["Eω"][:, :, :, end]
        dropdims(sum(abs2.(Eωr); dims=(1, 2)); dims=(1, 2))
    end
    inner = q.r .< 0.5R
    Ia = Maths.gauss.(q.r, wL/2)
    err(I) = norm(I[inner]/I[1] - Ia[inner]/Ia[1])/norm(Ia[inner]/Ia[1])
    er = err(rundiverge(:rate))
    en = err(rundiverge(:none))
    @info "transverse collar: inner-profile error with collar $er, without $en"
    @test er < en
    @test er < 0.05
end
