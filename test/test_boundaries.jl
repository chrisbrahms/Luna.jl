using Luna
import Test: @test, @testset, @test_throws
import FFTW
import LinearAlgebra: mul!, ldiv!
import Luna: Grid, Boundaries, Maths, RK45, Output, Processing, PhysData
import Luna.Interface: prop_capillary_args
import Logging
import Random: MersenneTwister

#= These tests cover the rate-based absorbing boundaries. The property they exist to
   protect is that the absorption depends on the propagation *distance* and not on how many
   steps the adaptive solver took to cover it -- the old scheme multiplied the solution by a
   fixed profile once per accepted step, so the answer depended on rtol. =#

@testset "profiles and rates" begin
grid = Grid.RealGrid(800e-9, (150e-9, 4e-6), 1e-12)
zmax = 0.3
ℓ = zmax/20

α = Boundaries.rate(grid.ωwin, ℓ)
@test all(isfinite, α)
@test all(α .>= 0)
@test all(α[grid.ωwin .== 1] .== 0)
@test maximum(α) <= Boundaries.MAX_αℓ/ℓ*(1 + 1e-12)
#= α is a power coefficient, so the field over one reference length goes as exp(-αℓ/2) and
   that is what must reproduce the profile, where it is not clamped. =#
unclamped = grid.ωwin .> exp(-Boundaries.MAX_αℓ/2)
@test isapprox(exp.(-α[unclamped].*ℓ./2), grid.ωwin[unclamped]; rtol=1e-12)

#= Telescoping: this is the whole point. Applying exp(-α dz) over any partition of [0, L]
   gives exactly exp(-α L), so the absorption cannot depend on the step layout. =#
i = findfirst(x -> 0 < x < 1, grid.ωwin)
for n in (3, 17, 200)
    dzs = rand(n)
    dzs .*= 0.3/sum(dzs)
    @test isapprox(prod(exp.(-α[i].*dzs./2)), exp(-α[i]*0.3/2); rtol=1e-12)
end

# ωwin is exactly zero outside the simulation band, so the hard mask is grid.sidx
@test (grid.ωwin .> 0) == grid.sidx
@test all(Boundaries.spectral_rate(grid, zmax)[.!grid.sidx] .== 0)
end

@testset "temporal collar is never degenerate" begin
#= grid.twin's collar is only the slack between the requested trange and the realised
   power-of-two window, so for some tranges it is a single sample wide -- those grids have
   no temporal absorber at all. tprofile guarantees a collar while never being weaker than
   grid.twin. =#
for (zmax, trange) in ((0.3, 1e-12), (0.3, 0.6825e-12), (1.0, 250e-15), (4e-5, 1e-12))
    for grid in (Grid.RealGrid(800e-9, (150e-9, 4e-6), trange),
                 Grid.EnvGrid(800e-9, (150e-9, 4e-6), trange))
        W = Boundaries.tprofile(grid)
        @test all(W .<= grid.twin .+ 1e-15) # never weaker than the historical profile
        @test count(0 .< W .< 1) > 0.03*length(grid.t) # and never degenerate
        @test all(isfinite, Boundaries.temporal_rate(grid, zmax))
    end
end
# a wider natural collar is kept as-is
grid = Grid.RealGrid(800e-9, (150e-9, 4e-6), 1e-12)
@test Boundaries.tprofile(grid; collar=0) == grid.twin
end

@testset "addloss" begin
grid = Grid.RealGrid(800e-9, (150e-9, 4e-6), 1e-12)
α = Boundaries.spectral_rate(grid, 0.3)
Nω = length(grid.ω)
# ω is axis 1 for every linop shape Luna uses
for sz in ((Nω,), (Nω, 4), (Nω, 8), (Nω, 8, 6))
    L = randn(ComplexF64, sz...)
    La = Boundaries.addloss(L, α)
    @test La isa AbstractArray # so linoptype stays "constant"
    @test size(La) == size(L)
    @test imag(La) == imag(L) # dispersion untouched
    trailing = ntuple(_ -> 1, length(sz) - 1)
    @test real(La)[:, trailing...] ≈ real(L)[:, trailing...] .- α./2 # linops carry -α/2
    # the closure form must agree with the array form
    Lc! = (out, z) -> (out .= L)
    Lac! = Boundaries.addloss(Lc!, α)
    out = similar(L)
    Lac!(out, 0.0)
    @test out ≈ La
    @test !isa(Lac!, AbstractArray) # so linoptype stays "variable"
end
end

@testset "propagator applies the rate exactly" begin
#= The spectral absorber rides the interaction-picture propagator, so exp(-α Δz) must be
   applied exactly for whatever sub-interval the stepper chooses. =#
grid = Grid.RealGrid(800e-9, (150e-9, 4e-6), 1e-12)
α = Boundaries.spectral_rate(grid, 0.3)
linop = Boundaries.addloss(zeros(ComplexF64, length(grid.ω)), α)
y = ones(ComplexF64, length(grid.ω))
prop! = RK45.make_prop!(linop, y)

onestep = copy(y)
prop!(onestep, 0.0, 0.3)
for n in (2, 7, 1000)
    many = copy(y)
    zs = range(0, 0.3, length=n+1)
    for j = 1:n
        prop!(many, zs[j], zs[j+1])
    end
    @test isapprox(many, onestep; rtol=1e-12)
end
@test isapprox(abs.(onestep), exp.(-α.*0.3./2); rtol=1e-12)
end

@testset "temporal absorber converges to the rate limit" begin
#= Zero RHS and zero linear operator, so the only thing acting on the field is the temporal
   absorber applied in Luna.run's stepfun, whose target is exp(-α_t L) applied once.

   The exponential factors telescope exactly among themselves, but the spectral absorber
   rides the propagator in between them and the two do not commute, so the result approaches
   exp(-α_t L/2) as the steps get smaller rather than matching it exactly. That is a
   splitting error of the same order as the temporal splitting itself, and — unlike the
   historical scheme — it converges. =#
grid = Grid.EnvGrid(800e-9, (400e-9, 2e-6), 1e-12)
zmax = 0.2
FT = FFTW.plan_fft(zeros(ComplexF64, length(grid.t)))
αt = Boundaries.temporal_rate(grid, zmax)
mask = float.(grid.sidx) # the band limit Luna.run applies once, at the start
# deliberately wide, so that the field actually overlaps the absorber collar
Et0 = Maths.gauss.(grid.t, fwhm=0.5*(maximum(grid.t) - minimum(grid.t))) .+ 0im
Eω0 = (FT*Et0) .* mask

function bare_run(; kwargs...)
    Eω = copy(Eω0)
    out = Output.MemoryOutput(0, zmax, 3)
    Logging.with_logger(Logging.NullLogger()) do
        Luna.run(Eω, grid, zeros(ComplexF64, length(grid.ω)),
                 (nl, Eω, z) -> fill!(nl, 0), FT, out; zmax, kwargs...)
    end
    out["Eω"][:, end]
end
reldiff(a, b) = sqrt(sum(abs2, a .- b))/sqrt(sum(abs2, b))

expected = (FT*((FT \ Eω0) .* exp.(-αt.*zmax./2))) .* mask
coarse = reldiff(bare_run(init_dz=1e-3), expected)                      # 20 steps
fine = reldiff(bare_run(init_dz=1e-5, max_dz=zmax/500), expected)  # 500 steps
@test coarse < 1e-3
@test fine < 10*coarse/25 # converges at least first order in the step size

# with nothing driving it and no absorber, :none must leave the field exactly alone
@test bare_run(boundary=:none) == Eω0
# :legacy applies the whole window once per step and does not approach the rate limit
legacy = Logging.with_logger(Logging.NullLogger()) do
    reldiff(bare_run(boundary=:legacy), expected)
end
@test legacy > 100*fine
end

@testset "legacy reproduces the historical scheme bit-identically" begin
#= boundary=:legacy must be an exact reproduction of the pre-change behaviour, so that
   published results stay reproducible. Compare against the historical stepfun driven
   straight through RK45, with no other part of `run` involved. =#
grid = Grid.EnvGrid(800e-9, (400e-9, 2e-6), 1e-12)
zmax = 0.2
FT = FFTW.plan_fft(zeros(ComplexF64, length(grid.t)))
Et0 = Maths.gauss.(grid.t, fwhm=0.5*(maximum(grid.t) - minimum(grid.t))) .+ 0im
Eω0 = FT*Et0
linop = -im.*(1e4.*(grid.ω .- grid.ω0).^2) # some dispersion, so the steps do something
transform = (nl, Eω, z) -> fill!(nl, 0)
maxdz = zmax/50 # force ~50 window applications, so bit-identity is a real statement

Eωa = copy(Eω0)
outa = Output.MemoryOutput(0, zmax, 3)
Logging.with_logger(Logging.NullLogger()) do
    Luna.run(Eωa, grid, linop, transform, FT, outa; zmax, boundary=:legacy, max_dz=maxdz)
end

Eωb = copy(Eω0)
outb = Output.MemoryOutput(0, zmax, 3)
Et = FT \ Eωb
function historical_stepfun(Eω, z, dz, interpolant)
    Eω .*= grid.ωwin
    ldiv!(Et, FT, Eω)
    Et .*= grid.twin
    mul!(Eω, FT, Et)
    outb(Eω, z, dz, interpolant)
end
Logging.with_logger(Logging.NullLogger()) do
    RK45.solve_precon(transform, linop, Eωb, 0.0, 1e-4, zmax;
                      stepfun=historical_stepfun, max_dt=maxdz, min_dt=0,
                      rtol=1e-6, atol=1e-10, safety=0.9, norm=RK45.weaknorm,
                      status_period=1)
end

@test outa["Eω"] == outb["Eω"] # bit-identical, not merely approximately equal
end

@testset "tolerance independence" begin
#= The measured defect: with the historical per-step windowing, tightening rtol changes the
   answer because it changes the number of times the window is applied. Use a case with
   real energy in the taper collars -- a strongly self-broadening capillary run. =#
function arm(boundary, rtol)
    Eω, grid, linop, transform, FT, output = prop_capillary_args(
        125e-6, 0.15, :He, 3.0; λ0=800e-9, energy=300e-6, τfwhm=10e-15,
        λlims=(120e-9, 4e-6), trange=1e-12, saveN=3, plasma=false, boundary,
        PPT_options=Dict(:cache => false))
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=0.15, rtol, boundary)
    #= Compare only inside the flat part of the window: interpolated saves reconstruct the
       out-of-band part from the stepper's stages, which is rtol-dependent by construction
       and is not what this test is about. =#
    output["Eω"][grid.ωwin .== 1, end]
end

reldiff(a, b) = sqrt(sum(abs2, a .- b))/sqrt(sum(abs2, a))

rate, legacy = Logging.with_logger(Logging.NullLogger()) do
    (reldiff(arm(:rate, 1e-6), arm(:rate, 1e-9)),
     reldiff(arm(:legacy, 1e-6), arm(:legacy, 1e-9)))
end
@info "tolerance sensitivity: :rate $rate, :legacy $legacy"
@test rate < 1e-5
@test legacy > 20*rate # the mechanism is genuinely exercised
end

@testset "temporal walk-off" begin
#= The temporal absorber's job is to destroy light that walks off the end of the time
   window *before* it wraps around and aliases onto the rest of the field. Nothing above
   exercises that: the temporal tests use a zero linear operator, so nothing moves.

   Use the case where this actually bites in practice -- a deep-UV resonant dispersive wave
   in a deliberately small time window. The DW is emitted at ~230 nm around z = 0.2 m and
   then walks away from the pump at the group-delay difference of the fibre (~100 fs/m
   here), so in a 171 fs window it leaves the flat region at ~0.7 m and reaches the window
   edge at ~1.0 m, i.e. halfway along the 2 m fibre. =#
RDWBAND = (200e-9, 260e-9)
ZMAX = 2.0

function walkoff_arm(boundary; max_dz=ZMAX/2, rtol=1e-6)
    Eω, grid, linop, transform, FT, output = prop_capillary_args(
        50e-6, ZMAX, :Ar, 0.4; λ0=800e-9, energy=29e-6, τfwhm=10e-15,
        λlims=(150e-9, 4e-6), trange=100e-15, saveN=51, plasma=false, boundary,
        rng=MersenneTwister(1234))
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=ZMAX, rtol, boundary, max_dz)
    output
end

#= Energy in the DW band, via a hard 0/1 spectral mask. Processing.energy's tuple bandpass
   picks its taper width from the frequency spacing, which differs between grids; a hard
   mask is exact by Parseval and comparable across everything here. =#
function bandE(output)
    grid = Processing.makegrid(output)
    mask = float.(RDWBAND[1] .< PhysData.wlfreq.(grid.ω) .< RDWBAND[2])
    Processing.energy(grid, output["Eω"]; bandpass=mask)
end

# power in the leading half of the window, where anything that wraps around must appear
function leadingpower(output)
    t, Et = Processing.getEt(output, ZMAX; bandpass=(RDWBAND..., 1e14))
    maximum(abs2.(Et[t .< -30e-15, 1]))
end

rate, legacy, none, ratefine, legacyfine = Logging.with_logger(Logging.NullLogger()) do
    (walkoff_arm(:rate), walkoff_arm(:legacy), walkoff_arm(:none),
     walkoff_arm(:rate; max_dz=ZMAX/4000), walkoff_arm(:legacy; max_dz=ZMAX/4000))
end

z = rate["z"]
i0 = argmin(abs.(z .- 0.5)) # before the DW has reached the taper
Erate, Elegacy, Enone = bandE(rate), bandE(legacy), bandE(none)

# the absorber does nothing at all until the pulse gets there
@test isapprox(Erate[i0], Enone[i0]; rtol=1e-2)

# ... and then removes the DW: two orders of magnitude, where :none keeps most of it
@test Erate[end]/Erate[i0] < 1e-2
@test Elegacy[end]/Elegacy[i0] < 1e-2
@test Enone[end]/Enone[i0] > 0.5

#= What :none keeps is not physics, it is wraparound: with no absorber the DW reappears at
   the *front* of the time window. That is the aliasing the boundary exists to prevent. =#
@test leadingpower(none) > 100*leadingpower(rate)
@test leadingpower(none) > 100*leadingpower(legacy)

#= The point of rate semantics: how much is absorbed depends on the distance travelled, not
   on how the solver chose to subdivide it. Same rtol in both arms -- so the same local
   error control and the same physics -- with max_dz forcing ~3.5x as many steps, hence
   3.5x as many applications of the historical window. =#
dr = abs(bandE(ratefine)[end] - Erate[end])/Erate[end]
dl = abs(bandE(legacyfine)[end] - Elegacy[end])/Elegacy[end]
@info "walk-off, step-count sensitivity of the surviving DW: :rate $dr, :legacy $dl"
@test dr < 0.05
@test dl > 0.1        # the mechanism is genuinely exercised
@test dl > 10*dr

#= The absorber reports what it removed. This has to fire here -- the DW really is being
   eaten -- and it has to stay quiet when nothing reaches the boundary. The old diagnostic
   predicted from the group delay across the band instead of measuring, which on a grid
   running to 4 um is dominated by band-edge components carrying nothing but shot noise, so
   it warned on runs where the absorber never touched anything real. =#
warnings = String[]
mutable struct CollectWarnings <: Logging.AbstractLogger; msgs::Vector{String}; end
Logging.min_enabled_level(::CollectWarnings) = Logging.Info
Logging.shouldlog(::CollectWarnings, args...) = true
Logging.catch_exceptions(::CollectWarnings) = false
function Logging.handle_message(l::CollectWarnings, lvl, msg, _m, g, id, file, line; kw...)
    occursin("Temporal absorbing boundary", string(msg)) && push!(l.msgs, string(msg))
    nothing
end

function warnings_from(; trange, flength)
    lg = CollectWarnings(String[])
    Logging.with_logger(lg) do
        Eω, grid, linop, transform, FT, output = prop_capillary_args(
            50e-6, flength, :Ar, 0.4; λ0=800e-9, energy=29e-6, τfwhm=10e-15,
            λlims=(150e-9, 4e-6), trange, saveN=11, plasma=false,
            rng=MersenneTwister(1234))
        Luna.run(Eω, grid, linop, transform, FT, output; zmax=flength)
    end
    lg.msgs
end

# the DW is well into the collar by 2 m, so this must warn -- exactly once
warned = warnings_from(trange=100e-15, flength=ZMAX)
@test length(warned) == 1
@test occursin("removed", warned[1])
# stop before the DW gets there and there is nothing to report
@test isempty(warnings_from(trange=100e-15, flength=0.3))
end

@testset "run interface" begin
grid = Grid.EnvGrid(800e-9, (400e-9, 2e-6), 1e-12)
zmax = 0.1

@test_throws ErrorException Logging.with_logger(Logging.NullLogger()) do
    FT = FFTW.plan_fft(zeros(ComplexF64, length(grid.t)))
    Luna.run(zeros(ComplexF64, length(grid.ω)), grid,
             zeros(ComplexF64, length(grid.ω)), (nl, Eω, z) -> fill!(nl, 0), FT,
             Output.MemoryOutput(0, zmax, 3); zmax, boundary=:nonsense)
end

#= A non-positive reference length is silently destructive rather than merely useless -- it
   NaNs the whole grid or turns the absorber into a gain -- so reject it at the one place
   every call site goes through. =#
@test Boundaries.reflength(zmax, 20, nothing) == zmax/20
@test Boundaries.reflength(zmax, 20, 0.05) == 0.05
for bad in (0, -1, Inf, NaN)
    @test_throws ErrorException Boundaries.reflength(zmax, bad, nothing)
    @test_throws ErrorException Boundaries.reflength(zmax, 20, bad)
end

# the chosen mode is recorded, so a saved run can be reproduced
out = Logging.with_logger(Logging.NullLogger()) do
    prop_capillary(125e-6, 0.05, :He, 1.0; λ0=800e-9, energy=1e-6, τfwhm=10e-15,
                   λlims=(150e-9, 4e-6), trange=1e-12, saveN=3, plasma=false,
                   boundary_N=5, PPT_options=Dict(:cache => false))
end
@test out["simulation_type"]["boundary"] == "rate"
@test out["prop_capillary_args"]["boundary_N"] == "5"
end

import Luna: LinearOps, NonlinearRHS, Fields, Nonlinear
import Luna.PhysData: wlfreq
Luna.set_fftw_mode(:estimate)

@testset "free space" begin
Rs = 50e-6
grid = Grid.RealGrid(800e-9, (400e-9, 4000e-9), 0.2e-12)
zmax = 1e-3
q = Grid.RadialGrid(Rs, 32)
xgrid = Grid.Free2DGrid(Rs, 32)
xygrid = Grid.FreeGrid(Rs, 16, Rs, 16)
nfunλ = PhysData.ref_index_fun(:Ar, 1)
nfun = (λ; z=0.0) -> nfunλ(λ)
nfunω = (ω; z) -> nfun(wlfreq(ω); z)
ℓ = zmax/20
ratemax = Boundaries.MAX_αℓ/(2ℓ)

# profiles: sized like the k axes, 1 in the interior, tapering to 0 at the edge
for sg in (q, xgrid, xygrid)
    kperp2, _ = LinearOps.transverse_k2(sg)
    Wk = Boundaries.kprofile(sg, 0.1)
    @test size(Wk) == size(kperp2)
    @test all(0 .<= Wk .<= 1)
    @test maximum(Wk) == 1
    @test minimum(Wk) < 0.5
    @test all(Boundaries.kprofile(sg, 0.0) .== 1) # no collar
    Wr = Boundaries.rprofile(sg, 0.1)
    @test all(0 .<= Wr .<= 1)
    @test maximum(Wr) == 1
    @test minimum(Wr) < 0.5
end
@test all(Boundaries.rprofile(q, 0.1)[q.r .< 0.85Rs] .== 1)
@test all(Boundaries.rprofile(q, 0.1)[q.r .> 0.99Rs] .< 0.1)
@test Boundaries.rprofile(xgrid, 0.1) == xgrid.xwin

# addloss_k broadcasts over the trailing k axes; clampdecay caps the evanescent decay
linop = LinearOps.make_const_linop(grid, q, nfun, true)
α = Boundaries.rate(Boundaries.kprofile(q, 0.1), ℓ)
l2 = Boundaries.addloss_k(linop, α)
@test size(l2) == size(linop)
@test l2[:, 1, end] ≈ linop[:, 1, end] .- α[end]/2
@test l2[:, 1, 1] == linop[:, 1, 1]
@test minimum(real(linop)) < -ratemax # the raw operator decays faster than the clamp
lc = Boundaries.clampdecay(linop, ratemax)
@test minimum(real(lc)) ≈ -ratemax
@test imag(lc) == imag(linop)
keep = real(linop) .> -ratemax
@test real(lc)[keep] == real(linop)[keep]
linopf = LinearOps.make_linop(grid, q, nfunω, true)
out = similar(linop)
Boundaries.addloss_k(linopf, α)(out, 0.0)
@test out ≈ l2
Boundaries.clampdecay(linopf, ratemax)(out, 0.0)
@test out ≈ lc

# FreeSpaceNorm: physical prefactor βz/(μ0 ω) on both sides of cutoff, then the taper
nf = NonlinearRHS.const_norm_radial(grid, q, nfunλ)
n0 = copy(nf(0.0))
@test eltype(n0) == ComplexF64
ωs = grid.ω[grid.sidx]
k = [real(nfunλ(wlfreq(ω)))*ω/PhysData.c for ω in ωs]
βsq = k.^2 .- (q.k.^2)'
ωa = ωs .* ones(1, q.N)
evan = βsq .< 0
nn = n0[grid.sidx, 1, :]
@test count(evan) > 0
@test all(nn[.!evan] .≈ (sqrt.(max.(βsq, 0))./(PhysData.μ_0 .* ωa))[.!evan])
@test all(nn[evan] .≈ (-im .* sqrt.(max.(-βsq, 0))./(PhysData.μ_0 .* ωa))[evan])
@test all(n0[.!grid.sidx, :, :] .== 1)
Wk = ones(q.N)
Wk[end] = 0.5
NonlinearRHS.reflength!(nf, ℓ; κmax=ratemax, kwin=Wk)
n1 = copy(nf(0.0))
κ = min.(sqrt.(max.(-βsq, 0)), ratemax)
W = exp.(-κ.*ℓ) .* Wk'
@test n1[grid.sidx, 1, :] ≈ n0[grid.sidx, 1, :] ./ W
@test !(n1 ≈ n0)
@test maximum(abs.(n1[grid.sidx, 1, :] ./ n0[grid.sidx, 1, :])) ≈ 2exp(Boundaries.MAX_αℓ/2) # clamped channels
NonlinearRHS.reflength!(nf, 0.0)
@test nf(0.0) ≈ n0
# a z-dependent normalisation re-tapers on every call
nfz = NonlinearRHS.norm_radial(grid, q, nfunω)
@test nfz(0.0) ≈ n0
NonlinearRHS.reflength!(nfz, ℓ; κmax=ratemax, kwin=Wk)
@test nfz(0.5e-3) ≈ n1
# and the crystal-optics form gives the same physics for an isotropic index
nfb = NonlinearRHS.const_norm_free2D(grid, xgrid, ((λ, δθ) -> real(nfunλ(λ)), λ -> real(nfunλ(λ))))
nfi = NonlinearRHS.const_norm_free2D(grid, xgrid, (λ -> (real(nfunλ(λ)), real(nfunλ(λ)))))
@test nfb(0.0) ≈ nfi(0.0)

# setup end-to-end: what Luna.run gets back for a free-space transform
dummy(args...; kwargs...) = nothing
resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:Ar)),)
inputs = Fields.GaussGaussField(;λ0=800e-9, τfwhm=20e-15, energy=1e-9, w0=20e-6)
for (sg, nfs) in ((q, NonlinearRHS.const_norm_radial(grid, q, nfunλ)),
                  (xgrid, NonlinearRHS.const_norm_free2D(grid, xgrid, nfunλ)),
                  (xygrid, NonlinearRHS.const_norm_free(grid, xygrid, nfunλ)))
    Eω, transform, FT = Luna.setup(grid, sg, z -> 1e25, nfs, resp, inputs)
    @test Boundaries.spacegrid(transform) === sg
    Et = FT \ Eω
    lin = LinearOps.make_const_linop(grid, sg, nfun, true)
    b = Boundaries.setup(:rate, grid, transform, lin, Et, FT, dummy, 0.0, zmax, zmax/2, 1e-4)
    @test b.max_dz == ℓ
    @test nfs.ℓ == ℓ && nfs.κmax == ratemax
    @test minimum(nfs.kwin) >= exp(-Boundaries.MAX_αℓ/2) # never zero: the norm divides by it
    @test minimum(real(b.linop)) >= -ratemax - 1.5Boundaries.MAX_αℓ/ℓ - 1e-9
    @test b.stepfun.spatial isa (sg isa Grid.RadialGrid ? Boundaries.RadialCollar : Boundaries.CartesianCollar)
    # one accepted step through the absorber: finite, and only ever removes energy
    E1 = copy(Eω)
    e0 = sum(abs2, E1)
    b.stepfun(E1, 1e-4, 1e-4, nothing)
    @test all(isfinite, E1)
    @test sum(abs2, E1) <= e0
    # without absorbers the evanescent clamp and taper are still applied, over max_dz
    b2 = Boundaries.setup(:none, grid, transform, lin, Et, FT, dummy, 0.0, zmax, zmax/2, 1e-4)
    @test b2.stepfun isa Boundaries.NoAbsorber
    @test nfs.ℓ == zmax/2
    @test minimum(real(b2.linop)) ≈ max(minimum(real(lin)), -Boundaries.MAX_αℓ/zmax)
    @test all(nfs.kwin .== 1)
end
# a modal transform is untouched
cgrid = Grid.RealGrid(800e-9, (150e-9, 4e-6), 1e-12)
Eω, cgrid, linop, transform, FT, output = prop_capillary_args(125e-6, 0.3, :He, 1.0;
    λ0=800e-9, energy=1e-9, τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12, saveN=11)
@test isnothing(Boundaries.spacegrid(transform))
b = Boundaries.setup(:rate, cgrid, transform, linop, FT \ Eω, FT, dummy, 0.0, 0.3, 0.15, 1e-4)
@test isnothing(b.stepfun.spatial)
@test b.linop == Boundaries.addloss(linop, Boundaries.spectral_rate(cgrid, 0.3))
end

#= gpu/11: RateAbsorber and LegacyAbsorber are broadcasts and reductions over mirrored
   arrays now, so they should run on any array type. Skipped without JLArrays, which is
   not part of the main test environment (see test_device.jl for the same gate). =#
have_jlarrays = try
    @eval import JLArrays
    true
catch
    false
end

if !have_jlarrays
    @warn "JLArrays is not available; the JLArray boundaries tests are skipped. Run "*
          "through `Pkg.test()` or add JLArrays to the environment."
else

import AbstractFFTs
import GPUArraysCore

#= The same host-backed AbstractFFTs shim test_device.jl defines: JLArrays has no FFT
   plans of its own. Duplicated rather than shared, since the two test files are
   `include`d independently and neither depends on the other. =#
mutable struct BoundaryJLPlan{T, N, P} <: AbstractFFTs.Plan{T}
    hp::P
    sz::NTuple{N, Int}
    dims::Any
    pinv::AbstractFFTs.ScaledPlan
    BoundaryJLPlan{T, N, P}(hp, sz, dims) where {T, N, P} = new{T, N, P}(hp, sz, dims)
end
BoundaryJLPlan(hp, sz::NTuple{N, Int}, dims, T=ComplexF64) where {N} =
    BoundaryJLPlan{T, N, typeof(hp)}(hp, sz, dims)
Base.size(p::BoundaryJLPlan) = p.sz
Base.eltype(::BoundaryJLPlan{T}) where {T} = T
AbstractFFTs.plan_fft(x::JLArrays.JLArray{ComplexF64}, dims) =
    BoundaryJLPlan(FFTW.plan_fft(Array(x), dims), size(x), dims)
AbstractFFTs.plan_rfft(x::JLArrays.JLArray{Float64}, dims) =
    BoundaryJLPlan(FFTW.plan_rfft(Array(x), dims), size(x), dims, Float64)
AbstractFFTs.plan_inv(p::BoundaryJLPlan) =
    AbstractFFTs.ScaledPlan(BoundaryJLPlan(inv(p.hp).p, p.sz, p.dims, ComplexF64),
                            AbstractFFTs.normalization(Float64, p.sz, p.dims))
Base.:*(p::BoundaryJLPlan, x::JLArrays.JLArray) = JLArrays.JLArray(p.hp * Array(x))
mul!(y::JLArrays.JLArray, p::BoundaryJLPlan, x::JLArrays.JLArray) =
    (copyto!(y, p.hp * Array(x)); y)

const BJLArray = JLArrays.JLArray

_bnoop(args...; kwargs...) = nothing

@testset "RateAbsorber and LegacyAbsorber on JLArray" begin
    GPUArraysCore.allowscalar(false)
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    zmax = 1e-2
    dz = zmax/20
    αt = Boundaries.temporal_rate(grid, zmax)
    Nt = length(grid.t)
    Eω0 = randn(ComplexF64, length(grid.ω))

    Et_h = zeros(Float64, Nt)
    FT_h = FFTW.plan_rfft(Et_h, 1)
    Et_d = BJLArray(zeros(Float64, Nt))
    FT_d = AbstractFFTs.plan_rfft(Et_d, 1)

    for (absorber_h, absorber_d) in (
            (Boundaries.RateAbsorber(αt, Et_h, FT_h, _bnoop, 0.0),
             Boundaries.RateAbsorber(αt, Et_d, FT_d, _bnoop, 0.0)),
            (Boundaries.LegacyAbsorber(grid, Et_h, FT_h, _bnoop),
             Boundaries.LegacyAbsorber(grid, Et_d, FT_d, _bnoop)))
        Eω_h = copy(Eω0)
        absorber_h(Eω_h, dz, dz, nothing)
        Eω_d = BJLArray(copy(Eω0))
        absorber_d(Eω_d, dz, dz, nothing)
        @test maximum(abs, Array(Eω_d) .- Eω_h)/maximum(abs, Eω_h) < 1e-10
    end
end

end # have_jlarrays
