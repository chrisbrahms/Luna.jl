import Luna
import Luna: PhysData, Grid, LinearOps, Modes, Capillary, RK45, Boundaries
import Test: @testset, @test, @test_logs, @test_throws
import Luna.PhysData: wlfreq

R = 5e-3
Nr = 256
Nx = 128
Ny = 64
gas = :Ar
pressure = 1

@testset "free space" begin
    rgrid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    egrid = Grid.EnvGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    q = Grid.RadialGrid(R, Nr)
    xygrid = Grid.FreeGrid(R, Nx, R, Ny)
    xgrid = Grid.Free2DGrid(R, Nx)

    getshape(grid, q::Grid.RadialGrid, pol) = (length(grid.ω), pol ? 2 : 1, q.N)
    getshape(grid, sg::Grid.Free2DGrid, pol) = (length(grid.ω), pol ? 2 : 1, length(sg.x))
    getshape(grid, sg::Grid.FreeGrid, pol) = (length(grid.ω), pol ? 2 : 1, length(sg.x), length(sg.y))

    @testset "$(typeof(grid)), $(typeof(sg)), pol = $pol, thg = $thg" for sg in (q, xgrid, xygrid),
                                                              pol in (false, true),
                                                              thg in (false, true),
                                                              grid in (rgrid, egrid)
        if grid isa Grid.RealGrid && ~thg
            continue
        end
        nfunλ = PhysData.ref_index_fun(gas, pressure)
        if pol
            nfun = (λ; z=0.0) -> (nfunλ(λ), nfunλ(λ))
        else
            nfun = (λ; z=0.0) -> nfunλ(λ)
        end
        nfunω = (ω; z) -> nfun(wlfreq(ω); z)

        linop = LinearOps.make_const_linop(grid, sg, nfun, thg)
        linopf = LinearOps.make_linop(grid, sg, nfunω, thg)
        out = similar(linop)

        @test size(linop) == getshape(grid, sg, pol)

        linopf(out, 0.0)
        @test all(imag(out) .≈ imag(linop))
        @test all(real(out) .≈ real(linop))
        linopf(out, 0.5)
        @test all(imag(out) .≈ imag(linop))
        @test all(real(out) .≈ real(linop))
    end
end

@testset "free space birefringent (tuple nfuns)" begin
    rgrid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    egrid = Grid.EnvGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    xygrid = Grid.FreeGrid(R, Nx, R, Ny)
    xgrid = Grid.Free2DGrid(R, Nx)

    #=
    Isotropic index: the tuple (crystal) path must match the generic vector-nfun path.
    The tuple path always subtracts β1*ω (a pure time shift, transparent to carrier-mixing
    nonlinearities), whereas the generic path subtracts β1*(ω - ω0) + β0—so for EnvGrid the
    two differ by the constant frame phase β1*ω0 - β0.
    =#
    nfunλ = PhysData.ref_index_fun(gas, pressure)
    nfunx = (λ, δθ=0.0) -> real(nfunλ(λ))
    nfuny = λ -> real(nfunλ(λ))
    nfun = (λ; z=0.0) -> (real(nfunλ(λ)), real(nfunλ(λ)))

    @testset "isotropic: $(typeof(grid)), $(typeof(sg)), thg = $thg" for sg in (xgrid, xygrid),
                                                              thg in (false, true),
                                                              grid in (rgrid, egrid)
        if grid isa Grid.RealGrid && ~thg
            continue
        end
        linop = LinearOps.make_const_linop(grid, sg, (nfunx, nfuny))
        linopv = LinearOps.make_const_linop(grid, sg, nfun, thg)
        @test size(linop) == size(linopv)
        β1 = PhysData.dispersion_func(1, nfuny)(grid.referenceλ)
        ω0 = LinearOps.getω0(grid, thg)
        β0 = LinearOps.getβ0_n(grid, λ -> nfun(λ), thg)
        offset = im*(β1*ω0 - β0) # 0 for RealGrid and for EnvGrid with thg=true
        # the tuple path only fills frequencies within grid.sidx
        sel = ntuple(_ -> Colon(), ndims(linop) - 1)
        @test all(isapprox.(linop[grid.sidx, sel...], linopv[grid.sidx, sel...] .+ offset;
                            atol=1e-3, rtol=0))
    end

    # for an envelope grid the phase subtracted is β1*ω: check the referencing directly
    # for the y polarisation at kperp = 0
    linop = LinearOps.make_const_linop(egrid, xgrid, (nfunx, nfuny))
    ik0 = argmin(abs.(xgrid.kx))
    β1 = PhysData.dispersion_func(1, nfuny)(egrid.referenceλ)
    for iω in (argmin(abs.(egrid.ω .- egrid.ω0)), findfirst(egrid.sidx))
        ωi = egrid.ω[iω]
        expected = -(nfuny(PhysData.wlfreq(ωi))*ωi/PhysData.c - β1*ωi)
        @test isapprox(imag(linop[iω, 2, ik0]), expected; atol=1e-6, rtol=0)
        @test real(linop[iω, 2, ik0]) == 0
    end

    # real birefringent crystal on an envelope grid
    θ = deg2rad(29.2)
    bbogrid = Grid.EnvGrid(800e-9, (250e-9, 2e-6), 120e-15; thg=true)
    bboxgrid = Grid.Free2DGrid(80e-6, 32)
    nfuns = PhysData.ref_index_fun_xy(:BBO, θ)
    linop = LinearOps.make_const_linop(bbogrid, bboxgrid, nfuns)
    @test size(linop) == (length(bbogrid.ω), 2, length(bboxgrid.kx))
    @test all(isfinite, linop)
    # birefringence: the two polarisations see different indices
    @test any(linop[bbogrid.sidx, 1, :] .!= linop[bbogrid.sidx, 2, :])
end

@testset "equivalence for fast z-dependent linops" begin
a = 125e-6
L = 1
grid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.5e-12)
coren, densityfun = Capillary.gradient(gas, L, pressure, 0)
m = Capillary.MarcatiliMode(a, coren)
dm = Modes.delegated(m) # delegated mode tricks make_linop into using the generic version

lom!, βm! = LinearOps.make_linop(grid, m, 800e-9)
lodm!, βdm! = LinearOps.make_linop(grid, dm, 800e-9)
@assert typeof(lom!) != typeof(lodm!) # ...but best to check

outm = complex(similar(grid.ω))
outdm = complex(similar(grid.ω))
for zi in range(0, L, length=10)
    lom!(outm, zi)
    lodm!(outdm, zi)
    @test outm == outdm
    βm!(outm, zi)
    βdm!(outdm, zi)
    @test outm == outdm
end

a = 125e-6
L = 1
# NO THG
thg = false
grid = Grid.EnvGrid(800e-9, (400e-9, 2000e-9), 0.5e-12; thg=thg)
coren, densityfun = Capillary.gradient(gas, L, pressure, 0)
m = Capillary.MarcatiliMode(a, coren)
dm = Modes.delegated(m) # delegated mode tricks make_linop into using the generic version...

lom!, βm! = LinearOps.make_linop(grid, m, 800e-9; thg=thg)
lodm!, βdm! = LinearOps.make_linop(grid, dm, 800e-9; thg=thg)
@assert typeof(lom!) != typeof(lodm!) # ...but best to check

outm = complex(similar(grid.ω))
outdm = complex(similar(grid.ω))
for zi in range(0, L, length=10)
    lom!(outm, zi)
    lodm!(outdm, zi)
    @test outm == outdm
    βm!(outm, zi)
    βdm!(outdm, zi)
    @test outm == outdm
end
# WITH THG
thg = true
grid = Grid.EnvGrid(800e-9, (400e-9, 2000e-9), 0.5e-12; thg=thg)
coren, densityfun = Capillary.gradient(gas, L, pressure, 0)
m = Capillary.MarcatiliMode(a, coren)
dm = Modes.delegated(m) # delegated mode tricks make_linop into using the generic version...

lom!, βm! = LinearOps.make_linop(grid, m, 800e-9; thg=thg)
lodm!, βdm! = LinearOps.make_linop(grid, dm, 800e-9; thg=thg)
@assert typeof(lom!) != typeof(lodm!) # ...but best to check

outm = complex(similar(grid.ω))
outdm = complex(similar(grid.ω))
for zi in range(0, L, length=10)
    lom!(outm, zi)
    lodm!(outdm, zi)
    @test outm == outdm
    βm!(outm, zi)
    βdm!(outdm, zi)
    @test outm == outdm
end
end

#= Evanescent region: components with k⊥ > k(ω) decay at exactly sqrt(k⊥² - k²), with no
   cap, and carry no propagation phase; out of band the operator is exactly zero. A small
   aperture makes k⊥,max exceed k(ω) at the long-wavelength end of the band. =#
@testset "evanescent region" begin
    @test LinearOps.βz(4.0) == 2
    @test LinearOps.βz(-4.0) == -2im
    @test LinearOps.βz(0.0) == 0

    Rs = 50e-6
    grid = Grid.RealGrid(800e-9, (400e-9, 4000e-9), 0.2e-12)
    qs = Grid.RadialGrid(Rs, 64)
    xygrids = Grid.FreeGrid(Rs, 32, Rs, 32)
    xgrids = Grid.Free2DGrid(Rs, 64)
    nfunλ = PhysData.ref_index_fun(gas, pressure)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    nfunω = (ω; z) -> nfun(wlfreq(ω); z)
    β1 = PhysData.dispersion_func(1, λ -> nfun(λ)[end])(grid.referenceλ)
    ωs = grid.ω[grid.sidx]
    k = [real(nfunλ(wlfreq(ω)))*ω/PhysData.c for ω in ωs]
    for sg in (qs, xgrids, xygrids)
        linop = LinearOps.make_const_linop(grid, sg, nfun, true)
        kperp2, idcs = LinearOps.transverse_k2(sg)
        nd = ndims(kperp2)
        cols = ntuple(_ -> :, nd)
        lin = linop[grid.sidx, 1, cols...]
        ωa = reshape(ωs, :, ntuple(_ -> 1, nd)...)
        βsq = reshape(k.^2, :, ntuple(_ -> 1, nd)...) .- reshape(kperp2, 1, size(kperp2)...)
        evan = βsq .< 0
        @test count(evan) > 0
        @test all(real(lin[evan]) .== -sqrt.(-βsq[evan]))
        @test all(imag(lin[evan]) .≈ (β1 .* ωa .* ones(size(βsq)))[evan])
        @test all(real(lin[.!evan]) .== 0)
        @test all(imag(lin[.!evan]) .≈ -(sqrt.(max.(βsq, 0)) .- β1 .* ωa)[.!evan])
        @test minimum(real(lin)) < -1e5 # far beyond the 200/m cap this replaces
        @test all(linop[.!grid.sidx, :, cols...] .== 0)
        linopf = LinearOps.make_linop(grid, sg, nfunω, true)
        out = similar(linop)
        linopf(out, 0.0)
        @test out ≈ linop
    end
end

#= Adaptive tabulation of a z-dependent operator (GPU_PLAN.md section 4.5 layer 2). What
   is checked here is the table and the propagator built from it, against an operator whose
   integral is known in closed form and against Luna's own z-dependent operators; what
   tabulation does to a whole propagation is in `test_device.jl` and `test_metal.jl`. =#
@testset "tabulated linear operator" begin
    L = 0.1
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    nω = length(grid.ω)
    proto = zeros(ComplexF64, nω)
    zend = 1.05L # zmax + max_dz, which is the span Luna.run tabulates over

    #= An operator with the awkwardness Luna's own have and an integral in closed form: a
       √z term, whose z derivative is infinite at z = 0, exactly as the density of a
       pressure gradient filled from zero behaves. The constants are of the size a
       capillary's are -- a propagation constant of ~1e7 rad/m and a loss of ~1/m -- so the
       absolute tolerance means the same thing here as it does there. =#
    k0 = @. 1e7*(1 + 0.1*sin(1:nω))
    k1 = @. 1e4*(1 + 0.5*cos(1:nω))
    k2 = @. 1.0*(1 + 0.5*sin(2*(1:nω)))
    cusp!(out, z) = (@. out = -im*(k0 + k1*sqrt(z)) - k2; out)
    cuspΦ(z) = @. -im*(k0*z + k1*2/3*z^1.5) - k2*z

    @testset "against a known integral, tol = $tol" for tol in (1e-4, 1e-6, 1e-8)
        tab = LinearOps.TabulatedLinop(cusp!, proto, 0.0, zend; tol, quiet=true)
        @test tab.z[1] == 0.0
        @test tab.z[end] == zend # the table covers zmax + max_dz
        @test issorted(tab.z)
        @test tab.err <= tol # the bisection reached the tolerance it reports
        @test size(tab.Φ) == (nω, length(tab.z))
        @test size(tab.dΦ) == size(tab.Φ)
        out = similar(proto)
        err = 0.0
        for z in (0.0, 1e-9, 1e-5, 1e-3, L/7, L/2, 0.83L, L, 1.04L)
            LinearOps.integrated!(out, tab, z)
            err = max(err, maximum(abs, out .- cuspΦ(z)))
        end
        #= The criterion is on one interval at a time and the differences accumulate over
           the intervals up to z, so the total is allowed a few of them. =#
        @test err < 20tol
    end

    @testset "tighter tolerance, more nodes, smaller error" begin
        coarse = LinearOps.TabulatedLinop(cusp!, proto, 0.0, zend; tol=1e-4, quiet=true)
        fine = LinearOps.TabulatedLinop(cusp!, proto, 0.0, zend; tol=1e-8, quiet=true)
        @test length(fine.z) > 3*length(coarse.z)
        out = similar(proto)
        errs = map((coarse, fine)) do tab
            maximum((1e-5, L/3, 0.71L)) do z
                LinearOps.integrated!(out, tab, z)
                maximum(abs, out .- cuspΦ(z))
            end
        end
        @test errs[2] < 1e-3*errs[1]
    end

    #= What is stored is the deviation from the secant, which the propagator adds back.
       Both ends of the table are on the secant by construction, and the deviation is what
       has to fit in Float32 on a device. =#
    @testset "secant subtraction" begin
        tab = LinearOps.TabulatedLinop(cusp!, proto, 0.0, L; tol=1e-8, quiet=true)
        out = similar(proto)
        LinearOps.phase!(out, tab, 0.0)
        @test all(iszero, out)
        LinearOps.phase!(out, tab, L)
        @test maximum(abs, out) < 1e-8
        LinearOps.integrated!(out, tab, L)
        @test out ≈ tab.secant .* L
        full = similar(proto)
        LinearOps.integrated!(full, tab, L/2)
        @test tab.scale < 0.1*maximum(abs, full)
        @test tab.scale ≈ maximum(abs, tab.Φ)
    end

    #= The propagator: `exp(Φ(t2) - Φ(t1))` against the closed-form integral, the
       backward propagation as its inverse, and the caches not changing the answer when
       the arguments repeat (six stages share one `t1`, and `t2` repeats between the
       forward and backward propagation of each stage). =#
    @testset "propagator" begin
        tab = LinearOps.TabulatedLinop(cusp!, proto, 0.0, L; tol=1e-10, quiet=true)
        prop! = RK45.make_prop!(tab, proto)
        y0 = ones(ComplexF64, nω)
        t1, t2 = 0.021L, 0.023L
        y = copy(y0)
        prop!(y, t1, t2)
        ref = exp.(cuspΦ(t2) .- cuspΦ(t1))
        @test maximum(abs, y .- ref)/maximum(abs, ref) < 1e-8
        prop!(y, t1, t2, true)
        @test maximum(abs, y .- y0) < 1e-12
        #= The cached readback is the same number as the fresh one: `prop!` above has
           already been called at these arguments and now takes the cache, while a
           propagator built here has not. =#
        y2 = copy(y0)
        prop!(y2, t1, t2)
        y3 = copy(y0)
        RK45.make_prop!(tab, proto)(y3, t1, t2)
        @test y2 == y3
        # ... including after an intervening step at other arguments evicted it
        prop!(y3, 0.5L, 0.6L)
        y4 = copy(y0)
        prop!(y4, t1, t2)
        @test y4 == y2
    end

    #= Luna's own z-dependent operators. No closed form, so the table is checked against
       one built to a much tighter tolerance, and the interesting property -- where the
       nodes went -- directly. =#
    function gradient_linop(p0, p1)
        coren, _ = Capillary.gradient(:Ar, L, p0, p1)
        m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
        (LinearOps.make_linop(grid, m, 800e-9)..., m)
    end
    function taper_linop()
        afun = z -> 75e-6 + (50e-6 - 75e-6)*z/L
        m = Capillary.MarcatiliMode(afun, :Ar, 1.0, loss=false, model=:full)
        (LinearOps.make_linop(grid, m, 800e-9)..., m)
    end

    @testset "$name" for (name, linop!) in (
            ("gradient 1 -> 0 bar", gradient_linop(1.0, 0.0)[1]),
            ("gradient 0 -> 1 bar", gradient_linop(0.0, 1.0)[1]),
            ("taper", taper_linop()[1]))
        tol = 1e-6
        tab = LinearOps.TabulatedLinop(linop!, proto, 0.0, zend; tol, quiet=true)
        ref = LinearOps.TabulatedLinop(linop!, proto, 0.0, zend; tol=1e-10, quiet=true)
        @test tab.z[end] == zend
        @test tab.err <= tol
        @test length(ref.z) > length(tab.z)
        out = similar(proto)
        rout = similar(proto)
        err = 0.0
        for z in (1e-7, 1e-4, 1e-3, L/7, L/2, 0.83L, L, 1.04L)
            LinearOps.integrated!(out, tab, z)
            LinearOps.integrated!(rout, ref, z)
            err = max(err, maximum(abs, out .- rout))
        end
        @test err < 20tol
    end

    #= The nodes go where the operator is hard and not where it is not: a p₀ = 0 entrance
       has a cusp at z = 0 and gets a bisection's worth of nodes there; a linear taper is
       smooth and gets none. =#
    @testset "node placement" begin
        tol = 1e-6
        cusp = LinearOps.TabulatedLinop(gradient_linop(0.0, 1.0)[1], proto, 0.0, zend;
                                        tol, quiet=true)
        smooth = LinearOps.TabulatedLinop(taper_linop()[1], proto, 0.0, zend;
                                          tol, quiet=true)
        @test cusp.z[2] - cusp.z[1] < 1e-3*zend
        @test smooth.z[2] - smooth.z[1] > 1e-2*zend
        @test count(<(0.01zend), cusp.z) > 5
        @test count(<(0.01zend), smooth.z) <= 2
    end

    #= A z-independent operator written as a closure: two nodes, nothing stored beyond the
       secant, and a propagator which reproduces the constant one. =#
    @testset "z-independent closure" begin
        m = Capillary.MarcatiliMode(75e-6, :Ar, 1.0, loss=false)
        linop!, _ = LinearOps.make_linop(grid, m, 800e-9)
        tab = LinearOps.TabulatedLinop(linop!, proto, 0.0, L; tol=1e-6, quiet=true)
        @test length(tab.z) == 2
        @test tab.scale < 1e-6
        const_linop = similar(proto)
        linop!(const_linop, 0.0)
        @test maximum(abs, Array(tab.secant) .- const_linop) < 1e-6
        y = ones(ComplexF64, nω)
        yref = copy(y)
        RK45.make_prop!(tab, proto)(y, 0.01, 0.013)
        RK45.make_prop!(const_linop, proto)(yref, 0.01, 0.013)
        @test maximum(abs, y .- yref) < 1e-10
    end

    #= The shape of the operator is not special-cased anywhere -- `_stack`, `selectdim`
       and `phase!` work on any number of axes -- so a multimode operator tabulates the
       same way. Four `MarcatiliMode`s on the same gradient, checked against a midpoint
       quadrature fine enough that its own error is well below the tolerance. =#
    @testset "a multimode operator" begin
        coren, _ = Capillary.gradient(:Ar, L, 0.0, 1.0)
        ms = [Capillary.MarcatiliMode(75e-6, coren, n=1, m=m, kind=:HE, loss=false)
              for m = 1:4]
        linop! = LinearOps.make_linop(grid, ms, 800e-9)
        mproto = zeros(ComplexF64, nω, length(ms))
        tab = LinearOps.TabulatedLinop(linop!, mproto, 0.0, zend; tol=1e-6, quiet=true)
        @test size(tab.Φ) == (nω, length(ms), length(tab.z))
        @test size(tab.secant) == (nω, length(ms))
        @test tab.err <= 1e-6
        out = similar(mproto)
        buf = similar(mproto)
        npanel = 20000
        for z in (L/3, L)
            LinearOps.integrated!(out, tab, z)
            #= Midpoint over `npanel` panels. The integrand has a √z cusp at 0, whose
               contribution to the midpoint error is O(h^{3/2}); at this panel count that
               is below the tolerance being checked. =#
            ref = zeros(ComplexF64, size(mproto))
            h = z/npanel
            for i = 1:npanel
                linop!(buf, (i - 0.5)*h)
                ref .+= buf .* h
            end
            @test maximum(abs, out .- ref) < 1e-4
        end
        y = ones(ComplexF64, size(mproto))
        RK45.make_prop!(tab, mproto)(y, 0.3L, 0.31L)
        @test all(isfinite, y)
    end

    #= Read outside the table, the value at the nearest end is used -- and the operator
       readback says so, because the propagator adds the secant term whatever it returns,
       so a step outside the table would propagate with the mean operator over the whole
       of it. `Luna.run` builds the table over every z the stepper can reach, so this
       cannot happen from there. =#
    @testset "reading outside the table" begin
        tab = LinearOps.TabulatedLinop(cusp!, proto, 0.0, L; tol=1e-6, quiet=true)
        out = similar(proto)
        ref = similar(proto)
        @test_logs (:warn, r"outside") LinearOps.phase!(out, tab, 1.5L)
        LinearOps.phase!(ref, tab, L)
        @test out == ref
        #= ... and a value table calls its source outside its span rather than holding its
           end value, so that a diagnostic recorded a fraction of a step past the end of
           the fibre is still right. =#
        atab = LinearOps.TabulatedScalar(z -> 1 + z^2, 0.0, L; tol=1e-6)
        @test atab(1.5L) == 1 + (1.5L)^2
        @test atab(-0.1L) == 1 + (-0.1L)^2
        @test atab(L) != atab(1.5L)
    end

    #= The value tables: β and Aeff are interpolated rather than integrated, to a
       tolerance relative to the largest value in the table. =#
    @testset "value tables" begin
        _, βfun!, m = gradient_linop(0.0, 1.0)
        βtab = LinearOps.TabulatedVector(βfun!, proto, nω, 0.0, zend; tol=1e-8)
        βref = zeros(Float64, nω)
        err = 0.0
        for z in (1e-7, 1e-3, L/3, 0.77L, L, 1.04L)
            βfun!(βref, z)
            err = max(err, maximum(abs, βtab(z) .- βref)/maximum(abs, βref))
        end
        @test err < 1e-7
        @test βtab(0.31L) === βtab.buf # the buffer is reused
        @test size(βtab.f) == (nω, length(βtab.z))

        _, _, mt = taper_linop()
        aeff = z -> Modes.Aeff(mt, z=z)
        atab = LinearOps.TabulatedScalar(aeff, 0.0, zend; tol=1e-6)
        for z in (0.0, 1e-4, L/3, L, 1.04L)
            @test isapprox(atab(z), aeff(z); rtol=1e-5)
        end
        #= Aeff of a fixed-radius mode does not depend on z, so this is the two-node table
           a uniform fibre gets, which is what makes tabulation one code path. =#
        flat = LinearOps.TabulatedScalar(z -> Modes.Aeff(m, z=z), 0.0, L; tol=1e-6)
        @test length(flat.z) == 2
        @test flat(0.5L) ≈ Modes.Aeff(m)
        #= The callable the table was built from is kept, so a table can be rebuilt over a
           wider span without going through its own interpolant. =#
        @test atab.src === aeff
        wider = LinearOps.TabulatedScalar(atab.src, 0.0, 1.5zend; tol=1e-6)
        @test wider.z[end] == 1.5zend
        @test isapprox(wider(L/3), aeff(L/3); rtol=1e-5)
    end
end

#= A caller-supplied analytic Φ, the third source of an integrated operator and the reason
   the interface is a type rather than a keyword: `-im*(k0 + k1*√z) - k2` is the shape a
   pressure gradient's operator has (a capillary filled by `Capillary.gradient` from
   p₀ = 0 has a density ∝ √z, and the propagation constant of a dilute gas is linear in
   density), and its integral is elementary. `phase!` returns Φ itself -- no secant is
   subtracted, so `secant` keeps its default of `nothing` -- and `derivative!` returns the
   operator, which is what the diagnostics need. Every scalar goes through `Luna.scalar`
   so that the broadcasts hold nothing but the state's own element type, which is what
   makes the same type run on a device (`test_device.jl`). =#
struct AnalyticLinop{T} <: LinearOps.AbstractIntegratedLinop
    k0::T
    k1::T
    k2::T
end

function LinearOps.phase!(out, op::AnalyticLinop, z)
    zs = Luna.scalar(out, z)
    zh = Luna.scalar(out, z^1.5*2/3)
    k0, k1, k2 = op.k0, op.k1, op.k2
    @. out = -im*(k0*zs + k1*zh) - k2*zs
    out
end

function LinearOps.derivative!(out, op::AnalyticLinop, z)
    sq = Luna.scalar(out, sqrt(z))
    k0, k1, k2 = op.k0, op.k1, op.k2
    @. out = -im*(k0 + k1*sq) - k2
    out
end

#= gpu/27: the interface the stepper actually propagates with. A z-dependent operator is
   given to `RK45.make_prop!` as its integral `Φ(z) = ∫linop dz'`, and there are three
   ways to get one: the table above, adaptive quadrature over each step
   (`QuadratureLinop`), and a closed form written by the caller. They are checked against
   each other and against the closed form here. =#
@testset "integrated linear operators" begin
    L = 0.1
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    nω = length(grid.ω)
    proto = zeros(ComplexF64, nω)

    #= The same operator the tabulation testset uses: a √z term on top of a constant, of
       the size a capillary's operator is. It is not an arbitrary choice -- a capillary
       filled by `Capillary.gradient` from p₀ = 0 has a density ∝ √z, and the propagation
       constant of a dilute gas is linear in density, so `a + b√z` is the shape of a
       pressure gradient's operator. Its integral is known in closed form, which is what
       makes it the reference for all three sources. =#
    k0 = @. 1e7*(1 + 0.1*sin(1:nω))
    k1 = @. 1e4*(1 + 0.5*cos(1:nω))
    k2 = @. 1.0*(1 + 0.5*sin(2*(1:nω)))
    cusp!(out, z) = (@. out = -im*(k0 + k1*sqrt(z)) - k2; out)
    cuspΦ(z) = @. -im*(k0*z + k1*2/3*z^1.5) - k2*z

    @testset "the three sources agree with the closed form" begin
        tab = LinearOps.TabulatedLinop(cusp!, proto, 0.0, L; tol=1e-8, quiet=true)
        quad = LinearOps.QuadratureLinop(cusp!, proto; tol=1e-8)
        ana = AnalyticLinop(k0, k1, k2)
        @test LinearOps.PhaseStyle(tab) === LinearOps.AbsolutePhase()
        @test LinearOps.PhaseStyle(quad) === LinearOps.IncrementalPhase()
        @test LinearOps.PhaseStyle(ana) === LinearOps.AbsolutePhase()
        @test LinearOps.secant(ana) === nothing
        out = similar(proto)
        for (z1, z2) in ((0.0, 1e-4), (0.0, L/20), (0.3L, 0.35L), (0.5L, L))
            ref = cuspΦ(z2) .- cuspΦ(z1)
            for op in (tab, quad, ana)
                LinearOps.phasediff!(out, op, z1, z2)
                @test maximum(abs, out .- ref) < 2e-7
            end
        end
        # the operator itself, which is what `derivative!` is for
        for z in (1e-6, L/3, 0.9L)
            cusp!(proto, z)
            for (op, atol) in ((tab, 1e-4), (quad, 1e-12), (ana, 1e-12))
                LinearOps.derivative!(out, op, z)
                @test maximum(abs, out .- proto)/maximum(abs, proto) < atol
            end
        end
    end

    #= The propagator is one method for all three, and `exp(Φ(t2) − Φ(t1))` is exact for
       the linear part: propagating forwards and back returns the state. =#
    @testset "the propagator" begin
        ops = (("tabulated", LinearOps.TabulatedLinop(cusp!, proto, 0.0, L;
                                                      tol=1e-10, quiet=true)),
               ("quadrature", LinearOps.QuadratureLinop(cusp!, proto; tol=1e-10)),
               ("analytic", AnalyticLinop(k0, k1, k2)))
        y0 = ones(ComplexF64, nω)
        t1, t2 = 0.31L, 0.33L
        ref = exp.(cuspΦ(t2) .- cuspΦ(t1))
        for (name, op) in ops
            y = copy(y0)
            prop! = RK45.make_prop!(op, y0)
            prop!(y, t1, t2)
            @test maximum(abs, y .- ref)/maximum(abs, ref) < 1e-8
            prop!(y, t1, t2, true)
            @test maximum(abs, y .- y0) < 1e-12
            # the caches do not change the answer when the arguments repeat
            y2 = copy(y0)
            prop!(y2, t1, t2)
            y3 = copy(y0)
            RK45.make_prop!(op, y0)(y3, t1, t2)
            @test y2 == y3
        end
    end

    #= Luna's own operators, where there is no closed form: the table and the quadrature
       have to agree with each other to their tolerances. =#
    function gradient_linop(p0, p1)
        coren, _ = Capillary.gradient(:Ar, L, p0, p1)
        LinearOps.make_linop(grid, Capillary.MarcatiliMode(75e-6, coren, loss=false),
                             800e-9)[1]
    end
    taper_linop() = LinearOps.make_linop(
        grid, Capillary.MarcatiliMode(z -> 75e-6 + (50e-6 - 75e-6)*z/L, :Ar, 1.0,
                                      loss=false, model=:full), 800e-9)[1]

    @testset "table and quadrature agree: $name" for (name, linop!) in (
            ("gradient 0 -> 1 bar", gradient_linop(0.0, 1.0)),
            ("taper", taper_linop()))
        tol = 1e-8
        tab = LinearOps.TabulatedLinop(linop!, proto, 0.0, L; tol, quiet=true)
        quad = LinearOps.QuadratureLinop(linop!, proto; tol)
        a = similar(proto)
        b = similar(proto)
        for (z1, z2) in ((0.0, L/20), (0.2L, 0.21L), (0.5L, L))
            LinearOps.phasediff!(a, tab, z1, z2)
            LinearOps.phasediff!(b, quad, z1, z2)
            # a few interval tolerances, as the interpolation error accumulates over them
            @test maximum(abs, a .- b) < 20tol
        end
        # and so do the operators they report
        linop!(proto, 0.43L)
        LinearOps.derivative!(a, tab, 0.43L)
        LinearOps.derivative!(b, quad, 0.43L)
        @test maximum(abs, b .- proto) == 0 # the quadrature just calls the closure
        @test maximum(abs, a .- proto)/maximum(abs, proto) < 1e-6
    end

    #= The quadrature over a step is one 15-point rule for a smooth operator. Over the
       whole of a gradient filled from zero it is not, because the √z cusp is at the
       entrance -- which is what `maxevals` bounds and what the table's bisection handles
       instead. =#
    @testset "quadrature cost" begin
        linop! = gradient_linop(1.0, 0.0) # smooth: no cusp at either end
        quad = LinearOps.QuadratureLinop(linop!, proto; tol=1e-6)
        out = similar(proto)
        LinearOps.phasediff!(out, quad, 0.5L, 0.5L + L/20)
        @test quad.nevals[] == 15
        @test quad.ncalls[] == 1
        # a zero-length step costs nothing at all
        LinearOps.phasediff!(out, quad, 0.3L, 0.3L)
        @test all(iszero, out)
        @test quad.nevals[] == 15
    end

    #= A constant added to an integrated operator -- which is how an absorbing boundary
       reaches one -- goes into the integral exactly, because it is linear in z. =#
    @testset "OffsetLinop" begin
        δ = @. -0.5*(1 + 0.2*cos(1:nω))
        tab = LinearOps.TabulatedLinop(cusp!, proto, 0.0, L; tol=1e-10, quiet=true)
        quad = LinearOps.QuadratureLinop(cusp!, proto; tol=1e-10)
        out = similar(proto)
        ref = similar(proto)
        for op in (tab, quad, AnalyticLinop(k0, k1, k2))
            off = LinearOps.OffsetLinop(op, δ)
            @test LinearOps.PhaseStyle(off) === LinearOps.PhaseStyle(op)
            LinearOps.phasediff!(out, off, 0.2L, 0.25L)
            LinearOps.phasediff!(ref, op, 0.2L, 0.25L)
            @test maximum(abs, out .- ref .- δ.*(0.05L)) < 1e-9
            LinearOps.derivative!(out, off, 0.2L)
            LinearOps.derivative!(ref, op, 0.2L)
            @test out ≈ ref .+ δ
        end
        # and the same through `Boundaries.addloss`, which is where it comes from
        α = @. 1.0*(1 + 0.5*sin(1:nω))
        off = Boundaries.addloss(AnalyticLinop(k0, k1, k2), α)
        @test off isa LinearOps.OffsetLinop
        LinearOps.derivative!(out, off, 0.4L)
        cusp!(ref, 0.4L)
        @test out ≈ ref .- α./2
        # the free-space clamp cannot be applied to an operator which is already integrated
        @test_throws ErrorException Boundaries.clampdecay(
            AnalyticLinop(k0, k1, k2), 100.0)
    end

    #= A bare `linop!(out, z)` callable is not a propagator any more: the one-point rule
       `exp(linop(t2)·(t2 − t1))` has been removed. =#
    @testset "a callable is refused" begin
        @test_throws ArgumentError RK45.make_prop!(cusp!, proto)
    end

    @testset "linop_integral" begin
        @test Luna._linop_integral(:tabulated, nothing) === :tabulated
        @test Luna._linop_integral(:quadrature, nothing) === :quadrature
        @test_throws ErrorException Luna._linop_integral(:hermite, nothing)
        # the deprecated keyword warns and does not override an explicit request
        @test_logs (:warn, r"deprecated") Luna._linop_integral(:tabulated, true)
        @test Luna._linop_integral(:quadrature, false) === :quadrature
    end
end
