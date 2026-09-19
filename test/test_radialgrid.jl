#= Tests for Grid.RadialGrid, the Luna-owned transverse grid for radially symmetric
   propagation. This is the only file in test/ which is allowed to build a Hankel.QDHT:
   everywhere else Luna is used through Grid.RadialGrid. =#
import Test: @test, @testset, @test_throws, @test_logs, TestLogger
import Luna
import Luna: Grid, Maths, Fields, LinearOps, NonlinearRHS, Nonlinear, PhysData
import Luna: Hankel
import Logging
import LinearAlgebra: mul!, ldiv!

R = 12e-3
N = 64

@testset "RadialGrid" begin

@testset "construction" begin
    rg = Grid.RadialGrid(R, N)
    q = Hankel.QDHT(R, N)
    @test rg isa Grid.SpaceGrid
    @test rg.R == R
    @test rg.N == N
    @test rg.order == 0
    @test rg.K == q.K
    @test rg.r == q.r
    @test rg.k == q.k
    @test rg.wr == q.scaleR
    @test rg.wk == q.scaleK
    @test size(rg.Tfwd) == (N, N)
    @test size(rg.Tbwd) == (N, N)
    @test Grid.kperp2(rg) == q.k.^2
    Grid.validate(rg)

    # higher orders are constructed but only the transform is defined for them
    rg1 = Grid.RadialGrid(R, N; order=1)
    q1 = Hankel.QDHT{1}(R, N)
    @test rg1.order == 1
    @test rg1.r == q1.r
    @test rg1.K == q1.K
    @test_throws DomainError Grid.symmetric(rg1, ones(N))
    @test_throws DomainError Grid.onaxis(rg1, ones(N))
end

@testset "conversion from QDHT" begin
    q = Hankel.QDHT(R, N; dim=3)
    #= `Logging.@warn(..., maxlog=1)` is counted per logger, so a plain `@test_logs` would
       depend on whether anything earlier in the session already converted a QDHT.
       `respect_maxlog=false` makes the test independent of that. =#
    logger = TestLogger(respect_maxlog=false)
    rg = Logging.with_logger(logger) do
        Grid.RadialGrid(q)
    end
    @test any(r -> r.level == Logging.Warn && occursin("deprecated", r.message),
              logger.logs)
    @test rg.R == q.R
    @test rg.K == q.K
    @test rg.N == q.N
    @test rg.r == q.r
    @test rg.k == q.k
    # dim is ignored: a RadialGrid always transforms along the last dimension
    q2 = Hankel.QDHT(R, N; dim=1)
    rg2 = Grid.RadialGrid(q2)
    @test rg2.Tfwd == rg.Tfwd
    @test rg2.Tbwd == rg.Tbwd
end

#= The right GEMM on a (:, N) reshape must reproduce Hankel's left multiplication along
   the radial axis, for arrays of any rank. =#
@testset "transform agreement with Hankel: $(length(sz))-D" for sz in ((N,), (16, N), (16, 2, N))
    rg = Grid.RadialGrid(R, N)
    q = Hankel.QDHT(R, N; dim=length(sz))
    for T in (Float64, ComplexF64)
        A = T <: Complex ? randn(ComplexF64, sz) : randn(sz)
        Ak_h = similar(A)
        mul!(Ak_h, q, A)
        Ak = Grid.to_kspace(rg, A)
        @test Ak ≈ Ak_h
        @test maximum(abs, Ak - Ak_h) < 1e-12*maximum(abs, Ak_h)

        Ar_h = similar(A)
        ldiv!(Ar_h, q, A)
        Ar = Grid.to_rspace(rg, A)
        @test Ar ≈ Ar_h
        @test maximum(abs, Ar - Ar_h) < 1e-12*maximum(abs, Ar_h)

        # round trip
        @test Grid.to_rspace(rg, Grid.to_kspace(rg, A)) ≈ A

        # in-place forms, including aliased input and output
        out = similar(A)
        Grid.to_kspace!(out, rg, A)
        @test out == Ak
        Ac = copy(A)
        Grid.to_kspace!(Ac, rg, Ac)
        @test Ac == Ak
        Grid.to_rspace!(out, rg, A)
        @test out == Ar
        Ac = copy(A)
        Grid.to_rspace!(Ac, rg, Ac)
        @test Ac == Ar
    end
end

@testset "transform along a leading dimension" begin
    rg = Grid.RadialGrid(R, N)
    q = Hankel.QDHT(R, N; dim=2)
    A = randn(8, N, 3)
    @test Grid.to_kspace(rg, A; dim=2) ≈ Hankel.mul!(similar(A), q, A)
    @test Grid.to_rspace(rg, A; dim=2) ≈ Hankel.ldiv!(similar(A), q, A)
    # rank 4, every dimension: the result keeps the axes of the input, in order
    for dim in 1:4
        sz = [3, 2, 5, 4]
        sz[dim] = N
        A4 = randn(sz...)
        q4 = Hankel.QDHT(R, N; dim=dim)
        Ak = Grid.to_kspace(rg, A4; dim=dim)
        @test size(Ak) == size(A4)
        @test Ak ≈ Hankel.mul!(similar(A4), q4, A4)
        @test Grid.to_rspace(rg, A4; dim=dim) ≈ Hankel.ldiv!(similar(A4), q4, A4)
    end
end

@testset "shape and dimension errors" begin
    rg = Grid.RadialGrid(R, N)
    @test_throws DimensionMismatch Grid.to_kspace(rg, randn(N+1))
    @test_throws DimensionMismatch Grid.to_kspace!(randn(N+1), rg, randn(N))
    @test_throws DimensionMismatch Grid.integrate_r(rg, randn(N); dim=2)
    # a dim outside the input's rank is an error, not a BoundsError from the permutation
    @test_throws DimensionMismatch Grid.to_kspace(rg, randn(N); dim=2)
    @test_throws DimensionMismatch Grid.to_rspace(rg, randn(4, N); dim=3)
    @test_throws DimensionMismatch Grid.integrate_k(rg, randn(4, N); dim=0)
end

@testset "integrals" begin
    rg = Grid.RadialGrid(R, N)
    # 1-D: matches Hankel and the analytical result
    w0 = 1e-3
    A = Maths.gauss.(rg.r, w0/2)
    q = Hankel.QDHT(R, N)
    @test Grid.integrate_r(rg, abs2.(A)) ≈ Hankel.integrateR(abs2.(A), q)
    @test Grid.integrate_r(rg, abs2.(A)) isa Float64
    @test Grid.integrate_r(rg, abs2.(A)) ≈ w0^2/8 rtol=1e-9 # ∫exp(-2r²/(w0/2)²·2)r dr
    Ak = Grid.to_kspace(rg, A)
    @test Grid.integrate_k(rg, abs2.(Ak)) ≈ Hankel.integrateK(abs2.(Ak), q) # Parseval
    @test Grid.integrate_k(rg, abs2.(Ak)) ≈ Grid.integrate_r(rg, abs2.(A))

    # 2-D (t, r), as Fields.energyfuncs sees it before the polarisation axis exists
    A2 = randn(16, N)
    q2 = Hankel.QDHT(R, N; dim=2)
    I2 = Grid.integrate_r(rg, A2)
    @test size(I2) == (16,)
    @test I2 ≈ dropdims(Hankel.integrateR(A2, q2; dim=2); dims=2)
    # the same numbers reached by integrating each row on its own
    @test I2 ≈ [Grid.integrate_r(rg, A2[ii, :]) for ii in axes(A2, 1)]

    # 3-D (t, pol, r)
    A3 = randn(16, 2, N)
    q3 = Hankel.QDHT(R, N; dim=3)
    I3 = Grid.integrate_k(rg, A3)
    @test size(I3) == (16, 2)
    @test I3 ≈ dropdims(Hankel.integrateK(A3, q3; dim=3); dims=3)

    #= Integrating along a leading dimension: the result must have the axes of the input
       with `dim` dropped, in order, not merely the right numbers in some order. =#
    A4 = randn(8, N, 3)
    q4 = Hankel.QDHT(R, N; dim=2)
    @test Grid.integrate_r(rg, A4; dim=2) ≈ dropdims(Hankel.integrateR(A4, q4; dim=2); dims=2)

    # rank 4, with two axes after the integrated one
    A5 = randn(N, 3, 2, 5)
    q5 = Hankel.QDHT(R, N; dim=1)
    I5 = Grid.integrate_r(rg, A5; dim=1)
    @test size(I5) == (3, 2, 5)
    @test I5 ≈ dropdims(Hankel.integrateR(A5, q5; dim=1); dims=1)
    A6 = randn(3, N, 2, 5)
    q6 = Hankel.QDHT(R, N; dim=2)
    I6 = Grid.integrate_k(rg, A6; dim=2)
    @test size(I6) == (3, 2, 5)
    @test I6 ≈ dropdims(Hankel.integrateK(A6, q6; dim=2); dims=2)
    O6 = Grid.onaxis(rg, A6; dim=2)
    @test size(O6) == (3, 2, 5)
    @test O6 ≈ dropdims(Hankel.onaxis(A6, q6; dim=2); dims=2)
end

@testset "onaxis and symmetric" begin
    rg = Grid.RadialGrid(R, N)
    q = Hankel.QDHT(R, N)
    w0 = 1e-3
    A = Maths.gauss.(rg.r, w0/2)
    @test Grid.onaxis(rg, Grid.to_kspace(rg, A)) ≈ 1 rtol=1e-9
    @test Grid.onaxis(rg, Grid.to_kspace(rg, A)) ≈ Hankel.onaxis(q*A, q)

    As = Grid.symmetric(rg, A)
    @test As ≈ Hankel.symmetric(A, q)
    rs = Grid.rsymmetric(rg)
    @test rs == Hankel.Rsymmetric(q)
    @test length(rs) == length(As) == 2N+1
    @test As[1:N] == A[N:-1:1]
    @test As[N+1] ≈ 1 rtol=1e-9
    @test As[N+2:end] == A

    # 2-D, symmetric along the last dimension
    A2 = Maths.gauss.(randn(4), 1.0) .* Maths.gauss.(rg.r, w0/2)'
    q2 = Hankel.QDHT(R, N; dim=2)
    @test Grid.symmetric(rg, A2) ≈ Hankel.symmetric(A2, q2; dim=2)
end

#= The deprecated entry points: a Hankel.QDHT where a RadialGrid is expected. `Luna.setup`
   is the one every pre-existing radial script calls, with six positional arguments, which
   is also the arity of the mode-averaged `setup`; nothing else in the test suite or the
   examples exercises it any more. =#
@testset "deprecated Hankel.QDHT entry points" begin
    Rq = 100e-6
    Nq = 16
    λ0 = 800e-9
    q = Hankel.QDHT(Rq, Nq; dim=3)
    rg = Grid.RadialGrid(Rq, Nq)
    nfunλ = PhysData.ref_index_fun(:Ar, 1)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    nfunω = (ω; z) -> nfun(PhysData.wlfreq(ω); z)
    dens = z -> PhysData.density(:Ar, 1)
    inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy=1e-9, w0=20e-6)
    grids = (Grid.RealGrid(λ0, (400e-9, 2000e-9), 0.2e-12),
             Grid.EnvGrid(λ0, (400e-9, 2000e-9), 0.2e-12))
    for grid in grids
        resp = grid isa Grid.RealGrid ?
            (Nonlinear.Kerr_field(PhysData.γ3_gas(:Ar)),) :
            (Nonlinear.Kerr_env(PhysData.γ3_gas(:Ar)),)
        # the legacy call: six positional arguments with a QDHT in second place
        Eωq, tq, _ = Luna.setup(grid, q, dens,
                                NonlinearRHS.const_norm_radial(grid, q, nfun), resp, inputs)
        Eωr, tr, _ = Luna.setup(grid, rg, dens,
                                NonlinearRHS.const_norm_radial(grid, rg, nfun), resp, inputs)
        @test Eωq == Eωr
        @test tq.rgrid.Tfwd == tr.rgrid.Tfwd
        @test tq.rgrid.N == Nq

        lq = LinearOps.make_const_linop(grid, q, nfun, true)
        lr = LinearOps.make_const_linop(grid, rg, nfun, true)
        @test lq == lr
        oq, or = similar(lr), similar(lr)
        LinearOps.make_linop(grid, q, nfunω, true)(oq, 0.0)
        LinearOps.make_linop(grid, rg, nfunω, true)(or, 0.0)
        @test oq == or

        Et = randn(length(grid.t), Nq)
        @test Fields.energyfuncs(grid, q)[1](Et) == Fields.energyfuncs(grid, rg)[1](Et)
    end
end

@testset "to_dict/from_dict" begin
    for order in (0, 1)
        rg = Grid.RadialGrid(R, N; order=order)
        d = Grid.to_dict(rg)
        @test sort(collect(keys(d))) == ["N", "R", "order"]
        rg2 = Grid.from_dict(Grid.RadialGrid, d)
        @test rg2.R == rg.R
        @test rg2.N == rg.N
        @test rg2.order == rg.order
        @test rg2.r == rg.r
        @test rg2.k == rg.k
        @test rg2.Tfwd == rg.Tfwd
        @test rg2.Tbwd == rg.Tbwd
        @test rg2.wr == rg.wr
        @test rg2.wk == rg.wk
        @test Grid.RadialGrid(d).R == rg.R
    end
end

end
