#= Tests for Grid.RadialGrid, the Luna-owned transverse grid for radially symmetric
   propagation. This is the only file in test/ which is allowed to build a Hankel.QDHT:
   everywhere else Luna is used through Grid.RadialGrid. =#
import Test: @test, @testset, @test_throws, @test_logs
import Luna: Grid, Maths
import Luna: Hankel
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
    rg = @test_logs (:warn,) match_mode=:any Grid.RadialGrid(q)
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
end

@testset "shape and dimension errors" begin
    rg = Grid.RadialGrid(R, N)
    @test_throws DimensionMismatch Grid.to_kspace(rg, randn(N+1))
    @test_throws DimensionMismatch Grid.to_kspace!(randn(N+1), rg, randn(N))
    @test_throws DimensionMismatch Grid.integrate_r(rg, randn(N); dim=2)
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

    # integrating along a leading dimension
    A4 = randn(8, N, 3)
    q4 = Hankel.QDHT(R, N; dim=2)
    @test Grid.integrate_r(rg, A4; dim=2) ≈ dropdims(Hankel.integrateR(A4, q4; dim=2); dims=2)
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
