#= Checks for `Utils.threaded`, included by test_utils.jl: in process when Julia has
   several threads, otherwise from a child process started with two. =#
import Luna
import Luna: Utils
import Test: @test, @testset

@testset "threaded broadcasts, $(Threads.nthreads()) threads" begin
    z = ComplexF64[cis(0.1i) * (1 + i/7) for i in 1:5000]
    a = z .* 0.3; y = copy(z)
    ref = @. y * exp(a - z/2)
    @. $(Utils.threaded(y; minlen=16)) = y * exp(a - z/2)
    @test y == ref                                  # in place, dest is an argument
    E = rand(3000, 2, 5); d = similar(E)
    v = rand(3000)                                  # broadcasts along the other axes
    @. $(Utils.threaded(d; minlen=16)) = E^2 * v + 1
    @test d == @. E^2 * v + 1
    @. $(Utils.threaded(d; minlen=16)) = 0          # scalar right-hand side
    @test all(iszero, d)
    V = view(E, :, 1:1, :); W = similar(V)          # view destination
    @. $(Utils.threaded(W; minlen=16)) = sqrt(V)
    @test W == sqrt.(V)
    #= An argument aliasing the destination without being it is copied first, as Base
       does, so the shifted read sees the old values. =#
    x = collect(1.0:6000.0); xs = view(x, [2:6000; 1])
    refx = x .+ xs
    @. $(Utils.threaded(x; minlen=16)) = x + xs
    @test x == refx
    # nested and switched-off regions take the plain path and agree too
    y2 = copy(z)
    Utils.serial_region() do
        @. $(Utils.threaded(y2; minlen=16)) = y2 * exp(a - z/2)
    end
    @test y2 == ref
    was = Luna.settings["threaded_broadcasts"]
    Luna.set_threaded_broadcasts(false)
    y3 = copy(z)
    @. $(Utils.threaded(y3; minlen=16)) = y3 * exp(a - z/2)
    Luna.set_threaded_broadcasts(was)
    @test y3 == ref
    @test Utils._THREADED_DEPTH[] == 0
    # the automatic FFTW count gives large plans the Julia threads, small ones one
    fftw0 = Luna.settings["fftw_threads"]
    try
        Luna.set_fftw_threads(0)
        @test Utils.FFTWthreads(Utils.FFTW_THREAD_MINLEN) == Threads.nthreads()
        @test Utils.FFTWthreads(Utils.FFTW_THREAD_MINLEN - 1) == 1
        @test Utils.serial_fftw(() -> Utils.FFTWthreads(Utils.FFTW_THREAD_MINLEN)) == 1
    finally
        Luna.set_fftw_threads(fftw0)
    end
    #= The modified shot noise's advance (LupoLab/Luna.jl#440): the phase readback and the
       `cis` broadcast are threaded above their thresholds, and agree with the serial path
       bit for bit, inside the table, in the extrapolation past it, and for a constant
       operator. =#
    n = Utils.THREAD_MINLEN
    dφ = [1e3*sin(0.37i) for i in 1:n]
    E0 = ComplexF64[cis(0.2i)*(1 + i/n) for i in 1:n]
    Φ = hcat(zeros(n), 0.1 .* dφ, 0.3 .* dφ .+ 1e-3 .* cos.(1:n))
    dΦ = hcat(dφ, dφ, 1.1 .* dφ)
    for p in (Luna.LinearOps.ConstUnitaryPhase(dφ),
              Luna.LinearOps.TabulatedUnitaryPhase([0.0, 0.1, 0.3], Φ, dΦ, 1e-6, 0.0))
        mk() = Luna.NonlinearRHS.ModifiedNoise(copy(E0), p, zeros(n), similar(E0), Ref(NaN))
        for z in (0.17, 0.45)
            nt, ns = mk(), mk()
            a = copy(Luna.NonlinearRHS.noise_at!(nt, z))
            b = Utils.serial_region(() -> copy(Luna.NonlinearRHS.noise_at!(ns, z)))
            @test a == b
            @test a == E0 .* cis.(Utils.serial_region(() -> p(zeros(n), z)))
        end
    end
    @test Utils._THREADED_DEPTH[] == 0
end
