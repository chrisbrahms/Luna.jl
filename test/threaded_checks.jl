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
    Luna.set_fftw_threads(0)
    @test Utils.FFTWthreads(Utils.FFTW_THREAD_MINLEN) == Threads.nthreads()
    @test Utils.FFTWthreads(Utils.FFTW_THREAD_MINLEN - 1) == 1
    @test Utils.serial_fftw(() -> Utils.FFTWthreads(Utils.FFTW_THREAD_MINLEN)) == 1
    Luna.set_fftw_threads(fftw0)
end
