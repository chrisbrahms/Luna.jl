import Test: @test, @testset, @test_throws
import Luna
import Luna: Utils, Grid, PhysData, NonlinearRHS, Nonlinear, Fields
import LinearAlgebra
import HDF5
import Dates

@testset "Utils" begin

@testset "dict->HDF5" begin
d = Dict{String, Any}()
d["float"] = 1.0
d["float[]"] = [1.0, 2.0, 3.0]
d["string"] = "foo"
d["nothing"] = nothing
d["dict"] = Dict("foo"=>5, "bar"=>[1, 2, 3], "baz"=>Dict("complex"=>5.0+2.0im))
fpath = joinpath(Utils.cachedir(), "output_test", "test.h5")
isfile(fpath) && rm(fpath)
isdir(dirname(fpath)) || mkpath(dirname(fpath))
Utils.save_dict_h5(fpath, d)
@test_throws ErrorException Utils.save_dict_h5(fpath, d, force=false)
HDF5.h5open(fpath) do file
    for k in ["float", "float[]", "string"]
        @test d[k] == read(file[k])
    end
    @test d["dict"]["foo"] == read(file["dict"]["foo"])
    @test d["dict"]["bar"] == read(file["dict"]["bar"])
    @test d["dict"]["baz"]["complex"] == read(file["dict"]["baz"]["complex"])
end
rm(fpath)

delete!(d, "nothing") # nothing values are converted to empty arrays and not re-converted
Utils.save_dict_h5(fpath, d)
dd = Utils.load_dict_h5(fpath)
@test d == dd

rm(fpath)
rm(dirname(fpath), force=true)
end

@testset "super/subscripts" begin
@test Utils.subscript(0) == Utils.subscript("0") == Utils.subscript('0') == "₀"
@test Utils.subscript(1) == Utils.subscript("1") == Utils.subscript('1') == "₁"
@test Utils.subscript(2) == Utils.subscript("2") == Utils.subscript('2') == "₂"

@test Utils.subscript(123456789) == Utils.subscript("123456789") == "₁₂₃₄₅₆₇₈₉"
@test Utils.subscript("0123456789") == "₀₁₂₃₄₅₆₇₈₉"
end


@testset "date formatting" begin
    start = Dates.now()
    sleep(2)
    finish = Dates.now()
    st = Utils.format_elapsed(finish-start)
    @test ~isnothing(match(r"2.[0-9]{3} seconds", st))

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2022, 01, 01, 0, 1, 0)
    @test Utils.format_elapsed(finish-start) == "1 minute, 0.000 seconds"

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2022, 01, 01, 1, 0, 0)
    @test Utils.format_elapsed(finish-start) == "1 hour, 0 minutes, 0.000 seconds"

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2022, 01, 01, 2, 0, 0)
    @test Utils.format_elapsed(finish-start) == "2 hours, 0 minutes, 0.000 seconds"

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2022, 01, 01, 1, 15, 32)
    @test Utils.format_elapsed(finish-start) == "1 hour, 15 minutes, 32.000 seconds"

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2022, 01, 01, 1, 15, 32, 123)
    @test Utils.format_elapsed(finish-start) == "1 hour, 15 minutes, 32.123 seconds"

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2022, 02, 01, 0, 15, 32, 123)
    @test Utils.format_elapsed(finish-start) == "744 hours, 15 minutes, 32.123 seconds"

    start = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    finish = Dates.DateTime(2023, 01, 01, 0, 0, 0)
    @test Utils.format_elapsed(finish-start) == "8760 hours, 0 minutes, 0.000 seconds"

    finish = Dates.DateTime(2022, 01, 01, 0, 0, 0)
    start = Dates.DateTime(2022, 01, 01, 1, 15, 32)
    @test Utils.format_elapsed(finish-start) == "-1 hour, -15 minutes, -32.000 seconds"

end


#= The FFTW wisdom switch. The regression gate (test/test_regression.jl) depends on it, so
   it is worth a test even though nothing else in the test suite turns it off. =#
@testset "FFTW wisdom switch" begin
    was = Luna.settings["fftw_wisdom"]
    @test was == true # the default, and nothing before this point may have changed it
    fpath = joinpath(Utils.cachedir(), "FFTWcache_$(Utils.FFTWthreads())threads")
    lockpath = joinpath(Utils.cachedir(), "FFTWlock")
    try
        @test Luna.set_fftw_wisdom(false) == false
        @test Luna.settings["fftw_wisdom"] == false
        #= Both must become no-ops. `mtime` of a missing file is 0.0, which is a fine
           before/after comparison either way: what matters is that neither call touches
           the cache the other worktrees share. =#
        before = mtime(fpath)
        @test Utils.loadFFTwisdom() === nothing
        @test Utils.saveFFTwisdom() === nothing
        @test mtime(fpath) == before
        @test !isfile(lockpath) # no pidlock taken either
        # The thread count is still re-asserted (see the loadFFTwisdom docstring).
        @test Utils.FFTWthreads() > 0

        @test Luna.set_fftw_wisdom(true) == true
        @test Luna.settings["fftw_wisdom"] == true
    finally
        Luna.set_fftw_wisdom(was)
    end
    @test Luna.settings["fftw_wisdom"] == was
end

#= Thread counts: the automatic FFTW and BLAS rules and the explicit overrides. The rules
   are those measured in benchmark/threads (THREADS_REPORT.md). =#
@testset "thread counts" begin
    fftw0, blas0 = Luna.settings["fftw_threads"], Luna.settings["blas_threads"]
    try
        Luna.set_fftw_threads(0)
        J = Threads.nthreads()
        big, small = Utils.FFTW_THREAD_MINLEN, Utils.FFTW_THREAD_MINLEN - 1
        @test Utils.FFTWthreads(small) == 1
        @test Utils.FFTWthreads(big) == (J > 1 ? J : 1)
        @test Utils.serial_fftw(() -> Utils.FFTWthreads(big)) == 1
        @test_throws ErrorException Utils.serial_fftw(() -> error("boom"))
        @test Utils._FFTW_SERIAL[] == false # restored after the throw
        Luna.set_fftw_threads(3) # explicit: every plan, also with one Julia thread
        @test Utils.FFTWthreads(small) == 3
        @test Utils.serial_fftw(() -> Utils.FFTWthreads(big)) == 3
        Luna.set_fftw_threads(0)

        nb = LinearAlgebra.BLAS.get_num_threads()
        other = nb == 1 ? 2 : 1
        @test Utils.with_BLAS_threads(() -> LinearAlgebra.BLAS.get_num_threads(), other) == other
        @test LinearAlgebra.BLAS.get_num_threads() == nb
        @test_throws ErrorException Utils.with_BLAS_threads(() -> error("boom"), other)
        @test LinearAlgebra.BLAS.get_num_threads() == nb # restored after the throw
        @test Utils.with_BLAS_threads(() -> LinearAlgebra.BLAS.get_num_threads(), nothing) == nb

        # the automatic BLAS rule, for J Julia threads and B₀ BLAS's own default
        @test Luna._auto_blas(:radial, 4, 8) == 4
        @test Luna._auto_blas(:radial, 8, 8) == 4
        @test Luna._auto_blas(:modal, 4, 8) == 4
        @test Luna._auto_blas(:modal, 8, 8) == 1
        @test Luna._auto_blas(:modal, 10, 8) == 1
        @test isnothing(Luna._auto_blas(:none, 4, 8))
        E = zeros(ComplexF64, 4)
        Luna.set_blas_threads(3)
        @test Luna.blas_threads(nothing, E) == 3
        Luna.set_blas_threads(0)
        @test isnothing(Luna.blas_threads(nothing, E)) # no matrix products
        # the dispatch on a real radial transform
        grid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
        rg = Grid.RadialGrid(1e-3, 8)
        nfun = PhysData.ref_index_fun(:Ar, 1.0)
        _, transform, _ = Luna.setup(grid, rg, z -> PhysData.density(:Ar, 1.0),
                                     NonlinearRHS.const_norm_radial(grid, rg, nfun),
                                     (Nonlinear.Kerr_field(PhysData.γ3_gas(:Ar)),),
                                     Fields.GaussGaussField(λ0=800e-9, τfwhm=30e-15,
                                                            energy=1e-9, w0=200e-6))
        @test Luna._gemmkind(transform) === :radial
    finally
        Luna.set_fftw_threads(fftw0)
        Luna.set_blas_threads(blas0)
    end
end

#= The threaded path needs several Julia threads; with one, the same checks run in a child
   process started with two, so the suite exercises it either way. =#
if Threads.nthreads() > 1
    include(joinpath(@__DIR__, "threaded_checks.jl"))
else
    @testset "threaded broadcasts (child process, 2 threads)" begin
        cmd = `$(Base.julia_cmd()) --startup-file=no -t 2 --project=$(Base.active_project())
               $(joinpath(@__DIR__, "threaded_checks.jl"))`
        @test success(pipeline(cmd; stdout=devnull, stderr=stderr))
    end
end

end
