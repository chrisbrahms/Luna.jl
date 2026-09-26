import Test: @test, @testset, @test_throws
import Luna
import Luna: Utils
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
