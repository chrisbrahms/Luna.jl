import Test: @test, @test_throws, @test_logs, @testset, TestLogger
import Logging
import FFTW
import Luna: Grid, Output, Luna

@testset "grid constructors" begin
λ0 = 800e-9
λlims = (160e-9, 3000e-9)
trange = 1e-12

rgrid = Logging.with_logger(Logging.NullLogger()) do
    Grid.RealGrid(λ0, λlims, trange)
end
egrid = Logging.with_logger(Logging.NullLogger()) do
    Grid.EnvGrid(λ0, λlims, trange)
end

# the propagation length is not part of the grid any more
@test !(:zmax in fieldnames(Grid.RealGrid))
@test !(:zmax in fieldnames(Grid.EnvGrid))
@test !hasproperty(rgrid, :zmax)
@test !hasproperty(egrid, :zmax)

# λ_lims as a vector works as well as a tuple
rgridv = Logging.with_logger(Logging.NullLogger()) do
    Grid.RealGrid(λ0, [λlims...], trange)
end
@test rgridv.ω == rgrid.ω

# δt is a keyword argument now; it sets the spacing of the oversampled grid
rgridδ = Logging.with_logger(Logging.NullLogger()) do
    Grid.RealGrid(λ0, λlims, trange; δt=1e-17)
end
@test rgridδ.to[2] - rgridδ.to[1] < rgrid.to[2] - rgrid.to[1]

egridthg = Logging.with_logger(Logging.NullLogger()) do
    Grid.EnvGrid(λ0, λlims, trange; thg=true)
end
@test maximum(egridthg.ω) > maximum(egrid.ω)
end

@testset "deprecated grid constructors" begin
λ0 = 800e-9
λlims = (160e-9, 3000e-9)
trange = 1e-12

#= Four positional numbers, the third a tuple, must still reach the deprecated method and
   produce the same grid as the new one, minus the discarded zmax. =#
rgrid_old = @test_logs((:warn, r"deprecated"), match_mode=:any, min_level=Logging.Warn,
                       Grid.RealGrid(1.0, λ0, λlims, trange))
rgrid_new = Logging.with_logger(Logging.NullLogger()) do
    Grid.RealGrid(λ0, λlims, trange)
end
for f in fieldnames(Grid.RealGrid)
    @test getfield(rgrid_old, f) == getfield(rgrid_new, f)
end

egrid_old = @test_logs((:warn, r"deprecated"), match_mode=:any, min_level=Logging.Warn,
                       Grid.EnvGrid(1.0, λ0, λlims, trange; thg=true))
egrid_new = Logging.with_logger(Logging.NullLogger()) do
    Grid.EnvGrid(λ0, λlims, trange; thg=true)
end
for f in fieldnames(typeof(egrid_new))
    @test getfield(egrid_old, f) == getfield(egrid_new, f)
end

# the deprecated RealGrid method also takes δt positionally, as it always did
rgrid_oldδ = Logging.with_logger(Logging.NullLogger()) do
    Grid.RealGrid(1.0, λ0, λlims, trange, 1e-17)
end
rgrid_newδ = Logging.with_logger(Logging.NullLogger()) do
    Grid.RealGrid(λ0, λlims, trange; δt=1e-17)
end
@test rgrid_oldδ.to == rgrid_newδ.to
end

@testset "grid serialisation" begin
λ0 = 800e-9
λlims = (160e-9, 3000e-9)
trange = 1e-12

for (T, grid) in ((Grid.RealGrid, Logging.with_logger(Logging.NullLogger()) do
                       Grid.RealGrid(λ0, λlims, trange)
                   end),
                  (Grid.EnvGrid, Logging.with_logger(Logging.NullLogger()) do
                       Grid.EnvGrid(λ0, λlims, trange)
                   end))
    d = Grid.to_dict(grid)
    @test !haskey(d, "zmax")
    rt = Grid.from_dict(T, d)
    for f in fieldnames(typeof(grid))
        @test getfield(rt, f) == getfield(grid, f)
    end
    #= Output files written before zmax left the grid carry it in the "grid" group, because
       to_dict iterates over the fields. They must still load. =#
    dlegacy = copy(d)
    dlegacy["zmax"] = 0.5
    legacy = Grid.from_dict(T, dlegacy)
    for f in fieldnames(typeof(grid))
        @test getfield(legacy, f) == getfield(grid, f)
    end
end
end

@testset "zmax in Luna.run" begin
grid = Logging.with_logger(Logging.NullLogger()) do
    Grid.EnvGrid(800e-9, (400e-9, 2e-6), 1e-12)
end
zmax = 0.1
FT = FFTW.plan_fft(zeros(ComplexF64, length(grid.t)))
linop = zeros(ComplexF64, length(grid.ω))
transform = (nl, Eω, z) -> fill!(nl, 0)
mkEω() = zeros(ComplexF64, length(grid.ω))

# zmax is required: the grid no longer carries it
@test_throws ErrorException Logging.with_logger(Logging.NullLogger()) do
    Luna.run(mkEω(), grid, linop, transform, FT, Output.MemoryOutput(0, zmax, 3))
end

# zmax must agree with the end of the output's save grid
@test_throws ErrorException Logging.with_logger(Logging.NullLogger()) do
    Luna.run(mkEω(), grid, linop, transform, FT, Output.MemoryOutput(0, 2zmax, 3); zmax)
end

# and the propagation length is recorded in the output
out = Output.MemoryOutput(0, zmax, 3)
Logging.with_logger(Logging.NullLogger()) do
    Luna.run(mkEω(), grid, linop, transform, FT, out; zmax)
end
@test out["zmax"] == zmax
@test out["z"][end] ≈ zmax

# whatever the caller wrote, the length reaches the output as a Float64, as it did when the
# grid carried it through float(zmax)
outi = Output.MemoryOutput(0, 1, 3)
Logging.with_logger(Logging.NullLogger()) do
    Luna.run(mkEω(), grid, linop, transform, FT, outi; zmax=1)
end
@test outi["zmax"] === 1.0
end

#= A wrapper which forwards everything, including `haskey`: several outputs can be driven at
   once this way and the result is still queryable, so the guards in `Luna.run` still work.
   Defined at top level because a struct cannot be defined inside a testset. =#
struct ForwardingOutput{oT} <: Output.AbstractOutput
    out::oT
end
(w::ForwardingOutput)(args...; kwargs...) = w.out(args...; kwargs...)
Base.getindex(w::ForwardingOutput, k) = w.out[k]
Base.haskey(w::ForwardingOutput, k) = haskey(w.out, k)
Output.check_cache(w::ForwardingOutput, y, t, dt) = Output.check_cache(w.out, y, t, dt)

@testset "zmax in an HDF5Output" begin
grid = Logging.with_logger(Logging.NullLogger()) do
    Grid.EnvGrid(800e-9, (400e-9, 2e-6), 1e-12)
end
zmax = 0.1
FT = FFTW.plan_fft(zeros(ComplexF64, length(grid.t)))
linop = zeros(ComplexF64, length(grid.ω))
transform = (nl, Eω, z) -> fill!(nl, 0)
mkEω() = zeros(ComplexF64, length(grid.ω))

dirpath = joinpath(homedir(), ".luna", "zmax_test")
isdir(dirpath) && rm(dirpath; recursive=true)
mkpath(dirpath)

warned(logger) = any(logger.logs) do r
    occursin("already has dataset zmax", string(r.message))
end

function run_to(output)
    logger = TestLogger(min_level=Logging.Warn)
    Logging.with_logger(logger) do
        Luna.run(mkEω(), grid, linop, transform, FT, output; zmax)
    end
    logger
end

try
    # the propagation length is in the file
    fpath = joinpath(dirpath, "direct.h5")
    out = Output.HDF5Output(fpath, 0, zmax, 3)
    @test !warned(run_to(out))
    @test out["zmax"] == zmax
    @test out["z"][end] ≈ zmax

    #= Resuming: the first run left a completed cache in the file, so this one picks it up
       and must not write zmax a second time. =#
    out2 = Output.HDF5Output(fpath, 0, zmax, 3)
    @test Output.hasdata(out2, "zmax")
    @test !warned(run_to(out2))
    @test out2["zmax"] == zmax

    # the same through a wrapper which forwards haskey
    fpathw = joinpath(dirpath, "wrapped.h5")
    w = ForwardingOutput(Output.HDF5Output(fpathw, 0, zmax, 3))
    @test !warned(run_to(w))
    @test w["zmax"] == zmax
    w2 = ForwardingOutput(Output.HDF5Output(fpathw, 0, zmax, 3))
    @test Output.hasdata(w2, "zmax")
    @test !warned(run_to(w2))

    #= A bare closure cannot be queried, so the guard cannot fire and the second run writes
       again and warns. Documented behaviour, asserted so that it stays documented. =#
    fpathc = joinpath(dirpath, "closure.h5")
    oc = Output.HDF5Output(fpathc, 0, zmax, 3)
    closure(args...; kwargs...) = oc(args...; kwargs...)
    @test !Output.hasdata(closure, "zmax")
    @test !warned(run_to(closure))
    oc2 = Output.HDF5Output(fpathc, 0, zmax, 3)
    closure2(args...; kwargs...) = oc2(args...; kwargs...)
    @test warned(run_to(closure2))
    @test oc2["zmax"] == zmax # written again, with the same value
finally
    rm(dirpath; recursive=true, force=true)
end
end
