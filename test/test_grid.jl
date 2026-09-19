import Test: @test, @test_throws, @test_logs, @testset
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
end
