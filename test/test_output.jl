import Test: @test, @testset, @test_throws
import Luna: Output, Processing
using EllipsisNotation

@testset "HDF5" begin
    import HDF5
    import Luna: Utils
    fpath = joinpath(tempname(), "test.h5")
    isfile(fpath) && rm(fpath)
    shape = (1024, 4, 2)
    n = 11
    stat = randn()
    statsfun(y, t, dt) = Dict("stat" => stat, "stat2d" => [stat, stat])
    t0 = 0
    t1 = 10
    t = collect(range(t0, stop=t1, length=n))
    ω = randn((1024,))
    wd = dirname(@__FILE__)
    gitc = Utils.git_commit()
    o = Output.HDF5Output(fpath, t0, t1, n, statsfun; yname="y", tname="t")
    extra = Dict()
    extra["ω"] = ω
    extra["git_commit"] = gitc
    o(extra)
    meta = Dict()
    meta["meta1"] = 100
    meta["meta2"] = "src"
    o(meta, meta=true)
    y0 = randn(ComplexF64, shape)
    y(t) = y0
    for (ii, ti) in enumerate(t)
        o(y0, ti, 0, y)
    end
    @test o(extra, force=true) === nothing
    HDF5.h5open(fpath, "r") do file
        @test all(read(file["t"]) == t)
        global yr = read(file["y"])
        @test all([all(yr[:, :, :, ii] == y0) for ii=1:n])
        @test all(ω == read(file["ω"]))
        @test gitc == read(file["git_commit"])
        @test all(read(file["stats"]["stat"]) .== stat)
        @test all(read(file["stats"]["stat2d"]) .== stat)
        @test Utils.git_commit() == read(file["meta"]["git_commit"])
        # Need to strip out date from sourcecode to compare
        src = read(file["meta"]["sourcecode"])
        @test split(Utils.sourcecode(), '\n')[2:end] == split(src, '\n')[2:end]
        @test 100 == read(file["meta"]["meta1"])
        @test "src" == read(file["meta"]["meta2"])
    end
    @test all(o["stats", "stat"] .== stat)
    @test all(yr .== o["y"])
    @test 100 == o["meta"]["meta1"]
    rm(fpath)
    rm(splitdir(fpath)[1], force=true)
end

@testset "Memory" begin
    import Luna: Utils
    shape = (1024, 4, 2)
    n = 11
    stat = randn()
    statsfun(y, t, dt) = Dict("stat" => stat, "stat2d" => [stat, stat])
    t0 = 0
    t1 = 10
    t = collect(range(t0, stop=t1, length=n))
    ω = randn((1024,))
    wd = dirname(@__FILE__)
    gitc = Utils.git_commit()
    o = Output.MemoryOutput(t0, t1, n, statsfun, yname="y", tname="t")
    extra = Dict()
    extra["ω"] = ω
    extra["git_commit"] = gitc
    o(extra)
    meta = Dict()
    meta["meta1"] = 100
    meta["meta2"] = "src"
    o(meta, meta=true)
    y0 = randn(ComplexF64, shape)
    y(t) = y0
    for (ii, ti) in enumerate(t)
        o(y0, ti, 0, y)
    end
    @test all(o.data["t"] == t)
    @test all([all(o.data["y"][:, :, :, ii] == y0) for ii=1:n])
    @test all(ω == o.data["ω"])
    @test gitc == o.data["git_commit"]
    @test_throws ErrorException o(extra)
    @test o(extra, force=true) === nothing
    @test all(o.data["stats"]["stat"] .== stat)
    @test all(o.data["stats"]["stat2d"] .== stat)
    @test all(o["stats", "stat"] .== stat)
    @test Utils.git_commit() == o.data["meta"]["git_commit"]
    # Need to strip out date from sourcecode to compare
    src = o.data["meta"]["sourcecode"]
    @test split(Utils.sourcecode(), '\n')[2:end] == split(src, '\n')[2:end]
    @test 100 == o.data["meta"]["meta1"]
    @test "src" == o.data["meta"]["meta2"]
end

dirpath = tempname()
fpath = joinpath(dirpath, "test.h5")
fpath_comp = joinpath(dirpath, "test_comp.h5")
@testset "HDF5 vs Memory" begin
    using Luna
    import FFTW
    import HDF5
    import LinearAlgebra: norm

    a = 13e-6
    gas = :Ar
    pres = 5
    τ = 30e-15
    λ0 = 800e-9
    grid = Grid.RealGrid(800e-9, (160e-9, 3000e-9), 1e-12)
    m = Capillary.MarcatiliMode(a, gas, pres, loss=false)
    aeff = let m=m
        z -> Modes.Aeff(m, z=z)
    end
    energyfun, energyfunω = Fields.energyfuncs(grid)
    densityfun = let dens0=PhysData.density(gas, pres)
        z -> dens0
    end
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun!, frame_vel, αfun = LinearOps.make_const_linop(grid, m, λ0)

    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                                   Stats.ω0(grid),
                                   Stats.energy(grid, energyfunω))
    hdf5 = Output.HDF5Output(fpath, 0, 5e-2, 51, statsfun)
    hdf5c = Output.HDF5Output(fpath_comp, 0, 5e-2, 51, statsfun,
                              compression=true)
    mem = Output.MemoryOutput(0, 5e-2, 51, statsfun)
    function outfun(args...; kwargs...)
        hdf5(args...; kwargs...)
        hdf5c(args...; kwargs...)
        mem(args...; kwargs...)
    end
    for o in (hdf5, hdf5c, mem)
        o(Dict("λ0" => λ0))
        o("τ", τ)
    end
    Luna.run(Eω, grid, linop, transform, FT, outfun, status_period=10, zmax=5e-2)
    HDF5.h5open(hdf5.fpath, "r") do file
        @test read(file["λ0"]) == mem.data["λ0"]
        Eω = reinterpret(ComplexF64, read(file["Eω"]))
        @test Eω == mem.data["Eω"]
        @test read(file["stats"]["ω0"]) == mem.data["stats"]["ω0"]
        @test read(file["z"]) == mem.data["z"]
        @test read(file["grid"]) == Grid.to_dict(grid)
        @test read(file["simulation_type"]["field"]) == "field-resolved"
        @test read(file["simulation_type"]["transform"]) == string(transform)
    end
    HDF5.h5open(hdf5c.fpath, "r") do file
        @test read(file["λ0"]) == mem.data["λ0"]
        Eω = reinterpret(ComplexF64, read(file["Eω"]))
        @test Eω == mem.data["Eω"]
        @test read(file["stats"]["ω0"]) == mem.data["stats"]["ω0"]
        @test read(file["z"]) == mem.data["z"]
        @test read(file["grid"]) == Grid.to_dict(grid)
        @test read(file["simulation_type"]["field"]) == "field-resolved"
        @test read(file["simulation_type"]["transform"]) == string(transform)
    end
    @test stat(hdf5.fpath).size >= stat(hdf5c.fpath).size
    # Test read-only
    o = Output.HDF5Output(fpath)
    HDF5.h5open(o.fpath, "r") do file
        @test read(file["λ0"]) == mem.data["λ0"]
        Eω = reinterpret(ComplexF64, read(file["Eω"]))
        @test Eω == mem.data["Eω"]
        @test read(file["stats"]["ω0"]) == mem.data["stats"]["ω0"]
        @test read(file["z"]) == mem.data["z"]
        @test read(file["grid"]) == Grid.to_dict(grid)
        @test read(file["simulation_type"]["field"]) == "field-resolved"
        @test read(file["simulation_type"]["transform"]) == string(transform)
    end
    # test slice reading
    @test o["Eω", :, 1] == mem["Eω"][:, 1]
    @test o["Eω", 1, :] == mem["Eω"][1, :]
    @test o["Eω", :, 1:5] == mem["Eω"][:, 1:5]
    @test o["Eω", :, [1, 2, 50]] == mem["Eω"][:, [1, 2, 50]]
    @test o["Eω", .., 1] == mem["Eω"][:, 1]
    @test o["Eω", .., [1, 2, 50]] == mem["Eω"][:, [1, 2, 50]]
    @test mem["Eω", :, 1] == mem["Eω"][:, 1]
    @test mem["Eω", 1, :] == mem["Eω"][1, :]
    @test mem["Eω", :, 1:5] == mem["Eω"][:, 1:5]
    @test mem["Eω", :, [1, 2, 50]] == mem["Eω"][:, [1, 2, 50]]
    @test mem["Eω", .., 1] == mem["Eω"][:, 1]
    @test mem["Eω", .., [1, 2, 50]] == mem["Eω"][:, [1, 2, 50]]
    # test slice reading in processing functions
    ω, Eω, zac = Processing.getEω(o, 5e-2)
    @test (ω, Eω[:, 1]) == Processing.getEω(grid, mem["Eω"][:, 51])
    @test zac[1] == 5e-2
end
rm(fpath)
rm(fpath_comp)
rm(splitdir(fpath)[1], force=true)

##
fpath = joinpath(homedir(), ".luna", "output_test", "test.h5")
@testset "Continuing" begin
    using Luna
    import FFTW
    import HDF5

    a = 13e-6
    gas = :Ar
    pres = 5
    τ = 30e-15
    λ0 = 800e-9
    grid = Grid.RealGrid(800e-9, (160e-9, 3000e-9), 1e-12)
    m = Capillary.MarcatiliMode(a, gas, pres, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    energyfun, energyfunω = Fields.energyfuncs(grid)
    dens0 = PhysData.density(gas, pres)
    densityfun(z) = dens0
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
    linop, βfun!, frame_vel, αfun = LinearOps.make_const_linop(grid, m, λ0)

    # Run with arbitrary error at 3 cm
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                                   Stats.ω0(grid),
                                   Stats.energy(grid, energyfunω))
    output = Output.HDF5Output(fpath, 0, 5e-2, 51, statsfun)
    function stepfun(Eω, z, dz, interpolant)
        output(Eω, z, dz, interpolant)
        if z > 3e-2
            error("Oh no!")
        end
    end
    stepfun(args...; kwargs...) = output(args...; kwargs...)
    try
        Luna.run(Eω, grid, linop, transform, FT, stepfun, status_period=10, z0=0.0, zmax=5e-2)
    catch
    end

    # Run again, starting from 3 cm
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                                   Stats.ω0(grid),
                                   Stats.energy(grid, energyfunω))
    output = Output.HDF5Output(fpath, 0, 5e-2, 51, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output, status_period=5, zmax=5e-2)

    # Run from scratch with MemoryOutput
    inputs = Fields.GaussField(λ0=λ0, τfwhm=τ, energy=1e-6)
    Eω, transform, FT = Luna.setup(
        grid, densityfun, responses, inputs, βfun!, aeff)
    statsfun = Stats.collect_stats(grid, Eω,
                                   Stats.ω0(grid),
                                   Stats.energy(grid, energyfunω))
    mem = Output.MemoryOutput(0, 5e-2, 51, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, mem, status_period=5, zmax=5e-2)

    idx1 = findfirst(grid.ωwin .!= 1)
    idx2 = findlast(grid.ωwin .== 1)
    Eωm = mem["Eω"][idx1:idx2, :]
    Eω = output["Eω"][idx1:idx2, :]
    Iω = abs2.(Eω)
    Iωm = abs2.(Eωm)
    @test norm(Iω - Iωm)/norm(Iω) < 1e-7
    @test all(isapprox.(Eωm, Eω, atol=1e-4*maximum(abs.(Eωm))))
    @test all(isapprox.(output["stats"]["ω0"], mem.data["stats"]["ω0"], rtol=1e-6))
    @test all(output["stats"]["energy"] .≈ mem.data["stats"]["energy"])
    @test output["z"] == mem.data["z"]
    @test output["grid"] == Grid.to_dict(grid)
    @test output["simulation_type"]["field"] == "field-resolved"
    @test output["simulation_type"]["transform"] == string(transform)
end
rm(fpath, force=true)
rm(splitdir(fpath)[1], force=true)

@testset "willsave" begin
    # GridCondition's own points: 0, 0.25, 0.5, 0.75, 1.0. `willsave` does not itself
    # save, so it has to be interleaved with real saves to track `o.saved` correctly,
    # exactly as `Luna.ScaledOutput` interleaves it with the real per-step call.
    o = Output.MemoryOutput(0, 1.0, 5, Output.nostats)
    y = randn(ComplexF64, 8)
    @test Output.willsave(o, y, 0.0, 0.1) == true # the first point, t = 0, is reached
    o(y, 0.0, 0.1, _ -> y)
    @test o.saved == 1
    @test Output.willsave(o, y, 0.1, 0.1) == false # next point (0.25) not reached yet
    @test Output.willsave(o, y, 0.25, 0.1) == true # now it is

    h = Output.HDF5Output(joinpath(tempname(), "willsave.h5"), 0, 1.0, 5, Output.nostats)
    @test Output.willsave(h, y, 0.0, 0.1) == true
    h(y, 0.0, 0.1, _ -> y)
    @test Output.willsave(h, y, 0.1, 0.1) == false
    @test Output.willsave(h, y, 0.25, 0.1) == true
    rm(dirname(h.fpath), recursive=true, force=true)

    # A save condition willsave does not know how to inspect falls back to `true`
    @test Output.willsave(Output.MemoryOutput(Output.always, "Eω", "z"), y, 0.3, 0.1) == true
    # Any other kind of output (a bare function, say) is conservatively `true`
    @test Output.willsave((args...; kwargs...) -> nothing, y, 0.3, 0.1) == true
end

@testset "PeriodicStats" begin
    calls = Ref(0)
    f(y, t, dt) = (calls[] += 1; Dict("n" => calls[]))
    p = Output.PeriodicStats(f, 3)
    # first call always fires, then every 3rd
    results = [p(nothing, 0.0, 0.0) for _ in 1:7]
    @test [isnothing(r) for r in results] == [false, true, true, false, true, true, false]
    # `calls` only increments when `f` actually runs, i.e. on the 1st, 4th and 7th call
    @test [r["n"] for r in results if !isnothing(r)] == [1, 2, 3]

    # willfire predicts the next call correctly, without mutating or running it
    p2 = Output.PeriodicStats(f, 3)
    fires = Bool[]
    for _ in 1:7
        push!(fires, Output.willfire(p2, 0.0))
        p2(nothing, 0.0, 0.0)
    end
    @test fires == [true, false, false, true, false, false, true]

    # validation: an integer step count must be >= 1
    @test_throws ArgumentError Output.PeriodicStats(f, 0)
    @test_throws ArgumentError Output.PeriodicStats(f, -3)
    # a non-integer (distance) period must be > 0
    @test_throws ArgumentError Output.PeriodicStats(f, 0.0)
    @test_throws ArgumentError Output.PeriodicStats(f, -2.5)

    # the "every Δz" form: fires when the propagation coordinate has advanced by period
    calls2 = Ref(0)
    g(y, t, dt) = (calls2[] += 1; Dict("z" => t))
    pd = Output.PeriodicStats(g, 0.5) # non-integer -> distance mode, every 0.5 m
    zs = [0.0, 0.2, 0.5, 0.6, 1.0, 1.1, 1.4]
    rd = [pd(nothing, z, 0.0) for z in zs]
    # fires at z=0.0 (first call), 0.5 (>= 0.5 since last fire) and 1.0 (>= 0.5 since 0.5);
    # 1.4 is only 0.4 past the fire at 1.0, so it does not fire
    @test [isnothing(r) for r in rd] == [false, true, false, true, false, true, true]
    @test [r["z"] for r in rd if !isnothing(r)] == [0.0, 0.5, 1.0]
    @test Output.willfire(pd, 1.49) == false
    @test Output.willfire(pd, 1.5) == true

    # maybe_periodic: the trivial integer 1 is not wrapped at all (zero overhead), a
    # non-trivial value is, and an invalid value is refused either way
    @test Output.maybe_periodic(f, 1) === f
    @test Output.maybe_periodic(f, 1.0) === f
    @test Output.maybe_periodic(f, 2) isa Output.PeriodicStats
    @test Output.maybe_periodic(f, 0.3) isa Output.PeriodicStats
    @test_throws ArgumentError Output.maybe_periodic(f, 0)
    @test_throws ArgumentError Output.maybe_periodic(f, -1.0)
    @test_throws ArgumentError Output.maybe_periodic(f, "3")

    # MemoryOutput/HDF5Output skip a `nothing` statistics result instead of erroring
    o = Output.MemoryOutput(0, 1.0, 4, Output.PeriodicStats((y, t, dt) -> Dict("s" => t), 2))
    y0 = randn(ComplexF64, 4)
    for t in (0.0, 1/3, 2/3, 1.0)
        o(y0, t, 0.1, _ -> y0)
    end
    @test o.data["stats"]["s"] == [0.0, 2/3] # only the 1st and 3rd calls collected stats
end

@testset "Float32 output eltype" begin
    using Luna
    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, :He, 1.0, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    dens = z -> PhysData.density(:He, 1.0)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, 800e-9)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6)

    Eω64, transform64, FT64 = Luna.setup(grid, dens, resp, inputs, βfun!, aeff; constβ=true)
    out64 = Output.MemoryOutput(0, 1e-2, 3, Output.nostats)
    Luna.run(Eω64, grid, linop, transform64, FT64, out64;
             zmax=1e-2, boundary=:rate, init_dz=5e-4, rtol=1e-8)

    Eω32, transform32, FT32 = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                         constβ=true, precision=Float32)
    out32 = Output.MemoryOutput(0, 1e-2, 3, Output.nostats)
    Luna.run(Eω32, grid, linop, transform32, FT32, out32;
             zmax=1e-2, boundary=:rate, init_dz=5e-4, rtol=1e-8)

    # ScaledOutput unscales and Output allocates with eltype(y): Float32 saves Float32
    @test eltype(out32["Eω"]) === ComplexF32
    @test eltype(out64["Eω"]) === ComplexF64
    for idx in axes(out64["Eω"], 2)
        h = out64["Eω"][:, idx]
        d = ComplexF64.(out32["Eω"][:, idx])
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-5
    end
end

@testset "HDF5 resume of a Float32 run" begin
    using Luna
    import HDF5
    fdir = tempname()
    fpath = joinpath(fdir, "resume32.h5")

    grid = Grid.RealGrid(800e-9, (300e-9, 2000e-9), 400e-15)
    m = Capillary.MarcatiliMode(75e-6, :He, 1.0, loss=false)
    aeff(z) = Modes.Aeff(m, z=z)
    dens = z -> PhysData.density(:He, 1.0)
    resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:He)),)
    linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, 800e-9)
    inputs = Fields.GaussField(λ0=800e-9, τfwhm=20e-15, energy=1e-6)

    # Run with an arbitrary failure partway through, so the file is left with a cache
    Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                   constβ=true, precision=Float32)
    output = Output.HDF5Output(fpath, 0, 1e-2, 11, Output.nostats)
    function stepfun(Eω, z, dz, interpolant)
        output(Eω, z, dz, interpolant)
        z > 6e-3 && error("interrupted")
    end
    # Luna.run treats its 6th argument as the output object for metadata, willsave and
    # check_cache too, so this stand-in must forward everything else to the real one.
    stepfun(args...; kwargs...) = output(args...; kwargs...)
    try
        Luna.run(Eω, grid, linop, transform, FT, stepfun;
                 zmax=1e-2, boundary=:rate, init_dz=5e-4, rtol=1e-8)
    catch e
        e isa ErrorException || rethrow()
    end

    # Resume: check_cache reads the cache back (host, physical units, ComplexF32) and
    # Luna.run rescales and reuploads it -- on the CPU that upload is the identity, but
    # the rescale (dividing by the *new* run's E_ref, which is deterministic from the
    # same input and therefore the same as the interrupted run's) is not.
    Eω2, transform2, FT2 = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                      constβ=true, precision=Float32)
    output2 = Output.HDF5Output(fpath, 0, 1e-2, 11, Output.nostats)
    Luna.run(Eω2, grid, linop, transform2, FT2, output2;
             zmax=1e-2, boundary=:rate, init_dz=5e-4, rtol=1e-8)

    # An uninterrupted run started fresh, for comparison
    Eω3, transform3, FT3 = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                      constβ=true, precision=Float32)
    output3 = Output.MemoryOutput(0, 1e-2, 11, Output.nostats)
    Luna.run(Eω3, grid, linop, transform3, FT3, output3;
             zmax=1e-2, boundary=:rate, init_dz=5e-4, rtol=1e-8)

    @test eltype(output2["Eω"]) === ComplexF32
    @test output2["z"] ≈ output3["z"]
    for idx in axes(output3["Eω"], 2)
        h = output3["Eω"][:, idx]
        d = output2["Eω"][:, idx]
        @test maximum(abs, d .- h)/maximum(abs, h) < 1e-5
    end
    rm(fdir, recursive=true, force=true)
end
