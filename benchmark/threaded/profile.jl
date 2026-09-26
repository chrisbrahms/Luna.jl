# Where does per-step time go? (julia --project=. -t 1 benchmark/threaded/profile.jl)
# Samples classified by the innermost frame that matches a
# category; broadcasts/reductions also attributed to the Luna source file that called them.
import LinearAlgebra
using Luna, Profile, Printf
Luna.set_fftw_mode(:estimate); Luna.set_fftw_threads(parse(Int, get(ENV, "FFTWT", "1")))
Luna.set_fftw_wisdom(false)
LinearAlgebra.BLAS.set_num_threads(parse(Int, get(ENV, "BLAST", "1")))

const CATS = [
    ("fft",       fr -> occursin("unsafe_execute", string(fr.func)) || occursin("fftw", lowercase(string(fr.file)))),
    ("gemm",      fr -> occursin("gemm", string(fr.func)) || occursin("gemv", string(fr.func))),
    ("scan",      fr -> occursin("accumulate", string(fr.func))),
    ("reduce",    fr -> string(fr.func) in ("mapreduce", "_mapreduce", "mapreduce_impl", "_mapreduce_dim", "_zipreduce", "foldl_impl")),
    ("broadcast", fr -> string(fr.func) in ("copyto!", "materialize!", "materialize") && occursin("broadcast.jl", string(fr.file))),
    ("gc",        fr -> occursin("gc", string(fr.func)) && fr.from_c),
]
lunafile(fr) = (f = string(fr.file); occursin("/src/", f) && occursin("Luna", f)) ? basename(f) : nothing

function classify(f)
    f()                                   # compile
    Profile.init(n=10^8, delay=0.0005)
    Profile.clear()
    t = @elapsed @profile f()
    data = Profile.fetch(include_meta=false)
    cache = Dict{UInt64, Vector{Base.StackTraces.StackFrame}}()
    counts = Dict{String, Int}(); callers = Dict{String, Int}(); total = 0
    i = 1
    while i <= length(data)
        j = findnext(==(0), data, i); j === nothing && break
        frames = Base.StackTraces.StackFrame[]
        for ip in data[i:j-1]
            append!(frames, get!(() -> Profile.lookup(ip), cache, ip))
        end
        i = j + 1
        isempty(frames) && continue
        total += 1
        cat = "other"; k0 = 0
        for (k, fr) in enumerate(frames)
            for (c, test) in CATS
                if test(fr); cat = c; k0 = k; break; end
            end
            k0 > 0 && break
        end
        if cat == "other"   # attribute "other" to the innermost Luna file
            for fr in frames
                lf = lunafile(fr); lf === nothing || (cat = "other:" * lf; break)
            end
        end
        counts[cat] = get(counts, cat, 0) + 1
        if cat in ("broadcast", "reduce", "scan")
            for fr in frames[k0:end]
                lf = lunafile(fr)
                if lf !== nothing
                    key = cat * " <- " * lf * ":" * string(fr.func)
                    callers[key] = get(callers, key, 0) + 1
                    break
                end
            end
        end
    end
    @printf("wall %.2f s, %d samples\n", t, total)
    for (k, v) in sort(collect(counts), by=last, rev=true)
        v/total > 0.005 && @printf("  %-40s %5.1f %%\n", k, 100v/total)
    end
    println("  -- broadcast/reduce/scan by Luna caller")
    for (k, v) in sort(collect(callers), by=last, rev=true)[1:min(end, 15)]
        v/total > 0.005 && @printf("  %-60s %5.1f %%\n", k, 100v/total)
    end
    flush(stdout)
end

# 1. mode-averaged field, Kerr + PPT plasma (default prop_capillary physics for Ar)
modeavg() = prop_capillary(125e-6, 0.15, :Ar, 1.0; λ0=800e-9, τfwhm=30e-15, energy=150e-6,
                           λlims=(150e-9, 4e-6), trange=1e-12, shotnoise=false,
                           status_period=1e9)
# 2. four-mode fixed-rule multimode, Kerr + plasma
modal() = prop_capillary(125e-6, 0.05, :Ar, 1.0; λ0=800e-9, τfwhm=30e-15, energy=150e-6,
                         λlims=(150e-9, 4e-6), trange=0.5e-12, shotnoise=false, modes=4,
                         modal_integral=:fixed, status_period=1e9)
# 3. radial free space, Kerr + plasma
function radial()
    gas = :Ar; pres = 1.2; λ0 = 800e-9; L = 0.05
    grid = Grid.RealGrid(800e-9, (400e-9, 2000e-9), 0.2e-12)
    q = Grid.RadialGrid(4e-3, 256)
    densityfun = let d = PhysData.density(gas, pres); z -> d; end
    ionrate = Ionisation.IonRatePPTCached(gas, λ0)
    responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
                 Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, PhysData.ionisation_potential(gas)))
    linop = LinearOps.make_const_linop(grid, q, PhysData.ref_index_fun(gas, pres))
    normfun = NonlinearRHS.const_norm_radial(grid, q, PhysData.ref_index_fun(gas, pres))
    inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=20e-15, energy=20e-6, w0=200e-6, propz=-0.02)
    Eω, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)
    statsfun = Stats.default(grid, Eω, linop, transform; gas=gas)
    output = Output.MemoryOutput(0, L, 11, statsfun)
    Luna.run(Eω, grid, linop, transform, FT, output; zmax=L, status_period=1e9)
end

for name in split(get(ENV, "CASES", "modeavg,modal,radial"), ",")
    println("==== ", name, "  (julia threads ", Threads.nthreads(), ")")
    classify(getfield(Main, Symbol(name)))
end
println("PROFILE_DONE")
