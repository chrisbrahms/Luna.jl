#= Whole-propagation wall time against FFTW and BLAS thread counts, for one Julia thread
   count per process. Appends to benchmark/threads/results/runs.csv.

       for J in 1 4 8; do julia --project=. -t $J benchmark/threads/runs.jl; done

   Environment: CASES=a,b (subset), FLAG=estimate|patient (FFTW planning; default
   estimate; with patient the time includes planning, since wisdom is off), DEVICE=metal
   (run the device-capable cases on Metal; needs an environment with Metal, e.g.
   `--project=<env>` with Luna developed into it), REPS (default 2, best of).

   Each configuration replans (`Luna.setup` reads the FFTW thread count before planning).
   Outputs are compared with the first configuration of the same case: threaded FFTW plans
   can differ at rounding level. =#
import LinearAlgebra
import LinearAlgebra: BLAS
using Luna, Printf
const DEVICE = get(ENV, "DEVICE", "cpu")
#= Loading Metal sets settings["device"] = :auto, so every case below runs on the GPU
   (Float32) without further changes. =#
DEVICE == "metal" && @eval using Metal

const J = Threads.nthreads()
const FLAG = get(ENV, "FLAG", "estimate")
const REPS = parse(Int, get(ENV, "REPS", "2"))
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results", DEVICE == "metal" ? "runs_metal.csv" : "runs.csv"))
Luna.set_fftw_mode(Symbol(FLAG)); Luna.set_fftw_wisdom(false)
#= J1_FFTW=1: let FFTW use its own pthreads with one Julia thread, which Luna's
   Utils.FFTWthreads forbids (it returns 1 when nthreads() == 1). Benchmark-only override. =#
const J1_FFTW = get(ENV, "J1_FFTW", "0") == "1"
J1_FFTW && @eval Luna.Utils FFTWthreads() = settings["fftw_threads"] == 0 ? 1 : settings["fftw_threads"]

include(joinpath(@__DIR__, "..", "threaded", "cases.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "regression", "cases.jl"))
const RC = RegressionCases
quiet(f) = RC.quiet(f)

# a 3-D free-space envelope case on a 64 × 64 transverse grid (the gate's is 16 × 8)
function free3d_big()
    grid = RC.makegrid(Grid.EnvGrid, RC.L_FREE, RC.Λ0, RC.ΛLIMS_FREE, RC.TRANGE_FREE)
    sg = Grid.FreeGrid(RC.R_FREE, 64, RC.R_FREE, 64)
    nfunλ = PhysData.ref_index_fun(RC.GAS_FREE, RC.P_FREE)
    nfun = (λ; z=0.0) -> nfunλ(λ)
    Eω, grid, linop, transform, FT, output = RC.setup_free(grid, RC.L_FREE, sg,
        NonlinearRHS.const_norm_free(grid, sg, nfun), (Nonlinear.Kerr_env(PhysData.γ3_gas(RC.GAS_FREE)),))
    quiet() do
        Luna.run(Eω, grid, linop, transform, FT, output; zmax=RC.L_FREE, status_period=1e9)
    end
    output
end

gate(name) = () -> RC.runcase(RC.getcase(name), :adaptive)

# name => (function, uses GEMM)
const ALLCASES = [
    "modeavg"      => (() -> quiet(() -> modeavg()), false),
    "modeavg_long" => (() -> quiet(() -> modeavg_long()), false),
    "modeavg_env_raman" => (gate("modeavg_env_raman"), false),
    "gnlse"        => (gate("gnlse_raman_shock"), false),
    "modal_fixed"  => (() -> quiet(() -> modal()), true),
    "modal_adaptive" => (gate("multimode_field_plasma"), true),
    "radial"       => (() -> quiet(() -> radial()), true),
    "radial_big"   => (() -> quiet(() -> radial_big()), true),
    "free2d_chi2"  => (gate("free2d_field_chi2"), false),
    "free3d"       => (gate("free3d_env_kerr"), false),
    "free3d_big"   => (free3d_big, false),
]
const METALCASES = ("modeavg", "modeavg_long", "modal_fixed", "radial", "radial_big")

field(out) = out isa Output.MemoryOutput ? out.data["Eω"] : out["Eω"]

fftw_settings() = J == 1 ? (J1_FFTW ? [1, 2, 4, 8] : [1]) : sort(unique([1, J, 2J, 4J]))
blas_settings(gemm) = gemm ? sort(unique([1, J, 8])) : [1]

function main()
    sel = haskey(ENV, "CASES") ? split(ENV["CASES"], ",") :
          DEVICE == "metal" ? collect(METALCASES) : first.(ALLCASES)
    isfile(OUT) || open(io -> println(io,
        "device,flag,julia_threads,fftw_threads,blas_threads,case,time_s,maxreldiff"), OUT, "w")
    for (name, (f, gemm)) in ALLCASES
        name in sel || continue
        ref = nothing
        for nf in fftw_settings(), nb in blas_settings(gemm)
            Luna.set_fftw_threads(nf); BLAS.set_num_threads(nb)
            ref === nothing && f() # compile once per case
            # best of REPS, and of as many more (up to 20) as fit in one second for fast cases
            t = Inf; E = nothing; total = 0.0; k = 0
            while k < REPS || (total < 1.0 && k < 20)
                t1 = @elapsed out = f()
                t = min(t, t1); E = Array(field(out)); total += t1; k += 1
            end
            ref === nothing && (ref = E)
            d = maximum(abs, E .- ref)/maximum(abs, ref)
            open(OUT, "a") do io
                @printf(io, "%s,%s,%d,%d,%d,%s,%.4f,%.2e\n", DEVICE * (J1_FFTW ? "_pthreads" : ""), FLAG, J, nf, nb, name, t, d)
            end
            @printf("J=%d fftw=%2d blas=%d %-16s %8.3f s  diff %.1e\n", J, nf, nb, name, t, d)
            flush(stdout)
        end
    end
end

main()
println("RUNS_DONE J=$J")
