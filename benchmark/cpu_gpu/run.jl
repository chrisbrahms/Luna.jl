#= Whole-propagation wall time on the CPU (with Luna's default threading) and on a GPU, for
   the cases in `cases.jl`. One process per setting (device and precision); each process
   appends one row per case to `results/<MACHINE>.csv` and saves the final field of each
   case to `fields/<MACHINE>/` for `summarise.jl`, which compares every run with the
   single-threaded CPU Float64 run (`cpu64_t1`, the baseline) of the same case on the same
   machine. A CPU setting's name carries its Julia thread count (`cpu64_t8`).

   From the repository root, CPU (the benchmark environment; see benchmark/Project.toml):

       julia --project=benchmark -t 1 benchmark/cpu_gpu/run.jl                # baseline
       julia --project=benchmark -t 8 benchmark/cpu_gpu/run.jl                # Float64
       PRECISION=32 julia --project=benchmark -t 8 benchmark/cpu_gpu/run.jl   # Float32

   GPU, from an environment which has Luna (developed from this checkout) and the GPU
   package, since neither Metal nor CUDA is a dependency of Luna, e.g. for CUDA

       julia --project=gpuenv -e 'using Pkg; Pkg.develop(path="."); Pkg.add("CUDA")'

   and then

       DEVICE=metal julia --project=gpuenv -t 8 benchmark/cpu_gpu/run.jl  # Float32 only
       DEVICE=cuda  julia --project=gpuenv -t 8 benchmark/cpu_gpu/run.jl  # Float64
       DEVICE=cuda PRECISION=32 julia --project=gpuenv -t 8 benchmark/cpu_gpu/run.jl

   and `julia --project=benchmark benchmark/cpu_gpu/summarise.jl`. Use the same `-t` as the
   number of performance cores for the CPU runs; `run.jl` records it.

   Environment:
   - `DEVICE` = cpu (default) | metal | cuda; `PRECISION` = 64 (default; 32 for metal) | 32
   - `CASES`: comma-separated subset of the case names below (default: all); cases the
     setting cannot run are skipped (the adaptive multimode rule is host Float64 only)
   - `STATS_N`: statistics points per case (default 0: every step; see `cases.jl`); a
     non-zero value is added to the setting's name, e.g. `metal32_s100`
   - `SCALE`: multiplies every propagation length (default 1)
   - `MACHINE`: label for the results (default: the host name)

   The CPU runs use Luna's defaults (`:measure` planning, automatic FFTW and BLAS threads,
   threaded broadcasts), except that FFTW wisdom is off so that no run depends on what an
   earlier one planned. Setup (including planning) and propagation are timed separately;
   each case is first run over a short length to compile it, and the timed run is a single
   full propagation. The comparison is of the time a user waits, so the step count is not
   fixed: it is reported, and a run in Float32 may take a different number of steps. =#
import LinearAlgebra
using Luna, Printf, Serialization
const DEVICE = Symbol(get(ENV, "DEVICE", "cpu"))
DEVICE == :metal && @eval using Metal
DEVICE == :cuda && @eval using CUDA
const PRECISION = get(ENV, "PRECISION", DEVICE == :metal ? "32" : "64") == "32" ? Float32 : Float64
const SCALE = parse(Float64, get(ENV, "SCALE", "1"))
const MACHINE = get(ENV, "MACHINE", replace(gethostname(), r"\.local$" => ""))
const OUT = joinpath(@__DIR__, "results", MACHINE * ".csv")
const FIELDS = joinpath(@__DIR__, "fields", MACHINE)
Luna.set_fftw_wisdom(false)

include(joinpath(@__DIR__, "cases.jl"))

const ALLCASES = [
    "modeavg"        => (; kw...) -> modeavg(; kw...),
    "modeavg_long"   => (; kw...) -> modeavg_long(; kw...),
    "modal_adaptive" => (; kw...) -> modal(:adaptive, 0; kw...),
    "modal_fixed64"  => (; kw...) -> modal(:fixed, 64; kw...),
    "modal_fixed128" => (; kw...) -> modal(:fixed, 128; kw...),
    "radial256"      => (; kw...) -> radial(256; kw...),
    "radial512"      => (; kw...) -> radial(512; kw...),
    "radial1024"     => (; kw...) -> radial(1024; kw...),
    "free3d64"       => (; kw...) -> free3d(64; kw...),
    "free3d128"      => (; kw...) -> free3d(128; kw...),
    "free3d256"      => (; kw...) -> free3d(256; kw...),
]
"Cases which run only on the host in Float64."
const HOSTONLY = ("modal_adaptive",)

field(out) = Array(out.data["Eω"][ntuple(_ -> Colon(), ndims(out.data["Eω"]) - 1)..., end])

function main()
    setting = string(DEVICE, PRECISION == Float32 ? "32" : "64",
                     DEVICE == :cpu ? "_t$(Threads.nthreads())" : "",
                     STATSN > 0 ? "_s$STATSN" : "")
    sel = haskey(ENV, "CASES") ? split(ENV["CASES"], ",") : first.(ALLCASES)
    hostf64 = DEVICE == :cpu && PRECISION == Float64
    sel = [c for c in sel if hostf64 || !(c in HOSTONLY)]
    mkpath(dirname(OUT)); mkpath(FIELDS)
    isfile(OUT) || open(io -> println(io,
        "machine,setting,julia_threads,case,scale,setup_s,run_s,steps,ms_per_step"), OUT, "w")
    for (name, f) in ALLCASES
        name in sel || continue
        # one failing case (e.g. out of device memory) must not lose the rest of the setting
        ts, tr, n, out = try
            f(; scale=0.02SCALE, device=DEVICE, precision=PRECISION) # compile
            GC.gc()
            f(; scale=SCALE, device=DEVICE, precision=PRECISION)
        catch err
            @error "case $name failed" exception=(err, catch_backtrace())
            GC.gc()
            continue
        end
        serialize(joinpath(FIELDS, "$(setting)_$(name).jls"), (; E=field(out), steps=n))
        open(OUT, "a") do io
            @printf(io, "%s,%s,%d,%s,%g,%.2f,%.2f,%d,%.2f\n", MACHINE, setting,
                    Threads.nthreads(), name, SCALE, ts, tr, n, 1e3tr/n)
        end
        @printf("%-8s %-15s setup %7.2f s  run %8.2f s  %5d steps  %8.2f ms/step\n",
                setting, name, ts, tr, n, 1e3tr/n)
        flush(stdout)
    end
end

main()
println("CPU_GPU_DONE")
