#= Runs of the fixed-versus-adaptive study, serial on the CPU: one Julia thread, one BLAS
   thread, one FFTW thread, `:measure` planning, no wisdom.

       RUNS="vortex/A1e-3m512/1e-6 vortex/F64x16/1e-6" \
           julia --project=benchmark -t 1 benchmark/modal_rules/run.jl

   A run is `case/rule/prtol` (`cases.jl`: `parserule`); `SCALE` (default 1) multiplies
   every length. Each run saves the field at the `NSAVE` saved z, the statistics' z and
   `transverse_points`, the attempted step count and the setup and propagation times to
   `fields/<case>/<rule>_p<prtol>_s<scale>.jls`, and appends a row to `results/runs.csv`
   (`results/<OUT>.csv` with `OUT`; `FIELDS=x` puts the files in `fields_x/`). A run whose file exists is skipped, so an interrupted
   queue can be restarted as it is. =#
import LinearAlgebra
using Luna, Printf, Serialization
Luna.set_fftw_wisdom(false)
Luna.set_fftw_threads(1)
Luna.set_blas_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Threads.nthreads() == 1 || @warn "run.jl is meant for serial runs (-t 1)"

include(joinpath(@__DIR__, "cases.jl"))

const SCALE = parse(Float64, get(ENV, "SCALE", "1"))
const CSV = joinpath(@__DIR__, "results", get(ENV, "OUT", "runs") * ".csv")
"Where the run files go: `fields/`, or `fields_<FIELDS>/` (the probes use their own)."
const FDIR = joinpath(@__DIR__, haskey(ENV, "FIELDS") ? "fields_" * ENV["FIELDS"] : "fields")

runfile(case, rule, prtol) =
    joinpath(FDIR, case, "$(rule)_p$(prtol)_s$(SCALE).jls")

function dorun(spec)
    case, rule, prtol = split(spec, "/")
    file = runfile(case, rule, prtol)
    isfile(file) && (println("skip ", spec); flush(stdout); return)
    mkpath(dirname(file)); mkpath(dirname(CSV))
    f = CASES[case]
    r = parserule(rule)
    p = parse(Float64, prtol)
    f(r; scale=0.01SCALE, prtol=p) # compile
    GC.gc()
    res = f(r; scale=SCALE, prtol=p)
    st = res.output.data["stats"]
    tp = haskey(st, "transverse_points") ? Array(st["transverse_points"]) : Float64[]
    serialize(file, (; E=Array(res.output.data["Eω"]), z=Array(res.output.data["z"]),
                       statz=Array(st["z"]), transverse_points=tp, steps=res.steps,
                       setup_s=res.setup_s, run_s=res.run_s))
    isfile(CSV) || open(io -> println(io,
        "case,rule,prtol,scale,setup_s,run_s,steps,points_median,points_max"), CSV, "w")
    tpv = filter(!isnan, tp)
    med = isempty(tpv) ? NaN : sort(tpv)[cld(length(tpv), 2)]
    mx = isempty(tpv) ? NaN : maximum(tpv)
    open(CSV, "a") do io
        @printf(io, "%s,%s,%s,%g,%.2f,%.2f,%d,%g,%g\n", case, rule, prtol, SCALE,
                res.setup_s, res.run_s, res.steps, med, mx)
    end
    @printf("%-12s %-12s p=%s  setup %6.1f s  run %8.1f s  %6d steps  points %g / %g\n",
            case, rule, prtol, res.setup_s, res.run_s, res.steps, med, mx)
    flush(stdout)
end

for spec in split(get(ENV, "RUNS", ""))
    try
        dorun(spec)
    catch err
        @error "run $spec failed" exception=(err, catch_backtrace())
    end
end
println("MODAL_RULES_DONE")
