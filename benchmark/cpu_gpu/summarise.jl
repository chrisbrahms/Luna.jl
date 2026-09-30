#= Tables from the CPU-versus-GPU runs (`run.jl`): for each machine in `results/`, per case,
   the propagation time of every setting, its speed-up over the baseline -- the
   serial CPU Float64 run, `cpu64_t1_b1` (one Julia thread, one BLAS thread) -- and how far its final field is from the
   baseline's.

       julia --project=benchmark benchmark/cpu_gpu/summarise.jl [MACHINE ...]

   `× base` and `× cpu64_t8` are the baseline's and the threaded CPU Float64 run's
   propagation time divided by this row's (above 1 is faster).

   Accuracy columns, against the baseline run of the same case on the same machine, at
   the last z (`fields/<MACHINE>/*.jls`):
   - `L2`: ‖E − E₆₄‖ / ‖E₆₄‖ over the whole final field
   - `dB`: max |10 log₁₀(S/S₆₄)| of the transverse-summed spectrum S = Σ |E|², where S₆₄
     is above −40 dB of its peak
   - `ΔU_m` (multimode only): max over modes of the relative difference of the mode
     energy Σ_ω |E_m|², however weak the mode
   A case the baseline did not run is compared with the threaded CPU Float64 run.
   For the multimode rows the reference is the baseline run with the same transverse
   rule; how the rules compare with each other is in `benchmark/threads/modal_accuracy.jl`.
   The last row of each machine's table is the sum over the cases every setting ran. =#
using Printf, Serialization

const DIR = @__DIR__

function readcsv(file)
    lines = readlines(file)
    hdr = split(lines[1], ",")
    [Dict(zip(hdr, split(l, ","))) for l in lines[2:end] if !isempty(l)]
end

function metrics(a, r; modal)
    l2 = sqrt(sum(abs2, a.E .- r.E)/sum(abs2, r.E))
    sumdims = Tuple(2:ndims(r.E))
    S = vec(sum(abs2, a.E; dims=sumdims)); Sr = vec(sum(abs2, r.E; dims=sumdims))
    keep = Sr .> 1e-4*maximum(Sr)
    db = maximum(abs.(10 .* log10.(S[keep] ./ Sr[keep])))
    du = modal ? maximum(abs(sum(abs2, a.E[:, m]) - sum(abs2, r.E[:, m])) /
                                   sum(abs2, r.E[:, m]) for m in axes(r.E, 2)) : NaN
    l2, db, du
end

"First of `names` that is among `settings`, or `nothing`."
pick(settings, names) = (i = findfirst(in(settings), names); isnothing(i) ? nothing : names[i])

"The baseline: fully serial CPU Float64 (`cpu64_t1_b1`), else one Julia thread."
baseline(settings) = pick(settings, ["cpu64_t1_b1", "cpu64_t1"])

"The threaded CPU reference: Luna's defaults at 8 Julia threads (one per performance core)."
threaded(settings) = pick(settings, ["cpu64_t8", "cpu64_t10", "cpu64_t4"])

function summarise(machine)
    rows = readcsv(joinpath(DIR, "results", machine * ".csv"))
    #= repeated runs of a (setting, case): the fastest is reported, `spread` is the slowest
       over the fastest; the step count and field are the last run's =#
    times = Dict{Tuple{String, String}, Vector{Float64}}()
    runs = Dict{Tuple{String, String}, Dict}()
    for r in rows
        k = (r["setting"], r["case"]); t = parse(Float64, r["run_s"])
        push!(get!(times, k, Float64[]), t)
        runs[k] = merge(r, Dict("run_s" => @sprintf("%.2f", minimum(times[k])),
                                "ms_per_step" => @sprintf("%.2f", 1e3minimum(times[k])/parse(Int, r["steps"]))))
    end
    spread = Dict(k => maximum(v)/minimum(v) for (k, v) in times)
    settings = unique(r["setting"] for r in rows)
    base = baseline(settings)
    thr = threaded(settings)
    settings = sort(settings; by=s -> (s != base, !startswith(s, "cpu"), s))
    cases = unique(r["case"] for r in rows)
    println("## ", machine, " (baseline ", base, ")\n")
    println("| case | setting | setup s | run s | steps | ms/step | n | spread | × base | × $thr | L2 | dB | ΔU_m |")
    println("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    load(s, c) = (f = joinpath(DIR, "fields", machine, "$(s)_$(c).jls"); isfile(f) ? deserialize(f) : nothing)
    for c in cases, s in settings
        haskey(runs, (s, c)) || continue
        r = runs[(s, c)]
        ratio(ref) = isnothing(ref) ? "" :
            @sprintf("%.2f", parse(Float64, ref["run_s"])/parse(Float64, r["run_s"]))
        sp = ratio(get(runs, (base, c), nothing)), ratio(get(runs, (thr, c), nothing))
        acc = ("", "", "")
        refset = haskey(runs, (base, c)) ? base : thr
        a, b = load(s, c), load(refset, c)
        if s != refset && !isnothing(a) && !isnothing(b)
            acc = map(x -> isnan(x) ? "" : @sprintf("%.1e", x), metrics(a, b; modal=startswith(c, "modal")))
        end
        nrep = count(x -> (x["setting"], x["case"]) == (s, c), rows)
        @printf("| %s | %s | %s | %s | %s | %s | %d | %.2f | %s | %s | %s | %s | %s |\n", c, s,
                r["setup_s"], r["run_s"], r["steps"], r["ms_per_step"], nrep, spread[(s, c)],
                sp..., acc...)
    end
    common = [c for c in cases if all(s -> haskey(runs, (s, c)), settings)]
    if !isempty(common)
        tot = Dict(s => sum(parse(Float64, runs[(s, c)]["run_s"]) for c in common) for s in settings)
        println("\nTotal run time over ", join(common, ", "), ": ",
                join((@sprintf("%s %.1f s", s, tot[s]) for s in settings), "; "))
    end
    println()
end

machines = isempty(ARGS) ?
    [splitext(f)[1] for f in readdir(joinpath(DIR, "results")) if endswith(f, ".csv")] : ARGS
foreach(summarise, machines)
