#= Markdown tables from benchmark/threads/results/*.csv.

       julia --project=. benchmark/threads/summarise.jl [fft|blas|runs|runs_metal] =#
using DelimitedFiles, Printf

const RES = joinpath(@__DIR__, "results")

function load(name)
    d, h = readdlm(joinpath(RES, name * ".csv"), ','; header=true)
    cols = Dict(string(h[i]) => d[:, i] for i in eachindex(h))
    [Dict(k => v[r] for (k, v) in cols) for r in 1:size(d, 1)]
end

# FFT: execution time relative to one FFTW thread, per Julia thread count and shape
function fft()
    rows = load("fft")
    for flag in ("estimate", "patient")
        R = filter(r -> r["flag"] == flag, rows)
        isempty(R) && continue
        println("\n### FFT execution, `$flag` plans: t(best) / t(1 FFTW thread) and best count\n")
        Js = sort(unique(r["julia_threads"] for r in R))
        shapes = unique((r["kind"], string(r["dims"]), r["elements"]) for r in R)
        println("| kind | dims | t(1 thread) | ", join(("J=$j" for j in Js), " | "), " |")
        println("|---|---|---:|", repeat("---:|", length(Js)))
        for (k, dims, n) in shapes
            S = filter(r -> r["kind"] == k && string(r["dims"]) == dims, R)
            base = filter(r -> r["fftw_threads"] == 1 && r["julia_threads"] == 1, S)
            t1 = isempty(base) ? NaN : base[1]["exec_s"]
            cells = String[]
            for j in Js
                Sj = filter(r -> r["julia_threads"] == j, S)
                isempty(Sj) && (push!(cells, "–"); continue)
                ref = filter(r -> r["fftw_threads"] == 1, Sj)[1]["exec_s"]
                b = argmin(r -> r["exec_s"], Sj)
                push!(cells, @sprintf("%.2f (%d)", b["exec_s"]/ref, b["fftw_threads"]))
            end
            @printf("| %s | %s | %.2e s | %s |\n", k, dims, t1, join(cells, " | "))
        end
    end
    # planning time, patient
    P = filter(r -> r["flag"] == "patient", rows)
    isempty(P) && return
    println("\n### PATIENT planning time (s), by FFTW threads, J = max\n")
    jm = maximum(r["julia_threads"] for r in P)
    Pm = filter(r -> r["julia_threads"] == jm, P)
    nfs = sort(unique(r["fftw_threads"] for r in Pm))
    println("| kind | dims | ", join(("$n" for n in nfs), " | "), " |")
    println("|---|---|", repeat("---:|", length(nfs)))
    for (k, dims) in unique((r["kind"], string(r["dims"])) for r in Pm)
        S = filter(r -> r["kind"] == k && string(r["dims"]) == dims, Pm)
        cells = [(x = filter(r -> r["fftw_threads"] == n, S); isempty(x) ? "–" : @sprintf("%.2f", x[1]["plan_s"])) for n in nfs]
        println("| $k | $dims | ", join(cells, " | "), " |")
    end
end

function blas()
    rows = load("blas")
    println("\n### GEMM: time relative to 1 BLAS thread (best count), per Julia threads\n")
    Js = sort(unique(r["julia_threads"] for r in rows))
    println("| case | m×k×n | t(1,1) | ", join(("J=$j" for j in Js), " | "), " |")
    println("|---|---|---:|", repeat("---:|", length(Js)))
    for key in unique((r["case"], r["m"], r["k"], r["n"]) for r in rows)
        S = filter(r -> (r["case"], r["m"], r["k"], r["n"]) == key, rows)
        t11 = filter(r -> r["julia_threads"] == 1 && r["blas_threads"] == 1, S)[1]["exec_s"]
        cells = String[]
        for j in Js
            Sj = filter(r -> r["julia_threads"] == j, S)
            ref = filter(r -> r["blas_threads"] == 1, Sj)[1]["exec_s"]
            b = argmin(r -> r["exec_s"], Sj)
            push!(cells, @sprintf("%.2f (%d)", b["exec_s"]/ref, b["blas_threads"]))
        end
        @printf("| %s | %d×%d×%d | %.2e s | %s |\n", key..., t11, join(cells, " | "))
    end
end

function runs(name="runs")
    rows = load(name)
    for flag in unique(r["flag"] for r in rows)
        R = filter(r -> r["flag"] == flag, rows)
        println("\n### Whole runs ($name, `$flag`): wall time (s) per (Julia, FFTW, BLAS) threads\n")
        cfgs = sort(unique((r["julia_threads"], r["fftw_threads"], r["blas_threads"]) for r in R))
        cases = unique(r["case"] for r in R)
        for c in cases
            S = filter(r -> r["case"] == c, R)
            b = argmin(r -> r["time_s"], S)
            cells = [@sprintf("J%d/F%d/B%d %.2f", r["julia_threads"], r["fftw_threads"], r["blas_threads"], r["time_s"]) for r in sort(S, by=r -> (r["julia_threads"], r["fftw_threads"], r["blas_threads"]))]
            maxd = maximum(r["maxreldiff"] for r in S)
            @printf("- **%s**: best J%d/F%d/B%d %.2f s; max diff vs first config %.1e\n  %s\n",
                    c, b["julia_threads"], b["fftw_threads"], b["blas_threads"], b["time_s"], maxd,
                    join(cells, ", "))
        end
    end
end

what = isempty(ARGS) ? ["fft", "blas", "runs"] : ARGS
for w in what
    isfile(joinpath(RES, w * ".csv")) || continue
    w == "fft" ? fft() : w == "blas" ? blas() : runs(w)
end
