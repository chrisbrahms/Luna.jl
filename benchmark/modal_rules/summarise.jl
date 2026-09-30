#= Accuracy and cost of every run of a case against its reference.

       julia --project=benchmark benchmark/modal_rules/summarise.jl [case ...]

   The reference of a case is `REF_<case>` if set (a run file's base name), otherwise the
   fixed run with the most transverse nodes, and of those the one at the tightest
   propagation tolerance. `permode>1e-3`/`energy>1e-3` are the same metrics over the modes
   above a 1e-3 energy share. A second table compares each rule run at several propagation
   tolerances with its own tightest run: the step control's floor. `τ` (FWHM) is
   supplementary: it jumps between sub-pulses of a structured, compressed pulse. Metrics, at every saved z after the first, reported as the maximum over z
   (and `end` at the last z):
   - `global`: max |ΔE| over (ω, mode) / max |E_ref|
   - `permode`: max over modes m of max|ΔE_m| / max|E_ref,m|, over the modes whose energy
     share is above `MINSHARE` (a mode that stays empty by symmetry would divide by noise)
   - `energy`: max over the same modes of |ΔU_m| / U_m, U_m = Σ_ω |E_m|²
   - `dB`: max |10 log₁₀(S/S_ref)| of the mode-summed spectrum where S_ref is above −40 dB
   - `Ppk`, `τ`: relative error of the peak power and FWHM duration of the mode-summed
     power
   plus the attempted steps, the adaptive rule's median and maximum transverse point
   count, the propagation time (`run s`, concurrent with other runs, indicative) and the
   serial time from the timed lane where it ran (`serial s`). =#
using Luna, Printf, Serialization
import Luna: Processing
include(joinpath(@__DIR__, "cases.jl"))

const MINSHARE = 1e-8
const DIR = joinpath(@__DIR__, haskey(ENV, "FIELDS") ? "fields_" * ENV["FIELDS"] : "fields")

function pulse(grid, E)
    (first(Processing.peakpower(grid, E; sumdims=2)), first(Processing.fwhm_t(grid, E; sumdims=2)))
end

function metrics(grid, E, R; minshare=MINSHARE)
    U = [sum(abs2, R[:, m]) for m in axes(R, 2)]
    keepm = findall(u -> u > minshare*sum(U), U)
    g = maximum(abs, E .- R)/maximum(abs, R)
    pm = maximum(maximum(abs, E[:, m] .- R[:, m])/maximum(abs, R[:, m]) for m in keepm)
    en = maximum(abs(sum(abs2, E[:, m]) - U[m])/U[m] for m in keepm)
    S = vec(sum(abs2, E; dims=2)); Sr = vec(sum(abs2, R; dims=2))
    k = Sr .> 1e-4*maximum(Sr)
    db = maximum(abs.(10 .* log10.(S[k] ./ Sr[k])))
    (P, τ), (Pr, τr) = pulse(grid, E), pulse(grid, R)
    (; global_=g, permode=pm, energy=en, dB=db, Ppk=abs(P - Pr)/Pr, τ=abs(τ - τr)/τr)
end

nodes(name) = (m = match(r"^F(\d+)x(\d+)", name); isnothing(m) ? 0 : parse(Int, m[1])*parse(Int, m[2]))
prtol(name) = parse(Float64, match(r"_p([0-9.e-]+)_s", name)[1])

"Serial propagation times from the timed lane (`results/timed.csv`), by (case, run name)."
function timed()
    f = joinpath(@__DIR__, "results", "timed.csv")
    isfile(f) || return Dict{Tuple{String, String}, Float64}()
    rows = split.(readlines(f)[2:end], ",")
    Dict((r[1], "$(r[2])_p$(r[3])_s$(parse(Float64, r[4]))") => parse(Float64, r[6]) for r in rows)
end

function summarise(case)
    files = filter(f -> endswith(f, ".jls"), readdir(joinpath(DIR, case)))
    names = first.(splitext.(files))
    isempty(names) && return
    ref = get(ENV, "REF_" * case, "")
    if isempty(ref)
        # most nodes first, then the tightest propagation tolerance among those
        fixed = filter(n -> startswith(n, "F"), names)
        ref = argmax(n -> (nodes(n), -prtol(n)), fixed)
    end
    R = deserialize(joinpath(DIR, case, ref * ".jls"))
    grid = GRIDS[case]()
    println("## ", case, " (reference ", ref, ")\n")
    T = timed()
    println("| run | steps | points med/max | run s | serial s | global | permode | energy | permode>1e-3 | energy>1e-3 | dB | Ppk | τ | permode end |")
    println("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for n in sort(names; by=n -> (n[1], nodes(n), n))
        r = deserialize(joinpath(DIR, case, n * ".jls"))
        size(r.E) == size(R.E) || (println("| $n | shape differs |"); continue)
        ms = [metrics(grid, r.E[:, :, iz], R.E[:, :, iz]) for iz in 2:size(R.E, 3)]
        ms3 = [metrics(grid, r.E[:, :, iz], R.E[:, :, iz]; minshare=1e-3) for iz in 2:size(R.E, 3)]
        mx(f) = maximum(getfield.(ms, f))
        mx3(f) = maximum(getfield.(ms3, f))
        tp = filter(!isnan, r.transverse_points)
        pts = isempty(tp) ? "" : @sprintf("%d/%d", sort(tp)[cld(length(tp), 2)], maximum(tp))
        ser = haskey(T, (case, n)) ? @sprintf("%.1f", T[(case, n)]) : ""
        @printf("| %s | %d | %s | %.1f | %s | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |\n",
                n, r.steps, pts, r.run_s, ser, mx(:global_), mx(:permode), mx(:energy),
                mx3(:permode), mx3(:energy), mx(:dB), mx(:Ppk), mx(:τ), ms[end].permode)
    end
    println()
    # the propagation floor: each rule run at more than one tolerance, against its own
    # tightest run (the difference is the step control's, not the transverse rule's)
    rules = unique(first.(split.(names, "_p")))
    println("Propagation floor (same rule, against its tightest tolerance):\n")
    println("| run | against | global | permode>1e-3 | dB |")
    println("| --- | --- | ---: | ---: | ---: |")
    for rl in rules
        rs = sort(filter(n -> startswith(n, rl * "_p"), names); by=prtol)
        length(rs) > 1 || continue
        T0 = deserialize(joinpath(DIR, case, rs[1] * ".jls"))
        for n in rs[2:end]
            r = deserialize(joinpath(DIR, case, n * ".jls"))
            ms = [metrics(grid, r.E[:, :, iz], T0.E[:, :, iz]; minshare=1e-3) for iz in 2:size(T0.E, 3)]
            @printf("| %s | %s | %.1e | %.1e | %.1e |\n", n, rs[1], maximum(getfield.(ms, :global_)),
                    maximum(getfield.(ms, :permode)), maximum(getfield.(ms, :dB)))
        end
    end
    println()
end

foreach(summarise, isempty(ARGS) ? filter(d -> isdir(joinpath(DIR, d)), readdir(DIR)) : ARGS)
