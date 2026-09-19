#= Storage and comparison helpers shared by the baseline generator
   (`test/regression/run_cases.jl`), the sensitivity study
   (`test/regression/sensitivity.jl`) and the gate (`test/test_regression.jl`).

   Like `cases.jl`, this file is copied into a worktree of an older commit by
   `generate.jl`, so it must only use API that exists there.
=#
module RegressionCompare

using Luna
import Luna: Output, Utils

export saverun, loadcase, rundict, compare, skipstats, Difference

"""
    rundict(output)

Extract the comparable content of a `MemoryOutput`: the saved field `Eω`, the save
positions `z`, and every statistic. Metadata (source code, git commit, the saved argument
lists) is deliberately left out.
"""
function rundict(output)
    Dict{String, Any}("Eω" => output["Eω"],
                      "z" => output["z"],
                      "stats" => Dict{String, Any}(output["stats"]))
end

"""
    saverun(dir, name, mode, output)

Save the result of case `name` in mode `mode` to `joinpath(dir, name * ".h5")`, under a
group named after the mode.
"""
function saverun(dir, name, mode, output)
    isdir(dir) || mkpath(dir)
    fpath = joinpath(dir, name * ".h5")
    Utils.save_dict_h5(fpath, Dict(string(mode) => rundict(output)); force=true)
    fpath
end

"""
    loadcase(dir, name)

Load every mode of case `name` from `joinpath(dir, name * ".h5")`.
"""
function loadcase(dir, name)
    fpath = joinpath(dir, name * ".h5")
    isfile(fpath) || error("No baseline for case $name at $fpath")
    Utils.load_dict_h5(fpath)
end

"""
    Difference(what, value, note)

One comparison result: `what` names the quantity ("Eω", "z" or a statistic), `value` is
the normalised maximum difference and `note` is either an empty string or an explanation
of why the comparison could not be made (in which case `value` is `Inf`).
"""
struct Difference
    what::String
    value::Float64
    note::String
end

Difference(what, value) = Difference(what, value, "")

#= The metric is `maximum(abs, Δ) / maximum(abs, reference)`: elementwise relative
   differences are meaningless in the window tapers and outside `grid.sidx`, where the
   field is many orders of magnitude below its peak. NaNs (which some statistics, e.g. the
   zero-dispersion wavelength, legitimately produce) compare equal to NaNs in the same
   places and unequal to anything else. =#
"""
    normdiff(ref, new)

Normalised maximum difference between `ref` and `new`: `maximum(abs, new - ref)` divided
by `maximum(abs, ref)`. Returns `(value, note)`; `value` is `Inf` and `note` says why when
the two cannot be compared.
"""
function normdiff(ref, new)
    size(ref) == size(new) || return (Inf, "size $(size(new)) != baseline $(size(ref))")
    nanr = isnan.(ref)
    nann = isnan.(new)
    nanr == nann || return (Inf, "NaN pattern differs")
    ok = .!nanr
    any(ok) || return (0.0, "")
    num = maximum(abs, new[ok] .- ref[ok])
    den = maximum(abs, ref[ok])
    if den == 0
        return (num == 0 ? 0.0 : Inf, num == 0 ? "" : "baseline is all zero")
    end
    (num/den, "")
end

normdiff(ref::Number, new::Number) = normdiff([ref], [new])

"""
    STEP_STATS

The statistics that record the step sequence rather than the field: `dz`, the size of each
accepted step, and `z`, its running sum. `Stats.zdz!` writes both at every accepted step.

They are not `rundict`'s top-level `"z"`, which is the save grid and is fixed by
`Output.GridCondition` — that one is compared in every mode and must be exact.
"""
const STEP_STATS = ("z", "dz")

"""
    skipstats(mode)

The statistics left out of the comparison in run `mode`.

In the `:fixed` mode nothing is left out: the step sequence is imposed, so `z` and `dz` must
match exactly and are a useful check that they were in fact imposed.

In the `:adaptive` mode [`STEP_STATS`](@ref) is left out. The step-size controller's
accept/reject decision and its PI update respond to a change in the last bits of the error
estimate far more strongly than the field does, and the response compounds over the steps it
takes to ramp `init_dz` up to `max_dz`; including `dz` would set the case's tolerance three
to six orders of magnitude above what `Eω` needs and make the adaptive run no constraint at
all. The field and the physical statistics are still compared, so a real change still shows
up — through its effect on the result rather than on the bookkeeping.
"""
skipstats(mode::Symbol) = mode === :adaptive ? STEP_STATS : ()

"""
    compare(ref::AbstractDict, new::AbstractDict; skip=())

Compare two run dictionaries (as produced by [`rundict`](@ref) or read back from HDF5) and
return a `Vector{Difference}`, one entry for `Eω`, one for `z` and one per statistic.

`skip` names statistics to leave out; pass [`skipstats(mode)`](@ref).
"""
function compare(ref::AbstractDict, new::AbstractDict; skip=())
    out = Difference[]
    for key in ("Eω", "z")
        v, note = normdiff(ref[key], new[key])
        push!(out, Difference(key, v, note))
    end
    rstats = ref["stats"]
    nstats = new["stats"]
    for key in sort(collect(keys(rstats)))
        key in skip && continue
        if !haskey(nstats, key)
            push!(out, Difference("stats/" * key, Inf, "missing from the new run"))
            continue
        end
        v, note = normdiff(rstats[key], nstats[key])
        push!(out, Difference("stats/" * key, v, note))
    end
    for key in sort(collect(keys(nstats)))
        key in skip && continue
        haskey(rstats, key) || push!(out, Difference("stats/" * key, Inf,
                                                     "not in the baseline"))
    end
    out
end

"""
    worst(diffs)

The largest difference in `diffs`.
"""
worst(diffs::AbstractVector{Difference}) = isempty(diffs) ? Difference("", 0.0) :
    diffs[argmax([d.value for d in diffs])]

end
