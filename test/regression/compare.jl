#= Storage and comparison helpers shared by the baseline generator
   (`test/regression/run_cases.jl`), the sensitivity study
   (`test/regression/sensitivity.jl`) and the gate (`test/test_regression.jl`).

   Like `cases.jl`, this file is copied into a worktree of an older commit by
   `generate.jl`, so it must only use API that exists there.
=#
module RegressionCompare

using Luna
import Luna: Output, Utils

export saverun, loadcase, rundict, compare, skipstats, classof, Difference

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
    Difference(what, class, value, where, note)

One comparison result.

- `what` names the quantity: `"Eω"`, `"z"`, `"stats/<name>"` or `"step count"`.
- `class` is the tolerance class the quantity belongs to, `:Eω` or `:stats`; see
  [`classof`](@ref).
- `value` is the normalised maximum difference.
- `where` locates it — for `Eω`, the component and save at which the maximum was reached —
  and is empty when there is nothing to point at.
- `note` is empty, or says why the comparison could not be made, in which case `value` is
  `Inf`.
"""
struct Difference
    what::String
    class::Symbol
    value::Float64
    where::String
    note::String
end

Difference(what, class, value) = Difference(what, class, value, "", "")

"""
    classof(what)

The tolerance class of a quantity: `:Eω` for the field itself, `:stats` for the save grid
`z`, the statistics and the step-count check.

Two classes rather than one tolerance for everything: `Eω` is the quantity every later
branch must not move, and sharing a tolerance with the statistics gave it several orders of
magnitude of slack in the adaptive mode, because the statistics are recorded once per step
and the step-size controller responds to a change in the last bits far more strongly than
the field does.
"""
classof(what::AbstractString) = what == "Eω" ? :Eω : :stats

"The two tolerance classes, in the order they are reported."
const CLASSES = (:Eω, :stats)

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
to six orders of magnitude above what the field needs.

The *number* of accepted steps is checked in both modes regardless — see [`compare`](@ref).
Excluding `dz` removes the graded part of the step sequence's sensitivity, not all of it.
"""
skipstats(mode::Symbol) = mode === :adaptive ? STEP_STATS : ()

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

Used for the statistics, which are one number — or one number per mode — per step. `Eω`
uses [`fielddiff`](@ref) instead.
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
    fieldaxes(n)

The axes of an `n`-dimensional `Eω` array that [`fielddiff`](@ref) reduces over: the
frequency axis, and the transverse axes when there are any.

`Eω` is `(Nω, Nz)` mode-averaged, `(Nω, Nm, Nz)` multimode or vector, `(Nω, Npol, Nk, Nz)`
in 2-D free space and `(Nω, Npol, Nkx, Nky, Nz)` in 3-D. Axis 1 is always frequency and the
last axis is always the save; axis 2 is the mode or polarisation whenever there are more
than two axes, and axes 3 … n-1 are transverse.
"""
fieldaxes(n::Integer) = n >= 4 ? (1, (3:n-1)...) : (1,)

"""
    whereidx(idx, sz)

Human-readable form of an index into the reduced `(component, save)` array, or `(save,)`
when the field has no component axis.
"""
function whereidx(idx::CartesianIndex, sz)
    if length(sz) == 1
        "save $(idx[1])/$(sz[1])"
    else
        "component $(idx[1])/$(sz[1]), save $(idx[2])/$(sz[2])"
    end
end

"""
    fielddiff(ref, new)

Normalised maximum difference between two `Eω` arrays, normalised **per component and per
save**: the frequency and transverse axes are reduced away and the ratio
`maximum(abs, Δ)/maximum(abs, ref)` is formed separately for every (mode or polarisation,
save) pair. Returns `(value, where, note)`, where `value` is the largest such ratio and
`where` says which pair produced it.

One global normalisation, as in [`normdiff`](@ref), hides weak components: in
`multimode_field_plasma` the strongest mode is 4e6 times the weakest, so a change that
rewrote the weakest mode entirely would sit far below any tolerance the strongest mode makes
sensible. GPU_PLAN.md §4.11 also specifies the metric per save.

A slice whose baseline is identically zero is skipped when the difference is zero there too,
and reported as `Inf` otherwise.
"""
function fielddiff(ref, new)
    size(ref) == size(new) || return (Inf, "", "size $(size(new)) != baseline $(size(ref))")
    nanr = isnan.(ref)
    nanr == isnan.(new) || return (Inf, "", "NaN pattern differs")
    Δ = abs.(new .- ref)
    R = abs.(ref)
    if any(nanr)
        Δ[nanr] .= 0
        R[nanr] .= 0
    end
    ax = fieldaxes(ndims(ref))
    num = dropdims(maximum(Δ; dims=ax); dims=ax)
    den = dropdims(maximum(R; dims=ax); dims=ax)
    best = 0.0
    bestidx = nothing
    for idx in CartesianIndices(num)
        if den[idx] == 0
            num[idx] == 0 && continue # nothing there in either run
            return (Inf, whereidx(idx, size(num)), "baseline slice is all zero")
        end
        r = num[idx]/den[idx]
        if r > best
            best = r
            bestidx = idx
        end
    end
    (best, isnothing(bestidx) ? "" : whereidx(bestidx, size(num)), "")
end

"""
    stepcount(d)

The number of accepted steps recorded in a run dictionary, or `nothing` when the case
records no statistics.
"""
stepcount(d::AbstractDict) = haskey(d["stats"], "z") ? length(d["stats"]["z"]) : nothing

"""
    compare(ref::AbstractDict, new::AbstractDict; skip=())

Compare two run dictionaries (as produced by [`rundict`](@ref) or read back from HDF5) and
return a `Vector{Difference}`: one entry for `Eω`, one for the save grid `z`, and one per
compared statistic.

`skip` names statistics to leave out; pass [`skipstats(mode)`](@ref).

If the number of accepted steps differs, every statistic would fail its size check with an
uninformative `Inf`. Instead this returns `Eω`, `z` and a single `"step count"` entry saying
what the two counts were, and compares no statistics. That is a hard failure in both modes:
a different number of steps is a different propagation, whether or not the controller was
free to choose it.
"""
function compare(ref::AbstractDict, new::AbstractDict; skip=())
    out = Difference[]
    v, w, note = fielddiff(ref["Eω"], new["Eω"])
    push!(out, Difference("Eω", classof("Eω"), v, w, note))
    v, note = normdiff(ref["z"], new["z"])
    push!(out, Difference("z", classof("z"), v, "", note))

    rsteps, nsteps = stepcount(ref), stepcount(new)
    if !isnothing(rsteps) && !isnothing(nsteps) && rsteps != nsteps
        push!(out, Difference("step count", :stats, Inf, "",
                              "changed: $nsteps steps, baseline $rsteps"))
        return out
    end

    rstats = ref["stats"]
    nstats = new["stats"]
    for key in sort(collect(keys(rstats)))
        key in skip && continue
        what = "stats/" * key
        if !haskey(nstats, key)
            push!(out, Difference(what, :stats, Inf, "", "missing from the new run"))
            continue
        end
        v, note = normdiff(rstats[key], nstats[key])
        push!(out, Difference(what, :stats, v, "", note))
    end
    for key in sort(collect(keys(nstats)))
        key in skip && continue
        haskey(rstats, key) || push!(out, Difference("stats/" * key, :stats, Inf, "",
                                                     "not in the baseline"))
    end
    out
end

"""
    worst(diffs)
    worst(diffs, class)

The largest difference in `diffs`, optionally restricted to one tolerance `class`. Returns a
placeholder named `"-"` when every difference is exactly zero: `argmax` over an all-zero
vector returns the first entry, which reads like a ranking and is not one.
"""
function worst(diffs::AbstractVector{Difference})
    isempty(diffs) && return Difference("-", :stats, 0.0)
    i = argmax([d.value for d in diffs])
    diffs[i].value == 0 ? Difference("-", :stats, 0.0) : diffs[i]
end

worst(diffs::AbstractVector{Difference}, class::Symbol) =
    worst([d for d in diffs if d.class === class])

"""
    describe(d::Difference)

`what`, with the location and the note appended when there are any.
"""
function describe(d::Difference)
    s = d.what
    isempty(d.where) || (s *= " [" * d.where * "]")
    isempty(d.note) || (s *= " (" * d.note * ")")
    s
end

end
