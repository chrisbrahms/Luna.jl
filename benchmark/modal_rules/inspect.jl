#= Physics check of one run: per saved z, the energy share of each mode, the peak power
   and FWHM duration of the total (mode-summed) power, and the adaptive rule's transverse
   point count over the steps up to that z.

       julia --project=benchmark benchmark/modal_rules/inspect.jl <case> <run file> =#
using Luna, Printf, Serialization
import Luna: Processing
include(joinpath(@__DIR__, "cases.jl"))

function inspect(case, file)
    r = deserialize(file)
    grid = GRIDS[case]()
    nz = size(r.E, 3)
    println(case, "  ", basename(file), "  steps ", r.steps, "  run ", round(r.run_s; digits=1), " s")
    tp = r.transverse_points
    @printf("%7s %9s %8s %9s  %s\n", "z (cm)", "Ppk (GW)", "τ (fs)", "points", "mode energy shares")
    for iz in 1:nz
        E = r.E[:, :, iz]
        U = [sum(abs2, E[:, m]) for m in axes(E, 2)]
        P = first(Processing.peakpower(grid, E; sumdims=2))
        τ = first(Processing.fwhm_t(grid, E; sumdims=2))
        sel = findall(z -> iz == 1 ? z <= r.z[1] : r.z[iz-1] < z <= r.z[iz], r.statz)
        pts = isempty(tp) || isempty(sel) ? "" : @sprintf("%d-%d", minimum(tp[sel]), maximum(tp[sel]))
        @printf("%7.2f %9.2f %8.2f %9s  %s\n", 100r.z[iz], 1e-9P, 1e15τ, pts,
                join((@sprintf("%.1e", u/sum(U)) for u in U), " "))
    end
end

inspect(ARGS[1], ARGS[2])
