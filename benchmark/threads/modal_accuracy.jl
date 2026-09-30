#= Accuracy of the multimode transverse rules at 0.95 P_cr (the modal_rules.jl case),
   with metrics that see weak modes and spectral wings, against two independent
   references. Writes results/modal_accuracy.csv.

       julia --project=. -t 8 benchmark/threads/modal_accuracy.jl

   References, both at a propagation rtol of 1e-7 (10× tighter than Luna.run's default):
   `:fixed` nr = 256 ("refF") and `:adaptive` radial_integral_rtol = 1e-5 ("refA").
   Metrics for each run against each reference, at the last z:
   - global: max |ΔE| over (ω, mode) / max |E_ref| over (ω, mode)
   - permode: max over modes m of  max|ΔE_m| / max|E_ref,m|  (every mode, however weak)
   - energy: max over modes of |ΔU_m| / U_m, U_m = Σ_ω |E_m|²
   - dB: max |10 log10(S/S_ref)| of the total spectrum S = Σ_m |E_m|², over the region
     where S_ref is above -40 dB of its peak
   plus the accepted step count (length of stats["z"]) and the wall time (BLAS 1). =#
import LinearAlgebra
import LinearAlgebra: BLAS
using Luna, Printf, Logging, Statistics
Luna.set_fftw_mode(:estimate); Luna.set_fftw_wisdom(false); Luna.set_fftw_threads(1)
# BLAS counts are set explicitly below; stop Luna.run from choosing its own
setblas(n) = (isdefined(Luna, :set_blas_threads) && Luna.set_blas_threads(n); BLAS.set_num_threads(n))

setblas(1)

const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results", "modal_accuracy.csv"))
const λ0, τ, NM, L, FRAC = 800e-9, 30e-15, 8, 0.3, 0.95

function pcrit()
    ω = 2π*PhysData.c/λ0
    _, n0, n2 = Tools.getN0n0n2(ω, :Ar; P=1.0)
    Tools.Pcr(ω, n0, n2)
end
const E = FRAC*pcrit()*τ/0.94

function run(rule; rtol=1e-6, kw...)
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        Eω, grid, linop, transform, FT, output = Luna.Interface.prop_capillary_args(
            125e-6, L, :Ar, 1.0; λ0, τfwhm=τ, energy=E, λlims=(150e-9, 4e-6),
            trange=0.5e-12, shotnoise=false, modes=NM, modal_integral=rule, kw...)
        Luna.run(Eω, grid, linop, transform, FT, output; zmax=L, rtol, status_period=1e9)
        output
    end
end

last(out) = out["Eω"][:, :, end]
nsteps(out) = length(out["stats"]["z"])

function metrics(E, R)
    g = maximum(abs, E .- R)/maximum(abs, R)
    pm = maximum(maximum(abs, E[:, m] .- R[:, m])/maximum(abs, R[:, m]) for m in axes(R, 2))
    U(X, m) = sum(abs2, X[:, m])
    en = maximum(abs(U(E, m) - U(R, m))/U(R, m) for m in axes(R, 2))
    S = vec(sum(abs2, E; dims=2)); Sr = vec(sum(abs2, R; dims=2))
    keep = Sr .> 1e-4*maximum(Sr)
    db = maximum(abs.(10 .* log10.(S[keep] ./ Sr[keep])))
    (g, pm, en, db)
end

println("references..."); flush(stdout)
tF = @elapsed refF = run(:fixed; rtol=1e-7, modal_nr=256)
tA = @elapsed refA = run(:adaptive; rtol=1e-7, radial_integral_rtol=1e-5)
RF, RA = last(refF), last(refA)
cases = [("fixed nr=256 rtol=1e-7 (refF)", refF, tF), ("adaptive 1e-5 rtol=1e-7 (refA)", refA, tA)]
for (name, rule, kw) in (("adaptive 1e-3", :adaptive, (radial_integral_rtol=1e-3,)),
                         ("adaptive 1e-4", :adaptive, (radial_integral_rtol=1e-4,)),
                         ("fixed nr=32", :fixed, (modal_nr=32,)),
                         ("fixed nr=64", :fixed, (modal_nr=64,)),
                         ("fixed nr=128", :fixed, (modal_nr=128,)),
                         ("fixed nr=256", :fixed, (modal_nr=256,)),
                         ("fixed nr=64 rtol=1e-7", :fixed, (modal_nr=64, rtol=1e-7)))
    t = @elapsed out = run(rule; kw...)
    push!(cases, (name, out, t))
    println("  done ", name); flush(stdout)
end
open(OUT, "w") do io
    println(io, "run,steps,time_s,ref,global,permode,energy,dB")
    for (name, out, t) in cases, (rname, R) in (("refF", RF), ("refA", RA))
        m = metrics(last(out), R)
        @printf(io, "%s,%d,%.1f,%s,%.2e,%.2e,%.2e,%.2e\n", name, nsteps(out), t, rname, m...)
    end
end
# mode energy shares of the reference, to read the per-mode metric against
U = [sum(abs2, RF[:, m]) for m in axes(RF, 2)]
println("mode energy shares (refF): ", join((@sprintf("%.1e", u/sum(U)) for u in U), " "))
println(read(OUT, String))
println("ACCURACY_DONE")
