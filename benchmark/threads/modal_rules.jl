#= Multimode on the CPU: the adaptive transverse integral against the fixed rule, at
   matched accuracy, against Julia and BLAS threads. Appends to results/modal_rules.csv.

       for J in 1 8; do julia --project=. -t $J benchmark/threads/modal_rules.jl; done

   The case has to make the transverse integral hard, or every rule is trivially exact:
   argon at 1 bar in a 125 µm capillary, 8 HE₁ₘ modes, 30 fs at 800 nm, 30 cm, with the
   peak power at 0.8 and 0.95 of the critical power for self-focusing (`Tools.Pcr`,
   11.3 GW). There the beam self-focuses inside the capillary until the plasma arrests
   it; at the default tolerance the adaptive rule refines from 31 to 127-255 points and
   27-37 % of the energy leaves HE₁₁ (at 0.5 P_cr it stays at 31-63 points and HE₁₁ keeps
   98.7 %). Accuracy: final Eω against `:fixed` with `nr=256`, as the largest difference
   over all modes divided by the peak over all modes. =#
import LinearAlgebra
import LinearAlgebra: BLAS
using Luna, Printf, Logging, Statistics
Luna.set_fftw_mode(:estimate); Luna.set_fftw_wisdom(false); Luna.set_fftw_threads(1)
# BLAS counts are set explicitly below; stop Luna.run from choosing its own
setblas(n) = (isdefined(Luna, :set_blas_threads) && Luna.set_blas_threads(n); BLAS.set_num_threads(n))

const J = Threads.nthreads()
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results", "modal_rules.csv"))
const FRACS = parse.(Float64, split(get(ENV, "FRACS", "0.8,0.95"), ","))
const λ0, τ, NM, L = 800e-9, 30e-15, 8, 0.3

function pcrit()
    ω = 2π*PhysData.c/λ0
    _, n0, n2 = Tools.getN0n0n2(ω, :Ar; P=1.0)
    Tools.Pcr(ω, n0, n2)
end

run(E, rule; kw...) = with_logger(ConsoleLogger(stderr, Logging.Warn)) do
    prop_capillary(125e-6, L, :Ar, 1.0; λ0, τfwhm=τ, energy=E, λlims=(150e-9, 4e-6),
                   trange=0.5e-12, shotnoise=false, modes=NM, modal_integral=rule,
                   status_period=1e9, kw...)
end

relerr(E, R) = maximum(abs, E[:, :, end] .- R[:, :, end])/maximum(abs, R[:, :, end])

const RULES = [(:adaptive, 0, 1e-3), (:adaptive, 0, 1e-4),
               (:fixed, 32, 0.0), (:fixed, 64, 0.0), (:fixed, 128, 0.0)]

isfile(OUT) || open(io -> println(io,
    "julia_threads,blas_threads,p_over_pcr,rule,nr,rtol,time_s,relerr,points_median,points_max"), OUT, "w")
Pcr = pcrit()
for frac in FRACS
    E = frac*Pcr*τ/0.94 # peak power of a Gaussian: 0.94 E/τ
    setblas(1)
    ref = run(E, :fixed; modal_nr=256)["Eω"]
    first = true
    for (rule, nr, rtol) in RULES, nb in (1, 8)
        setblas(nb)
        kw = rule === :fixed ? (modal_nr=nr,) : (radial_integral_rtol=rtol,)
        first && (run(E, rule; kw...); first = false) # compile
        t = @elapsed out = run(E, rule; kw...)
        tp = out["stats"]["transverse_points"]
        e = relerr(out["Eω"], ref)
        open(OUT, "a") do io
            @printf(io, "%d,%d,%.2f,%s,%d,%.0e,%.2f,%.2e,%d,%d\n", J, nb, frac, rule, nr,
                    rtol, t, e, round(Int, median(tp)), maximum(tp))
        end
        @printf("J=%d blas=%d P/Pcr=%.2f %-8s nr=%3d rtol=%.0e  %7.2f s  err %.1e  points %d/%d\n",
                J, nb, frac, rule, nr, rtol, t, e, round(Int, median(tp)), maximum(tp))
        flush(stdout)
    end
end
println("MODAL_DONE J=$J")
