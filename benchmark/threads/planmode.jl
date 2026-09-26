#= Execution time of FFTW plans made with ESTIMATE, MEASURE and PATIENT, one FFTW thread,
   for the batched and multi-dimensional shapes where fft.jl found PATIENT plans slower.
   Each plan is made on a copy of the input (MEASURE/PATIENT planning overwrites the
   array it is given) and executed on the original data. Appends to
   results/planmode.csv; run it in several processes to see the spread:

       for i in 1 2 3; do julia --project=. -t 1 benchmark/threads/planmode.jl; done =#
import FFTW, LinearAlgebra
using Printf

const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results", "planmode.csv"))
FFTW.set_num_threads(1)

function exectime(f)
    f()
    t1 = @elapsed f()
    reps = max(1, min(10_000, ceil(Int, 0.1/max(t1, 1e-9))))
    best = Inf
    for _ in 1:5
        t = @elapsed for _ in 1:reps; f(); end
        best = min(best, t/reps)
    end
    best
end

const SHAPES = [(:c, (1024, 32, 32), (1, 2, 3)), (:r, (4096, 256), 1), (:r, (16384, 32), 1),
                (:r, (131072,), 1), (:c, (16384,), 1)]
const FLAGS = [("estimate", FFTW.ESTIMATE), ("measure", FFTW.MEASURE), ("patient", FFTW.PATIENT)]

isfile(OUT) || open(io -> println(io, "pid,kind,dims,flag,plan_s,exec_s"), OUT, "w")
for (kind, dims, region) in SHAPES, (fname, flag) in FLAGS
    x = kind === :r ? rand(dims...) : rand(ComplexF64, dims...)
    xc = copy(x)
    FFTW.forget_wisdom()
    tp = @elapsed p = kind === :r ? FFTW.plan_rfft(xc, region; flags=flag) :
                                    FFTW.plan_fft(xc, region; flags=flag)
    y = p * x
    te = exectime(() -> LinearAlgebra.mul!(y, p, x))
    open(OUT, "a") do io
        @printf(io, "%d,%s,%s,%s,%.3e,%.4e\n", getpid(), kind, join(dims, "x"), fname, tp, te)
    end
end
println("PLANMODE_DONE")
