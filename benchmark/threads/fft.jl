#= FFTW execution and planning time for Luna-shaped transforms against the FFTW thread
   count, in one Julia process (so for one Julia thread count). Appends to
   benchmark/threads/results/fft.csv.

       for J in 1 2 4 8 10; do
           julia --project=. -t $J benchmark/threads/fft.jl
       done

   With `Threads.nthreads() > 1`, FFTW.jl runs FFTW's threads as Julia tasks
   (`Threads.@spawn`); with one Julia thread FFTW uses its own pthreads. Both are measured.
   `FFTW.set_num_threads` applies to plans made after it, so every setting replans.
   Set QUICK=1 for a short smoke run. =#
import FFTW, LinearAlgebra
using Printf

const J = Threads.nthreads()
const QUICK = get(ENV, "QUICK", "0") == "1"
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results", "fft.csv"))
const FLAGS = QUICK ? (FFTW.ESTIMATE,) : (FFTW.ESTIMATE, FFTW.PATIENT)
const FLAGNAME = Dict(FFTW.ESTIMATE => "estimate", FFTW.PATIENT => "patient")

thread_settings() = sort(unique(filter(n -> n >= 1,
    J == 1 ? [1, 2, 4, 8] : [1, 2, 4, 8, J, 2J, 4J])))

# (kind, dims, region): kind is :r (real forward) or :c (complex forward)
function shapes()
    s = Tuple{Symbol, Tuple, Any}[]
    for p in (QUICK ? (12, 16) : 10:18)
        push!(s, (:r, (2^p,), 1))
        push!(s, (:c, (2^p,), 1))
    end
    for nt in (QUICK ? (2^13,) : (2^12, 2^14)), nc in (QUICK ? (32,) : (2, 8, 32, 256, 1024))
        push!(s, (:r, (nt, nc), 1))
    end
    for d in (QUICK ? ((2^10, 32, 32),) : ((2^10, 32, 32), (2^10, 64, 64), (2^12, 64, 64)))
        push!(s, (:c, d, (1, 2, 3)))
    end
    s
end

# mean time per call over batches of at least 0.1 s, minimum over 3 batches
function exectime(f)
    f()
    t1 = @elapsed f()
    reps = max(1, min(10_000, ceil(Int, 0.1/max(t1, 1e-9))))
    best = Inf
    for _ in 1:3
        t = @elapsed for _ in 1:reps; f(); end
        best = min(best, t/reps)
    end
    best
end

function run()
    isfile(OUT) || open(io -> println(io, "julia_threads,fftw_threads,flag,kind,dims,elements,plan_s,exec_s"), OUT, "w")
    for flag in FLAGS, (kind, dims, region) in shapes()
        n = prod(dims)
        # PATIENT on the largest shapes takes minutes per plan; keep it to <= 2^20 elements
        flag == FFTW.PATIENT && n > 2^20 && continue
        x = kind === :r ? rand(dims...) : rand(ComplexF64, dims...)
        for nf in thread_settings()
            FFTW.set_num_threads(nf)
            tp = @elapsed p = kind === :r ? FFTW.plan_rfft(x, region; flags=flag) :
                                            FFTW.plan_fft(x, region; flags=flag)
            y = p * x
            te = exectime(() -> LinearAlgebra.mul!(y, p, x))
            open(OUT, "a") do io
                @printf(io, "%d,%d,%s,%s,%s,%d,%.6e,%.6e\n", J, nf, FLAGNAME[flag], kind,
                        join(dims, "x"), n, tp, te)
            end
        end
        FFTW.forget_wisdom() # each plan measured from scratch
        @printf("J=%d %-8s %s %-14s done\n", J, FLAGNAME[flag], kind, join(dims, "x")); flush(stdout)
    end
end

run()
println("FFT_DONE J=$J")
