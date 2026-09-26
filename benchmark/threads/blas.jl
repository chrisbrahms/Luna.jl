#= BLAS (OpenBLAS) GEMM time for the matrix products Luna's transforms do, against the
   BLAS thread count, in one Julia process. Appends to benchmark/threads/results/blas.csv.

       for J in 1 4 8; do julia --project=. -t $J benchmark/threads/blas.jl; done

   - radial Hankel step: `mul!(out, reshape(E, nt*npol, nr), T)`, `T` nr × nr, E real
     (field-resolved time domain) or complex (frequency domain);
   - modal synthesis: `mul!(Et, Emt, M)`, `Emt` (nt, nmodes), `M` (nmodes, npts), and the
     projection back with `M'`.
   OpenBLAS runs its own pthreads, not Julia's, so Julia's thread count only matters through
   contention; it is recorded to check that. =#
import LinearAlgebra
import LinearAlgebra: mul!, BLAS
using Printf

const J = Threads.nthreads()
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results", "blas.csv"))

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

function cases()
    c = Tuple{String, Any, Any, Any}[]
    for nr in (64, 256, 1024), T in (Float64, ComplexF64)
        A = rand(T, 2^12, nr); M = rand(T, nr, nr); C = similar(A)
        push!(c, ("radial_$(T == Float64 ? "real" : "complex")", A, M, C))
    end
    for nm in (4, 16), np in (16, 64, 256)
        A = rand(2^13, nm); M = rand(nm, np); C = zeros(2^13, np)
        push!(c, ("modal_synth_$(nm)modes", A, M, C))
        B = rand(2^13, np); P = rand(np, nm); D = zeros(2^13, nm)
        push!(c, ("modal_proj_$(nm)modes", B, P, D))
    end
    c
end

isfile(OUT) || open(io -> println(io, "julia_threads,blas_threads,case,m,k,n,exec_s"), OUT, "w")
for (name, A, M, C) in cases(), nb in (1, 2, 4, 8)
    BLAS.set_num_threads(nb)
    t = exectime(() -> mul!(C, A, M))
    open(OUT, "a") do io
        @printf(io, "%d,%d,%s,%d,%d,%d,%.6e\n", J, nb, name, size(A, 1), size(A, 2),
                size(M, 2), t)
    end
end
println("BLAS_DONE J=$J")
