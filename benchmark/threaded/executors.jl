# Serial Base broadcast vs a KernelAbstractions CPU kernel vs Threads.@threads/@spawn
# chunks, over the same Broadcasted. Needs KernelAbstractions, which Luna does not depend
# on: run it in an environment that has it, e.g.
#   julia --project=<env with KernelAbstractions> -t 4 benchmark/threaded/executors.jl
using KernelAbstractions, Printf
import Base.Broadcast: broadcasted, instantiate, Broadcasted
@kernel function bck!(dest, bc)
    I = @index(Global, Linear)
    @inbounds dest[I] = bc[I]
end
ka!(dest, bc) = (bck!(CPU(), 1024)(dest, bc; ndrange=length(dest)); dest)
function th!(dest, bc)
    n = length(dest); nt = Threads.nthreads(); ch = cld(n, nt)
    Threads.@threads :static for k in 1:nt
        @inbounds for i in (k-1)*ch+1:min(k*ch, n)
            dest[i] = bc[i]
        end
    end
    dest
end
function sp!(dest, bc)
    n = length(dest); nt = Threads.nthreads(); ch = cld(n, nt)
    @sync for k in 1:nt
        Threads.@spawn @inbounds for i in (k-1)*ch+1:min(k*ch, n)
            dest[i] = bc[i]
        end
    end
    dest
end
function tm(f, n=200)
    f(); t = time_ns(); for _ in 1:n; f(); end; (time_ns()-t)/n/1e3
end
for n in (2^12, 2^14, 2^16, 2^18, 2^20)
    y = rand(ComplexF64, n); a = rand(ComplexF64, n); b = rand(ComplexF64, n); d = similar(y)
    # propagator-like: y * exp(a - b); plasma-like real
    bc = instantiate(broadcasted((y, a, b) -> y*exp(a - b), y, a, b))
    x = rand(n); r = similar(x)
    bcr = instantiate(broadcasted((x) -> x^2*sqrt(abs(x)) + 1.3x, x))
    t0 = tm(() -> copyto!(d, bc)); t1 = tm(() -> ka!(d, bc)); t2 = tm(() -> th!(d, bc)); t3 = tm(() -> sp!(d, bc))
    s0 = tm(() -> copyto!(r, bcr)); s1 = tm(() -> ka!(r, bcr)); s2 = tm(() -> th!(r, bcr)); s3 = tm(() -> sp!(r, bcr))
    @printf("n=%7d  cexp: serial %8.1f us  KA %8.1f  @threads %8.1f  @spawn %8.1f | cheap real: %7.1f %7.1f %7.1f %7.1f\n",
            n, t0, t1, t2, t3, s0, s1, s2, s3)
end
