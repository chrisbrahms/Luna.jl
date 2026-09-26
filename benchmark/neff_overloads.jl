import LinearAlgebra
using Luna, Printf
import Luna: LinearOps, Capillary, Grid
Luna.set_fftw_mode(:estimate); Luna.set_fftw_threads(1); Luna.set_fftw_wisdom(false)
LinearAlgebra.BLAS.set_num_threads(1)

args(; kw...) = (125e-6, 1.0, :He, (0.0, 1.0))
kw(; k...) = (λ0=800e-9, τfwhm=10e-15, energy=150e-6, λlims=(150e-9, 4e-6), trange=400e-15,
              shotnoise=false, status_period=1e9, k...)

function timeit(f, n)
    f(); GC.gc(); t = time_ns(); for _ in 1:n; f(); end; (time_ns() - t)/n/1e9
end

# per-call cost of the z-dependent linop! closure
grid = Grid.RealGrid(800e-9, (150e-9, 4e-6), 400e-15)
coren, dens = Capillary.gradient(:He, 1.0, 0.0, 1.0)
m = Capillary.MarcatiliMode(125e-6, coren)
ms = Tuple(Capillary.MarcatiliMode(125e-6, coren; n=1, m=i) for i in 1:4)
function linopcost()
    l!, _ = LinearOps.make_linop(grid, m, 800e-9)
    out = zeros(ComplexF64, length(grid.ω))
    lm! = LinearOps.make_linop(grid, ms, 800e-9)
    outm = zeros(ComplexF64, length(grid.ω), 4)
    (timeit(() -> l!(out, 0.37), 20), timeit(() -> lm!(outm, 0.37), 5), copy((l!(out, 0.37); out)), copy((lm!(outm, 0.37); outm)))
end

function runs(tag)
    res = Dict{String,Any}()
    for (name, extra) in (("modeavg_auto", (;)), ("modeavg_quad", (linop_integral=:quadrature,)),
                          ("modal4_auto", (modes=4,)))
        f() = prop_capillary(args()...; kw(; extra...)...)
        f()  # compile
        t = @elapsed out = f()
        res[name] = (t, copy(out["Eω"]))
        @printf("%-6s %-14s %8.2f s\n", tag, name, t)
    end
    res
end

c1, cm1, l1, lm1 = linopcost()
@printf("with overload:    linop! %.3e s/call (N=%d), 4-mode linop! %.3e s/call\n", c1, length(grid.ω), cm1)
r1 = runs("with")

for f in (LinearOps.neff_β_grid, LinearOps.neff_grid)
    for mt in methods(f)
        mt.module === Capillary && Base.delete_method(mt)
    end
end
c2, cm2, l2, lm2 = linopcost()
@printf("generic neff:     linop! %.3e s/call, 4-mode linop! %.3e s/call\n", c2, cm2)
@printf("ratio: mode-avg %.2f, 4-mode %.2f\n", c2/c1, cm2/cm1)
@printf("linop value diff: mode-avg %.3e, 4-mode %.3e (max abs rel)\n",
        maximum(abs, l2 .- l1)/maximum(abs, l1), maximum(abs, lm2 .- lm1)/maximum(abs, lm1))
r2 = runs("generic")
for k in sort(collect(keys(r1)))
    d = maximum(abs, r2[k][2] .- r1[k][2]) / maximum(abs, r1[k][2])
    @printf("%-14s  time %.2f -> %.2f s  Eω diff %.3e\n", k, r1[k][1], r2[k][1], d)
end
println("BENCH_DONE")
