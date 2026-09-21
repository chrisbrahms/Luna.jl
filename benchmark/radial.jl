#= Radially symmetric free-space Kerr propagation timed on each device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/radial.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/radial.jl")'

   `benchmark/device.jl` times the mode-averaged transform, which has one column and is
   launch-bound on any GPU. This one times the radial transform, which has `N` of them:
   two Hankel GEMMs of `(nto*npol, N) x (N, N)` and a batched FFT over `N` columns per
   right-hand side. That is where a GPU is supposed to win, so the sweep is over `N` and
   the number to read off is the crossover -- the `N` at which the device beats the host.

   As for `run.jl` and `device.jl`: one Julia thread, one FFTW thread, one BLAS thread,
   `:estimate` planning and no wisdom, or the numbers are not comparable between runs.
=#
using Luna
import Luna: RK45, Grid, Modes, Fields, LinearOps, Nonlinear, NonlinearRHS, Output,
             PhysData, Utils, Boundaries, DeviceSpec, HostSpec
import LinearAlgebra
import BenchmarkTools: @benchmarkable, run as brun, minimum as bminimum
import Printf: @printf, @sprintf
import Logging

Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Luna.set_fftw_wisdom(false)

"Time budget in seconds for each `rhs` and `step` benchmark."
const BUDGET = 2.0

"Number of fixed steps in the timed propagation."
const NSTEPS = 20

"""
Numbers of radial points to sweep. `LUNA_BENCH_NRADIAL` (comma-separated) overrides it,
which is how the crossover between the host and the device was narrowed down.
"""
const NRADIAL = haskey(ENV, "LUNA_BENCH_NRADIAL") ?
    Tuple(parse.(Int, split(ENV["LUNA_BENCH_NRADIAL"], ','))) : (64, 256, 1024)

const GAS = :Ar
const PRES = 1.0
const λ0 = 800e-9
const TRANGE = 100e-15
const ΛLIMS = (400e-9, 2000e-9)
const R = 1e-3
const W0 = 200e-6
const FLENGTH = 1e-2
const ENERGY = 1e-6

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

"""
    prepare(spec, N)

Set up the radial Kerr propagation with `N` radial points on `spec` and return everything
the timings need. `boundary=:none`, so that what is timed is the transform and the
propagator rather than the absorbers (`benchmark/boundaries.jl` times those).
"""
function prepare(spec, N)
    quiet() do
        grid = Grid.RealGrid(λ0, ΛLIMS, TRANGE)
        rg = Grid.RadialGrid(R, N)
        nfunλ = PhysData.ref_index_fun(GAS, PRES)
        nfun = (λ; z=0.0) -> nfunλ(λ)
        linop = LinearOps.make_const_linop(grid, rg, nfun)
        ρ = PhysData.density(GAS, PRES)
        dens = z -> ρ
        resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS)),)
        normfun = NonlinearRHS.const_norm_radial(grid, rg, nfun)
        inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy=ENERGY, w0=W0,
                                         propz=-FLENGTH)
        Eω, transform, FT = Luna.setup(grid, rg, dens, normfun, resp, inputs; device=spec)
        #= Both forms of the operator. `Luna.run` needs the host `Float64` one, because
           `Boundaries.setup` wraps it there (the evanescent clamp and the two absorbers
           are host scalar code) and `Luna.run` uploads the result itself; a stepper built
           by hand, as `steptime` does, needs it already on the device.

           The device copy is made from the *clamped* operator, not the raw one:
           `boundary=:none` still puts the free-space operator through
           `Boundaries.evanescent`, and without it the evanescent `κ` on this grid makes
           the interaction picture's `exp(κ Δz)` overflow, so a hand-built stepper would
           be timed on `Inf`/`NaN` data. `evanescent` also tapers the transform's
           normalisation to match, which is what `Luna.run` does. =#
        clamped = Boundaries.evanescent(linop, transform, FLENGTH/NSTEPS)
        (Eω, linop, Luna.upload_like(Eω, clamped), transform, grid, FT)
    end
end

"An output which copies the saved field to the host, as `Luna.ScaledOutput` does."
struct ToHost{O}
    o::O
end
(h::ToHost)(y, t, dt, yfun) = h.o(Array(y), t, dt, ti -> Array(yfun(ti)))
(h::ToHost)(args...; kwargs...) = h.o(args...; kwargs...)

#= Device work is queued asynchronously, so a timing which does not synchronise measures
   the launch, not the kernel. `device_synchronize` is a no-op on the CPU. =#
sync(spec) = Luna.device_synchronize(spec)

function rhstime(spec, N)
    Eω, _, _, transform, _, _ = prepare(spec, N)
    nl = similar(Eω)
    b = @benchmarkable (($transform)($nl, $Eω, 0.0); sync($spec)) seconds=BUDGET
    (length(Eω), bminimum(brun(b)).time/1e9)
end

function hankeltime(spec, N)
    _, _, _, transform, _, _ = prepare(spec, N)
    b = @benchmarkable (Grid.radial_matmul!($(transform.Eto_r), $(transform.Eto_k),
                                            $(transform.Tbwd)); sync($spec)) seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

function steptime(spec, N)
    Eω, _, linop, transform, _, _ = prepare(spec, N)
    dz = FLENGTH/NSTEPS
    b = @benchmarkable (RK45.step!(s); sync($spec)) setup=(
            s = RK45.PreconStepper($transform, $linop, $Eω, 0.0, $dz;
                                   max_dt=$dz, min_dt=$dz)
        ) evals=1 seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

function proptime(spec, N; samples=3)
    minimum(1:samples) do _
        Eω, linop, _, transform, grid, FT = prepare(spec, N)
        out = Output.MemoryOutput(0, FLENGTH, 3, Output.nostats)
        output = Utils.isdevice(Eω) ? ToHost(out) : out
        dz = FLENGTH/NSTEPS
        t = @elapsed quiet() do
            Luna.run(Eω, grid, linop, transform, FT, output;
                     zmax=FLENGTH, boundary=:none,
                     init_dz=dz, min_dz=dz, max_dz=dz)
            sync(spec)
        end
        t
    end
end

fmttime(t) = t >= 1    ? @sprintf("%7.3f s ", t)    :
             t >= 1e-3 ? @sprintf("%7.3f ms", 1e3t) :
             t >= 1e-6 ? @sprintf("%7.3f µs", 1e6t) :
                         @sprintf("%7.3f ns", 1e9t)

specs = Any[("CPU Float64", HostSpec()), ("CPU Float32", DeviceSpec(Array, Float32))]
for name in Luna.devicenames()
    spec = Luna.DEVICES[name].spec
    Luna.device_functional(name) &&
        push!(specs, (string(name)*" "*string(Luna.realtype(spec)), spec))
end

@printf("radial Kerr: %s at %g bar, R = %g m, %g m, %d fixed steps, boundary=:none\n",
        GAS, PRES, R, FLENGTH, NSTEPS)
@printf("%d Julia threads, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n",
        Threads.nthreads())
@printf("\n%-14s %8s %10s %12s %12s %12s %12s\n",
        "device", "N", "state", "hankel", "rhs", "step", "prop")
@printf("%s\n", "-"^86)
for N in NRADIAL, (name, spec) in specs
    n, tr = rhstime(spec, N)
    th = hankeltime(spec, N)
    ts = steptime(spec, N)
    tp = proptime(spec, N)
    @printf("%-14s %8d %10d %12s %12s %12s %12s\n",
            name, N, n, fmttime(th), fmttime(tr), fmttime(ts), fmttime(tp))
    flush(stdout)
end
