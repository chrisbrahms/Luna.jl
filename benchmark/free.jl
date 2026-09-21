#= 3-D Cartesian free-space envelope Kerr propagation timed on each device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/free.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/free.jl")'

   `benchmark/radial.jl` times the radial transform, whose transverse step is a GEMM.
   This one times `TransFree`, whose transverse step is part of the joint FFT: one
   region-(1, 3, 4) transform each way over `Nx*Ny` columns per right-hand side, with no
   matrix multiplication anywhere. The sweep is over the transverse grid and the number
   to read off is the crossover -- the size at which the device beats the host.

   As for `run.jl`, `device.jl` and `radial.jl`: one Julia thread, one FFTW thread, one
   BLAS thread, `:estimate` planning and no wisdom, or the numbers are not comparable
   between runs.
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

"Time budget in seconds for each `fft`, `rhs` and `step` benchmark."
const BUDGET = 2.0

"Number of fixed steps in the timed propagation."
const NSTEPS = 10

"""
Timed runs of the whole propagation per row, after one discarded warm-up run; the column
reports the minimum. `LUNA_BENCH_PROPSAMPLES` overrides it. See [`proptime`](@ref) for why
this is not 1.
"""
const PROPSAMPLES = parse(Int, get(ENV, "LUNA_BENCH_PROPSAMPLES", "5"))

"""
Transverse grid sizes to sweep; each is used for both `Nx` and `Ny`.
`LUNA_BENCH_NFREE` (comma-separated) overrides it.
"""
const NFREE = haskey(ENV, "LUNA_BENCH_NFREE") ?
    Tuple(parse.(Int, split(ENV["LUNA_BENCH_NFREE"], ','))) : (32, 64, 128)

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

Set up the 3-D free-space envelope Kerr propagation on an `N x N` transverse grid on
`spec` and return everything the timings need. `boundary=:none`, so that what is timed is
the transform and the propagator rather than the absorbers.
"""
function prepare(spec, N)
    quiet() do
        grid = Grid.EnvGrid(λ0, ΛLIMS, TRANGE)
        xygrid = Grid.FreeGrid(R, N, R, N)
        nfunλ = PhysData.ref_index_fun(GAS, PRES)
        nfun = (λ; z=0.0) -> nfunλ(λ)
        linop = LinearOps.make_const_linop(grid, xygrid, nfun)
        ρ = PhysData.density(GAS, PRES)
        dens = z -> ρ
        resp = (Nonlinear.Kerr_env(PhysData.γ3_gas(GAS)),)
        normfun = NonlinearRHS.const_norm_free(grid, xygrid, nfun)
        inputs = Fields.GaussGaussField(;λ0, τfwhm=20e-15, energy=ENERGY, w0=W0,
                                         propz=-FLENGTH)
        Eω, transform, FT = Luna.setup(grid, xygrid, dens, normfun, resp, inputs;
                                       device=spec)
        #= The device copy is made from the *clamped* operator: `boundary=:none` still
           puts a free-space operator through `Boundaries.evanescent`, and without it the
           evanescent κ on this grid makes the interaction picture's exp(κ Δz) overflow,
           so a hand-built stepper would be timed on Inf/NaN data. `evanescent` also
           tapers the transform's normalisation to match, which is what `Luna.run` does. =#
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

"One joint inverse transform: (ω, kx, ky) -> (t, x, y), which is half the FFT cost of a
right-hand side."
function ffttime(spec, N)
    Eω, _, _, t, _, _ = prepare(spec, N)
    b = @benchmarkable (NonlinearRHS.to_time!($(t.Eto), $Eω, $(t.Eωo), $(t.IFT));
                        sync($spec)) seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

function rhstime(spec, N)
    Eω, _, _, transform, _, _ = prepare(spec, N)
    nl = similar(Eω)
    b = @benchmarkable (($transform)($nl, $Eω, 0.0); sync($spec)) seconds=BUDGET
    (length(Eω), bminimum(brun(b)).time/1e9)
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

"""
    proptime(spec, N; samples=PROPSAMPLES)

Wall time of the whole fixed-step propagation, as the minimum of `samples` runs after one
discarded warm-up run.

The warm-up and the sample count matter. A propagation is timed with `@elapsed` rather
than by `BenchmarkTools`, so the first one in a process pays for compiling the stepper for
this combination of types, and on a device for building and caching the FFT graphs and
kernels; each run also allocates a fresh set of field-sized buffers, so the garbage
collector can land in the middle of one. Review 1 of `gpu/21-free-device` measured a
factor of 2.2 between two runs of this column at 64 x 64 with `samples=2` and no warm-up,
while the `fft`, `rhs` and `step` columns -- which `BenchmarkTools` runs many times --
reproduced to under 2 %. With the warm-up and `samples = $(PROPSAMPLES)` the column is
stable, but the ratios worth quoting are still the ones from `step`.
"""
function proptime(spec, N; samples=PROPSAMPLES)
    runonce() = begin
        Eω, linop, _, transform, grid, FT = prepare(spec, N)
        out = Output.MemoryOutput(0, FLENGTH, 3, Output.nostats)
        output = Utils.isdevice(Eω) ? ToHost(out) : out
        dz = FLENGTH/NSTEPS
        @elapsed quiet() do
            Luna.run(Eω, grid, linop, transform, FT, output;
                     zmax=FLENGTH, boundary=:none,
                     init_dz=dz, min_dz=dz, max_dz=dz)
            sync(spec)
        end
    end
    runonce() # warm-up: compilation, device graph caching, first-touch allocation
    minimum(_ -> runonce(), 1:samples)
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

@printf("3-D free-space envelope Kerr: %s at %g bar, R = %g m, %g m, %d fixed steps, \
         boundary=:none\n", GAS, PRES, R, FLENGTH, NSTEPS)
@printf("%d Julia threads, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n",
        Threads.nthreads())
@printf("\n%-14s %10s %12s %12s %12s %12s %12s\n",
        "device", "NxN", "state", "fft", "rhs", "step", "prop")
@printf("%s\n", "-"^88)
for N in NFREE, (name, spec) in specs
    tf = ffttime(spec, N)
    n, tr = rhstime(spec, N)
    ts = steptime(spec, N)
    tp = proptime(spec, N)
    @printf("%-14s %10s %12d %12s %12s %12s %12s\n",
            name, "$(N)x$(N)", n, fmttime(tf), fmttime(tr), fmttime(ts), fmttime(tp))
    flush(stdout)
end
