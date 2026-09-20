#= Mode-averaged Kerr propagation timed on each device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/device.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/device.jl")'

   Metal is never a dependency of Luna or of this environment; the script uses whatever
   is already loaded and skips the device it cannot see.

   `benchmark/run.jl` times the regression case matrix on the CPU in double precision.
   This one times the *same* propagation on `Array{Float64}`, `Array{Float32}` and the
   GPU, which is the comparison the device model is for. It sweeps the time-grid size,
   because a small single-column case is launch-bound on any GPU and the scaling with
   problem size is the number that matters.

   As for `run.jl`: one Julia thread, one FFTW thread, one BLAS thread, `:estimate`
   planning and no wisdom, or the numbers are not comparable between runs.
=#
using Luna
import Luna: RK45, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear, Output,
             PhysData, Utils, DeviceSpec, HostSpec
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

"Time windows to sweep, in seconds. The grid size follows from the frequency limits."
const TRANGES = (400e-15, 1600e-15, 6400e-15)

const GAS = :He
const PRES = 1.0
const λ0 = 800e-9
const FLENGTH = 1e-2

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

"""
    prepare(spec, trange)

Set up the mode-averaged Kerr propagation on `spec` and return everything the timings
need. `boundary=:none`, since the absorbing boundaries are host code until `gpu/11`.
"""
function prepare(spec, trange)
    quiet() do
        grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), trange)
        m = Capillary.MarcatiliMode(75e-6, GAS, PRES, loss=false)
        aeff(z) = Modes.Aeff(m, z=z)
        #= The density is constant here, and `PhysData.density` goes through CoolProp,
           which costs ~85 us a call -- more than the whole right-hand side. The simple
           interface precomputes it the same way. =#
        ρ = PhysData.density(GAS, PRES)
        dens = z -> ρ
        resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS)),)
        linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0)
        inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=1e-6)
        Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                       constβ=true, device=spec)
        (Eω, Luna.upload_like(Eω, linop), transform, grid, FT)
    end
end

"An output which copies the saved field to the host; `ScaledOutput` proper is gpu/11's."
struct ToHost{O}
    o::O
end
(h::ToHost)(y, t, dt, yfun) = h.o(Array(y), t, dt, ti -> Array(yfun(ti)))
(h::ToHost)(args...; kwargs...) = h.o(args...; kwargs...)

#= Device work is queued asynchronously, so a timing which does not synchronise measures
   the launch, not the kernel. `device_synchronize` is a no-op on the CPU. =#
sync(spec) = Luna.device_synchronize(spec)

function rhstime(spec, trange)
    Eω, _, transform, _, _ = prepare(spec, trange)
    nl = similar(Eω)
    b = @benchmarkable (($transform)($nl, $Eω, 0.0); sync($spec)) seconds=BUDGET
    (length(Eω), bminimum(brun(b)).time/1e9)
end

function steptime(spec, trange)
    Eω, linop, transform, _, _ = prepare(spec, trange)
    dz = FLENGTH/NSTEPS
    b = @benchmarkable (RK45.step!(s); sync($spec)) setup=(
            s = RK45.PreconStepper($transform, $linop, $Eω, 0.0, $dz;
                                   max_dt=$dz, min_dt=$dz)
        ) evals=1 seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

function proptime(spec, trange; samples=3)
    minimum(1:samples) do _
        Eω, linop, transform, grid, FT = prepare(spec, trange)
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
    Luna.device_functional(name) && push!(specs, (string(name)*" "*string(Luna.realtype(spec)), spec))
end

@printf("Mode-averaged Kerr, %s at %g bar, %g m, %d fixed steps, boundary=:none\n",
        GAS, PRES, FLENGTH, NSTEPS)
@printf("1 Julia thread, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n\n")
@printf("%-14s %8s %10s %12s %12s %12s\n",
        "device", "trange", "state", "rhs", "step", "prop")
@printf("%s\n", "-"^72)
for trange in TRANGES, (name, spec) in specs
    n, tr = rhstime(spec, trange)
    ts = steptime(spec, trange)
    tp = proptime(spec, trange)
    @printf("%-14s %6.0f fs %10d %12s %12s %12s\n",
            name, trange*1e15, n, fmttime(tr), fmttime(ts), fmttime(tp))
    flush(stdout)
end
