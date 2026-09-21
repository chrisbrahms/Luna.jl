#= A pressure-graded and a tapered mode-averaged propagation timed with the two ways of
   integrating a z-dependent linear operator, `linop_integral=:tabulated` and
   `:quadrature`, on each device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/tabulated.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/tabulated.jl")'

   These are the cases GPU_PLAN.md section 4.5 is about. With `:quadrature` a z-dependent
   operator is a host scalar loop over `Modes.neff` run fifteen times per Gauss--Kronrod
   rule at every stage, plus a host evaluation and upload of β at every right-hand side;
   with `:tabulated`, both are interpolations in a table which lives where the state does.
   On a device that is the difference between a propagation which synchronises with the
   host at every stage and one which does not.

   The two compute the same integral to their tolerances (agreement measured at 1e-7 on
   these cases), so unlike the one-point rule they replaced this is a like-for-like
   comparison of cost. Two numbers are reported per row: the setup cost of the tables,
   paid once, and the time for a fixed-step propagation.

   As for the other benchmark scripts: one Julia thread, one FFTW thread, one BLAS thread,
   `:estimate` planning and no wisdom, or the numbers are not comparable between runs.
=#
using Luna
import Luna: Grid, Modes, Capillary, Fields, LinearOps, Nonlinear, Output, PhysData,
             Utils, DeviceSpec, HostSpec
import LinearAlgebra
import Printf: @printf, @sprintf
import Logging

Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Luna.set_fftw_wisdom(false)

const λ0 = 800e-9
const FLENGTH = 0.1
const NSTEPS = 20
const TRANGES = (400e-15, 1600e-15)

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

"An output which copies the saved field to the host, for the timing runs."
struct ToHost{O}
    o::O
end
(h::ToHost)(y, t, dt, yfun) = h.o(Array(y), t, dt, ti -> Array(yfun(ti)))
(h::ToHost)(args...; kwargs...) = h.o(args...; kwargs...)

sync(spec) = Luna.device_synchronize(spec)

function prepare(spec, trange; kind=:gradient)
    quiet() do
        grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), trange)
        if kind === :gradient
            coren, densityfun = Capillary.gradient(:Ar, FLENGTH, 1.0, 0.0)
            m = Capillary.MarcatiliMode(75e-6, coren, loss=false)
        else
            afun = z -> 75e-6 + (50e-6 - 75e-6)*z/FLENGTH
            m = Capillary.MarcatiliMode(afun, :Ar, 1.0, loss=false, model=:full)
            ρ = PhysData.density(:Ar, 1.0)
            densityfun = z -> ρ
        end
        aeff(z) = Modes.Aeff(m, z=z)
        resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(:Ar)),)
        linop, βfun! = LinearOps.make_linop(grid, m, λ0)
        inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=1e-6)
        Eω, transform, FT = Luna.setup(grid, densityfun, resp, inputs, βfun!, aeff;
                                       device=spec)
        (Eω, linop, transform, grid, FT)
    end
end

#= The propagation, and separately the setup of the tables. `Luna.run` builds them, so the
   setup cost is measured as a run of one step: everything before the stepping is the same
   and the difference between the two is the table build. =#
function proptime(spec, trange; kind=:gradient, linop_integral=:tabulated, samples=3)
    minimum(1:samples) do _
        Eω, linop, transform, grid, FT = prepare(spec, trange; kind)
        out = Output.MemoryOutput(0, FLENGTH, 3, Output.nostats)
        output = Utils.isdevice(Eω) ? ToHost(out) : out
        dz = FLENGTH/NSTEPS
        t = @elapsed quiet() do
            Luna.run(Eω, grid, linop, transform, FT, output;
                     zmax=FLENGTH, boundary=:none, init_dz=dz, min_dz=dz, max_dz=dz,
                     linop_integral)
            sync(spec)
        end
        t
    end
end

function tabletime(spec, trange; kind=:gradient)
    Eω, linop, _, _, _ = prepare(spec, trange; kind)
    minimum(1:3) do _
        @elapsed LinearOps.TabulatedLinop(linop, Eω, 0.0, FLENGTH*1.05;
                                          tol=LinearOps.DEFAULT_LINOP_TOL, quiet=true)
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

@printf("Ar, %g m, %d fixed steps, boundary=:none, linop_tol = %g\n",
        FLENGTH, NSTEPS, LinearOps.DEFAULT_LINOP_TOL)
@printf("%d Julia threads, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n",
        Threads.nthreads())
for kind in (:gradient, :taper)
    @printf("\n%s\n\n", kind === :gradient ? "pressure gradient 1 -> 0 bar" :
                                             "taper 75 -> 50 um")
    @printf("%-14s %8s %10s %12s %12s %8s %12s\n",
            "device", "trange", "state", "quadrature", "tabulated", "speedup", "table setup")
    @printf("%s\n", "-"^82)
    for trange in TRANGES, (name, spec) in specs
        n = length(prepare(spec, trange; kind)[1])
        t0 = proptime(spec, trange; kind, linop_integral=:quadrature)
        t1 = proptime(spec, trange; kind, linop_integral=:tabulated)
        ts = tabletime(spec, trange; kind)
        @printf("%-14s %6.0f fs %10d %12s %12s %7.2fx %12s\n",
                name, trange*1e15, n, fmttime(t0), fmttime(t1), t0/t1, fmttime(ts))
        flush(stdout)
    end
end
