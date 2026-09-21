#= Multimode propagation timed on each transverse integral, device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/modal.jl
       julia --project=benchmark -t 8 benchmark/modal.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/modal.jl")'

   Metal is never a dependency of Luna or of this environment; the script uses whatever is
   already loaded and skips the device it cannot see.

   What the rows are:

   - `adaptive nb=N`: `NonlinearRHS.TransModal`, the default. The cubature driver hands
     over a round of transverse points at a time and the transform evaluates them in
     blocks of at most `maxbatch = N` points. `nb=1` is one block per point, which is the
     shape of the per-point implementation this replaced -- but not the same code, since
     the modal field still goes to the time domain once per right-hand side rather than
     once per point. The difference between `nb=1` and `nb=16` is what batching the
     points buys on its own.
   - `fixed nr=N`: `NonlinearRHS.TransModalFixed`, `modal_integral=:fixed`, an N-node
     Gauss rule in r. Every right-hand side costs the same, and the number of transform
     columns does not grow with N -- only the two matrix products and the responses do.

   As for `run.jl` and `device.jl`: one FFTW thread, one BLAS thread, `:estimate` planning
   and no wisdom, or the numbers are not comparable between runs. The Julia thread count
   is what it is; the plasma response shares a block's columns out over threads, which is
   the one part of this which is threaded.
=#
using Luna
import Luna: RK45, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear, Output,
             NonlinearRHS, PhysData, Utils, DeviceSpec, HostSpec, Ionisation
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
const NSTEPS = 10

"Time windows to sweep, in seconds. The grid size follows from the frequency limits."
const TRANGES = (400e-15, 1600e-15)

const NMODES = 4
const GAS = :Ar
const PRES = 0.1
const λ0 = 800e-9
const FLENGTH = 2e-3
const ENERGY = 50e-6

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

#= The tabulated rate a `plasma=true` run uses, built here rather than pre-calculated, as
   `benchmark/device.jl` does: what is timed is the response, not the PPT series. =#
function tablerate(gas)
    Ebs = Ionisation.barrier_suppression(PhysData.ionisation_potential(gas), 1.0)
    E = collect(range(2Ebs/5000, 2Ebs, length=1<<16))
    Ionisation.IonRatePPTAccel(E, Ionisation.IonRateADK(gas).(E))
end

"""
    prepare(spec, trange; plasma, modal_integral, nr, maxbatch)

Set up the `NMODES`-mode propagation on `spec` and return everything the timings need.
`boundary=:none`, so that what is measured is the transform and the stepper.
"""
function prepare(spec, trange; plasma=false, modal_integral=:adaptive, nr=64,
                 maxbatch=NonlinearRHS.MODAL_MAXBATCH)
    quiet() do
        grid = Grid.RealGrid(λ0, (200e-9, 3000e-9), trange)
        modes = Tuple(Capillary.MarcatiliMode(75e-6, GAS, PRES; n=1, m=mi, loss=false)
                      for mi in 1:NMODES)
        #= The density is constant here, and `PhysData.density` goes through CoolProp,
           which costs more per call than the whole right-hand side. =#
        ρ = PhysData.density(GAS, PRES)
        dens = z -> ρ
        resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS)),)
        if plasma
            resp = (resp...,
                    Nonlinear.PlasmaCumtrapz(grid.to, zeros(length(grid.to)),
                                             tablerate(GAS),
                                             PhysData.ionisation_potential(GAS)))
        end
        linop = LinearOps.make_const_linop(grid, modes, λ0)
        inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=ENERGY)
        Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, modes, :y;
                                       modal_integral, nr, maxbatch, device=spec)
        (Eω, Luna.upload_like(Eω, linop), transform, grid, FT)
    end
end

#= Device work is queued asynchronously, so a timing which does not synchronise measures
   the launch, not the kernel. `device_synchronize` is a no-op on the CPU. =#
sync(spec) = Luna.device_synchronize(spec)

function rhstime(spec, trange; kwargs...)
    Eω, _, transform, _, _ = prepare(spec, trange; kwargs...)
    nl = similar(Eω)
    #= The first evaluation of an adaptive transform allocates the blocks and plans the
       transforms of every round width the driver asks for, so it is run once before the
       benchmark rather than being the benchmark's first sample. =#
    transform(nl, Eω, 0.0)
    b = @benchmarkable (($transform)($nl, $Eω, 0.0); sync($spec)) seconds=BUDGET
    (length(Eω), bminimum(brun(b)).time/1e9)
end

function steptime(spec, trange; kwargs...)
    Eω, linop, transform, _, _ = prepare(spec, trange; kwargs...)
    dz = FLENGTH/NSTEPS
    nl = similar(Eω)
    transform(nl, Eω, 0.0)
    b = @benchmarkable (RK45.step!(s); sync($spec)) setup=(
            s = RK45.PreconStepper($transform, $linop, $Eω, 0.0, $dz;
                                   max_dt=$dz, min_dt=$dz)
        ) evals=1 seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

function proptime(spec, trange; samples=3, kwargs...)
    minimum(1:samples) do _
        Eω, linop, transform, grid, FT = prepare(spec, trange; kwargs...)
        out = Output.MemoryOutput(0, FLENGTH, 3, Output.nostats)
        dz = FLENGTH/NSTEPS
        t = @elapsed quiet() do
            Luna.run(Eω, grid, linop, transform, FT, out;
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

#= The rules to time. The adaptive rows are host-Float64 only: `Cubature` is host scalar
   code returning `Vector{Float64}`. =#
const HOSTONLY = (("adaptive nb=1", (; modal_integral=:adaptive, maxbatch=1)),
                  ("adaptive nb=4", (; modal_integral=:adaptive, maxbatch=4)),
                  ("adaptive nb=16", (; modal_integral=:adaptive, maxbatch=16)),
                  ("adaptive nb=32", (; modal_integral=:adaptive, maxbatch=32)))
const EVERYWHERE = (("fixed nr=32", (; modal_integral=:fixed, nr=32)),
                    ("fixed nr=64", (; modal_integral=:fixed, nr=64)))

specs = Any[("CPU Float64", HostSpec()), ("CPU Float32", DeviceSpec(Array, Float32))]
for name in Luna.devicenames()
    spec = Luna.DEVICES[name].spec
    Luna.device_functional(name) &&
        push!(specs, (string(name)*" "*string(Luna.realtype(spec)), spec))
end

@printf("%d %s modes, %s at %g bar, %g m, %d fixed steps, boundary=:none\n",
        NMODES, "HE1m", GAS, PRES, FLENGTH, NSTEPS)
@printf("%d Julia threads, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n",
        Threads.nthreads())

for (label, rkw) in (("Kerr", (;)),
                     ("Kerr + plasma (tabulated rate)", (; plasma=true)))
    @printf("\nresponses: %s\n\n", label)
    @printf("%-14s %-14s %8s %10s %12s %12s %12s\n",
            "device", "integral", "trange", "state", "rhs", "step", "prop")
    @printf("%s\n", "-"^88)
    for trange in TRANGES
        for (rule, rulekw) in HOSTONLY
            n, tr = rhstime(HostSpec(), trange; rkw..., rulekw...)
            ts = steptime(HostSpec(), trange; rkw..., rulekw...)
            tp = proptime(HostSpec(), trange; rkw..., rulekw...)
            @printf("%-14s %-14s %6.0f fs %10d %12s %12s %12s\n",
                    "CPU Float64", rule, trange*1e15, n,
                    fmttime(tr), fmttime(ts), fmttime(tp))
            flush(stdout)
        end
        for (name, spec) in specs, (rule, rulekw) in EVERYWHERE
            n, tr = rhstime(spec, trange; rkw..., rulekw...)
            ts = steptime(spec, trange; rkw..., rulekw...)
            tp = proptime(spec, trange; rkw..., rulekw...)
            @printf("%-14s %-14s %6.0f fs %10d %12s %12s %12s\n",
                    name, rule, trange*1e15, n,
                    fmttime(tr), fmttime(ts), fmttime(tp))
            flush(stdout)
        end
    end
end
