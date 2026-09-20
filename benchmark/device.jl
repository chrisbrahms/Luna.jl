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
             PhysData, Utils, DeviceSpec, HostSpec, Ionisation, Raman
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

"The gas the Raman rows use; `GAS` has no Raman response."
const RAMANGAS = :N2
const PRES = 1.0
const λ0 = 800e-9
const FLENGTH = 1e-2

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

#= The tabulated rate a `plasma=true` run uses, built here rather than pre-calculated:
   `IonRatePPTAccel(E, rate)` is the constructor the PPT cache calls, the axis is uniform
   (so the kernel is the spline lookup a real cached rate uses), and it costs
   milliseconds. What is timed is the response, not the PPT series. =#
function tablerate(gas)
    Ebs = Ionisation.barrier_suppression(PhysData.ionisation_potential(gas), 1.0)
    E = collect(range(2Ebs/5000, 2Ebs, length=1<<16))
    Ionisation.IonRatePPTAccel(E, Ionisation.IonRateADK(gas).(E))
end

"""
    prepare(spec, trange; plasma=false, raman=false, gas=GAS)

Set up the mode-averaged propagation on `spec` and return everything the timings need.
`boundary=:none`. With `plasma=true` the response set is Kerr *and*
`Nonlinear.PlasmaCumtrapz` on a tabulated rate, which is what `prop_capillary` builds by
default for a field-resolved run in a non-Raman gas; with `raman=true` it is Kerr and
`Nonlinear.RamanPolarField`, which is what it builds for a molecular one.
"""
function prepare(spec, trange; plasma=false, raman=false, gas=GAS)
    quiet() do
        grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), trange)
        m = Capillary.MarcatiliMode(75e-6, gas, PRES, loss=false)
        aeff(z) = Modes.Aeff(m, z=z)
        #= The density is constant here, and `PhysData.density` goes through CoolProp,
           which costs ~85 us a call -- more than the whole right-hand side. The simple
           interface precomputes it the same way. =#
        ρ = PhysData.density(gas, PRES)
        dens = z -> ρ
        resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),)
        if plasma
            resp = (resp...,
                    Nonlinear.PlasmaCumtrapz(grid.to, zeros(length(grid.to)),
                                             tablerate(gas),
                                             PhysData.ionisation_potential(gas)))
        end
        if raman
            resp = (resp..., Nonlinear.RamanPolarField(grid.to,
                                                       Raman.raman_response(grid.to, gas)))
        end
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

function rhstime(spec, trange; kwargs...)
    Eω, _, transform, _, _ = prepare(spec, trange; kwargs...)
    nl = similar(Eω)
    b = @benchmarkable (($transform)($nl, $Eω, 0.0); sync($spec)) seconds=BUDGET
    (length(Eω), bminimum(brun(b)).time/1e9)
end

function steptime(spec, trange; kwargs...)
    Eω, linop, transform, _, _ = prepare(spec, trange; kwargs...)
    dz = FLENGTH/NSTEPS
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

@printf("%s at %g bar, %g m, %d fixed steps, boundary=:none\n",
        GAS, PRES, FLENGTH, NSTEPS)
@printf("%d Julia threads, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n",
        Threads.nthreads())
for (label, rkw) in (("Kerr", (;)),
                     ("Kerr + plasma (tabulated rate)", (; plasma=true)),
                     ("Kerr + Raman ($RAMANGAS)", (; raman=true, gas=RAMANGAS)))
    @printf("\nresponses: %s\n\n", label)
    @printf("%-14s %8s %10s %12s %12s %12s\n",
            "device", "trange", "state", "rhs", "step", "prop")
    @printf("%s\n", "-"^72)
    for trange in TRANGES, (name, spec) in specs
        n, tr = rhstime(spec, trange; rkw...)
        ts = steptime(spec, trange; rkw...)
        tp = proptime(spec, trange; rkw...)
        @printf("%-14s %6.0f fs %10d %12s %12s %12s\n",
                name, trange*1e15, n, fmttime(tr), fmttime(ts), fmttime(tp))
        flush(stdout)
    end
end

#= The plasma response on its own, over a block of transverse columns -- the shape a
   radial or free-space transform passes, and the only case where the host path threads
   the columns. Run the script with `-t 1` and with `-t N` to see the threading; the
   thread count is fixed for the process. =#
@printf("\nplasma response alone, %d-sample block of columns, %d threads\n\n",
        1<<11, Threads.nthreads())
@printf("%-14s %10s %12s\n", "device", "columns", "batched!")
@printf("%s\n", "-"^40)
let nt = 1<<11, ionpot = PhysData.ionisation_potential(GAS), ir = tablerate(GAS),
    ρ = PhysData.density(GAS, PRES)
    t = collect(range(-60e-15, 60e-15, length=nt))
    E = @. 6e10*exp(-t^2/(2*(10e-15/1.66)^2))*cos(2π*PhysData.c/800e-9*t)
    p0 = Nonlinear.PlasmaCumtrapz(t, E, ir, ionpot)
    for (name, spec) in specs, ncols in (1, 16, 128)
        Eh = repeat(reshape(E, nt, 1, 1), 1, 1, ncols)
        sc = Luna.realtype(spec) === Float64 ? Luna.UNIT_SCALING :
             Luna.unitscaling(Luna.realtype(spec), Eh, PhysData.ε_0)
        Ed = Luna.todevice(spec, Eh ./ sc.Eref)
        p = Nonlinear.rescale(p0, spec, sc, Ed)
        out = Luna.alloc(spec, Luna.realtype(spec), size(Ed))
        b = @benchmarkable (Nonlinear.batched!($p, $out, $Ed, $ρ, $sc); sync($spec)
                            ) seconds=BUDGET
        @printf("%-14s %10d %12s\n", name, ncols, fmttime(bminimum(brun(b)).time/1e9))
        flush(stdout)
    end
end

#= The Raman response on its own, over a block of transverse columns: what the batched
   FFT plans buy over the per-column ones this response used to make. The response
   function is evaluated once here, at the first call, and cached on the density, which
   is what a constant-pressure propagation sees too. =#
@printf("\nRaman response alone (%s), %d-sample block of columns\n\n",
        RAMANGAS, 1<<11)
@printf("%-14s %10s %12s\n", "device", "columns", "batched!")
@printf("%s\n", "-"^40)
let nt = 1<<11, ρ = PhysData.density(RAMANGAS, PRES)
    t = collect(range(-60e-15, 60e-15, length=nt))
    E = @. 1e10*exp(-t^2/(2*(10e-15/1.66)^2))*cos(2π*PhysData.c/800e-9*t)
    for (name, spec) in specs, ncols in (1, 16, 128)
        Eh = repeat(reshape(E, nt, 1, 1), 1, 1, ncols)
        sc = Luna.realtype(spec) === Float64 ? Luna.UNIT_SCALING :
             Luna.unitscaling(Luna.realtype(spec), Eh, PhysData.ε_0)
        Ed = Luna.todevice(spec, Eh ./ sc.Eref)
        R = quiet() do
            Nonlinear.rescale(
                Nonlinear.RamanPolarField(t, Raman.raman_response(t, RAMANGAS)),
                spec, sc, Ed)
        end
        out = Luna.alloc(spec, Luna.realtype(spec), size(Ed))
        b = @benchmarkable (Nonlinear.batched!($R, $out, $Ed, $ρ, $sc); sync($spec)
                            ) seconds=BUDGET
        @printf("%-14s %10d %12s\n", name, ncols, fmttime(bminimum(brun(b)).time/1e9))
        flush(stdout)
    end
end
