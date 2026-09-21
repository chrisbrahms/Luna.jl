#= The per-step cost of the default statistics, on each device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/stats.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/stats.jl")'

   Metal is never a dependency of Luna or of this environment; the script uses whatever is
   already loaded and skips the device it cannot see.

   What is timed is one call of the statistics function, which is what `Luna.run` does on
   every accepted step through the output handler. Two columns:

   - "host": the state copied to the host and unscaled first, then the host branches of
     `Stats.jl` -- what every device run did before `gpu/24`, and what a statistics set
     containing a user-written closure still does. The copy is included, because it is
     part of the cost.
   - "device": the statistics evaluated on the state where it is, which is what the
     default sets do now.

   The "step" column is one RK45 step of the same propagation, so that the overhead can
   be read as a fraction of a step rather than in the abstract.

   As for `run.jl` and `device.jl`: one Julia thread, one FFTW thread, one BLAS thread,
   `:estimate` planning and no wisdom, or the numbers are not comparable between runs.
=#
using Luna
import Luna: RK45, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear, Output,
             PhysData, Utils, Stats, Ionisation, DeviceSpec, HostSpec
import LinearAlgebra
import BenchmarkTools: @benchmarkable, run as brun, minimum as bminimum
import Printf: @printf, @sprintf
import Logging

Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Luna.set_fftw_wisdom(false)

"Time budget in seconds for each benchmark."
const BUDGET = 2.0

"Time windows to sweep, in seconds. The grid size follows from the frequency limits."
const TRANGES = (400e-15, 1600e-15, 6400e-15)

const GAS = :Ar
const PRES = 1.0
const λ0 = 800e-9
const FLENGTH = 1e-2

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

#= The tabulated rate a `plasma=true` run uses, as in `device.jl`: the axis is uniform, so
   the kernel is the spline lookup a real cached rate uses, and it costs milliseconds. =#
function tablerate(gas)
    Ebs = Ionisation.barrier_suppression(PhysData.ionisation_potential(gas), 1.0)
    E = collect(range(2Ebs/5000, 2Ebs, length=1<<16))
    Ionisation.IonRatePPTAccel(E, Ionisation.IonRateADK(gas).(E))
end

"""
    prepare(spec, trange; plasma=false)

A mode-averaged propagation on `spec` and the default statistics built for it. With
`plasma=true` the response set is Kerr *and* `Nonlinear.PlasmaCumtrapz`, which puts
`Stats.electrondensity` into the statistics set.
"""
function prepare(spec, trange; plasma=false)
    quiet() do
        grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), trange)
        m = Capillary.MarcatiliMode(75e-6, GAS, PRES, loss=false)
        aeff(z) = Modes.Aeff(m, z=z)
        ρ = PhysData.density(GAS, PRES)
        dens = z -> ρ
        resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS)),)
        if plasma
            resp = (resp...,
                    Nonlinear.PlasmaCumtrapz(grid.to, zeros(length(grid.to)),
                                             tablerate(GAS),
                                             PhysData.ionisation_potential(GAS)))
        end
        linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0)
        inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=150e-6)
        Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                       constβ=true, device=spec)
        sf = Stats.default(grid, Eω, m, linop, transform; gas=GAS)
        Eref = Luna.runscaling(transform).Eref
        # the host path: the same set built for a host, physical-unit state
        Eh = Array(Eω).*Eref
        sfh = Stats.default(grid, Eh, m, linop, transform; gas=GAS)
        (Eω, Luna.upload_like(Eω, linop), transform, sf, sfh, Eh, Eref)
    end
end

#= Device work is queued asynchronously, so a timing which does not synchronise measures
   the launch, not the kernel. `device_synchronize` is a no-op on the CPU. =#
sync(spec) = Luna.device_synchronize(spec)

"The statistics evaluated on the state where it lives."
function devicetime(spec, p)
    Eω, _, _, sf, _, _, _ = p
    b = @benchmarkable (($sf)($Eω, 0.0, 1e-4); sync($spec)) seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

"The state copied to the host and unscaled, then the host branches -- the `gpu/int-D` path."
function hosttime(spec, p)
    Eω, _, _, _, sfh, Eh, Eref = p
    f = let Eh=Eh, Eω=Eω, Eref=Eref, sfh=sfh
        function ()
            copyto!(Eh, Eω)
            Eref == 1.0 || (Eh .*= Eref)
            sfh(Eh, 0.0, 1e-4)
        end
    end
    b = @benchmarkable ($f(); sync($spec)) seconds=BUDGET
    bminimum(brun(b)).time/1e9
end

"One RK45 step of the same propagation, for scale."
function steptime(spec, p)
    Eω, linop, transform, _, _, _, _ = p
    dz = FLENGTH/20
    b = @benchmarkable (RK45.step!(s); sync($spec)) setup=(
            s = RK45.PreconStepper($transform, $linop, $Eω, 0.0, $dz;
                                   max_dt=$dz, min_dt=$dz)
        ) evals=1 seconds=BUDGET
    bminimum(brun(b)).time/1e9
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

@printf("default statistics, %s at %g bar, one call (one accepted step)\n", GAS, PRES)
@printf("%d Julia threads, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n",
        Threads.nthreads())
for (label, rkw) in (("Kerr", (;)), ("Kerr + plasma (tabulated rate)", (; plasma=true)))
    @printf("\nresponses: %s\n\n", label)
    @printf("%-14s %8s %10s %12s %12s %12s %8s\n",
            "device", "trange", "state", "host+copy", "device", "step", "dev/step")
    @printf("%s\n", "-"^82)
    for trange in TRANGES, (name, spec) in specs
        p = prepare(spec, trange; rkw...)
        th = hosttime(spec, p)
        td = devicetime(spec, p)
        ts = steptime(spec, p)
        @printf("%-14s %6.0f fs %10d %12s %12s %12s %7.1f%%\n",
                name, trange*1e15, length(p[1]), fmttime(th), fmttime(td), fmttime(ts),
                100*td/ts)
        flush(stdout)
    end
end
