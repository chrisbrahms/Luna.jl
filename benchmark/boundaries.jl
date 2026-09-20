#= Mode-averaged Kerr propagation with boundary=:rate and default statistics, timed on
   each device and precision.

   Usage (CPU only, from the repository root):

       julia --project=benchmark -t 1 benchmark/boundaries.jl

   To include Metal, run it from an environment which has Luna, BenchmarkTools and Metal,
   with Metal loaded:

       julia --project=<metalenv> -t 1 -e 'using Metal; include("benchmark/boundaries.jl")'

   `benchmark/device.jl` times the same propagation with `boundary=:none` and no
   statistics -- the RHS and step in isolation. This one adds what `gpu/11` made
   device-capable: `Boundaries.RateAbsorber` (a broadcast/reduction per accepted step) and
   the default statistics (`Stats.default`, still host code -- copied down every accepted
   step, gpu/24's to fix). The comparison that matters is the *overhead* boundary=:rate
   and the statistics add over the bare RHS/step numbers `device.jl` reports, on each
   backend.

   One Julia thread, one FFTW thread, one BLAS thread, `:estimate` planning, no wisdom.
=#
using Luna
import Luna: RK45, Grid, Modes, Capillary, Fields, LinearOps, Nonlinear, Output, Stats,
             PhysData, Utils, DeviceSpec, HostSpec
import LinearAlgebra
import Printf: @printf, @sprintf
import Logging

Luna.set_fftw_mode(:estimate)
Luna.set_fftw_threads(1)
LinearAlgebra.BLAS.set_num_threads(1)
Luna.set_fftw_wisdom(false)

"Number of accepted steps in the timed propagation (boundary=:rate needs at least this
many for max_dz to reach its reference-length cap; see Boundaries.jl)."
const NSTEPS = 20

"Time windows to sweep, in seconds. The grid size follows from the frequency limits."
const TRANGES = (400e-15, 1600e-15, 6400e-15)

const GAS = :He
const PRES = 1.0
const λ0 = 800e-9
const FLENGTH = 1e-2

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

"""
    prepare(spec, trange; boundary, stats)

Set up the mode-averaged Kerr propagation on `spec` and run it once, returning
`(n, elapsed)`: the state size and the wall-clock time of `Luna.run` (`boundary` and
`stats` as given; `stats=true` uses `Stats.default`, `false` uses `Output.nostats`).
"""
function proptime(spec, trange; boundary, stats, samples=3)
    quiet() do
        minimum(1:samples) do _
            grid = Grid.RealGrid(λ0, (300e-9, 2000e-9), trange)
            m = Capillary.MarcatiliMode(75e-6, GAS, PRES, loss=false)
            aeff(z) = Modes.Aeff(m, z=z)
            ρ = PhysData.density(GAS, PRES)
            dens = z -> ρ
            resp = (Nonlinear.Kerr_field(PhysData.γ3_gas(GAS)),)
            linop, βfun!, _, _ = LinearOps.make_const_linop(grid, m, λ0)
            inputs = Fields.GaussField(λ0=λ0, τfwhm=20e-15, energy=1e-6)
            Eω, transform, FT = Luna.setup(grid, dens, resp, inputs, βfun!, aeff;
                                           constβ=true, device=spec)
            n = length(Eω)
            statsfun = stats ? Stats.default(grid, Eω, m, linop, transform; gas=GAS) :
                               Output.nostats
            out = Output.MemoryOutput(0, FLENGTH, 3, statsfun)
            dz = FLENGTH/NSTEPS
            t = @elapsed Luna.run(Eω, grid, linop, transform, FT, out;
                                  zmax=FLENGTH, boundary, init_dz=dz, min_dz=dz, max_dz=dz)
            Luna.device_synchronize(spec)
            (n, t)
        end
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

@printf("Mode-averaged Kerr, %s at %g bar, %g m, %d fixed steps\n", GAS, PRES, FLENGTH, NSTEPS)
@printf("1 Julia thread, 1 FFTW thread, 1 BLAS thread, :estimate, no wisdom\n\n")
@printf("%-14s %8s %10s %14s %14s %14s\n",
        "device", "trange", "state", "none/nostats", "rate/nostats", "rate/stats")
@printf("%s\n", "-"^90)
for trange in TRANGES, (name, spec) in specs
    n, t_none = proptime(spec, trange; boundary=:none, stats=false)
    _, t_rate = proptime(spec, trange; boundary=:rate, stats=false)
    _, t_stats = proptime(spec, trange; boundary=:rate, stats=true)
    @printf("%-14s %6.0f fs %10d %14s %14s %14s\n",
            name, trange*1e15, n, fmttime(t_none), fmttime(t_rate), fmttime(t_stats))
    flush(stdout)
end
