#= Multimode cases for the fixed-versus-adaptive transverse-integral study (`run.jl`).

   Every case is built with the low-level interface, so that the adaptive rule's `mfcn`
   (the cap on transverse points per right-hand side, 512 by default, which
   `prop_capillary` does not expose) can be lifted for the references. All are
   field-resolved with Kerr and PPT plasma. Strengths are set against the critical power
   for self-focusing (`Tools.Pcr`).

   - `he1m_strong`: 8 HE₁ₘ modes, 125 µm core radius, 1 bar Ar, 800 nm, 30 fs at 0.95 P_cr,
     30 cm; radial integral. The hard case of `benchmark/threads/modal_accuracy.jl`.
   - `he1m_weak`: 4 HE₁ₘ modes at 0.04 P_cr, otherwise the same; nothing spatial happens.
   - `vortex`: after `LG_soliton_1030nm.jl`: HE₂₁ at ϕ = 0 and π/4 in quadrature (a
     circularly polarised vortex), 1030 nm, 12 fs, 250 µJ, 0.4 bar Ar, 125 µm, 60 cm
     (soliton self-compression); the mode set adds the higher-order modes of the same
     azimuthal order (`VORTEX_M`) that the nonlinearity couples into. Full polar integral.
   - `mixed`: HE₁₁ with an HE₂₁ seed (see [`mixed`](@ref)): the case that tests the θ rule.
   - `rect`: a silver-clad rectangular guide (half-widths 150 × 75 µm), 1 bar Ar, 800 nm,
     30 fs at `RECT_FRAC` P_cr, modes n ∈ (1, 3, 5), m ∈ (1, 3) polarised along x -- the
     modes of the input's symmetry; the others stay empty. Cartesian integral.

   A run is `case(rule; scale, prtol)`: `rule` is `(:adaptive, rtol, mfcn)` or
   `(:fixed, nr, nθ, kronrod)`, `prtol` the propagation tolerance of `Luna.run`, `scale`
   multiplies the length. Returns `(; setup_s, run_s, steps, output)`. =#
import Luna: Grid, Fields, LinearOps, Nonlinear, NonlinearRHS, Output, PhysData, Stats,
             Ionisation, Tools, Capillary, RectModes
import Logging

const NSAVE = parse(Int, get(ENV, "NSAVE", "31"))

function pcrit(gas, pres, λ0)
    ω = 2π*PhysData.c/λ0
    _, n0, n2 = Tools.getN0n0n2(ω, gas; P=pres)
    Tools.Pcr(ω, n0, n2)
end

"Energy of a Gaussian pulse of FWHM `τ` with peak power `frac` P_cr."
energy_at(frac, gas, pres, λ0, τ) = frac*pcrit(gas, pres, λ0)*τ/0.94

#= Forwards warnings to the console and picks the step count out of the integrator's
   "Propagation finished in ..., N steps" message (all attempted steps). =#
struct StepLogger <: Logging.AbstractLogger
    steps::Base.RefValue{Int}
    inner::Logging.ConsoleLogger
end
StepLogger() = StepLogger(Ref(0), Logging.ConsoleLogger(stderr, Logging.Warn))
Logging.min_enabled_level(::StepLogger) = Logging.Info
Logging.shouldlog(::StepLogger, args...) = true
Logging.catch_exceptions(::StepLogger) = false
function Logging.handle_message(l::StepLogger, level, msg, args...; kw...)
    m = match(r"^Propagation finished in .*, (\d+) steps", string(msg))
    isnothing(m) || (l.steps[] = parse(Int, m[1]))
    level >= Logging.Warn && Logging.handle_message(l.inner, level, msg, args...; kw...)
    nothing
end

function responses(grid, gas, pol)
    ionrate = Ionisation.IonRatePPTCached(gas, grid.referenceλ)
    Et = pol ? Array{Float64}(undef, length(grid.to), 2) : grid.to
    (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
     Nonlinear.PlasmaCumtrapz(grid.to, Et, ionrate, PhysData.ionisation_potential(gas)))
end

rulekw(r) = r[1] == :adaptive ? (modal_integral=:adaptive, rtol=r[2], mfcn=r[3]) :
            (modal_integral=:fixed, nr=r[2], nθ=r[3], kronrod=r[4])

function modalrun(grid, gas, pres, modes, components, full, inputs, L, rule; prtol)
    l = StepLogger()
    r = Logging.with_logger(l) do
        dens0 = PhysData.density(gas, pres)
        densityfun = z -> dens0
        resp = responses(grid, gas, components == :xy)
        ts = @elapsed begin
            Eω, transform, FT = Luna.setup(grid, densityfun, resp, inputs, modes,
                                           components; full, rulekw(rule)...)
            linop = LinearOps.make_const_linop(grid, modes, grid.referenceλ)
            statsfun = Stats.default(grid, Eω, modes, linop, transform; gas)
            output = Output.MemoryOutput(0, L, NSAVE, statsfun)
        end
        tr = @elapsed Luna.run(Eω, grid, linop, transform, FT, output; zmax=L, rtol=prtol,
                               status_period=1e9)
        (; setup_s=ts, run_s=tr, output)
    end
    (; r..., steps=l.steps[])
end

grid_he1m() = Grid.RealGrid(800e-9, (150e-9, 4e-6), 0.5e-12)
grid_vortex() = Grid.RealGrid(1030e-9, (160e-9, 4000e-9), 0.4e-12)
grid_rect() = Grid.RealGrid(800e-9, (200e-9, 4e-6), 0.5e-12)
"The time-frequency grid of each case, for the analysis scripts."
const GRIDS = Dict("he1m_strong" => grid_he1m, "he1m_weak" => grid_he1m,
                   "vortex" => grid_vortex, "mixed" => grid_vortex, "rect" => grid_rect)

gauss(λ0, τ, energy; kw...) = (Fields.GaussField(; λ0, τfwhm=τ, energy, kw...),)

function he1m(nm, frac, rule; scale=1.0, prtol=1e-6)
    gas, pres, λ0, τ, a, L = :Ar, 1.0, 800e-9, 30e-15, 125e-6, 0.3scale
    grid = grid_he1m()
    modes = Tuple(Capillary.MarcatiliMode(a, gas, pres; n=1, m) for m in 1:nm)
    inputs = ((mode=1, fields=gauss(λ0, τ, energy_at(frac, gas, pres, λ0, τ))),)
    modalrun(grid, gas, pres, modes, :y, false, inputs, L, rule; prtol)
end
he1m_strong(rule; kw...) = he1m(8, 0.95, rule; kw...)
he1m_weak(rule; kw...) = he1m(4, 0.04, rule; kw...)

"Radial orders of the HE₂ₘ pairs in the vortex case."
const VORTEX_M = parse.(Int, split(get(ENV, "VORTEX_M", "1,2,3"), ","))
const VORTEX_ENERGY = parse(Float64, get(ENV, "VORTEX_ENERGY", "250e-6"))

function vortex(rule; scale=1.0, prtol=1e-6)
    gas, pres, λ0, τ, a, L = :Ar, 0.4, 1030e-9, 12e-15, 125e-6, 0.6scale
    grid = grid_vortex()
    modes = Tuple(Capillary.MarcatiliMode(a, gas, pres; kind=:HE, n=2, m, ϕ)
                  for m in VORTEX_M for ϕ in (0.0, π/4))
    E = VORTEX_ENERGY
    inputs = ((mode=1, fields=gauss(λ0, τ, E/2; ϕ=[-π/2])),
              (mode=2, fields=gauss(λ0, τ, E/2)))
    modalrun(grid, gas, pres, modes, :xy, true, inputs, L, rule; prtol)
end

"""
Mixed azimuthal orders: HE₁₁ (ϕ = 0, polarised along y) carrying 90 % of the energy with an
in-phase HE₂₁ (ϕ = 0) seed carrying 10 %, in the vortex case's capillary, gas and pulse,
over 30 cm, with HE₁₂, both HE₂₁, TE₀₁, TM₀₁ and both HE₃₁ in the set (azimuthal orders
0, 1 and 2, so the θ rule is exact only for nθ ≥ 4·2+1 = 9 and only for the Kerr part).
The two input modes interfere into a θ-dependent intensity, unlike any combination of the
HE₂ₘ modes alone (whose |E|² is θ-independent, so the vortex case tests only r).
"""
function mixed(rule; scale=1.0, prtol=1e-6)
    gas, pres, λ0, τ, a, L = :Ar, 0.4, 1030e-9, 12e-15, 125e-6, 0.3scale
    grid = grid_vortex()
    M(kind, n, m, ϕ=0.0) = Capillary.MarcatiliMode(a, gas, pres; kind, n, m, ϕ)
    modes = (M(:HE, 1, 1), M(:HE, 1, 2), M(:HE, 2, 1), M(:HE, 2, 1, π/4), M(:TE, 0, 1),
             M(:TM, 0, 1), M(:HE, 3, 1), M(:HE, 3, 1, π/6))
    E = VORTEX_ENERGY
    inputs = ((mode=1, fields=gauss(λ0, τ, 0.9E)), (mode=3, fields=gauss(λ0, τ, 0.1E)))
    modalrun(grid, gas, pres, modes, :xy, true, inputs, L, rule; prtol)
end

const RECT_FRAC = parse(Float64, get(ENV, "RECT_FRAC", "0.9"))
const RECT_LENGTH = parse(Float64, get(ENV, "RECT_LENGTH", "0.1"))

function rect(rule; scale=1.0, prtol=1e-6)
    gas, pres, λ0, τ, L = :Ar, 1.0, 800e-9, 30e-15, RECT_LENGTH*scale
    grid = grid_rect()
    modes = Tuple(RectModes.RectMode(150e-6, 75e-6, gas, pres, :Ag; n, m, pol=:x)
                  for m in (1, 3) for n in (1, 3, 5))
    inputs = ((mode=1, fields=gauss(λ0, τ, energy_at(RECT_FRAC, gas, pres, λ0, τ))),)
    modalrun(grid, gas, pres, modes, :x, true, inputs, L, rule; prtol)
end

const CASES = Dict("he1m_strong" => he1m_strong, "he1m_weak" => he1m_weak,
                   "vortex" => vortex, "mixed" => mixed, "rect" => rect)

"""
    parserule(s)

`A1e-3m512` is the adaptive rule with rtol 1e-3 and mfcn 512; `F128x16` the fixed rule with
nr = 128 and nθ = 16, `F129x16k` with the Kronrod embedded estimate.
"""
function parserule(s)
    if s[1] == 'A'
        m = match(r"^A([0-9.e-]+)m(\d+)$", s)
        (:adaptive, parse(Float64, m[1]), parse(Int, m[2]))
    else
        m = match(r"^F(\d+)x(\d+)(k?)$", s)
        (:fixed, parse(Int, m[1]), parse(Int, m[2]), m[3] == "k")
    end
end
