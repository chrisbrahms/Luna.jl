module Ionisation
import SpecialFunctions: gamma, dawson
import HCubature: hquadrature
import HDF5
import FileWatching.Pidfile: mkpidlock
import HypergeometricFunctions: pFq
import Logging: @info
import Luna.PhysData: c, ħ, electron, m_e, au_energy, au_time, au_Efield, wlfreq, polarisability_difference, polarisability, au_polarisability
import Luna.PhysData: ionisation_potential, quantum_numbers
import Luna: Maths, Utils
import Luna
import Adapt
import Printf: @sprintf

abstract type AbstractIonRate end

"""
    IonRateADK(ionpot::Float64, threshold=true)
    IonRateADK(material::Symbol)

Ionisation rate based on the ADK formula. If `threshold` is true, use [`ADK_threshold`](@ref)
to avoid calculation below floating-point precision. If `cycle_average` is `true`, calculate
the cycle-averaged ADK ionisation rate instead.

The struct is parametric in the real type of its constants (`Float64` as built). Nothing
reachable from a device kernel may hold a `Float64` (GPU_PLAN.md §4.1) and the rate is
evaluated inside a broadcast which may be compiled for a GPU, so [`device_rate`](@ref)
makes the copy in the run's precision.
"""
struct IonRateADK{T} <: AbstractIonRate
    ionpot::T
    threshold::Bool
    cycle_average::Bool
    nstar::T
    cn_sq::T
    ω_p::T
    ω_t_prefac::T
    thr::T
    avfac::T
    occupancy::Int
end

function IonRateADK(material::Symbol; kwargs...)
    IonRateADK(ionisation_potential(material); kwargs...)
end

function IonRateADK(ionpot::Number; occupancy=2, threshold=true, cycle_average=false)
    nstar = sqrt(0.5/(ionpot/au_energy))
    cn_sq = 2^(2*nstar)/(nstar*gamma(nstar+1)*gamma(nstar))
    ω_p = ionpot/ħ
    ω_t_prefac = electron/sqrt(2*m_e*ionpot)

    if threshold
        thr = ADK_threshold(ionpot)
    else
        thr = 0.0
    end

    if cycle_average
        # Zenghu Chang: Fundamentals of Attosecond Optics (2011) p. 184
        # Section 4.2.3.1 Cycle-Averaged Rate
        # ̄w_ADK(Fₐ) = √(3/π) √(Fₐ/F₀) w_ADK(Fₐ) where Fₐ is the field amplitude
        Ip_au = ionpot / au_energy
        F0_au = (2Ip_au)^(3/2)
        F0 = F0_au*au_Efield
        avfac = sqrt.(3/(π*F0))
    else
        avfac = 1.0
    end

    #= One real type for every constant, so that the struct is `isbits` and a copy of it
       in Float32 (`device_rate`) has no Float64 field left. =#
    ip, ns, cn, wp, wt, th, av = promote(float(ionpot), nstar, cn_sq, ω_p, ω_t_prefac,
                                         thr, avfac)
    IonRateADK(ip, threshold, cycle_average, ns, cn, wp, wt, th, av, occupancy)
end

#= `-4/3` would be a Float64 literal and would promote a Float32 kernel; `-(4*one(T)/3)`
   is the same value in the type of the argument. Everything else is an `Int` times a
   field, which takes the field's type. The expression and its association are otherwise
   exactly as they were, so the Float64 values are unchanged. =#
function (ir::IonRateADK)(E)
    aE = abs(E)
    if aE >= ir.thr
        c43 = 4*one(aE)/3
        r = (ir.occupancy * ir.ω_p * ir.cn_sq *
             (4 * ir.ω_p / (ir.ω_t_prefac * aE))^(2 * ir.nstar - 1)
             * exp(-c43 * ir.ω_p / (ir.ω_t_prefac * aE)))
        if ir.avfac ≠ 1
            r *= ir.avfac*sqrt(aE)
        end
        return r
    else
        return zero(aE)
    end
end

#= The historical array form: a broadcast of the scalar rate, which is the same kernel
   body `ratekernel` wraps. It stays separate from `ionrate!` because `ionrate!`'s
   fallback calls this form for a rate with no kernel, and a rate whose array form
   delegated back to `ionrate!` would recurse. =#
function (ir::IonRateADK)(out::AbstractArray, E::AbstractArray)
    out .= ir.(E)
end

"""
    ionrate_ADK(IP_or_material, E::Number; kwargs...) -> Float64

Calculate the ionisation rate based on the ADK model.

# Arguments

- `IP_or_material`: Ionisation potential (`Float64`) or gas species (`Symbol`)
- `E::Number`: Electric field in SI units (V/M)

# Keywords

- `kwargs...`: See `IonRateADK`
"""
function ionrate_ADK(IP_or_material, E; kwargs...)
    IonRateADK(IP_or_material; kwargs...).(E)
end

"""
    ADK_threshold(ionpot)

Determine the lowest electric field strength at which the ADK ionisation rate for the
ionisation potential `ionpot` is non-zero to within 64-bit floating-point precision.
"""
function ADK_threshold(ionpot)
    ADKfun = IonRateADK(ionpot; threshold=false)
    E = 1e3
    out = ADKfun(E)
    while out == 0
        E *= 1.01
        out = ADKfun(E)
    end
    return E
end

"""
    IonRatePPT(ionpot::Float64, λ0, Z, l; kwargs...)

PPT ionisation rate for a state with ionisation potential `ionpot`, when driven at
wavelength `λ0`, also given charge state `Z` and angular momentum `l`

# Keyword arguments
- `sum_tol::Number`: Relative tolerance used to truncate the infinite sum. Defaults to 1e-6.
- `cycle_average::Bool`: If `true`, calculate the cycle-averaged rate. Defaults to `false`.
- `sum_integral::Bool`: whether to approximate the infinite sum in the PPT rate equation with
    an integral (this neglects the multiphoton thresholds).
- `Δα::Number`: polarisability difference between the ground state and the cation (in SI units)
    to calculate the Stark shift of the ground-state energy levels. Defaults to 0.
- `α_ion::Number`: polarisability of the cation (in SI units) to calculate the dipole correction
    to the rate. Defaults to 0.
- `msum::Bool`: for l ≠ 0, whether or not to sum over different m states. Defaults to `true`.
- `Cnl::Real` : Pre-calculated `Cₙₗ` constant. If not given, defaults to the approximate expression from
    the PPT papers.
- `occupancy`: Occupancy of the state(s) from which ionisation is considered. Defaults to 2 for
    a state with two electrons (spin up/down).

# References
[1] Ilkov, F. A., Decker, J. E. & Chin, S. L.
Ionization of atoms in the tunnelling regime with experimental evidence
using Hg atoms. Journal of Physics B: Atomic, Molecular and Optical
Physics 25, 4005–4020 (1992)

[2] Bergé, L., Skupin, S., Nuter, R., Kasparian, J. & Wolf, J.-P.
Ultrashort filaments of light in weakly ionized, optically transparent
media. Rep. Prog. Phys. 70, 1633–1713 (2007)
(Appendix A)

[3] A. Couairon and A. Mysyrowicz,
"Femtosecond filamentation in transparent media,"
Physics Reports 441(2–4), 47–189 (2007).

"""
struct IonRatePPT{oT, CT} <: AbstractIonRate
    Z::Float64 # Charge state
    l::Int # Orbital angular momentum quantum number
    Δα::Float64 # Polarisbility difference between ground state and ion
    α_ion_au::Float64 # Polarisbility of the ion in atomic units
    ω0_au::Float64 # central frequency in atomic units
    Cnl::CT # either missing (default) or a pre-defined value for Cₙₗ
    ionpot::Float64 # ionisation potential in SI units
    occupancy::oT # occupancy (integer) or function occ(m) returning occupancy in m level
    msum::Bool # sum over m levels?
    sum_integral::Bool # replace the infinite sum with an integral?
    sum_tol::Float64 # relative tolerance for convergence of the infinite sum
    cycle_average::Bool # average over one cycle?
end

"""
    IonRatePPT(material::Symbol, λ0; kwargs...)

PPT ionisation rate for the given `material` when driven at wavelength `λ0`.

# Keyword arguments
- `stark_shift::Bool`: whether to include the Stark shift
- `dipole_corr::Bool`: whether to include the dipole correction factor

Other keyword arguments are identical to `IonRatePPT(ionpot::Float64, λ0, Z, l; kwargs...)`
"""
function IonRatePPT(material::Symbol, λ0; stark_shift=true, dipole_corr=true, kwargs...)
    _, l, Z = quantum_numbers(material)
    Δα = stark_shift ? polarisability_difference(material) : 0.0
    α_ion = dipole_corr ? polarisability(material, true) : 0.0
    ip = ionisation_potential(material)
    IonRatePPT(ip, λ0, float(Z), l; Δα, α_ion, kwargs...)
end

function IonRatePPT(ip, λ0, Z, l; Δα=0, α_ion=0, sum_tol=1e-6,
    cycle_average=false, sum_integral=false, msum=true, Cnl=missing, occupancy=2)

    if ismissing(Δα)
        Δα = 0.0
    end

    if ismissing(α_ion)
        α_ion = 0.0
    end

    α_ion_au = α_ion/au_polarisability

    ω0 = 2π*c/λ0
    ω0_au = au_time*ω0

    IonRatePPT(float(Z), l, Δα, α_ion_au, ω0_au, Cnl, ip, occupancy,
        msum, sum_integral, sum_tol, cycle_average)
end

function (ir::IonRatePPT)(E)
    Ip_au = (ir.ionpot + ir.Δα/2 * E^2) / au_energy # Δα/2 * E^2 includes the Stark shift
    ns = ir.Z/sqrt(2Ip_au)
    ls = ns-1
    Cnl2 = ismissing(ir.Cnl) ? 2^(2ns)/(ns*gamma(ns + ls + 1)*gamma(ns - ls)) : ir.Cnl^2

    E0_au = (2*Ip_au)^(3/2)

    E_au = abs(E)/au_Efield
    γ = ir.ω0_au*sqrt(2Ip_au)/E_au
    γ2 = γ*γ
    β = 2γ/sqrt(1 + γ2)
    α = 2*(asinh(γ) - γ/sqrt(1+γ2))
    Up_au = E_au^2/(4*ir.ω0_au^2)
    Uit_au = Ip_au + Up_au
    v = Uit_au/ir.ω0_au
    ret = 0.0
    mrange = ir.msum ? (-ir.l:ir.l) : (0:0)
    for m in mrange
        mabs = abs(m)
        flm = ((2ir.l + 1)*factorial(ir.l + mabs)
            / (2^mabs*factorial(mabs)*factorial(ir.l - mabs)))
        # Following 5 lines are [1] eq. 8 and lead to identical results:
        # G = 3/(2γ)*((1 + 1/(2γ2))*asinh(γ) - sqrt(1 + γ2)/(2γ))
        # Am = 4/(sqrt(3π)*factorial(mabs))*γ2/(1 + γ2)
        # lret = sqrt(3/(2π))*Cnl2*flm*Ip_au
        # lret *= (2*E0_au/(E_au*sqrt(1 + γ2))) ^ (2ns - mabs - 3/2)
        # lret *= Am*exp(-2*E0_au*G/(3E_au))
        # [2] eq. (A14)
        lret = 4sqrt(2)/π*Cnl2
        lret *= (2*E0_au/(E_au*sqrt(1 + γ2))) ^ (2ns - mabs - 3/2)
        lret *= flm/factorial(mabs)
        lret *= exp(-2v*(asinh(γ) - γ*sqrt(1+γ2)/(1+2γ2)))
        lret *= Ip_au * γ2/(1+γ2)
        # Remove cycle average factor, see eq. (2) of [1]
        if !ir.cycle_average
            lret *= sqrt(π*E0_au/(3E_au))
        end
        n0 = ceil(v)
        if ir.sum_integral
            s = sqrt(π)*factorial(mabs)*β^mabs/(2*(α+β)^(mabs+1))*sqrt(β/α)
        else
            s, _, _ = Maths.converge_series(0, n0=n0, rtol=ir.sum_tol, maxiter=Inf) do x, n
                diff = n-v
                x + exp(-α*diff)*φ(m, sqrt(β*diff))
            end

        end
        lret *= s
        ret += occ(ir.occupancy, m)*lret
    end
    if ir.α_ion_au ≠ 0
        ret *= exp(-2*ir.α_ion_au*E_au)
    end
    return ret/au_time
end

occ(occupancy::Number, m) = occupancy
occ(occupancy, m) = occupancy(m)

"""
    φ(m, x)

Calculate the φ function for the PPT ionisation rate.

Note that w_m(x) in [1] and φ_m(x) in [2] look slightly different but
are in fact identical.
"""
function φ(m, x)
    #= second half of [3], eq. 81
        for m = 0, φ₀(x) is just the Dawson integral so we can get this directly.
        for m ≠ 0, we calculate it using the hypergeometric function where possible.
        for m ≠ 0 and large x, we need to do it brute force with BigFloats (slow)
    =#
    if m == 0
        return dawson(x)
    end

    if x <= 26
        mabs = abs(m)
        return (exp(-x^2)
            * sqrt(π)
            * x^(2mabs+1)
            * gamma(mabs+1)
            * pFq((1/2,), (3/2 + mabs,), x^2)
            / (2*gamma(3/2 + mabs)))
    else
        i, _ = hquadrature(0, x) do y
            y = BigFloat(y)
            x = BigFloat(x)
            (x^2 - y^2)^(abs(m))*exp(y^2)
        end
        return Float64(exp(-x^2) * i)
    end
end

function (ir::IonRatePPT)(out::AbstractArray, E::AbstractArray)
    out .= ir.(E)
end

function ionrate_PPT(ionpot, λ0, Z, l, E; kwargs...)
    return IonRatePPT(ionpot, λ0, Z, l; kwargs...).(E)
end

function ionrate_PPT(material::Symbol, λ0, E;
                     stark_shift=true, dipole_corr=true, kwargs...)
    _, l, Z = quantum_numbers(material)
    Δα = stark_shift ? polarisability_difference(material) : 0.0
    α_ion = dipole_corr ? polarisability(material, true) : 0.0
    ip = ionisation_potential(material)
    return ionrate_PPT(ip, λ0, Z, l, E; Δα, α_ion, kwargs...)
end

struct IonRatePPTAccel{ST, T} <: AbstractIonRate
    spline::ST # spline interpolant of log(rate)
    Emin::T # minimum electric field strength
    Emax::T # maximum electric field strength
end

"""
    IonRatePPTAccel(material::Symbol, λ0; kwargs...)
    IonRatePPTAccel(ionpot::Float64, λ0, Z, l; kwargs...)
    IonRatePPTAccel(E, rate)

Create a cached (saved) interpolated PPT ionisation rate function. If a saved lookup table
exists, load this rather than recalculate.

# Keyword arguments
- `N::Int`: Number of samples with which to create the `CSpline` interpolant.
- `Emax::Number`: Maximum field strength to include in the interpolant.
- `cache::Bool`: Whether to save the pre-calculated rate to a file
- `cachedir::String`: Path to the directory where the cache should be stored and loaded from.
    Defaults to \$HOME/.luna/pptcache

Other keyword arguments are passed on to [`IonRatePPT`](@ref)
"""
function IonRatePPTAccel(E, rate)
    # first remove points where the rate is zero within floating-point
    # precision to avoid NaNs in the CSpline
    idcs = rate .> 0
    E = E[idcs]
    rate = rate[idcs]
    # Interpolating the log and re-exponentiating makes the spline more accurate
    cspl = Maths.CSpline(E, log.(rate); bounds_error=true)
    Emin, Emax = promote(minimum(E), maximum(E))
    IonRatePPTAccel(cspl, Emin, Emax)
end

function IonRatePPTAccel(material::Symbol, λ0; stark_shift=true, dipole_corr=true, kwargs...)
    _, l, Z = quantum_numbers(material)
    Δα = stark_shift ? polarisability_difference(material) : 0.0
    α_ion = dipole_corr ? polarisability(material, true) : 0.0
    ip = ionisation_potential(material)
    IonRatePPTAccel(ip, λ0, Z, l; Δα, α_ion, kwargs...)
end

function IonRatePPTCached(args...; kwargs...)
    IonRatePPTAccel(args...; cache=true, kwargs...)
end

function IonRatePPTAccel(ionpot::Float64, λ0, Z, l;
    N=2^16, Emax=nothing, cache=true,
    cachedir=joinpath(Utils.cachedir(), "pptcache"),
    stale_age=60 * 10,
    kwargs...)
    h = hash((ionpot, λ0, Z, l, N, Emax, collect(kwargs)))
    fname = string(h, base=16) * ".h5"
    fpath = joinpath(cachedir, fname)
    if cache && isfile(fpath)
        lockpath = joinpath(cachedir, "pptlock")
        E, rate = mkpidlock(lockpath; stale_age) do
            @info @sprintf("Found cached PPT rate for %.2f eV, %.1f nm", ionpot / electron, 1e9λ0)
            HDF5.h5open(fpath, "r") do file
                (read(file["E"]), read(file["rate"]))
            end
        end
    else
        E, rate = makePPTcache(ionpot::Float64, λ0, Z, l;
            N, Emax, kwargs...)
    end

    if cache && ~isfile(fpath)
        lockpath = joinpath(cachedir, "pptlock")
        isdir(cachedir) || mkpath(cachedir)
        mkpidlock(lockpath; stale_age) do
            if ~isfile(fpath) # makePPTcache takes a while - has another process saved first?
                @info @sprintf(
                    "Saving PPT rate for %.2f eV, %.1f nm in %s",
                    ionpot / electron, 1e9λ0, fpath
                )
                HDF5.h5open(fpath, "cw") do file
                    file["E"] = E
                    file["rate"] = rate
                end
            end
        end
    end

    return IonRatePPTAccel(E, rate)
end

function (ir::IonRatePPTAccel)(E)
    aE = abs(E)
    if aE > ir.Emax
        error(
            "Field strength $aE V/m exceeds maximum for PPT ionisation rate ($(ir.Emax) V/m)."
            )
    end
    _pptaccel(ir, E)
end

#= The kernel: the same arithmetic without the error path, whose string interpolation
   cannot be compiled for a device and which a kernel could not report anyway. Above the
   table it returns the table's last value (`min` puts `aE` on the last knot, where the
   spline returns `y[end]` exactly); the host path errors instead, through the
   `maximum(abs, E)` check in `ionrate!`. GPU_PLAN.md §4.3. =#
@inline function _pptaccel(ir::IonRatePPTAccel, E)
    aE = abs(E)
    aE < ir.Emin && return zero(aE)
    exp(Maths.spline_eval(ir.spline, min(aE, ir.Emax)))
end

# See the note on `(::IonRateADK)(out, E)`.
function (ir::IonRatePPTAccel)(out::AbstractArray, E::AbstractArray)
    out .= ir.(E)
end

#=================================================#
#===========  ARRAY-LEVEL EVALUATION  ============#
#=================================================#

"""
    ratekernel(ir, Eref)

A callable `f(e)` giving the ionisation rate of the field `Eref*e`, in the element type
of `Eref`. This is the body of a broadcast which may be compiled for a GPU, so it
captures nothing but the rate object (which must be `isbits` in the run's precision, see
[`device_rate`](@ref)) and `Eref`.

`Eref` is the unit the field array is expressed in (see [`Luna.UnitScaling`](@ref)): a
scaled state holds `e = E/E_ref`, and the rate is not polynomial in the field, so the
physical field is reconstructed here rather than folded into a coefficient. `Eref == 1`
for every `Float64` run, where `1*e` is exact and the arithmetic is unchanged.
"""
function ratekernel end

ratekernel(ir::IonRateADK, Eref) = let ir=ir, Eref=Eref
    e -> ir(Eref*e)
end

ratekernel(ir::IonRatePPTAccel, Eref) = let ir=ir, Eref=Eref
    e -> _pptaccel(ir, Eref*e)
end

"""
    ionrate!(out, ir, E, Eref=1)

Ionisation rate of every element of the field array `E`, placed into `out`.

`E` holds `E_phys/Eref` (see [`Luna.UnitScaling`](@ref)); `Eref` defaults to `1`, i.e.
physical units. `out` and `E` must have the same shape.

**Which path it takes is decided by [`device_capable`](@ref), not by the type.** A rate
with a kernel is one broadcast of [`ratekernel`](@ref), on the host, on a GPU, at any
shape. Anything else — the direct [`IonRatePPT`](@ref), a cached rate on a non-uniform
table, a rate somebody wrote — is called as `ir(out, E)`, which is what it always was and
which needs host arrays in physical units; a device or a scaled run is refused with a
message naming the alternatives.

`check=false` skips the range check for a caller which has already made it on the whole
block (see [`check_field_range`](@ref)). The check is what raises the error for a field
above a cached rate's table, so skipping it without making it elsewhere would leave the
rate saturating silently.
"""
function ionrate!(out, ir, E, Eref=1; check=true)
    if device_capable(ir)
        check && check_field_range(ir, E, Eref)
        f = ratekernel(ir, Luna.scalar(out, Eref))
        out .= f.(E)
    else
        (Eref == 1 && !Utils.isdevice(E)) || error(
            "the ionisation rate $(nameof(typeof(ir))) has no device kernel, so it can "*
            "only be evaluated on host arrays in physical units. Use "*
            "`Ionisation.IonRateADK` or a cached PPT rate (`IonRatePPTCached`) on a "*
            "device or in reduced precision, or run on the CPU with `device=:cpu`.")
        ir(out, E)
    end
    out
end

"""
    check_field_range(ir, E, Eref)

Raise if any element of `Eref*E` is outside the range the rate `ir` can be evaluated
over. Only a cached PPT rate has one: its spline is built on a table which stops at twice
the barrier-suppression field.

This is one reduction over the whole array rather than a branch per element, so a caller
which splits an array into columns makes the check once, on the whole thing, and passes
`check=false` to [`ionrate!`](@ref) — both so that the reduction happens once and so that
the error is raised from the calling task rather than from inside a `@threads` loop.

A no-op on a device array: a device kernel cannot raise, and there the rate saturates at
the table's last value instead (see [`IonRatePPTAccel`](@ref)).
"""
check_field_range(ir, E, Eref) = nothing

function check_field_range(ir::IonRatePPTAccel, E, Eref)
    Utils.isdevice(E) && return nothing
    m = maximum(abs, E)*Eref
    m > ir.Emax && error(
        "Field strength $m V/m exceeds maximum for PPT ionisation rate ($(ir.Emax) V/m).")
    nothing
end

#=================================================#
#===========  DEVICE (GPU) EVALUATION  ===========#
#=================================================#

"""
    device_capable(ir) -> Bool

Whether the ionisation rate `ir` can be evaluated inside a device (GPU) kernel and in
reduced precision: the analytic ADK rate, and a cached PPT rate whose table is uniformly
spaced (which every table `makePPTcache` builds is).

`false` for the direct [`IonRatePPT`](@ref), whose series summation, `BigFloat` fallback
and `factorial`s cannot be compiled for a device, for a cached rate which ended up on a
`Maths.FastFinder`, and for a user-supplied callable.
"""
device_capable(ir) = false
device_capable(::IonRateADK) = true
device_capable(ir::IonRatePPTAccel) = ir.spline.ifun isa Maths.UniformIndex

"""
    device_rate(ir, spec)

The same ionisation rate with its constants in the precision of `spec` and its lookup
tables on its array type (see [`Luna.DeviceSpec`](@ref)), ready to be captured by a
broadcast kernel. Returns `ir` itself for the default host `Float64` run.

Errors, naming the alternatives, for a rate which is not
[`device_capable`](@ref).
"""
function device_rate(ir, spec)
    (Luna.isdevicespec(spec) || Luna.realtype(spec) !== Float64) || return ir
    device_capable(ir) || error(
        "the ionisation rate $(nameof(typeof(ir))) cannot be evaluated on "*
        "$(Luna.arraytype(spec)) in $(Luna.realtype(spec)): its rate function has no "*
        "device kernel. Use `Ionisation.IonRateADK`, or a cached PPT rate "*
        "(`IonRatePPTCached`/`IonRatePPTAccel`, whose table is uniformly spaced), or "*
        "run on the CPU with `device=:cpu`.")
    _device_rate(ir, spec)
end

_device_rate(ir::IonRateADK, spec) =
    IonRateADK{Luna.realtype(spec)}(ir.ionpot, ir.threshold, ir.cycle_average, ir.nstar,
                                    ir.cn_sq, ir.ω_p, ir.ω_t_prefac, ir.thr, ir.avfac,
                                    ir.occupancy)

function _device_rate(ir::IonRatePPTAccel, spec)
    T = Luna.realtype(spec)
    IonRatePPTAccel(Maths.todevice_spline(spec, ir.spline),
                    convert(T, ir.Emin), convert(T, ir.Emax))
end

"""
    resident_arrays(ir)

The arrays an ionisation rate carries which a device kernel indexes, for the residency
assertion of the response which holds it
(see [`Nonlinear.resident_arrays`](@ref Luna.Nonlinear.resident_arrays)).
"""
resident_arrays(ir) = ()
resident_arrays(ir::IonRatePPTAccel) = (ir.spline.x, ir.spline.y, ir.spline.D)

#= Structural moves for the kernel adaptor: when a broadcast kernel which captured one of
   these is compiled for a device, every array inside it has to become a device pointer.
   The scalars are already in the right precision (`device_rate`). =#
Adapt.adapt_structure(to, ir::IonRateADK) = ir
Adapt.adapt_structure(to, ir::IonRatePPTAccel) =
    IonRatePPTAccel(Adapt.adapt(to, ir.spline), ir.Emin, ir.Emax)

function makePPTcache(ionpot::Float64, λ0, Z, l;
                      N=2^16, Emax=nothing, kwargs...)
    Emax = isnothing(Emax) ? 2*barrier_suppression(ionpot, Z) : Emax

    # ω0 = 2π*c/λ0
    # Emin = ω0*sqrt(2m_e*ionpot)/electron/0.5 # Keldysh parameter of 0.5
    Emin = Emax/5000

    E = collect(range(Emin, stop=Emax, length=N));
    @info @sprintf("Pre-calculating PPT rate for %.2f eV, %.1f nm...", ionpot/electron, 1e9λ0)
    flush(stderr) # pre-calculating can take a while, so make sure this message is shown
    rate = ionrate_PPT(ionpot, λ0, Z, l, E; kwargs...)
    @info "...PPT pre-calcuation done"
    flush(stderr)
    return E, rate
end

"""
    barrier_suppression(ionpot, Z)

Calculate the barrier-suppresion **field strength** for the ionisation potential `ionpot`
and charge state `Z`.
"""
function barrier_suppression(ionpot, Z)
    Ip_au = ionpot / au_energy
    ns = Z/sqrt(2*Ip_au)
    Z^3/(16*ns^4) * au_Efield
end

"""
    keldysh(material, λ, E)

Calculate the Keldysh parameter for the given `material` at wavelength `λ` and electric field
strength `E`.
"""
function keldysh(material, λ, E)
    Ip_au = ionisation_potential(material)/au_energy
    E_au = E/au_Efield
    ω0_au = wlfreq(λ)*au_time
    ω0_au*sqrt(2Ip_au)/E_au
end

"""
    ionfrac(rate, E, δt)

Given an ionisation rate function `rate` and an electric field array `E` sampled with time
spacing `δt`, calculate the ionisation fraction as a function of time on the same time axis.

The function `rate` should have the signature `rate!(out, E)` and place its results into
`out`, like the functions returned by e.g. `IonRateADK` or `IonRatePPTCached`.
"""
function ionfrac(rate, E, δt)
    frac = similar(E)
    ionfrac!(frac, rate, E, δt)
end

function ionfrac!(frac, rate, E, δt)
    rate(frac, E)
    Maths.cumtrapz!(frac, δt)
    @. frac = 1 - exp(-frac)
end

end
