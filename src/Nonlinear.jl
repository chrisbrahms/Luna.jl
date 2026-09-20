module Nonlinear
import Luna
import Luna.PhysData: ε_0, e_ratio
import Luna: Maths, Utils
import Adapt
import FFTW
import LinearAlgebra: mul!, ldiv!
import Rotations: RotZY, RotYZ, RotMatrix, RotMatrix3
import StaticArrays: SMatrix, MArray

"""
    rescale(response, spec, scaling)

The same nonlinear response, with its coefficients expressed in the units of `scaling`
(see [`Luna.UnitScaling`](@ref)) and every array and constant it carries converted to
`spec`'s array type and precision (see [`Luna.DeviceSpec`](@ref)).

A transform calls this on each of its responses at construction, so that the powers of
`E_ref` are combined with the physical constants once, on the host, in `Float64`. There
is one kernel body per response; the units are a property of the constants it holds.

The fallback passes the response through unchanged for an unscaled `Float64` run -- which
is every run on the default CPU path, so an ad hoc response written as a closure keeps
working -- and errors otherwise. Responses gain their own methods as they are made
device-capable.
"""
function rescale(r, spec, scaling)
    (Luna.realtype(spec) === Float64 && Luna.isunity(scaling)) && return r
    error("the nonlinear response $(typeof(r)) has no `Nonlinear.rescale` method, so it "*
          "cannot be used in a reduced-precision or device run. On the default CPU path "*
          "(Float64, unscaled) any callable `resp!(out, E, ρ)` works.")
end

"""
    resident_arrays(response) -> Tuple

The arrays a nonlinear response carries which a kernel broadcasts against, and which
therefore have to live on the run's array type and precision. A transform passes them to
[`Luna.assert_resident`](@ref) at construction, together with its own buffers and grid
mirrors.

Empty by default, which is right for a response whose coefficients are all scalars and
for anything which does not run on a device at all.
"""
resident_arrays(r) = ()

#= The Kerr responses are structs rather than closures so that they can be parametric in
   the real element type (nothing reachable from a Metal kernel may hold a Float64), carry
   an `Adapt` rule for their arrays, and take a `rescale` method. The constructors
   `Kerr_field(γ3)` etc. keep their signatures and return the structs. =#

"""
    KerrField(γ3)

Kerr response for a real (field-resolved) field; built by `Kerr_field(γ3)`. In physical
units `γ3` is the third-order hyperpolarisability, e.g.
[`PhysData.γ3_gas`](@ref Luna.PhysData.γ3_gas); after [`rescale`](@ref) it carries the
powers of `E_ref` as well.
"""
struct KerrField{T}
    γ3::T
end

"Kerr response for real field"
Kerr_field(γ3) = KerrField(γ3)

function (k::KerrField)(out, E, ρ)
    fac = Luna.scalar(E, ρ*ε_0*k.γ3)
    if size(E, 2) == 1
        KerrScalar!(out, E, fac)
    else
        KerrVector!(out, E, fac)
    end
end

rescale(k::KerrField, spec, scaling) =
    KerrField(convert(Luna.realtype(spec), k.γ3*scaling.Eref^2/scaling.Pref))

function KerrScalar!(out, E, fac)
    @. out += fac*E^3
end

#= One broadcast per polarisation component, over views of the two columns. Each
   component's expression is the one the scalar loop evaluated, in the same order. =#
function KerrVector!(out, E, fac)
    Ex = view(E, :, 1)
    Ey = view(E, :, 2)
    ox = view(out, :, 1)
    oy = view(out, :, 2)
    @. ox += fac*(Ex^2 + Ey^2)*Ex
    @. oy += fac*(Ex^2 + Ey^2)*Ey
    out
end

"Kerr response for real field but without THG"
function Kerr_field_nothg(γ3, n)
    E = Array{Float64}(undef, n)
    hilbert = Maths.plan_hilbert(E)
    Kerr = let γ3 = γ3, hilbert = hilbert
        function Kerr(out, E, ρ)
            out .+= ρ*3/4*ε_0*γ3.*abs2.(hilbert(E)).*E
        end
    end
end

"""
    KerrEnv(γ3)

Kerr response for an envelope field without THG; built by `Kerr_env(γ3)`. See
[`KerrField`](@ref) for `γ3`.
"""
struct KerrEnv{T}
    γ3::T
end

"Kerr response for envelope"
Kerr_env(γ3) = KerrEnv(γ3)

function (k::KerrEnv)(out, E, ρ)
    #= The 3/4 is folded into the scalar factor here rather than left inside the
       broadcast, where its Float64 literal would promote a Float32 kernel. The value is
       unchanged: the broadcast evaluated the same product of scalars for every element. =#
    fac = Luna.scalar(E, 3/4*(ρ*ε_0*k.γ3))
    if size(E, 2) == 1
        KerrScalarEnv!(out, E, fac)
    else
        KerrVectorEnv!(out, E, fac)
    end
end

rescale(k::KerrEnv, spec, scaling) =
    KerrEnv(convert(Luna.realtype(spec), k.γ3*scaling.Eref^2/scaling.Pref))

"`fac` includes the factor 3/4; see [`KerrEnv`](@ref)."
function KerrScalarEnv!(out, E, fac)
    @. out += fac*abs2(E)*E
end

@doc (@doc KerrScalarEnv!)
function KerrVectorEnv!(out, E, fac)
    R = real(eltype(E))
    c23 = convert(R, 2/3)
    c13 = convert(R, 1/3)
    Ex = view(E, :, 1)
    Ey = view(E, :, 2)
    ox = view(out, :, 1)
    oy = view(out, :, 2)
    @. ox += fac*((abs2(Ex) + c23*abs2(Ey))*Ex + c13*conj(Ex)*Ey^2)
    @. oy += fac*((abs2(Ey) + c23*abs2(Ex))*Ey + c13*conj(Ey)*Ex^2)
    out
end

"""
    KerrEnvTHG(γ3, C)

Kerr response for an envelope field including THG; built by `Kerr_env_thg(γ3, ω0, t)`.
`C` is the carrier factor `exp(2iω₀t)` on the (oversampled) time grid. See Eq. 4, Genty
et al., Opt. Express 15 5382 (2007).
"""
struct KerrEnvTHG{T, V}
    γ3::T
    C::V
end

"Kerr response for envelope but with THG"
Kerr_env_thg(γ3, ω0, t) = KerrEnvTHG(γ3, exp.(2im*ω0.*t))

function (k::KerrEnvTHG)(out, E, ρ)
    fac = Luna.scalar(E, ρ*ε_0*k.γ3/4)
    C = k.C
    @. out += fac*(3*abs2(E) + C*E^2)*E
end

rescale(k::KerrEnvTHG, spec, scaling) = KerrEnvTHG(
    convert(Luna.realtype(spec), k.γ3*scaling.Eref^2/scaling.Pref),
    Luna.todevice(spec, k.C))

Adapt.adapt_structure(to, k::KerrEnvTHG) = KerrEnvTHG(k.γ3, Adapt.adapt(to, k.C))

resident_arrays(k::KerrEnvTHG) = (k.C,)

struct Chi2Field{χT}
    χ2::χT
    toCrystal::RotMatrix3{Float64}
    toLab::RotMatrix3{Float64}
    χ2_toLab::χT # combined matrix to multiply by toLab * χ2
    El::Vector{Float64} # field in the lab frame
    Ec::Vector{Float64} # field in the crystal frame
    Enl::Vector{Float64} # field products in the crystal frame
    Pl::Vector{Float64} # polarisation in the lab frame
end

"""
    Chi2Field(θ, ϕ, χ2)

Construct a second-order nonlinear polarisation response for real, two-component
electric fields in the lab frame.

`θ` and `ϕ` (radians) define the crystal orientation relative to the lab frame.
`χ2` must be a 3×6 second-order susceptibility tensor in contracted notation,
with column order `[xx, yy, zz, yz, xz, xy]` (where mixed terms are multiplied by 2
in [`field_products!`](@ref)).

The returned callable adds \$ε_0 P_{NL}\$ to `out` when invoked as
`response(out, E, ρ)`. Note that the density `ρ` is ignored.
"""
function Chi2Field(θ, ϕ, χ2)
    toCrystal = RotMatrix(RotZY(-ϕ, -θ)) # RotMatrix converts to static matrix
    toLab = RotMatrix(RotYZ(θ, ϕ))
    sv3 = zeros(3)
    sv6 = zeros(6)
    χ2 = SMatrix{3, 6}(χ2) # just χ2
    χ2_toLab = SMatrix{3, 6}(toLab * χ2) # χ2 and coordinate transform in one step
    Chi2Field(χ2, toCrystal, toLab, χ2_toLab, sv3, copy(sv3), sv6, copy(sv3))
end

function (c::Chi2Field)(out, E, ρ)
    for i in axes(E, 1)
        @inbounds c.El[1] = E[i, 1]
        @inbounds c.El[2] = E[i, 2]
        # note c.El[3] (Ez in the lab frame) is always zero here
        mul!(c.Ec, c.toCrystal, c.El) # transform to crystal frame
        @inbounds field_products!(c.Enl, c.Ec) # calculate nonlinear products
        mul!(c.Pl, c.χ2_toLab, c.Enl) # multiply by χ2 tensor and transform to lab frame
        @inbounds out[i, 1] += ε_0*c.Pl[1]
        @inbounds out[i, 2] += ε_0*c.Pl[2]
    end
end


"""
    field_products!(Enl, Ec)

Fill the contracted second-order field-product vector `Enl` from crystal-frame
field components `Ec`.

The output ordering is
`[Ex^2, Ey^2, Ez^2, 2EyEz, 2ExEz, 2ExEy]`, matching the 3x6 `χ2` tensor
column order `[xx, yy, zz, yz, xz, xy]` used by [`Chi2Field`](@ref).

Both `Enl` and `Ec` are mutated/read in place and are expected to have length 6
and 3, respectively.
"""
function field_products!(Enl, Ec)
    Enl[1] = Ec[1]^2
    Enl[2] = Ec[2]^2
    Enl[3] = Ec[3]^2
    Enl[4] = 2*Ec[2]*Ec[3]
    Enl[5] = 2*Ec[1]*Ec[3]
    Enl[6] = 2*Ec[1]*Ec[2]
end

struct Chi2Env{χT}
    χ2::χT
    toCrystal::RotMatrix3{Float64}
    toLab::RotMatrix3{Float64}
    χ2_toLab::χT # combined matrix to multiply by toLab * χ2
    C::Vector{ComplexF64} # carrier phase exp(iω0t) on the (oversampled) time grid
    Al::Vector{ComplexF64} # envelope in the lab frame
    Ac::Vector{ComplexF64} # envelope in the crystal frame
    Anl::Vector{ComplexF64} # combined SFG+DFG envelope products in the crystal frame
    Pl::Vector{ComplexF64} # polarisation envelope in the lab frame
end

"""
    Chi2Env(θ, ϕ, χ2, ω0, t)

Construct a second-order nonlinear polarisation response for complex envelope,
two-component electric fields in the lab frame. Envelope counterpart of
[`Chi2Field`](@ref).

`θ` and `ϕ` (radians) define the crystal orientation relative to the lab frame.
`χ2` must be a 3×6 second-order susceptibility tensor in contracted notation,
with column order `[xx, yy, zz, yz, xz, xy]`. `ω0` is the carrier frequency and `t`
the time axis on which the response is evaluated—for propagation simulations these
must be `grid.ω0` and `grid.to` (the oversampled time axis) of a `Grid.EnvGrid`.

The response includes both the sum-frequency term (``ω + ω → 2ω``, carrying the phase
factor ``e^{+iω_0t}``) and the difference-frequency term (``2ω - ω → ω``, back-conversion,
carrying ``e^{-iω_0t}``), see [`env_products!`](@ref). Optical-rectification content near
zero absolute frequency lies outside the frequency window and is removed by the grid
apodisation. Note that the grid must contain the second harmonic—use
`Grid.EnvGrid(...; thg=true)` or wavelength limits reaching below `λ0/2`.

!!! warning
    The linear operator must use a reference frame which is transparent to carrier-mixing
    nonlinearities, i.e. one whose subtracted phase is strictly linear in the absolute
    frequency (a pure time shift), like the crystal operators
    `LinearOps.make_const_linop(grid, xgrid, nfuns::Tuple)`. Envelope operators which
    subtract the carrier phase `β0` at `grid.ω0` (`thg=false`) introduce a spurious phase
    mismatch `β0 - β1ω_0` into χ⁽²⁾ processes.

The returned callable adds \$ε_0 P_{NL}\$ to `out` when invoked as
`response(out, E, ρ)`. Note that the density `ρ` is ignored.
"""
function Chi2Env(θ, ϕ, χ2, ω0, t)
    toCrystal = RotMatrix(RotZY(-ϕ, -θ)) # RotMatrix converts to static matrix
    toLab = RotMatrix(RotYZ(θ, ϕ))
    χ2 = SMatrix{3, 6}(χ2) # just χ2
    χ2_toLab = SMatrix{3, 6}(toLab * χ2) # χ2 and coordinate transform in one step
    C = exp.(1im*ω0.*t)
    sv3 = zeros(ComplexF64, 3)
    sv6 = zeros(ComplexF64, 6)
    Chi2Env(χ2, toCrystal, toLab, χ2_toLab, C, sv3, copy(sv3), sv6, copy(sv3))
end

function (c::Chi2Env)(out, E, ρ)
    size(E, 2) == 2 || error("Chi2Env requires a two-component (Nt×2) envelope field")
    length(c.C) == size(E, 1) || error(
        "Chi2Env carrier phase array does not match the field length. "
        * "The response must be constructed with the oversampled time axis `grid.to`.")
    for i in axes(E, 1)
        cp = 0.5*c.C[i] # ½exp(iω0t): sum-frequency (ω + ω → 2ω)
        cm = conj(c.C[i]) # exp(-iω0t): difference-frequency (2ω - ω → ω)
        @inbounds c.Al[1] = E[i, 1]
        @inbounds c.Al[2] = E[i, 2]
        # note c.Al[3] (Az in the lab frame) is always zero here
        mul!(c.Ac, c.toCrystal, c.Al) # transform to crystal frame
        @inbounds env_products!(c.Anl, c.Ac, cp, cm) # calculate nonlinear products
        mul!(c.Pl, c.χ2_toLab, c.Anl) # multiply by χ2 tensor and transform to lab frame
        @inbounds out[i, 1] += ε_0*c.Pl[1]
        @inbounds out[i, 2] += ε_0*c.Pl[2]
    end
end

"""
    env_products!(Anl, Ac, cp, cm)

Fill the contracted second-order envelope-product vector `Anl` from crystal-frame
envelope components `Ac`. Envelope counterpart of [`field_products!`](@ref).

With the real field given by ``E_j = \\mathrm{Re}[A_j e^{iω_0t}]``, the envelope of each
second-order product is

``(E_jE_k)_{env} = \\frac{1}{2}A_jA_k e^{+iω_0t} + \\frac{1}{2}(A_jA_k^* + A_j^*A_k)e^{-iω_0t}``

where the first (sum-frequency) and second (difference-frequency) terms enter here via
`cp` = ``\\frac{1}{2}e^{+iω_0t}`` and `cm` = ``e^{-iω_0t}`` at the current time sample.

The output ordering is `[xx, yy, zz, yz, xz, xy]` with mixed terms multiplied by 2,
matching the 3x6 `χ2` tensor column order used by [`Chi2Env`](@ref).

Both `Anl` and `Ac` are mutated/read in place and are expected to have length 6
and 3, respectively.
"""
function env_products!(Anl, Ac, cp, cm)
    Anl[1] = cp*Ac[1]^2 + cm*abs2(Ac[1])
    Anl[2] = cp*Ac[2]^2 + cm*abs2(Ac[2])
    Anl[3] = cp*Ac[3]^2 + cm*abs2(Ac[3])
    Anl[4] = 2*(cp*Ac[2]*Ac[3] + cm*real(Ac[2]*conj(Ac[3])))
    Anl[5] = 2*(cp*Ac[1]*Ac[3] + cm*real(Ac[1]*conj(Ac[3])))
    Anl[6] = 2*(cp*Ac[1]*Ac[2] + cm*real(Ac[1]*conj(Ac[2])))
end

"Response type for cumtrapz-based plasma polarisation, adapted from:
M. Geissler, G. Tempea, A. Scrinzi, M. Schnürer, F. Krausz, and T. Brabec, Physical Review Letters 83, 2930 (1999)."
struct PlasmaCumtrapz{R, EType, tType}
    ratefunc::R # the ionization rate function
    ionpot::Float64 # the ionization potential (for calculation of ionization loss)
    rate::tType # buffer to hold the rate
    fraction::tType # buffer to hold the ionization fraction
    phase::EType # buffer to hold the plasma induced (mostly) phase modulation
    J::EType # buffer to hold the plasma current
    P::EType # buffer to hold the plasma polarisation
    δt::Float64 # the time step
    preionfrac::Float64 # the pre-ionisation fraction
end

"""
    PlasmaCumtrapz(t, E, ratefunc, ionpot)

Construct the Plasma polarisation response for a field on time grid `t`
with example electric field like `E`, an ionization rate callable
`ratefunc` and ionization potential `ionpot`.
"""
function PlasmaCumtrapz(t, E, ratefunc, ionpot; preionfrac=0.0)
    rate = similar(t)
    fraction = similar(t)
    phase = similar(E)
    J = similar(E)
    P = similar(E)
    !(0.0 <= preionfrac <= 1.0) && throw(DomainError(preionfrac, "preionfrac must be between 0 and 1"))
    if preionfrac > 0.0
        @warn("Using preionfrac > 0.0 is not a well founded physical model. Use only after careful consideration.")
    end
    return PlasmaCumtrapz(ratefunc, ionpot, rate, fraction, phase, J, P, t[2]-t[1], preionfrac)
end

"The plasma response for a scalar electric field"
function PlasmaScalar!(Plas::PlasmaCumtrapz, E)
    Plas.ratefunc(Plas.rate, E)
    Maths.cumtrapz!(Plas.fraction, Plas.rate, Plas.δt)
    @. Plas.fraction = Plas.preionfrac + 1 - exp(-Plas.fraction)
    @. Plas.phase = Plas.fraction * e_ratio * E
    Maths.cumtrapz!(Plas.J, Plas.phase, Plas.δt)
    for ii in eachindex(E)
        if abs(E[ii]) > 0
            Plas.J[ii] += Plas.ionpot * Plas.rate[ii] * (1-Plas.fraction[ii])/E[ii]
        end
    end
    Maths.cumtrapz!(Plas.P, Plas.J, Plas.δt)
end

"""
The plasma response for a vector electric field.

We take the magnitude of the electric field to calculate the ionization
rate and fraction, and then solve the plasma polarisation component-wise
for the vector field.

A similar approach was used in: C Tailliez et al 2020 New J. Phys. 22 103038.
"""
function PlasmaVector!(Plas::PlasmaCumtrapz, E)
    Ex = E[:,1]
    Ey = E[:,2]
    Em = @. hypot.(Ex, Ey)
    Plas.ratefunc(Plas.rate, Em)
    Maths.cumtrapz!(Plas.fraction, Plas.rate, Plas.δt)
    @. Plas.fraction = Plas.preionfrac + 1 - exp(-Plas.fraction)
    @. Plas.phase = Plas.fraction * e_ratio * E
    Maths.cumtrapz!(Plas.J, Plas.phase, Plas.δt)
    for ii in eachindex(Em)
        if abs(Em[ii]) > 0
            pre = Plas.ionpot * Plas.rate[ii] * (1-Plas.fraction[ii])/Em[ii]^2
            Plas.J[ii,1] += pre*Ex[ii]
            Plas.J[ii,2] += pre*Ey[ii]
        end
    end
    Maths.cumtrapz!(Plas.P, Plas.J, Plas.δt)
end

"Handle plasma polarisation routing to `PlasmaVector` or `PlasmaScalar`."
function (Plas::PlasmaCumtrapz)(out, Et, ρ)
    if ndims(Et) > 1
        if size(Et, 2) == 1 # handle scalar case but within modal simulation
            PlasmaScalar!(Plas, reshape(Et, size(Et,1)))
            out .+= ρ .* reshape(Plas.P, size(Et))
        else
            PlasmaVector!(Plas, Et) # vector case
            out .+= ρ .* Plas.P
        end
    else
        PlasmaScalar!(Plas, Et) # straight scalar case
        out .+= ρ .* Plas.P
    end
end

"Raman polarisation response type"
abstract type RamanPolar end

"Raman polarisation response type for a carrier resolved field"
struct RamanPolarField{TR, Tt, Thv, Tω, Tv, FTt, HTt} <: RamanPolar
    r::TR # Raman response
    h::Tt # doubled buffer to hold response + padding
    ht::Thv # buffer to hold time domain response
    hω::Tω # the frequency domain Raman response function
    Eω2::Tω # buffer to hold the Fourier transform of E^2
    Pω::Tω # buffer to hold the frequency domain polarisation
    E2::Tt # buffer to hold E^2
    E2v::Tv # view into first half of E2
    P::Tt # buffer to hold the time domain polarisation
    Pout::Tt # buffer to hold the output portion of the time domain polarisation
    FT::FTt # Fourier transform plan
    HT::HTt # Hilbert transform
    thg::Bool # do we include third harmonic generation
    dt::Float64 # time step for scaling
end

"Raman polarisation response type for an envelope"
struct RamanPolarEnv{TR, Tt, Thv, Tω, Tv, FTt} <: RamanPolar
    r::TR # Raman response
    h::Tt # doubled buffer to hold response + padding
    ht::Thv # buffer to hold time domain response
    hω::Tω # the frequency domain Raman response function
    Eω2::Tω # buffer to hold the Fourier transform of E^2
    Pω::Tω # buffer to hold the frequency domain polarisation
    E2::Tω # buffer to hold E^2
    E2v::Tv # view into first half of E2
    P::Tω # buffer to hold the time domain polarisation
    Pout::Tω # buffer to hold the output portion of the time domain polarisation
    FT::FTt # Fourier transform plan
    dt::Float64 # time step for scaling
end

"""
    RamanPolarField(t, ht; thg=true)

Construct Raman polarisation response for a field on time grid `t`
using response function `r`. If `thg=false` then exclude the third
harmonic generation component of the response.
"""
function RamanPolarField(t, r; thg=true)
    h = zeros(length(t)*2) # note double grid size, see explanation below
    ht = view(h, 1:length(t))
    Utils.loadFFTwisdom()
    FT = FFTW.plan_rfft(h, 1, flags=Luna.settings["fftw_flag"])
    inv(FT)
    Utils.saveFFTwisdom()
    hω = FT * h
    Eω2 = similar(hω)
    Pω = similar(hω)
    E2 = similar(h)
    E2v = view(E2, 1:length(t))
    P = similar(h)
    Pout = similar(t)
    HT = Maths.plan_hilbert(Pout)
    fill!(E2, 0.0)
    RamanPolarField(r, h, ht, hω, Eω2, Pω, E2, E2v, P, Pout, FT, HT, thg, t[2] - t[1])
end

"""
    RamanPolarEnv(t, ht)

Construct Raman polarisation response for an envelope on time grid `t`
using response function `r`.
"""
function RamanPolarEnv(t, r)
    h = zeros(length(t)*2) # note double grid size, see explanation below
    ht = view(h, 1:length(t))
    Utils.loadFFTwisdom()
    FT = FFTW.plan_fft(h, 1, flags=Luna.settings["fftw_flag"])
    inv(FT)
    Utils.saveFFTwisdom()
    hω = FT * h
    Eω2 = similar(hω)
    Pω = similar(hω)
    E2 = similar(hω)
    P = similar(hω)
    Pout = Array{ComplexF64,}(undef,size(t))
    E2v = view(E2, 1:length(t))
    fill!(E2, 0.0)
    RamanPolarEnv(r, h, ht, hω, Eω2, Pω, E2, E2v, P, Pout, FT, t[2] - t[1])
end

"Square the field or envelope"
function sqr!(R::RamanPolarField, E)
    if !R.thg
        # see documentation for factor of 1/2 here
        R.E2v .= 1/2 .* abs2.(R.HT(E))
    else
        R.E2v .= E.^2
    end
end

function sqr!(R::RamanPolarEnv, E)
    # see documentation for factor of 1/2 here
    R.E2v .= 1/2 .* abs2.(E)
end

"Calculate Raman polarisation for field/envelope Et"
function (R::RamanPolar)(out, Et, ρ)
    # get the field as a 1D Array
    n = size(Et, 1)
    if ndims(Et) > 1
        if size(Et, 2) == 1 # handle scalar case but within modal simulation
            E = reshape(Et, n)
        else
            # handle vector case
            error("vector Raman not yet implemented")
        end
    else
        E = Et # handle straight scalar case
    end

    # square the field or envelope in first half
    # corresponding to the field/envelope grid size
    sqr!(R, E)

    # update frequency domain response function `hω`.
    # we fill only up to the first half of h (using the view ht)
    # i.e. only the part corresponding to the original time grid
    # note that the response function time 0 is put into the first element of the response array
    # this ensures that causality is maintained, and no artificial delay between the field and
    # the start of the response function occurs, at each convolution point.
    R.r(R.ht, ρ)
    R.hω .= R.FT * R.h

    # convolution by multiplication in frequency domain
    # The double grid gives us accurate full convolution between the full field grid
    # and full response function. It is unnecessary for highly damped responses, like
    # in glass. But for gases with very long decay times it prevents artefacts due to
    # truncation of the response function. There is likely a more efficient way. But
    # this is safe, until we come up with one.
    # we scale to correct for missing dt*dt*df from IFFT(FFT*FFT)
    # the ifft already scales by 1/n = dt*df, so we need an additional dt
    R.Eω2 .= R.FT * R.E2
    @. R.Pω = R.hω * R.Eω2 * R.dt
    R.P .= R.FT \ R.Pω

    # calculate full polarisation, extracting only the valid
    # grid region, which is the first length(E) part.
    for i = 1:length(E)
        R.Pout[i] = ρ*E[i]*R.P[i]
    end

    # copy to output in dimensions requested
    if ndims(Et) > 1
        out .+= reshape(R.Pout, size(Et))
    else
        out .+= R.Pout
    end
end

end
