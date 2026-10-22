module LinearOps
import FFTW
import Luna: Modes, Grid, PhysData, Maths
import Luna.PhysData: wlfreq, c, crystal_internal_angle

#=
Reference (carrier) frequency ω0 for the linear-operator frame, in which the total phase
subtracted from the propagation constant is β1*(ω - ω0) + β0:
- `thg=false` (envelope grids only): ω0 = grid.ω0 and β0 is the propagation constant at the
  carrier, so the envelope is phase-stationary at the carrier. Only appropriate when no
  carrier-mixing nonlinearity (THG, χ⁽²⁾) is present: the constant frame phase β0 - β1*ω0
  appears as a spurious phase mismatch in carrier-mixing processes.
- `thg=true`: there is no carrier reference (ω0 = 0, β0 = 0), so the subtracted phase is
  β1*ω—strictly linear in the absolute frequency, i.e. a pure time shift, which is
  transparent to all instantaneous nonlinear responses including carrier-mixing ones.
The single-argument forms give the historical defaults.
=#
getω0(grid::Grid.EnvGrid, thg=false) = thg ? 0.0 : grid.ω0
getω0(grid::Grid.RealGrid, thg=true) = 0.0

#=================================================#
#===============    FREE SPACE     ===============#
#=================================================#

"""
    βz(βsq)

Longitudinal wavevector ``β_z`` of a plane-wave component from its square
``β_z^2 = k^2 - k_\\perp^2``, as a complex number.

Above cutoff (`βsq ≥ 0`) it is real and the component propagates. Below cutoff (`βsq < 0`)
the component is evanescent: ``β_z = -iκ`` with ``κ = \\sqrt{k_\\perp^2 - k^2}``, the sign
chosen so that Luna's propagator `exp(-im*βz*z)` *decays* as `exp(-κz)`. This is the exact
one-way solution of the Helmholtz equation, and it is used as is. There is deliberately no
cap on `κ` here: a strongly damped component makes the interaction-picture stepper stiff,
but that is a property of the solver, not the physics, and is handled by
[`Luna.Boundaries`](@ref) (which clamps the decay rate of the operator and tapers the
nonlinear source to match) rather than by falsifying the operator.
"""
βz(βsq) = βsq < 0 ? -im*sqrt(-βsq) : complex(sqrt(βsq))

#= Out of band (`!grid.sidx`) the operator is exactly zero. Nothing propagates there: the
   field is band-limited at the start of `Luna.run` and `NonlinearRHS` zeroes the nonlinear
   polarisation there. A nonzero value would be harmless in principle, but any decay written
   there (as the evanescent `-κ` would be, since `k2 == 0` out of band) is amplified by the
   interaction-picture back-propagation, and `exp(k⊥ Δz)` overflows for the k⊥ of any fine
   transverse grid. Zero is the only value that is safe for every step size. =#
function fill_linop_matrix!(out, grid, β1::Number, βref::Number, ω0::Number, k2, kperp2, idcs)
    for ii in idcs
        for ip in axes(k2, 2)
            for iω in eachindex(grid.ω)
                if !grid.sidx[iω]
                    out[iω, ip, ii] = 0
                    continue
                end
                βsq = k2[iω, ip] - kperp2[ii]
                out[iω, ip, ii] = -im*(βz(βsq) - β1*(grid.ω[iω] - ω0) - βref)
            end
        end
    end
end

function transverse_k2(xygrid::Grid.FreeGrid)
    kperp2 = @. xygrid.kx^2 + (xygrid.ky^2)'
    idcs = CartesianIndices((length(xygrid.kx), length(xygrid.ky)))
    kperp2, idcs
end

function transverse_k2(xgrid::Grid.Free2DGrid)
    kperp2 = xgrid.kx.^2
    idcs = CartesianIndices(xgrid.kx)
    kperp2, idcs
end

function transverse_k2(rg::Grid.RadialGrid)
    kperp2 = Grid.kperp2(rg)
    idcs = CartesianIndices(rg.k)
    kperp2, idcs
end

# Hankel.QDHT is deprecated as a Luna transverse grid; convert it (Grid.RadialGrid warns)
transverse_k2(q::Grid.HankelTransform) = transverse_k2(Grid.RadialGrid(q))


"""
    make_const_linop(grid, xygrid, n, β1, β0, ω0=getω0(grid))

Low-level constructor for a constant (z-invariant) free-space linear operator.

Arguments:
- `grid`: `Grid.AbstractGrid` (`RealGrid` or `EnvGrid`)
- `xygrid`: transverse grid (`Grid.FreeGrid`, `Grid.Free2DGrid`, or `Grid.RadialGrid`)
- `n`: refractive-index table on `grid.ω`, with one column per polarisation
- `β1`: inverse reference-frame velocity
- `β0`: reference wavevector offset (typically zero for `RealGrid` and optional for `EnvGrid`)
- `ω0`: reference frequency of the frame; the total subtracted phase is `β1*(ω - ω0) + β0`.
    Defaults to `grid.ω0` for `EnvGrid` (co-rotating frame) and `0` for `RealGrid`; pass `0`
    (with `β0 = 0`) for the pure time-shift frame required by carrier-mixing nonlinearities
    (see `getω0`).

The output has shape `(Nω, Npol, N⊥...)`, where `N⊥...` matches the transverse grid.
"""
function make_const_linop(grid::Grid.AbstractGrid,
                          xygrid::Grid.TransverseGrid,
                          n::AbstractVecOrMat, β1::Number, β0::Number, ω0::Number=getω0(grid))
    kperp2, idcs = transverse_k2(xygrid)
    k2 = @. (n*grid.ω/c)^2
    out = zeros(ComplexF64, (length(grid.ω), size(n, 2), size(idcs)...))
    fill_linop_matrix!(out, grid, β1, β0, ω0, k2, kperp2, idcs)
    return out
end

thg_default(grid::Grid.RealGrid) = true
thg_default(grid::Grid.EnvGrid) = false

checkthg(grid::Grid.EnvGrid, thg) = nothing
checkthg(grid::Grid.RealGrid, thg) = thg || error("`thg` must be `true` for `RealGrid`s")

getβ0_n(grid::Grid.RealGrid, nfun, thg) = 0.0
getβ0_n(grid::Grid.EnvGrid, nfun, thg) = thg ? 0.0 : grid.ω0/c * nfun(wlfreq(grid.ω0))[end]

"""
    make_const_linop(grid, xygrid, nfun; thg=thg_default(grid))

Build a constant free-space operator from a refractive-index function `nfun`.

`nfun(λ)` must return either a scalar index (single polarisation) or a vector/tuple of indices
(e.g. both x and y polarisation). For `EnvGrid`, `thg=false` subtracts the reference propagation
constant at `grid.ω0` (envelope phase-stationary at the carrier, but **not** transparent to
carrier-mixing nonlinearities); `thg=true` subtracts only the group delay `β1*ω`, which is a
pure time shift and keeps the phase bookkeeping of carrier-mixing responses
(`Kerr_env_thg`, [`Luna.Nonlinear.Chi2Env`](@ref)) exact.
"""
function make_const_linop(grid::Grid.AbstractGrid,
                          xygrid::Grid.TransverseGrid,
                          nfun, thg::Bool=thg_default(grid))
    checkthg(grid, thg)
    ωfirst = grid.ω[findfirst(grid.sidx)]
    np = length(nfun(wlfreq(ωfirst))) # 1 if single ref index, 2 if nx, ny
    n = zeros(Float64, (length(grid.ω), np))
    for (ii, si) in enumerate(grid.sidx)
        if si
            n[ii, :] .= nfun(wlfreq(grid.ω[ii]))
        end
    end
    β1 = PhysData.dispersion_func(1, λ -> nfun(λ)[end])(grid.referenceλ)
    β0 = getβ0_n(grid, nfun, thg)
    make_const_linop(grid, xygrid, n, β1, β0, getω0(grid, thg))
end

"""
    make_const_linop(grid, xygrid::Grid.FreeGrid, nfuns::Tuple)

Constant full-3D free-space operator for crystal optics with two polarisation branches.

`nfuns = (nfunx, nfuny)`, where `nfunx(λ, δθ)` is the refractive index for x-polarisation
and `nfuny(λ)` is the refractive index for y-polarisation. The reference frame velocity
is calculated from `nfuny(λ)`.

For both `RealGrid` and `EnvGrid`, the phase subtracted is `β1*ω`—strictly linear in the
absolute frequency, i.e. a pure time shift. Note that for `EnvGrid` this differs from the
other envelope operators, which additionally subtract the residual carrier phase `β0 - β1*ω0`
at `grid.ω0` (`thg=false`): a constant frame phase is not transparent to carrier-mixing
nonlinearities such as [`Luna.Nonlinear.Chi2Env`](@ref)—it appears as a spurious phase mismatch
in e.g. second-harmonic generation—so these operators subtract the group delay only. The
envelope therefore acquires a residual global phase rotation `(β0 - β1*ω0)z`, which does
not affect intensities or spectra.
"""
function make_const_linop(grid::Grid.AbstractGrid, xygrid::Grid.FreeGrid, nfuns::Tuple)
    nfunx, nfuny = nfuns
    # here nfunx(λ, δθ) also takes the angle and returns n_x(λ, θ)
    # nfuny(λ) just takes wavelength
    out = zeros(ComplexF64, (length(grid.ω), 2, length(xygrid.kx), length(xygrid.ky)))
    β1 = PhysData.dispersion_func(1, nfuny)(grid.referenceλ)
    for (iω, si) in enumerate(grid.sidx)
        if si
            ny = nfuny(wlfreq(grid.ω[iω]))
            ksq_ypol = (ny*grid.ω[iω]/c)^2
            for (ikx, kxi) in enumerate(xygrid.kx)
                δθ = crystal_internal_angle(nfunx, grid.ω[iω], kxi)
                nx = nfunx(wlfreq(grid.ω[iω]), δθ)
                for (iky, kyi) in enumerate(xygrid.ky)
                    kperp2 = kxi^2 + kyi^2
                    k_xpol = nx*grid.ω[iω]/c
                    # evanescent components (βsq < 0) decay at their exact rate, see βz
                    out[iω, 1, ikx, iky] = -im*(βz(k_xpol^2 - kperp2) - β1*grid.ω[iω])
                    out[iω, 2, ikx, iky] = -im*(βz(ksq_ypol - kperp2) - β1*grid.ω[iω])
                end
            end
        end
    end
    out
end


"""
    make_const_linop(grid, xgrid::Grid.Free2DGrid, nfuns::Tuple)

Constant free-space operator for 2D (`x-z`) crystal propagation with two polarisations.

`nfuns = (nfunx, nfuny)` follows the same convention as the full-3D overload, as does
the reference frame (`β1*ω` from `nfuny`, a pure time shift, for both `RealGrid` and
`EnvGrid`—see the full-3D overload for why).
"""
function make_const_linop(grid::Grid.AbstractGrid, xgrid::Grid.Free2DGrid, nfuns::Tuple)
    nfunx, nfuny = nfuns
    # here nfunx(λ, δθ) also takes the angle and returns n_x(λ, θ)
    # nfuny(λ) just takes wavelength
    out = zeros(ComplexF64, (length(grid.ω), 2, length(xgrid.kx)))
    β1 = PhysData.dispersion_func(1, nfuny)(grid.referenceλ)
    for (iω, si) in enumerate(grid.sidx)
        if si
            ny = nfuny(wlfreq(grid.ω[iω]))
            ksq_ypol = (ny*grid.ω[iω]/c)^2
            for (ik, kxi) in enumerate(xgrid.kx)
                δθ = crystal_internal_angle(nfunx, grid.ω[iω], kxi)
                nx = nfunx(wlfreq(grid.ω[iω]), δθ)
                k_xpol = nx*grid.ω[iω]/c
                # evanescent components (βsq < 0) decay at their exact rate, see βz
                out[iω, 1, ik] = -im*(βz(k_xpol^2 - kxi^2) - β1*grid.ω[iω])
                out[iω, 2, ik] = -im*(βz(ksq_ypol - kxi^2) - β1*grid.ω[iω])
            end
        end
    end
    out
end

#= Deprecated entry points: a Hankel.QDHT in place of a Grid.RadialGrid. These repeat the
   concrete arities of the methods above rather than taking `args...`, so that they cannot
   be ambiguous with the modal or the βfun!/αfun! methods. =#
function make_const_linop(grid::Grid.AbstractGrid, q::Grid.HankelTransform,
                          n::AbstractVecOrMat, β1::Number, β0::Number,
                          ω0::Number=getω0(grid))
    make_const_linop(grid, Grid.RadialGrid(q), n, β1, β0, ω0)
end

function make_const_linop(grid::Grid.AbstractGrid, q::Grid.HankelTransform,
                          nfun, thg::Bool=thg_default(grid))
    make_const_linop(grid, Grid.RadialGrid(q), nfun, thg)
end

function make_linop(grid::Grid.AbstractGrid, q::Grid.HankelTransform,
                    nfun, thg::Bool=thg_default(grid))
    make_linop(grid, Grid.RadialGrid(q), nfun, thg)
end

"""
    make_linop(grid, xygrid, nfun)

Create a z-dependent free-space linear operator closure.

Applies to `xygrid::Grid.FreeGrid` (full 3D), `Grid.Free2DGrid` (x-z), and
`Grid.RadialGrid` (radial symmetry).

Returns `linop!(out, z)`, which fills `out` in-place for propagation distance `z`.
`nfun(ω; z)` may return one or multiple refractive indices; the last branch
defines the reference-frame velocity and, for `EnvGrid` with `thg=false`, the
reference phase subtraction at `grid.ω0`.
"""
function make_linop(grid::Grid.AbstractGrid,
                    xygrid::Grid.TransverseGrid,
                    nfun, thg::Bool=thg_default(grid))
    checkthg(grid, thg)
    kperp2, idcs = transverse_k2(xygrid)
    ωfirst = grid.ω[findfirst(grid.sidx)]
    np = length(nfun(ωfirst; z=0)) # 1 if single ref index, 2 if nx, ny
    k2 = zeros(Float64, (length(grid.ω), np))
    nfunλ(z) = λ -> nfun(wlfreq(λ); z)[end]
    ω0 = getω0(grid, thg)
    function linop!(out, z)
        β1 = PhysData.dispersion_func(1, nfunλ(z))(grid.referenceλ)
        β0 = getβ0_n(grid, nfunλ(z), thg)
        for (ii, si) in enumerate(grid.sidx)
            if si
                k2[ii, :] .= (nfun(grid.ω[ii]; z) .* grid.ω[ii] ./ c).^2
            end
        end
        fill_linop_matrix!(out, grid, β1, β0, ω0, k2, kperp2, idcs)
    end
end


#=================================================#
#===============   MODE AVERAGE   ================#
#=================================================#

"""
    αlim!(α)

Limit α so that we do not get overflow in exp(α*dz)
"""
function αlim!(α)
    # magic number: this is 130 dB/cm
    # a test script sensitive to this is test_main_rect_env.jl
    clamp!(α, 0.0, 3000.0)
end

"""
    conj_clamp(n, ω)

Simultaneously conjugate and clamp the effective index `n` to safe levels.

The real part is lower-bounded at 1e-3 and the imaginary part upper-bounded at an attenuation
coefficient `α` of 3000 (130 dB/cm). The limits are somewhat arbitrary and chosen empirically
from previous bugs. See https://github.com/LupoLab/Luna/pull/142.

See also [`αlim!`](@ref).
"""
conj_clamp(n, ω) = clamp(real(n), 1e-3, Inf) - im*clamp(imag(n), 0, 3000*c/ω)

function make_const_linop(grid::Grid.AbstractGrid, βfun!, αfun!, β1::Number, β0::Number,
                          ω0::Number=getω0(grid))
    β = similar(grid.ω)
    βfun!(β, 0)
    α = similar(grid.ω)
    αfun!(α, 0)
    αlim!(α)
    linop = @. -im*(β - β1*(grid.ω - ω0) - β0) - α/2
    linop[.!grid.sidx] .= 0
    return linop
end

# see getω0 for the reference-phase conventions
getβ0_mode(grid::Grid.RealGrid, mode, λ0, thg) = 0.0
getβ0_mode(grid::Grid.EnvGrid, mode, λ0, thg) = thg ? 0.0 : Modes.β(mode, wlfreq(λ0))

"""
    make_const_linop(grid, mode, λ0)

Make constant linear operator for mode-averaged propagation in mode `mode` with a reference
wavelength `λ0`. For the meaning of `thg` on an `EnvGrid` see
[`make_const_linop(grid, xygrid, nfun)`](@ref).
"""
function make_const_linop(grid::Grid.AbstractGrid, mode::Modes.AbstractMode, λ0;
                          thg::Bool=thg_default(grid))
    checkthg(grid, thg)
    β1 = Modes.dispersion(mode, 1, wlfreq(λ0))
    β0 = getβ0_mode(grid, mode, λ0, thg)
    βconst = zero(grid.ω)
    βconst[grid.sidx] = Modes.β.(mode, grid.ω[grid.sidx])
    βconst[.!grid.sidx] .= 1
    function βfun!(out, z)
        out .= βconst
    end
    αconst = zero(grid.ω)
    αconst[grid.sidx] = Modes.α.(mode, grid.ω[grid.sidx])
    function αfun!(out, z)
        out .= αconst
    end
    make_const_linop(grid, βfun!, αfun!, β1, β0, getω0(grid, thg)), βfun!, β1, αfun!
end


"""
    neff_β_grid(grid, mode, λ0)

Create closures which return the effective index and propagation constant
as a function of the frequency grid **index**, rather than the frequency itself.
Any [`Modes.AbstractMode`](@ref) may define its own method for `neff_β_grid` to
accelerate repeated calculation on the same frequency grid.
"""
function neff_β_grid(grid, mode, λ0)
    let grid=grid, mode=mode
        _neff(iω; z) = Modes.neff(mode, grid.ω[iω]; z=z)
        _β(iω; z) = Modes.β(mode, grid.ω[iω]; z=z)
        _neff, _β
    end
end

"""
    make_linop(grid, mode, λ0; thg=false)

Create z-dependent mode-averaged linear-operator closures.

For `RealGrid`, returns `(linop!, βfun!)`; for `EnvGrid`, `thg=false` additionally
subtracts the reference phase at `λ0` while `thg=true` subtracts the group delay `β1*ω`
(see [`make_const_linop(grid, xygrid, nfun)`](@ref)).
"""
function make_linop(grid::Grid.RealGrid, mode::Modes.AbstractMode, λ0)
    sidcs = (1:length(grid.ω))[grid.sidx]
    neff, β = neff_β_grid(grid, mode, λ0)
    linop! = let neff=neff, ω=grid.ω, mode=mode, λ0=λ0
        function linop!(out, z)
            fill!(out, 0.0)
            β1 = Modes.dispersion(mode, 1, wlfreq(λ0), z=z)::Float64
            for iω in sidcs
                nc = conj_clamp(neff(iω; z=z), ω[iω])
                out[iω] = -im*(ω[iω]/c*nc - ω[iω]*β1)
            end
        end
    end
    βfun! = let β=β, ω=grid.ω
        function βfun!(out, z)
            fill!(out, 1.0)
            for iω in sidcs
                out[iω] = β(iω; z=z)
            end
        end
    end
    return linop!, βfun!
end

function make_linop(grid::Grid.EnvGrid, mode::Modes.AbstractMode, λ0; thg=false)
    sidcs = (1:length(grid.ω))[grid.sidx]
    neff, β = neff_β_grid(grid, mode, λ0)
    # see getω0 for the reference-phase conventions
    linop! = let neff=neff, ω=grid.ω, mode=mode, λ0=λ0, ω0=getω0(grid, thg), sidcs=sidcs
        function linop!(out, z)
            fill!(out, 0.0)
            β1 = Modes.dispersion(mode, 1, wlfreq(λ0), z=z)::Float64
            βref = thg ? 0.0 : Modes.β(mode, wlfreq(λ0), z=z)
            for iω in sidcs
                nc = conj_clamp(neff(iω; z=z), ω[iω])
                out[iω] = -im*(ω[iω]/c*nc - (ω[iω] - ω0)*β1 - βref)
            end
        end
    end
    βfun! = let β=β, sidcs=sidcs
        function βfun!(out, z)
            fill!(out, 1.0)
            for iω in sidcs
                out[iω] = β(iω, z=z)
            end
        end
    end
    return linop!, βfun!
end

#=================================================#
#=================   MULTIMODE   =================#
#=================================================#
"""
    make_const_linop(grid, modes, λ0; ref_mode=1)

Make constant (z-invariant) linear operator for multimode propagation. The frame velocity is
taken as the group velocity at wavelength `λ0` in the mode given by `ref_mode` (which
indexes into `modes`).

Returns an array of shape `(Nω, Nmodes)`.
"""
function make_const_linop(grid::Grid.RealGrid, modes::Modes.ModeCollection, λ0; ref_mode=1)
    β1 = Modes.dispersion(modes[ref_mode], 1, wlfreq(λ0))
    nmodes = length(modes)
    linops = zeros(ComplexF64, length(grid.ω), nmodes)
    for i = 1:nmodes
        βconst = zero(grid.ω)
        βconst[grid.sidx] = Modes.β.(modes[i], grid.ω[grid.sidx])
        βconst[.!grid.sidx] .= 1
        α = zeros(length(grid.ω))
        α[grid.sidx] .= Modes.α.(modes[i], grid.ω[grid.sidx])
        αlim!(α)
        linops[:,i] = im.*(-βconst .+ grid.ω.*β1) .- α./2
    end
    linops
end

function make_const_linop(grid::Grid.EnvGrid, modes::Modes.ModeCollection, λ0; ref_mode=1, thg=false)
    β1 = Modes.dispersion(modes[ref_mode], 1, wlfreq(λ0))
    # see getω0 for the reference-phase conventions
    ω0 = getω0(grid, thg)
    βref = thg ? 0.0 : Modes.β(modes[ref_mode], wlfreq(λ0))
    nmodes = length(modes)
    linops = zeros(ComplexF64, length(grid.ω), nmodes)
    for i = 1:nmodes
        βconst = zero(grid.ω)
        βconst[grid.sidx] = Modes.β.(modes[i], grid.ω[grid.sidx])
        βconst[.!grid.sidx] .= 1
        α = Modes.α.(modes[i], grid.ω)
        αlim!(α)
        linops[:,i] = -im.*(βconst .- (grid.ω .- ω0).*β1 .- βref) .- α./2
    end
    linops
end

"""
    neff_grid(grid, modes, λ0; ref_mode=1)

Create a closure that returns the effective index as a function of the frequency grid and mode
**index**, rather than the mode and frequency themselves. Any [`Modes.AbstractMode`](@ref)
may define its own method for `neff_grid` to accelerate repeated calculation on the same
frequency grid.
"""
function neff_grid(grid, modes, λ0; ref_mode=1)
    _neff = let grid=grid, modes=modes
        _neff(iω, iim; z) = Modes.neff(modes[iim], grid.ω[iω]; z=z)
    end
    _neff
end

"""
    make_linop(grid, modes, λ0; ref_mode=1, thg=false)

Create a z-dependent multimode linear-operator closure `linop!(out, z)`.

The output is filled in-place with shape `(Nω, Nmodes)`. For `EnvGrid`, setting
`thg=false` subtracts the reference phase of `modes[ref_mode]` at `λ0`, while `thg=true`
subtracts the group delay `β1*ω` (see [`make_const_linop(grid, xygrid, nfun)`](@ref)).
"""
function make_linop(grid::Grid.RealGrid, modes::Modes.ModeCollection, λ0; ref_mode=1)
    sidcs = (1:length(grid.ω))[grid.sidx]
    neff = neff_grid(grid, modes, λ0; ref_mode=ref_mode)
    linop! = let neff=neff, ω=grid.ω, modes=modes, λ0=λ0, ref_mode=ref_mode
        function linop!(out, z)
            β1 = Modes.dispersion(modes[ref_mode], 1, wlfreq(λ0), z=z)::Float64
            fill!(out, 0.0)
            for i in eachindex(modes)
                for iω in sidcs
                    nc = conj_clamp(neff(iω, i; z=z), ω[iω])
                    out[iω, i] = -im*(ω[iω]/c*nc - ω[iω]*β1)
                end
            end
        end
    end
end

function make_linop(grid::Grid.EnvGrid, modes::Modes.ModeCollection, λ0; ref_mode=1, thg=false)
    sidcs = (1:length(grid.ω))[grid.sidx]
    neff = neff_grid(grid, modes, λ0; ref_mode=ref_mode)
    # see getω0 for the reference-phase conventions
    linop! = let neff=neff, ω=grid.ω, modes=modes, λ0=λ0, ω0=getω0(grid, thg), ref_mode=ref_mode
        function linop!(out, z)
            β1 = Modes.dispersion(modes[ref_mode], 1, wlfreq(λ0), z=z)::Float64
            fill!(out, 0.0)
            βref = thg ? 0.0 : Modes.β(modes[ref_mode], wlfreq(λ0), z=z)
            for i in eachindex(modes)
                for iω in sidcs
                    nc = conj_clamp(neff(iω, i; z=z), ω[iω])
                    out[iω, i] = -im*(ω[iω]/c*nc - (ω[iω] - ω0)*β1 - βref)
                end
            end
        end
    end
end



#=================================================#
#=========  UNITARY PHASE OF THE LINOP  ==========#
#=================================================#

"""
    AbstractUnitaryPhase

Callable which fills an array with the accumulated **unitary** phase of a linear operator,
`Φ(ω, z) = ∫₀ᶻ imag(linop(ω, z')) dz'`, i.e. the phase the linear operator would impart in
the absence of loss or gain. Subtypes implement `(p::AbstractUnitaryPhase)(out, z)`, which
**mutates** `out` and returns it.

This is the propagator required by the modified shot-noise model: the noise field must feel
dispersion (unitary, and the phase evolution of the medium's own vacuum modes) but neither
loss nor gain. See [`unitary_phase`](@ref) and
[`NonlinearRHS.ModifiedNoise`](@ref Luna.NonlinearRHS.ModifiedNoise).
"""
abstract type AbstractUnitaryPhase end

"""
    ConstUnitaryPhase(dφ)

[`AbstractUnitaryPhase`](@ref) for a `z`-invariant linear operator, for which
`Φ(ω, z) = imag(linop)·z` exactly.
"""
struct ConstUnitaryPhase{aT<:AbstractArray} <: AbstractUnitaryPhase
    dφ::aT # imag(linop), i.e. dΦ/dz, constant in z
end

#= Scalars enter the broadcasts in the element type of `out`: on a device a Float64
   kernel argument is refused (Metal) or slow, and on the host Float64 this is a no-op. =#
function (p::ConstUnitaryPhase)(out, z)
    zT = convert(eltype(out), z)
    @. out = p.dφ * zT
    out
end

"""
    TabulatedUnitaryPhase(z, Φ, dΦ, atol, err)

[`AbstractUnitaryPhase`](@ref) for a `z`-dependent linear operator. `Φ` and `dΦ` hold the
accumulated phase and its derivative `imag(linop)` at the nodes `z`, stacked along the last
dimension. The nodes are **not** uniformly spaced: they are placed by
[`unitary_phase`](@ref) so that the interpolation error stays below `atol`, and `err` records
the largest error actually measured during construction. Between nodes the phase is
interpolated with a cubic Hermite, which uses both the tabulated value and the tabulated
derivative. Outside `[z[1], z[end]]` it extrapolates linearly.
"""
struct TabulatedUnitaryPhase{aT<:AbstractArray} <: AbstractUnitaryPhase
    z::Vector{Float64}
    Φ::aT
    dΦ::aT
    atol::Float64 # tolerance the nodes were placed to satisfy, in radians
    err::Float64 # largest interpolation error measured while placing them, in radians
end

function (p::TabulatedUnitaryPhase)(out, z)
    zs = p.z
    d = ndims(p.Φ)
    if z >= zs[end] || z <= zs[1] # linear extrapolation off either end of the table
        k = z <= zs[1] ? 1 : length(zs)
        Φk = selectdim(p.Φ, d, k)
        dk = selectdim(p.dΦ, d, k)
        δ = convert(eltype(out), z - zs[k])
        @. out = Φk + dk*δ
        return out
    end
    k = searchsortedlast(zs, z)
    h = zs[k+1] - zs[k]
    s = (z - zs[k])/h
    # cubic Hermite basis
    s2 = s*s
    s3 = s2*s
    h00 = 2s3 - 3s2 + 1
    h10 = h*(s3 - 2s2 + s)
    h01 = -2s3 + 3s2
    h11 = h*(s3 - s2)
    h00, h10, h01, h11 = convert.(eltype(out), (h00, h10, h01, h11))
    Φk = selectdim(p.Φ, d, k)
    Φk1 = selectdim(p.Φ, d, k+1)
    dk = selectdim(p.dΦ, d, k)
    dk1 = selectdim(p.dΦ, d, k+1)
    @. out = h00*Φk + h10*dk + h01*Φk1 + h11*dk1
    out
end

"""
    unitary_phase(linop, grid, sz=size(linop); zmax, atol=1e-6, margin=0.05, maxdepth=24, maxnodes=4096)

Create an [`AbstractUnitaryPhase`](@ref) for the linear operator `linop`, which may be either
a constant array or a mutating function `linop!(out, z)`. `sz` is the size of the array
`linop!` fills (equal to `size(linop)` in the constant case, and to the size of the noise
field in general).

A constant `linop` needs no work: `Φ = imag(linop)·z` exactly.

For a `z`-dependent `linop!` the accumulated phase is tabulated by **adaptive bisection**.
Each candidate interval is integrated with Simpson's rule and then checked against the thing
that actually matters — the cubic Hermite interpolant that will be used to read the table
back. The interpolant's error peaks at the interval midpoint, where it is compared against
the directly integrated value there, and the interval is bisected until that error is at most
`atol` radians (across every element of `linop`). The largest error that survives is recorded
in the returned object and reported.

This matters because Luna's `z`-dependent operators are usually *not* smooth. A pressure
gradient built by [`Capillary.gradient`](@ref Luna.Capillary.gradient) goes as
`√(p₀² + z/L(p₁² − p₀²))`, which has a `1/√z` cusp in its derivative at the entrance whenever
`p₀ = 0`, and a multi-section fill has a derivative discontinuity at every junction; tapers
are whatever function the user supplies. Uniform nodes converge at second order or worse
across such a feature, and give no way to tell how far off they are. Bisection puts nodes
only where they are needed — a kink costs a handful of extra intervals rather than a finer
grid everywhere — so the table is usually *smaller* than a uniform one of equal accuracy.

`zmax` is the propagation length (the grid does not carry it); a constant `linop` does not
need it. The table covers `[0, (1 + margin)·zmax]` rather than `[0, zmax]`, because `RK45.solve` runs
`while tn <= tmax` and so overshoots the end of the fibre on its final step -- and that step's
stages are what the last saved plane is interpolated from. `linop!` has to be defined a little
past `zmax` for this, which it already does: `RK45.make_prop!` evaluates it there too.
Anything beyond the margin still falls back to linear extrapolation.

`maxdepth` and `maxnodes` bound the work if `linop!` is genuinely discontinuous in `z`; if
either bound stops refinement before `atol` is met, a warning reports the error achieved.
"""
unitary_phase(linop::AbstractArray, grid, sz=size(linop); kwargs...) =
    ConstUnitaryPhase(imag.(linop))

function unitary_phase(linop!, grid, sz; zmax, atol=1e-6, margin=0.05, maxdepth=24,
                       maxnodes=4096)
    N = length(sz)
    buf = Array{ComplexF64}(undef, sz)
    neval = Ref(0)
    dat = function (z)
        neval[] += 1
        linop!(buf, z)
        imag.(buf)
    end
    zend = float(zmax)*(1 + margin)
    znodes = Float64[0.0]
    dnodes = Array{Float64, N}[dat(0.0)]
    deltas = Array{Float64, N}[] # accumulated phase across each accepted interval
    worst = Ref(0.0)
    _refine!(znodes, dnodes, deltas, worst, dat, 0.0, zend, dnodes[1], dat(zend), nothing,
             atol, 0, maxdepth, maxnodes)

    n = length(znodes)
    Φ = zeros(Float64, (sz..., n))
    dΦ = zeros(Float64, (sz..., n))
    d = N + 1
    selectdim(dΦ, d, 1) .= dnodes[1]
    for k = 2:n
        selectdim(Φ, d, k) .= selectdim(Φ, d, k-1) .+ deltas[k-1]
        selectdim(dΦ, d, k) .= dnodes[k]
    end
    if worst[] > atol
        @warn("Noise phase table did not reach its tolerance: $n nodes, "*
              "max interpolation error $(worst[]) rad against a tolerance of $atol rad. "*
              "The linear operator may be discontinuous in z.")
    else
        @info("Noise phase table: $n nodes, $(neval[]) linop evaluations, "*
              "max interpolation error $(worst[]) rad.")
    end
    TabulatedUnitaryPhase(znodes, Φ, dΦ, float(atol), worst[])
end

#= Bisect [a, b] until the cubic Hermite interpolant built from the endpoint values and
   derivatives is within atol of the true phase at the midpoint, which is where its error
   peaks. `dmid` is the already-evaluated derivative at the midpoint, if the caller has it
   (a bisection's two children each inherit one of the parent's quarter points), so each
   call costs two new evaluations of linop rather than three. =#
function _refine!(znodes, dnodes, deltas, worst, dat, a, b, da, db, dmid,
                  atol, depth, maxdepth, maxnodes)
    h = b - a
    m = a + h/2
    dm = isnothing(dmid) ? dat(m) : dmid
    dq1 = dat(a + h/4)
    dq2 = dat(a + 3h/4)
    ΔΦ = @. h/12*(da + 4dq1 + 2dm + 4dq2 + db) # Φ(b) - Φ(a), Simpson on each half
    err = 0.0
    for i in eachindex(ΔΦ)
        hermite = ΔΦ[i]/2 + h*(da[i] - db[i])/8 # Hermite at the midpoint, minus Φ(a)
        exact = h/12*(da[i] + 4dq1[i] + dm[i]) # Simpson over [a, m]
        err = max(err, abs(hermite - exact))
    end
    if err <= atol || depth >= maxdepth || length(znodes) >= maxnodes
        push!(znodes, b)
        push!(dnodes, db)
        push!(deltas, ΔΦ)
        worst[] = max(worst[], err)
        return
    end
    _refine!(znodes, dnodes, deltas, worst, dat, a, m, da, dm, dq1,
             atol, depth+1, maxdepth, maxnodes)
    _refine!(znodes, dnodes, deltas, worst, dat, m, b, dm, db, dq2,
             atol, depth+1, maxdepth, maxnodes)
end

end
