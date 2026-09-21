module LinearOps
import FFTW
import Luna: Modes, Grid, PhysData, Maths, RK45
import Luna: upload_like, scalar
import Printf: @sprintf
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
#============  TABULATED OPERATORS  ==============#
#=================================================#

#= Adaptive tabulation of the z-dependent quantities a propagation needs at every stage:
   the integrated linear operator (`TabulatedLinop`), the propagation constant β
   (`TabulatedVector`) and the effective area (`TabulatedScalar`).

   All three are built by the same bisection: an interval is accepted when the
   interpolant that will actually be used to read the table back is within a tolerance of
   the directly computed value at the interval midpoint, which is where that interpolant's
   error peaks, and is bisected otherwise. Luna's z-dependent quantities are usually not
   smooth -- a pressure gradient built by `Capillary.gradient` goes as
   √(p₀² + z/L(p₁² − p₀²)), which has a 1/√z cusp in its derivative at the entrance when
   p₀ = 0, a multi-section fill has a derivative discontinuity at every junction, and a
   taper is whatever function the user wrote -- so uniform nodes converge at second order
   or worse across such a feature and give no way to tell how far off they are. Bisection
   puts nodes only where they are needed: a kink costs a number of intervals proportional
   to the depth it is resolved to, not to the resolution everywhere.

   The design is PR 440's `TabulatedUnitaryPhase` generalised from ∫imag(linop) dz to the
   full complex integrated operator, and to quantities which are tabulated rather than
   integrated. =#

"The default tolerance of the adaptive tabulation; see [`TabulatedLinop`](@ref)."
const DEFAULT_LINOP_TOL = 1e-6

"Largest number of nodes the adaptive tabulation will place before giving up and warning."
const DEFAULT_MAXNODES = 1024

"Deepest bisection the adaptive tabulation will go to before giving up and warning."
const DEFAULT_MAXDEPTH = 40

"""
Table size above which [`TabulatedLinop`](@ref) warns. A mode-averaged operator is a few
megabytes at any sensible tolerance; a multimode, radial or free-space one is the size of
the whole state, so the same node count costs `2·nnodes` times that.
"""
const TABLE_WARN_BYTES = 256*1024^2

_tabmax(x::Number) = abs(x)
_tabmax(x::AbstractArray) = isempty(x) ? 0.0 : maximum(abs, x)

# Stack the per-node values along a new trailing axis, which is the layout the readback
# broadcasts over (`selectdim` of the last dimension is a contiguous view on every backend).
_stack(vs::Vector{T}) where {T<:Number} = copy(vs)
function _stack(vs::Vector{<:AbstractArray})
    out = Array{eltype(vs[1])}(undef, (size(vs[1])..., length(vs)))
    for k in eachindex(vs)
        selectdim(out, ndims(out), k) .= vs[k]
    end
    out
end

#= Locate z in the node vector: the interval index, the normalised position in it and its
   width. Outside the table the nearest end node is returned exactly (s = 0 or 1), i.e. the
   quantity is held constant rather than extrapolated: a cubic Hermite run past its interval
   diverges, and holding the value is the safe thing to do.

   For the operator that is a guard and nothing more -- `Luna.run` builds the table over
   everything the stepper can ask for, and `phase!` warns if it is ever read outside it. For
   a value table it can be deliberate: `prop_capillary` tabulates `Aeff` over the fibre, and
   the statistics of the last step are recorded a fraction of a step past the end of it. =#
function _locate(zs::Vector{Float64}, z::Real)
    if z <= zs[1]
        return 1, 0.0, zs[2] - zs[1]
    elseif z >= zs[end]
        n = length(zs)
        return n-1, 1.0, zs[n] - zs[n-1]
    end
    k = searchsortedlast(zs, z)
    h = zs[k+1] - zs[k]
    k, (z - zs[k])/h, h
end

#=------------------------- the integrated linear operator -------------------------=#

"""
    TabulatedLinop(linop!, proto, z0, z1; tol, maxdepth, maxnodes)

The integrated linear operator `Φ(z) = ∫_{z0}^{z} linop(z') dz'` of the z-dependent
operator `linop!(out, z)`, tabulated on adaptively placed nodes covering `[z0, z1]` and
read back with a cubic Hermite interpolant. `proto` is the propagating field, whose array
type and element type the tables are built in (so they are device-resident for a device
run).

[`RK45.make_prop!`](@ref Luna.RK45.make_prop!) builds the interaction-picture propagator
`exp(Φ(t2) − Φ(t1))` from it. That is the exact propagator of the linear part over the
step, where the untabulated path uses `exp(linop(t2)·(t2 − t1))`, a one-point rule; it is
a different discretisation of the same equation, and it is opt-in
(`tabulate_linop=true` on [`Luna.run`](@ref)) for that reason.

# What is stored
`Φ` holds not the integral itself but its deviation from the straight line through the two
ends of the table,

```
Φ̃(z) = Φ(z) − L̄·(z − z0),   L̄ = Φ(z1)/(z1 − z0),
```

and `dΦ` holds `linop(z) − L̄`. The propagator adds `L̄·(t2 − t1)` back in the same
broadcast, so the result is unchanged in exact arithmetic while the tabulated numbers are
as small as the operator's *variation* along z rather than as large as its accumulated
phase. `make_linop` already works in a co-moving frame, so `max|Φ|` is tens to hundreds of
radians rather than thousands: the subtraction reduces the rounding error of
`Φ(t2) − Φ(t1)` in `Float32` by a measured factor of 4 to 11 over spans from 0.1 m to 10 m,
rather than being the difference between working and not working. It costs nothing -- a
cubic Hermite reproduces a linear function exactly, so subtracting the secant changes
neither the node placement nor the interpolation error -- and it is exact for a
z-independent operator, which then stores nothing at all.

# Fields
- `z`: the nodes, ascending, `z[1] == z0` and `z[end] == z1`
- `Φ`, `dΦ`: `(size(linop)..., length(z))`, as described above
- `secant`: `L̄`, the mean operator over `[z0, z1]`
- `z0`: the lower end of the table, the origin `Φ` is measured from
- `tol`: the tolerance the nodes were placed to satisfy
- `err`: the largest interpolation error measured while placing them
- `nevals`: how many times `linop!` was called to build the table
- `scale`: `maximum(abs, Φ̃)`, the size of the stored numbers

`maxdepth` and `maxnodes` bound the work if `linop!` is discontinuous in z; if either stops
the refinement before `tol` is met, a warning reports the error achieved. `quiet=true`
suppresses the summary the constructor otherwise logs, but not the warnings.
"""
struct TabulatedLinop{aT, sT}
    z::Vector{Float64}
    Φ::aT
    dΦ::aT
    secant::sT
    z0::Float64
    tol::Float64
    err::Float64
    nevals::Int
    scale::Float64
end

function TabulatedLinop(linop!, proto::AbstractArray, z0::Real, z1::Real;
                        tol=DEFAULT_LINOP_TOL, maxdepth=DEFAULT_MAXDEPTH,
                        maxnodes=DEFAULT_MAXNODES, quiet=false)
    z0, z1 = float(z0), float(z1)
    z1 > z0 || error("TabulatedLinop needs z1 > z0, got $z0 and $z1")
    sz = size(proto)
    buf = Array{ComplexF64}(undef, sz)
    nevals = Ref(0)
    dat = function (z)
        nevals[] += 1
        linop!(buf, z)
        copy(buf)
    end
    znodes = Float64[z0]
    dnodes = Array{ComplexF64, length(sz)}[dat(z0)]
    deltas = Array{ComplexF64, length(sz)}[] # ∫linop dz across each accepted interval
    worst = Ref(0.0)
    _refine_int!(znodes, dnodes, deltas, worst, dat, z0, z1, dnodes[1], dat(z1), nothing,
                 float(tol), 0, maxdepth, maxnodes)

    n = length(znodes)
    Φ = zeros(ComplexF64, (sz..., n))
    d = length(sz) + 1
    for k = 2:n
        selectdim(Φ, d, k) .= selectdim(Φ, d, k-1) .+ deltas[k-1]
    end
    dΦ = _stack(dnodes)
    # Subtract the secant (see the docstring): Φ̃(z0) = Φ̃(z1) = 0 by construction.
    secant = selectdim(Φ, d, n)./(z1 - z0)
    for k = 1:n
        selectdim(Φ, d, k) .-= secant.*(znodes[k] - z0)
        selectdim(dΦ, d, k) .-= secant
    end
    scale = _tabmax(Φ)
    #= Two tables of `(size(linop)..., nnodes)` in the state's element type. For a
       mode-averaged operator that is megabytes; for a multimode, radial or free-space one
       the operator is the size of the whole state and this is `2·nnodes` times it. =#
    bytes = 2*length(Φ)*sizeof(Complex{real(eltype(proto))})
    if worst[] > tol
        @warn("The tabulated linear operator did not reach its tolerance: $n nodes, "*
              "largest interpolation error $(worst[]) against a tolerance of $tol. The "*
              "operator may be discontinuous in z; raise `linop_tol` or check it.")
    elseif !quiet
        @info(@sprintf("Tabulated linear operator: %d nodes over [%.4g, %.4g] m, %d evaluations, largest interpolation error %.2e, largest stored value %.2e, %.1f MB.",
                       n, z0, z1, nevals[], worst[], scale, bytes/1024^2))
    end
    if bytes > TABLE_WARN_BYTES
        @warn(@sprintf("The tabulated linear operator needs %.1f MB: %d nodes of an operator of size %s. Tabulation stores two copies of the operator per node, which is cheap for a mode-averaged run and not for a multimode or free-space one. Raise `linop_tol` to place fewer nodes, or leave `tabulate_linop` off for this geometry.",
                       bytes/1024^2, n, string(sz)))
    end
    TabulatedLinop(znodes, upload_like(proto, Φ), upload_like(proto, dΦ),
                   upload_like(proto, secant), z0, float(tol), worst[], nevals[], scale)
end

#= Bisect [a, b] until the cubic Hermite interpolant built from the endpoint values and
   derivatives is within `tol` of the true integral at the midpoint, which is where its
   error peaks. `dmid` is the already-evaluated derivative at the midpoint when the caller
   has it -- a bisection's two children each inherit one of the parent's quarter points --
   so each call costs two new evaluations of `linop!` rather than three. =#
function _refine_int!(znodes, dnodes, deltas, worst, dat, a, b, da, db, dmid,
                      tol, depth, maxdepth, maxnodes)
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
    if err <= tol || depth >= maxdepth || length(znodes) >= maxnodes
        push!(znodes, b)
        push!(dnodes, db)
        push!(deltas, ΔΦ)
        worst[] = max(worst[], err)
        return
    end
    _refine_int!(znodes, dnodes, deltas, worst, dat, a, m, da, dm, dq1,
                 tol, depth+1, maxdepth, maxnodes)
    _refine_int!(znodes, dnodes, deltas, worst, dat, m, b, dm, db, dq2,
                 tol, depth+1, maxdepth, maxnodes)
end

"""
    phase!(out, tab::TabulatedLinop, z)

Fill `out` with the stored deviation `Φ̃(z) = Φ(z) − L̄·(z − z0)` (see
[`TabulatedLinop`](@ref)), as one broadcast over four node slices with four scalar
weights. Returns `out`.
"""
function phase!(out, tab::TabulatedLinop, z)
    #= Not an error, because it is recoverable; not silent, because the propagator adds the
       secant term whatever `phase!` returns, so a step taken outside the table would
       propagate with the mean operator over the whole table -- a plausible-looking wrong
       answer rather than a failure. `Luna.run` builds the table over every z the stepper
       can reach, so this cannot fire from there. =#
    if z < tab.z[1] || z > tab.z[end]
        @warn(@sprintf("The tabulated linear operator was read at z = %.6g m, outside the [%.6g, %.6g] m it was built for; the value at the nearest end is used and the propagator adds the mean operator over the table.",
                       z, tab.z[1], tab.z[end]), maxlog=1)
    end
    k, s, h = _locate(tab.z, z)
    s2 = s*s
    s3 = s2*s
    w00 = scalar(out, 2s3 - 3s2 + 1)
    w10 = scalar(out, h*(s3 - 2s2 + s))
    w01 = scalar(out, -2s3 + 3s2)
    w11 = scalar(out, h*(s3 - s2))
    d = ndims(tab.Φ)
    Φk = selectdim(tab.Φ, d, k)
    Φk1 = selectdim(tab.Φ, d, k+1)
    dk = selectdim(tab.dΦ, d, k)
    dk1 = selectdim(tab.dΦ, d, k+1)
    @. out = w00*Φk + w10*dk + w01*Φk1 + w11*dk1
    out
end

"""
    integrated!(out, tab::TabulatedLinop, z)

Fill `out` with the integrated operator `Φ(z) = ∫_{z0}^{z} linop dz'` itself. Only
differences of `Φ` enter the propagator, which forms them from [`phase!`](@ref) without
ever building this; this is for checking the table against a direct integration.
"""
function integrated!(out, tab::TabulatedLinop, z)
    phase!(out, tab, z)
    dz = scalar(out, z - tab.z0)
    sec = tab.secant
    @. out += sec*dz
    out
end

"""
    RK45.make_prop!(tab::TabulatedLinop, y0)

The interaction-picture propagator of a tabulated z-dependent operator,
`y *= exp(Φ(t2) − Φ(t1))`.

`Φ(t1)` and `Φ(t2)` are each read out of the table only when their argument changes: the
six stages of a step share one `t1`, and `t2` repeats between the forward and the backward
propagation of each stage, which is the same argument the untabulated propagator's
last-`t2` cache rests on. Both readbacks and the exponential are broadcasts over the
state's own array type, so nothing here touches the host on a device run.
"""
function RK45.make_prop!(tab::TabulatedLinop, y0)
    Φ1 = similar(y0)
    Φ2 = similar(y0)
    lastt1 = Ref(NaN)
    lastt2 = Ref(NaN)
    prop! = let tab=tab, Φ1=Φ1, Φ2=Φ2, lastt1=lastt1, lastt2=lastt2
        function prop!(y, t1, t2, bwd=false)
            if lastt1[] != t1
                phase!(Φ1, tab, t1)
                lastt1[] = t1
            end
            if lastt2[] != t2
                phase!(Φ2, tab, t2)
                lastt2[] = t2
            end
            #= The secant term is put back here rather than in the readback: it is the
               whole of a constant operator's contribution and is formed exactly as the
               constant propagator forms it, from the step length. =#
            dt = scalar(y, bwd ? (t1 - t2) : (t2 - t1))
            sec = tab.secant
            if bwd
                @. y *= exp(Φ1 - Φ2 + sec*dt)
            else
                @. y *= exp(Φ2 - Φ1 + sec*dt)
            end
        end
    end
    return prop!
end

#=--------------------------- tabulated values (not integrals) ---------------------------=#

#= β and Aeff are needed as values, not as integrals, and no derivative of either is
   available: the mode interface gives the quantity and nothing else. The interpolant is
   therefore linear rather than cubic Hermite, and the acceptance check compares the linear
   interpolant's midpoint value against the true one -- the same bisection, measuring the
   error of the interpolant that is actually used. Second order in the interval width
   instead of fourth, which for quantities that only scale the nonlinear polarisation costs
   a handful of extra nodes and nothing else.

   The tolerance is relative to the largest value in the table, because β is ~1e7 in SI
   units and Aeff ~1e-8, and one absolute tolerance cannot serve both. =#
function _refine_val!(znodes, vnodes, worst, dat, a, b, va, vb, tol, depth,
                      maxdepth, maxnodes)
    h = b - a
    m = a + h/2
    vm = dat(m)
    err = _tabmax(@. vm - (va + vb)/2)
    if err <= tol || depth >= maxdepth || length(znodes) >= maxnodes
        push!(znodes, b)
        push!(vnodes, vb)
        worst[] = max(worst[], err)
        return
    end
    _refine_val!(znodes, vnodes, worst, dat, a, m, va, vm, tol, depth+1, maxdepth, maxnodes)
    _refine_val!(znodes, vnodes, worst, dat, m, b, vm, vb, tol, depth+1, maxdepth, maxnodes)
end

#= Place the nodes for a value table. `f(z)` returns the quantity (a scalar or an array);
   `rtol` is relative to the largest value seen at the two ends. =#
function _tabulate_value(f, z0, z1, rtol, maxdepth, maxnodes)
    z0, z1 = float(z0), float(z1)
    z1 > z0 || error("a z table needs z1 > z0, got $z0 and $z1")
    nevals = Ref(0)
    dat = function (z)
        nevals[] += 1
        f(z)
    end
    #= The relative tolerance is taken against the two ends of the interval. That is enough
       for β and Aeff, which are monotonic in z over any fibre Luna describes, so neither is
       small at both ends and large in between; a quantity which was would be tabulated to
       an inappropriate absolute tolerance and would need the maximum over the nodes as
       they are placed instead. =#
    v0, v1 = dat(z0), dat(z1)
    scale = max(_tabmax(v0), _tabmax(v1))
    scale == 0 && (scale = 1.0)
    znodes = Float64[z0]
    vnodes = typeof(v0)[v0]
    worst = Ref(0.0)
    _refine_val!(znodes, vnodes, worst, dat, z0, z1, v0, v1, rtol*scale, 0,
                 maxdepth, maxnodes)
    if worst[] > rtol*scale
        @warn("A z table did not reach its tolerance: $(length(znodes)) nodes, largest "*
              "interpolation error $(worst[]/scale) relative against a tolerance of $rtol.")
    end
    znodes, vnodes, worst[]/scale, nevals[]
end

"""
    TabulatedScalar(f, z0, z1; tol, maxdepth, maxnodes)

A scalar function of `z` -- the effective area of a tapered or pressure-graded waveguide --
tabulated on adaptively placed nodes covering `[z0, z1]` and read back by linear
interpolation. Callable as `t(z)`.

`Modes.Aeff` is memoised on `(mode, z)`, so a z-dependent mode grows one `Dict` entry per
distinct `z` it is asked about, i.e. per stage of every step, for the whole propagation.
Tabulating it bounds that at the number of nodes and takes the quadrature out of the step.

`src` is the callable the table was built from, kept so that a table can be rebuilt over a
wider span without going through the interpolant. `NonlinearRHS.tabulate` does that when
[`Luna.prop_capillary`](@ref) has already tabulated `Aeff` over the fibre for the
statistics and the propagation needs it a little past the end.
"""
struct TabulatedScalar{F}
    z::Vector{Float64}
    f::Vector{Float64}
    src::F
    tol::Float64
    err::Float64
    nevals::Int
end

function TabulatedScalar(f, z0, z1; tol=DEFAULT_LINOP_TOL, maxdepth=DEFAULT_MAXDEPTH,
                         maxnodes=DEFAULT_MAXNODES)
    zs, vs, err, nevals = _tabulate_value(z -> float(f(z)), z0, z1, float(tol),
                                          maxdepth, maxnodes)
    TabulatedScalar(zs, _stack(vs), f, float(tol), err, nevals)
end

function (t::TabulatedScalar)(z)
    k, s, _ = _locate(t.z, z)
    (1 - s)*t.f[k] + s*t.f[k+1]
end

"""
    TabulatedVector(f!, proto, n, z0, z1; tol, maxdepth, maxnodes)

A vector-valued function of `z` -- the propagation constant `β(z)` the mode-averaged
normalisation divides by -- tabulated on adaptively placed nodes covering `[z0, z1]` and
read back by linear interpolation into a buffer on `proto`'s array type and precision.
`f!(out, z)` fills a length-`n` host buffer; `proto` is the propagating field.

Callable as `t(z)`, returning the buffer. The buffer is reused, and the readback is skipped
when `z` has not changed since the last call.

This replaces the per-evaluation host call and upload of
[`Luna.HostMirror`](@ref): with a table there is no host work left inside the right-hand
side of a tapered or pressure-graded run.
"""
struct TabulatedVector{aT, bT}
    z::Vector{Float64}
    f::aT
    buf::bT
    lastz::Base.RefValue{Float64}
    tol::Float64
    err::Float64
    nevals::Int
end

function TabulatedVector(f!, proto::AbstractArray, n::Integer, z0, z1;
                         tol=DEFAULT_LINOP_TOL, maxdepth=DEFAULT_MAXDEPTH,
                         maxnodes=DEFAULT_MAXNODES)
    host = zeros(Float64, n)
    f = function (z)
        f!(host, z)
        copy(host)
    end
    zs, vs, err, nevals = _tabulate_value(f, z0, z1, float(tol), maxdepth, maxnodes)
    tab = upload_like(proto, _stack(vs))
    buf = similar(proto, real(eltype(tab)), (n,))
    TabulatedVector(zs, tab, buf, Ref(NaN), float(tol), err, nevals)
end

function (t::TabulatedVector)(z)
    if t.lastz[] != z
        k, s, _ = _locate(t.z, z)
        w0 = scalar(t.buf, 1 - s)
        w1 = scalar(t.buf, s)
        fk = selectdim(t.f, 2, k)
        fk1 = selectdim(t.f, 2, k+1)
        out = t.buf
        @. out = w0*fk + w1*fk1
        t.lastz[] = z
    end
    t.buf
end

end
