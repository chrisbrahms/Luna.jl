"""
    Grid

The sampling grids on which `Luna` solves the propagation equation: the time and frequency
axes, the oversampled axes used for the nonlinear polarisation, and the apodisation windows
which go with them.

A grid also carries `sidx`, the indices of the frequency axis which are inside the
simulation band, and the profiles `ωwin`, `twin` and `towin` from which
[`Boundaries`](@ref Luna.Boundaries) derives the absorbing boundaries.
"""
module Grid
import Logging
import FFTW
import Hankel
import LinearAlgebra: mul!
import Printf: @sprintf
import Luna: PhysData, Maths

abstract type AbstractGrid end

abstract type TimeGrid <: AbstractGrid end
abstract type SpaceGrid <: AbstractGrid end

struct RealGrid <: TimeGrid
    zmax::Float64
    referenceλ::Float64
    t::Array{Float64, 1}
    ω::Array{Float64, 1}
    to::Array{Float64, 1}
    ωo::Array{Float64, 1}
    sidx::BitArray{1}
    ωwin::Array{Float64, 1}
    twin::Array{Float64, 1}
    towin::Array{Float64, 1}
end

"""
    RealGrid(zmax, referenceλ, λ_lims, trange; δt=1)

Time grid for simulations with real-valued (field-resolved) fields

# Arguments
- `zmax::Real` : Total distance to propagate
- `referenceλ::Real` : Reference wavelength (e.g. centre wavelength of input pulse)
- `λ_lims::Tuple{Real, Real}` : Wavelength limits of the frequency window
- `trange::Real` : Total extent of the time window required
- `δt::Real` : Sample spacing in time. The value actually used is either δt or the value
    required to satisfy `trange` and `λ_lims`, whichever is smaller.
"""
function RealGrid(zmax, referenceλ, λ_lims, trange, δt=1)
    f_lims = PhysData.c./λ_lims
    Logging.@info @sprintf("Freq limits %.2f - %.2f PHz", f_lims[2]*1e-15, f_lims[1]*1e-15)
    δto = min(1/(6*maximum(f_lims)), δt) # 6x maximum freq, or user-defined if finer
    samples = 2^(ceil(Int, log2(trange/δto))) # samples for fine grid (power of 2)
    trange_even = δto*samples # keep frequency window fixed, expand time window as necessary
    Logging.@info @sprintf("Samples needed: %.2f, samples: %d, δt = %.2f as",
                            trange/δto, samples, δto*1e18)
    Logging.@info @sprintf("Requested time window: %.1f fs, actual time window: %.1f fs", trange*1e15, trange_even*1e15)
    δωo = 2π/trange_even # frequency spacing for fine grid
    # Make fine grid
    Nto = collect(range(0, length=samples))
    to = @. (Nto-samples/2)*δto # centre on 0
    Nωo = collect(range(0, length=Int(samples/2 +1)))
    ωo = Nωo*δωo

    ωmin = 2π*minimum(f_lims)
    ωmax = 2π*maximum(f_lims)
    ωmax_win = 1.1*ωmax

    cropidx = findfirst(x -> x>ωmax_win, ωo)
    cropidx = 2^(ceil(Int, log2(cropidx))) + 1 # make coarse grid power of 2 as well
    ω = ωo[1:cropidx]
    δt = π/maximum(ω)
    tsamples = (cropidx-1)*2
    Nt = collect(range(0, length=tsamples))
    t = @. (Nt-tsamples/2)*δt

    # Make apodisation windows
    ωwindow = Maths.planck_taper(ω, ωmin/2, ωmin, ωmax, ωmax_win)

    twindow = Maths.planck_taper(t, minimum(t), -trange/2, trange/2, maximum(t))
    towindow = Maths.planck_taper(to, minimum(to), -trange/2, trange/2, maximum(to))

    # Indices to select real frequencies (for dispersion relation)
    sidx = (ω .> ωmin/2) .& (ω .< ωmax_win)

    @assert δt/δto ≈ length(to)/length(t)
    @assert δt/δto ≈ maximum(ωo)/maximum(ω)

    Logging.@info @sprintf("Grid: samples %d / %d, ωmax %.2e / %.2e",
                           length(t), length(to), maximum(ω), maximum(ωo))
    return RealGrid(float(zmax), referenceλ, t, ω, to, ωo, sidx, ωwindow, twindow, towindow)
end

function RealGrid(;zmax, referenceλ, t, ω, to, ωo, sidx, ωwin, twin, towin)
    RealGrid(zmax, referenceλ, t, ω, to, ωo, sidx, ωwin, twin, towin)
end

struct EnvGrid{T} <: TimeGrid
    zmax::Float64
    referenceλ::Float64
    ω0::Float64
    t::Array{Float64, 1}
    ω::Array{Float64, 1}
    to::Array{Float64, 1}
    ωo::Array{Float64, 1}
    sidx::T
    ωwin::Array{Float64, 1}
    twin::Array{Float64, 1}
    towin::Array{Float64, 1}
end

"""
    EnvGrid(zmax, referenceλ, λ_lims, trange; δt=1, thg=false)

Time grid for simulations with envelope (a.k.a. analytic) fields

# Arguments
- `zmax::Real` : Total distance to propagate
- `referenceλ::Real` : Reference wavelength (e.g. centre wavelength of input pulse)
- `λ_lims::Tuple{Real, Real}` : Wavelength limits of the frequency window
- `trange::Real` : Total extent of the time window required
- `δt::Real` : Sample spacing in time. The value actually used is either δt or the value
    required to satisfy `trange` and `λ_lims`, whichever is smaller.
- `thg::Bool` : Whether the grid should include space for the third hamonic (default: false)
"""
function EnvGrid(zmax, referenceλ, λ_lims, trange; δt=1, thg=false)
    fmin = PhysData.c/maximum(λ_lims)
    fmax = PhysData.c/minimum(λ_lims)
    fmax_win = 1.1*fmax # extended frequency window to accommodate apodisation
    ω0 = 2π*PhysData.c/referenceλ # central frequency
    Logging.@info @sprintf("Freq limits %.2f - %.2f PHz", fmin*1e-15, fmax*1e-15)
    if thg
        ffac = 6 # need to sample at 6x maximum desired frequency
        f_lims = (fmin, fmax)
    else
        ffac = 2 # need to sample at Nyquist limit for desired frequency only
        # Frequency window is not extended without THG, so need to extend additionally
        # to accommodate apodisation
        f_lims = (fmin, fmax_win)
    end
    rf_lims = f_lims .- PhysData.c/referenceλ # Relative frequency limits
    δt_f = 1/(ffac*maximum(rf_lims)) # time spacing as required by frequency window
    if δt_f <= δt
        # Demands of frequency window are more stringent
        δto = δt_f
        oversampling = thg # if no thg, then no oversampling
    else
        # User-defined time spacing is more stringent
        δto = δt
        oversampling = true
    end
    samples = 2^(ceil(Int, log2(trange/δto))) # samples for grid (power of 2)
    trange_even = δto*samples # keep frequency window fixed, expand time window as necessary
    Logging.@info @sprintf("Samples needed: %.2f, samples: %d, δt = %.2f as",
    trange/δto, samples, δto*1e18)
    δωo = 2π/trange_even # frequency spacing for grid
    # Make fine grid
    No = collect(range(0, length=samples))
    to = @. (No-samples/2)*δto # time grid, centre on 0
    vo = @. (No-samples/2)*δωo # freq grid relative to ω0
    vo = FFTW.fftshift(vo)
    ωo = vo .+ ω0

    ωmin = 2π*fmin
    ωmax = 2π*fmax
    ωmax_win = 2π*fmax_win

    # Find cropping area for coarse grid (contains frequencies of interest + apodisation)
    if oversampling
        cropidx = findfirst(x -> x>=(ωmax_win-δωo), ωo)
        cropidx = 2^(ceil(Int, log2(cropidx))) # make coarse grid power of 2 as well
        v = vcat(vo[1:cropidx], vo[end-cropidx+1:end])
    else
        v = vo
    end

    δt = -π/minimum(v)
    tsamples = length(v)
    Nt = collect(range(0, length=tsamples))
    t = @. (Nt - tsamples/2)*δt

    ω = v .+ ω0 # True frequency grid
    # Indices to select real frequencies (for dispersion relation)
    sidx = (ω .> ωmin/2) .& (ω .< ωmax_win)

    # Make apodisation windows
    ωwindow = Maths.planck_taper(ω, ωmin/2, ωmin, ωmax, ωmax_win)
    twindow = Maths.planck_taper(t, minimum(t), -trange/2, trange/2, maximum(t))
    towindow = Maths.planck_taper(to, minimum(to), -trange/2, trange/2, maximum(to))

    # Check that grids are correct
    @assert δt/δto ≈ length(to)/length(t)
    @assert δt/δto ≈ minimum(vo)/minimum(v) # FFT grid -> sample at -fs/2 but not +fs/2
    factor = Int(length(to)/length(t))
    zeroidx = findfirst(x -> x==0, to)
    # The time samples should be exactly the same (except fewer in t)
    @assert all(to[zeroidx:factor:end] .≈ t[t .>= 0])
    @assert all(to[zeroidx:-factor:1] .≈ t[t .<= 0][end:-1:1])

    return EnvGrid(float(zmax), referenceλ, ω0, t, ω, to, ωo, sidx, ωwindow, twindow, towindow)
end

function EnvGrid(;zmax, referenceλ, ω0, t, ω, to, ωo, sidx, ωwin, twin, towin)
    EnvGrid(zmax, referenceλ, ω0, t, ω, to, ωo, sidx, ωwin, twin, towin)
end

struct FreeGrid <: SpaceGrid
    x::Vector{Float64}
    y::Vector{Float64}
    kx::Vector{Float64}
    ky::Vector{Float64}
    r::Array{Float64, 3}
    xywin::Array{Float64, 3}
end

"""
    FreeGrid(Rx, Nx, Ry, Ny; window_factor=0.1)
    FreeGrid(R, N; window_factor=0.1)

Spatial grid for full 3D freespace propagation with `x`/`y` half-width `Rx`/`Ry` and
`Nx`/`Ny` samples. `window_factor` determines by how much the grid size is extended to fit
a filtering window.

If only `R` and `N` are given, it is assumed that `Rx = Ry = R` and `Nx = Ny = N`.
"""
function FreeGrid(Rx, Nx, Ry, Ny; window_factor=0.1)
    Rxw = Rx * (1 + window_factor)
    Ryw = Ry * (1 + window_factor)

    δx = 2Rxw/Nx
    nx = collect(range(0, length=Nx))
    x = @. (nx-Nx/2) * δx
    kx = 2π*FFTW.fftfreq(Nx, 1/δx)

    δy = 2Ryw/Ny
    ny = collect(range(0, length=Ny))
    y = @. (ny-Ny/2) * δy
    ky = 2π*FFTW.fftfreq(Ny, 1/δy)

    r = sqrt.(reshape(x, (1, Nx)).^2 .+ reshape(y, (1, 1, Ny)).^2)

    xwin = Maths.planck_taper(x, -Rxw, -Rx, Rx, Rxw)
    ywin = Maths.planck_taper(y, -Ryw, -Ry, Ry, Ryw)
    xywin = reshape(xwin, (1, length(xwin))) .* reshape(ywin, (1, 1, length(ywin)))

    FreeGrid(x, y, kx, ky, r, xywin)
end

FreeGrid(R, N) = FreeGrid(R, N, R, N)

struct Free2DGrid
    x::Vector{Float64}
    kx::Vector{Float64}
    xwin::Vector{Float64}
    r::Vector{Float64}
end

"""
    Free2DGrid(R, N; window_factor=0.1)

Spatial grid for 2D freespace propagation with `x` half-width `R` and
`N` samples. `window_factor` determines by how much the grid size is extended to fit
a filtering window.
"""
function Free2DGrid(R, N; window_factor=0.1)
    Rw = R * (1 + window_factor) # size including window

    δx = 2Rw/N
    n = collect(range(0, length=N))
    x = @. (n-N/2) * δx
    kx = 2π*FFTW.fftfreq(N, 1/δx)

    xwin = Maths.planck_taper(x, -Rw, -R, R, Rw)

    Free2DGrid(x, kx, xwin, copy(x))
end

#===================================================#
#================  RADIAL GRID  ====================#
#===================================================#

"""
    Grid.HankelTransform

Alias for `Hankel.QDHT`, the quasi-discrete Hankel transform of
[Hankel.jl](https://github.com/LupoLab/Hankel.jl). Luna's radial entry points accept one
for backwards compatibility and convert it to a [`RadialGrid`](@ref); `Hankel` is not used
anywhere else in `Luna`.
"""
const HankelTransform = Hankel.QDHT

"""
    RadialGrid(R, N; order=0)
    RadialGrid(q::HankelTransform)

Transverse grid for radially symmetric free-space propagation: a quasi-discrete Hankel
transform of order `order` over aperture radius `R` with `N` samples.

The grid owns the transform matrices and the integration weights, so that no per-step code
needs Hankel.jl. The matrices are laid out for multiplication *from the right* on an array
reshaped to `(:, N)`, i.e. they transform along the **last** dimension of an array of any
rank (see [`to_kspace!`](@ref), [`to_rspace!`](@ref)).

The second form converts an existing `Hankel.QDHT`. Its `dim` field is ignored, since a
`RadialGrid` always transforms along the last dimension.

# Fields
- `R`: aperture radius (the transform imposes `E(R) = 0`)
- `K`: largest transverse wavevector on the grid
- `N`: number of samples
- `order`: order of the Hankel transform (0 for a radially symmetric field)
- `r`, `k`: sample points in real and reciprocal space
- `Tfwd`, `Tbwd`: r→k and k→r transform matrices with the scale factors folded in
- `wr`, `wk`: integration weights in real and reciprocal space
    ([`integrate_r`](@ref), [`integrate_k`](@ref))
"""
struct RadialGrid <: SpaceGrid
    R::Float64
    K::Float64
    N::Int
    order::Int
    r::Vector{Float64}
    k::Vector{Float64}
    Tfwd::Matrix{Float64}
    Tbwd::Matrix{Float64}
    wr::Vector{Float64}
    wk::Vector{Float64}
end

#= The transform matrices are the transpose of Hankel's with its scalar scale factor folded
   in, exactly as NonlinearRHS.TransRadial built them for itself before. Hankel's matrix is
   symmetric, so right-multiplying a `(:, N)` reshape by these is the same transform as
   Hankel's left-multiplication along the radial axis, without the permutedims. =#
function _fromqdht(q::Hankel.QDHT)
    Hankel.sphericaldim(q) == 1 || error(
        "Only cylindrical (spherical dimension 1) Hankel transforms are supported, " *
        "got spherical dimension $(Hankel.sphericaldim(q))")
    Tfwd = Matrix{Float64}(transpose(q.T) .* q.scaleRK)
    Tbwd = Matrix{Float64}(transpose(q.T) ./ q.scaleRK)
    RadialGrid(float(q.R), float(q.K), q.N, Hankel.order(q),
               convert(Vector{Float64}, q.r), convert(Vector{Float64}, q.k),
               Tfwd, Tbwd,
               convert(Vector{Float64}, q.scaleR), convert(Vector{Float64}, q.scaleK))
end

RadialGrid(R, N; order::Int=0) = _fromqdht(Hankel.QDHT{order, 1}(R, N))

function RadialGrid(q::HankelTransform)
    Logging.@warn(
        "Hankel.QDHT is deprecated as a Luna transverse grid; use Grid.RadialGrid(R, N) " *
        "instead. The QDHT is being converted (its `dim` field is ignored).",
        maxlog=1)
    _fromqdht(q)
end

RadialGrid(;R, N, order=0) = RadialGrid(R, N; order=order)

RadialGrid(d::AbstractDict) = from_dict(RadialGrid, d)

"""
    radial_matmul!(out, A, T)

Apply the `N×N` matrix `T` along the last dimension of `A`, storing the result in `out`
(which must have the same size as `A`). This is one `mul!` on `reshape(A, :, N)`: no
permutation, no scalar loop, and no allocation unless `out === A`, in which case a
temporary copy of `A` is made.

`T` is one of a [`RadialGrid`](@ref)'s transform matrices, converted to the element type
being multiplied where that matters (see [`RadialGrid`](@ref)).
"""
function radial_matmul!(out, A, T)
    N = size(T, 1)
    size(T, 2) == N || throw(DimensionMismatch("transform matrix must be square"))
    size(A, ndims(A)) == N || throw(DimensionMismatch(
        "last dimension of the input is $(size(A, ndims(A))), expected $N"))
    size(out) == size(A) || throw(DimensionMismatch(
        "output size $(size(out)) does not match input size $(size(A))"))
    src = (out === A) ? copy(A) : A
    mul!(reshape(out, :, N), reshape(src, :, N), T)
    out
end

# permutation which swaps dimensions `dim` and `d`; it is its own inverse
function _swapperm(d, dim)
    perm = collect(1:d)
    perm[dim] = d
    perm[d] = dim
    Tuple(perm)
end

function _matmul_dim(rg::RadialGrid, A, T, dim)
    d = ndims(A)
    dim == d && return radial_matmul!(similar(A), A, T)
    perm = _swapperm(d, dim)
    Ap = permutedims(A, perm)
    permutedims(radial_matmul!(similar(Ap), Ap, T), perm)
end

"""
    to_kspace!(out, rg::RadialGrid, A)

Hankel transform `A` from real space to reciprocal space along its last dimension, storing
the result in `out`. `out === A` is allowed but then allocates a temporary.
"""
to_kspace!(out, rg::RadialGrid, A) = radial_matmul!(out, A, rg.Tfwd)

"""
    to_rspace!(out, rg::RadialGrid, A)

Inverse Hankel transform of `A` from reciprocal space to real space along its last
dimension, storing the result in `out`. `out === A` is allowed but then allocates a
temporary.
"""
to_rspace!(out, rg::RadialGrid, A) = radial_matmul!(out, A, rg.Tbwd)

"""
    to_kspace(rg::RadialGrid, A; dim=ndims(A))

Hankel transform `A` from real space to reciprocal space along dimension `dim`, which
defaults to the last. Allocates the output; for the per-step path use
[`to_kspace!`](@ref).
"""
to_kspace(rg::RadialGrid, A; dim=ndims(A)) = _matmul_dim(rg, A, rg.Tfwd, dim)

"""
    to_rspace(rg::RadialGrid, A; dim=ndims(A))

Inverse Hankel transform of `A` from reciprocal space to real space along dimension `dim`,
which defaults to the last. Allocates the output; for the per-step path use
[`to_rspace!`](@ref).
"""
to_rspace(rg::RadialGrid, A; dim=ndims(A)) = _matmul_dim(rg, A, rg.Tbwd, dim)

function _weighted_sum(A, w, dim)
    d = ndims(A)
    dim <= d || throw(DimensionMismatch(
        "cannot integrate along dimension $dim of a $d-dimensional array"))
    dim == d || (A = permutedims(A, _swapperm(d, dim)))
    N = length(w)
    size(A, d) == N || throw(DimensionMismatch(
        "dimension $dim of the input is $(size(A, dim)), expected $N"))
    out = reshape(A, :, N) * w
    d == 1 && return out[1]
    reshape(out, size(A)[1:d-1])
end

"""
    integrate_r(rg::RadialGrid, A; dim=ndims(A))

Radial integral of `A` over the aperture of `rg` in real space, along dimension `dim`
(the last by default), which is dropped from the result. A vector input gives a scalar.

Assuming `A` holds samples of ``f(r)`` at `rg.r`, this approximates ``\\int f(r) r dr``
from 0 to ∞. Together with [`integrate_k`](@ref) it fulfils Parseval's theorem:
`integrate_r(rg, abs2.(A))` equals `integrate_k(rg, abs2.(to_kspace(rg, A)))`.
"""
integrate_r(rg::RadialGrid, A; dim=ndims(A)) = _weighted_sum(A, rg.wr, dim)

"""
    integrate_k(rg::RadialGrid, A; dim=ndims(A))

Radial integral of `A` over the aperture of `rg` in reciprocal space, along dimension `dim`
(the last by default), which is dropped from the result. See [`integrate_r`](@ref).
"""
integrate_k(rg::RadialGrid, A; dim=ndims(A)) = _weighted_sum(A, rg.wk, dim)

"""
    onaxis(rg::RadialGrid, Ak; dim=ndims(Ak))

The on-axis (``r = 0``) sample of a field given in reciprocal space, obtained from
`Ak` by integration over `k` along dimension `dim`, which is dropped from the result.
Only defined for a 0-order grid.
"""
function onaxis(rg::RadialGrid, Ak; dim=ndims(Ak))
    rg.order == 0 || throw(DomainError(
        rg.order, "`onaxis` is only supported for 0-order RadialGrids"))
    integrate_k(rg, Ak; dim=dim)
end

"""
    kperp2(rg::RadialGrid)

The squared transverse wavevector `k⊥²` of every sample point of `rg`.
"""
kperp2(rg::RadialGrid) = rg.k .^ 2

"""
    symmetric(rg::RadialGrid, A)

Mirror `A`, sampled at `rg.r` along its last dimension, about the axis, inserting the
on-axis sample: the result is sampled at [`rsymmetric(rg)`](@ref rsymmetric), i.e.
`[...-r₂, -r₁, 0, r₁, r₂...]`, and is one longer than twice `rg.N` along that dimension.
Only defined for a 0-order grid.
"""
function symmetric(rg::RadialGrid, A)
    rg.order == 0 || throw(DomainError(
        rg.order, "`symmetric` is only supported for 0-order RadialGrids"))
    d = ndims(A)
    A0 = onaxis(rg, to_kspace(rg, A))
    A0 = reshape([A0;], (size(A)[1:d-1]..., 1))
    cat(reverse(A; dims=d), A0, A; dims=d)
end

"""
    rsymmetric(rg::RadialGrid)

The radial coordinate array which goes with [`symmetric(rg, A)`](@ref symmetric).
"""
rsymmetric(rg::RadialGrid) = vcat(-reverse(rg.r), 0.0, rg.r)

"""
    Grid.TransverseGrid

The transverse grids `Luna` accepts for free-space propagation: [`RadialGrid`](@ref),
[`FreeGrid`](@ref) and [`Free2DGrid`](@ref), plus [`HankelTransform`](@ref), which the
radial entry points convert to a `RadialGrid`.
"""
const TransverseGrid = Union{RadialGrid, FreeGrid, Free2DGrid, HankelTransform}


function to_dict(g::GT) where GT <: AbstractGrid
    d = Dict{String, Any}()
    for field in fieldnames(GT)
        d[string(field)] = getfield(g, field)
    end
    d
end

function from_dict(gridtype, d)
    kwargs = (Symbol(k) => v for (k, v) in pairs(d))
    grid = gridtype(;kwargs...)

    # Make sure the grid is valid
    validate(grid)
    return grid
end

RealGrid(d::AbstractDict) = from_dict(RealGrid, d)
EnvGrid(d::AbstractDict) = from_dict(EnvGrid, d)

#= The transverse grids are described by the few numbers they are built from rather than by
   their fields: the transform matrices and the Cartesian windows are large and are
   rebuilt exactly by the constructors. =#
to_dict(g::RadialGrid) = Dict{String, Any}("R" => g.R, "N" => g.N, "order" => g.order)
to_dict(g::FreeGrid) = Dict{String, Any}(
    "x" => g.x, "y" => g.y, "kx" => g.kx, "ky" => g.ky)
to_dict(g::Free2DGrid) = Dict{String, Any}("x" => g.x, "kx" => g.kx)

function validate(grid::TimeGrid)
    δt = grid.t[2] - grid.t[1]
    δto = grid.to[2] - grid.to[1]
    @assert δt/δto ≈ length(grid.to)/length(grid.t)
    @assert length(grid.towin) == length(grid.to)
    @assert length(grid.ωwin) == length(grid.ω)
    @assert length(grid.sidx) == length(grid.ω)
    if grid isa EnvGrid
        δω = grid.ω[2] - grid.ω[1]
        Δω = length(grid.ω)*δω
        @assert δt ≈ 2π/Δω
        @assert length(grid.t) == length(grid.ω)
        @assert length(grid.to) == length(grid.ωo)
    else
        Δω = maximum(grid.ω)
        @assert δt ≈ π/Δω
        @assert length(grid.t) == 2*(length(grid.ω)-1)
        @assert length(grid.to) == 2*(length(grid.ωo)-1)
    end
end

function validate(grid::RadialGrid)
    @assert length(grid.r) == grid.N
    @assert length(grid.k) == grid.N
    @assert length(grid.wr) == grid.N
    @assert length(grid.wk) == grid.N
    @assert size(grid.Tfwd) == (grid.N, grid.N)
    @assert size(grid.Tbwd) == (grid.N, grid.N)
    @assert issorted(grid.r)
    @assert 0 < grid.r[1] && grid.r[end] < grid.R
    # the real and reciprocal sample points and weights are the same numbers, rescaled
    @assert grid.k ≈ grid.r .* (grid.K/grid.R)
    @assert grid.wk ≈ grid.wr .* (grid.K/grid.R)^2
end

end
