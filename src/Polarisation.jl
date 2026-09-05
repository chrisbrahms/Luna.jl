module Polarisation

using StaticArrays
using LinearAlgebra

"""
    H()

Make horizontally polarised Jones vector.
"""
function H()
   SVector(1.0, 0.0)
end

"""
    LP()

Make horizontally orientated linear polariser Jones matrix.
"""
function LP()
    @SMatrix [ 1.0  0.0 ;
               0.0  0.0 ]
end

"""
    WP(ϕ)

Arbitrary waveplate with phase `ϕ`; fast axis horizontal.
"""
function WP(ϕ)
    @SMatrix [ exp(im*ϕ/2)  0.0 ;
               0.0          exp(-im*ϕ/2) ]
end

"""
    rot(θ)

Rotation operator.
"""
function rot(θ)
    @SMatrix [ cos(θ)  sin(θ) ;
              -sin(θ)  cos(θ) ]
end

"""
    rotate(J, θ)

Rotate Jones matrix `J` by `θ`.
"""
function rotate(J, θ)
    rot(-θ)*J*rot(θ)
end

"""
    Stokes(E; normalise=false)

Get Stokes parameters for input field `E = (Ex, Ey)`.
"""
function Stokes(E; normalise=false)
    Ex = E[1]
    Ey = E[2]
    I = abs(Ex)^2 + abs(Ey)^2
    Q = abs(Ex)^2 - abs(Ey)^2
    U = 2*real(Ex*conj(Ey))
    V = -2*imag(Ex*conj(Ey))
    S = SVector(I, Q, U, V)
    if normalise
        S = S./S[1]
    end
    S
end

"""
    cartesian(S)

Get normalised cartesian coordinates from Stokes parameters
(for Poincare sphere).
"""
function cartesian(S)
    S[2:end]./S[1]
end

"""
    ellipse(S)

Get polarisation ellipse parameters from Stokes parameters.
"""
function ellipse(S)
    I = S[1]
    Q = S[2]
    U = S[3]
    V = S[4]
    aL = sqrt(Q^2 + U^2)
    θ = angle(Q + 1im*U)/2
    A = sqrt((I + aL)/2)
    B = sqrt((abs(I - aL))/2)
    h = sign(V)
    A, B, θ, h
end

"""
    ellipticity(S)

Calculate ellipticity from Stokes parameters.
"""
function ellipticity(S)
    A, B, θ, h = ellipse(S)
    r = A/B
    r > 1 ? 1/r : r
end

end
