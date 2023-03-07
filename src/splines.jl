abstract type AbstractSplineComplex{T} end

"""
    SplineComplex{T}

1D B-Spline bases and collocation matrices

# Initialization 
By knot sequence or domain interval and number of points

* SplineComplex(knots::Vector{T}, order::Int=4)
* SplineComplex(domain::Tuple=(-1.0, 1.0), n::Int=11, order::Int=4)

# Fields

* `B`: Spline basis
* `R`: Recombined spline basis (satisfying Neumann boundary condition)
* `ξ`: Greville points 
* `I₁`: Collocation matrix of R
* `I₂`: Collocation matrix of B
"""
struct SplineComplex{T<:AbstractFloat} <: AbstractSplineComplex{T}
    B::AbstractBSplineBasis
    R::AbstractBSplineBasis
    ξ::Vector{T}
    I₁::SparseMatrixCSC{T}
    I₂::SparseMatrixCSC{T}
end

function SplineComplex(knots::Vector{T}, order::Int=4) where {T<:AbstractFloat}

    B = BSplineBasis(BSplineOrder(order), knots)
    R = RecombinedBSplineBasis(Derivative(1), B)

    ξ = collocation_points(B)
    I₁ = collocation_matrix(R, ξ, SparseMatrixCSC{T})
    I₂ = collocation_matrix(B, ξ, SparseMatrixCSC{T})

    return SplineComplex{T}(B, R, ξ, I₁, I₂)
end

function SplineComplex(domain::Tuple{T,T}=(-1.0, 1.0), n::Int=11, order::Int=4) where {T}
    knots = collect(range(domain[1], domain[2], length=n))
    return SplineComplex(knots, order)
end

"""
    approximate_nonlinear(a, f, SC::SplineComplex)

Return coefficients for approximation of nonlinear function f(x, Tₕ(a)) in spline space B.

The approximation is based on interpolation at the Greville points.
"""
function _approximate_nonlinear(a::Vector{T}, f::Function, SC::SplineComplex{T}) where {T}
    d = ones(T, length(SC.ξ))
    _approximate_nonlinear!(d, a, f, SC)
    return d
end

"""
    approximate_nonlinear!(d, a, f, SC::SplineComplex)

Fill preallocated coefficient vector, see also [`_approximate_nonlinear`](@ref).
"""
function _approximate_nonlinear!(d::Vector{T}, a::Vector{T}, f::Function, SC::SplineComplex{T}) where {T}
    d .= convert.(T, SC.I₂ \ f.(SC.ξ, SC.I₁ * a))
    return nothing
end

"""
    approximate_nonlinear(f, SC::SplineComplex)

Return coefficients for approximation of function f(x) in spline space B, see also [`_approximate_nonlinear`](@ref).
"""
function _approximate_nonlinear(f::Function, SC::SplineComplex{T}) where {T}
    d = convert.(T, SC.I₂ \ f.(SC.ξ))
    return d
end