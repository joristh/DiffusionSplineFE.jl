abstract type AbstractSplineComplex{k,T} end

"""
    SplineComplex{k, T}

1D B-Spline bases (order k with numeric type T) and collocation matrices

# Initialization 
By knot sequence or domain interval and number of points

* SplineComplex(knots::Vector{T}, order::Int=4)
* SplineComplex(domain::Tuple(T,T)=(-1.0, 1.0), n::Int=11, order::Int=4)

# Fields

* `B`: Spline basis
* `R`: Recombined spline basis (satisfying Neumann boundary condition)
* `ξ`: Greville points 
* `I₁`: Collocation matrix of R
* `I₂`: Collocation matrix of B
"""
struct SplineComplex{k,T<:AbstractFloat} <: AbstractSplineComplex{k,T}
    B::AbstractBSplineBasis{k,T}
    R::AbstractBSplineBasis{k,T}
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

    return SplineComplex{order,T}(B, R, ξ, I₁, I₂)
end

function SplineComplex(domain::Tuple{T,T}=(-1.0, 1.0), n::Int=11, order::Int=4) where {T}
    knots = collect(range(domain[1], domain[2], length=n))
    return SplineComplex(knots, order)
end

"""
    _approximate_nonlinear(SC::SplineComplex, f, a)

Return coefficients for approximation of nonlinear function f(x, Tₕ(a)) in spline space B.

The approximation is based on interpolation at the Greville points.
"""
function _approximate_nonlinear(SC::SplineComplex{k,T}, f, a::Vector{T}) where {k,T}
    d = ones(T, length(SC.ξ))
    _approximate_nonlinear!(d, SC, f, a)
    return d
end

"""
    _approximate_nonlinear!(d, SC::SplineComplex, f, a)

Fill preallocated coefficient vector, see also [`_approximate_nonlinear`](@ref).
"""
function _approximate_nonlinear!(d::Vector{T}, SC::SplineComplex{k,T}, f, a::Vector{T}) where {k,T}
    d .= convert.(T, SC.I₂ \ f.(SC.ξ, SC.I₁ * a))
    return nothing
end

"""
    _approximate_nonlinear(f, SC::SplineComplex)

Return coefficients for approximation of space-dependent function f(x) in spline space B.
"""
function _approximate_nonlinear(SC::SplineComplex{k,T}, f) where {k,T}
    d = convert.(T, SC.I₂ \ f.(SC.ξ))
    return d
end