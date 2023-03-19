"""
    _assemble_tensors(SC::SplineComplex)

Helper function to compute all needed Galerkin tensors and matrices.
"""
function _assemble_tensors(SC::SplineComplex{k,T}) where {k,T}
    A_tensor = galerkin_tensor((SC.R, SC.R, SC.B), (Derivative(1), Derivative(1), Derivative(0)), T)
    M_tensor = galerkin_tensor((SC.R, SC.R, SC.B), (Derivative(0), Derivative(0), Derivative(0)), T)
    S_matrix = galerkin_matrix((SC.R, SC.B), SparseMatrixCSC{T})
    return M_tensor, A_tensor, S_matrix
end

"""
    initial_coefficients(SC::SplineComplex, f::Function)

Compute coefficient vector of initial condition in recombined spline space R.
"""
function initial_coefficients(SC::SplineComplex{k,T}, f::Function) where {k,T}
    coefficients(approximate(f, SC.R))
end

"""
    DiffusionProblem(SC::SplineComplex, C::Function, D::Function, S::Function, u0::Vector, tspan::Tuple)

Set up OrdinaryDiffEq.ODEProblem for equation with only space-dependent capacity C, diffusion D, and source S.

Alternative initialization with u0::Function instead of u0::Vector.
"""
function DiffusionProblem(SC::SplineComplex{k,T}, C::Function, D::Function, S::Function, u0::Vector, tspan::Tuple) where {k,T}

    if length(u0) != SC.R.M.M
        throw(DimensionMismatch("Initial coefficient vector (n=$(length(u0))) and spline complex (n=$(SC.R.M.M)) are not compatible."))
    end

    u0 = convert.(T, u0)

    M_tensor, A_tensor, S_matrix = _assemble_tensors(SC)

    c = _approximate_nonlinear(SC, C)
    d = _approximate_nonlinear(SC, D)
    s = _approximate_nonlinear(SC, S)

    function f!(du, u, params, t)
        mul!(du, params.A, u)
        du .+= params.S
        ldiv!(du, params.M, du)
    end

    params = (
        M=factorize(M_tensor * c),
        A=A_tensor * (-d),
        S=S_matrix * s,
    )

    return ODEProblem(f!, u0, tspan, params)
end

function DiffusionProblem(SC::SplineComplex, C::Function, D::Function, S::Function, u0::Function, tspan::Tuple)
    u0_vector = initial_coefficients(SC, u0)
    DiffusionProblem(SC, C, D, S, u0_vector, tspan)
end

"""
    DiffusionProblem(SC::SplineComplex, C::Function, D::Function, S::Function, u0::Vector, tspan::Tuple)

Set up OrdinaryDiffEq.ODEProblem for equation with space- and value-dependent diffusion D, and source S.

Alternative initialization with u0::Function instead of u0::Vector.
"""
function GeneralDiffusionProblem(SC::SplineComplex{k,T}, C::Function, D::Function, S::Function, u0::Vector, tspan::Tuple) where {k,T}

    if length(u0) != SC.R.M.M
        throw(DimensionMismatch("Initial coefficient vector (n=$(length(u0))) and spline complex (n=$(SC.R.M.M)) are not compatible."))
    end

    u0 = convert.(T, u0)

    M_tensor, A_tensor, S_matrix = _assemble_tensors(SC)

    c = _approximate_nonlinear(SC, C)
    M = M_tensor * c

    function f!(du, u, params, t)
        d = _approximate_nonlinear(SC, D, u)
        mul!(du, params.A_tensor * -d, u)
        s = _approximate_nonlinear(SC, S, u)
        du .+= params.S_matrix * s
        ldiv!(du, params.M, du)
    end

    params = (
        M=factorize(M),
        A_tensor=A_tensor,
        S_matrix=S_matrix,
        D=D,
        S=S,
        SC=SC,
    )

    return ODEProblem(f!, u0, tspan, params)
end

function GeneralDiffusionProblem(SC::SplineComplex, C::Function, D::Function, S::Function, u0::Function, tspan::Tuple)
    u0_vector = initial_coefficients(SC, u0)
    GeneralDiffusionProblem(SC, C, D, S, u0_vector, tspan)
end

