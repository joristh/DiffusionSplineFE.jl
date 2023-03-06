function _assemble_tensors(SC::SplineComplex)
    A_tensor = galerkin_tensor((SC.R, SC.R, SC.B), (Derivative(1), Derivative(1), Derivative(0)))
    M_tensor = galerkin_tensor((SC.R, SC.R, SC.B), (Derivative(0), Derivative(0), Derivative(0)))
    S_matrix = galerkin_matrix((SC.R, SC.B));
    return M_tensor, A_tensor, S_matrix
end

function LinearDiffusionProblem(SC::SplineComplex, C::Function, D::Function, S::Function, u0::Vector, tspan::Tuple)
    M_tensor, A_tensor, S_matrix = _assemble_tensors(SC)
    
    c = approximate_nonlinear(C, SC)
    d = approximate_nonlinear(D, SC)
    s = approximate_nonlinear(S, SC)

    function f!(du, u, params, t)
        mul!(du, params.A, u)
        du .+= params.S
        ldiv!(du, params.M, du)
    end

    params = (
        M = factorize(M_tensor*c),
        A = A_tensor*(-d),
        S = S_matrix*s,
    )

    return ODEProblem(f!, u0, tspan, params)
end

function NonlinearDiffusionProblem(SC::SplineComplex, C::Function, D::Function, S::Function, u0::Vector, tspan::Tuple)
    M_tensor, A_tensor, S_matrix = _assemble_tensors(SC)
    
    c = approximate_nonlinear(C, SC)
    M = M_tensor*c

    function f!(du, u, params, t)
        d = approximate_nonlinear(u, D, SC)
        mul!(du, params.A_tensor*-d, u)
        s = approximate_nonlinear(u, S, SC)
        du .+= params.S_matrix*s
        ldiv!(du, params.M, du)
    end

    params = (
        M = factorize(M),
        A_tensor = A_tensor,
        S_matrix = S_matrix,
        D = D,
        S = S,
        SC = SC,
    )

    return ODEProblem(f!, u0, tspan, params)
end

