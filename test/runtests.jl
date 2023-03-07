using DiffusionSplineFE
using Test
using OrdinaryDiffEq

@testset verbose = true "DiffusionSplineFE.jl" begin

    @testset "SplineComplex" begin
        SC1 = SplineComplex()
        SC2 = SplineComplex(collect(range(-1, 1, length=11)), 4)
        SC3 = SplineComplex((-1., 1.), 11, 4)

        @test SC1.B == SC2.B == SC3.B
        @test SC1.I₁ == SC2.I₁ == SC3.I₁
    end

    @testset "DiffusionProblem" begin
        SC = SplineComplex()

        # test calculation of initial coefficients
        @test isapprox(initial_coefficients(SC, one), ones(11))

        # test steady state
        u0 = initial_coefficients(SC, x -> x^2)
        tspan = (0, 10)
        prob = DiffusionProblem(SC, one, one, zero, u0, tspan)
        sol = solve(prob, Trapezoid(autodiff=false))
        @test all(isapprox.(sol(5), sol(10), atol=1e-6))

        # test source
        u0 = initial_coefficients(SC, zero)
        prob = DiffusionProblem(SC, one, one, one, u0, tspan)
        sol = solve(prob, Trapezoid(autodiff=false))
        @test all(isapprox.(sol[end], 10))

        #test construction of problem via initial function
        prob = DiffusionProblem(SC, one, one, one, zero, tspan)
        @test isapprox(prob.u0, u0)
    end

    @testset "GeneralDiffusionProblem" begin
        SC = SplineComplex()

        # test T-dependent source
        C(x) = 1.0
        D1(x, T) = 1.0
        S1(x, T) = -T

        u0 = initial_coefficients(SC, x -> x^2)
        tspan = (0, 100)
        prob = GeneralDiffusionProblem(SC, C, D1, S1, u0, tspan)
        sol = solve(prob, Trapezoid(autodiff=false))
        @test all(isapprox.(sol[end], 0, atol=1e-10))

        # test T-dependent diffusion
        D2(x, T) = T^2 + 0.01
        S2(x, T) = sin(π * x)

        tspan = (0, 100)
        prob = GeneralDiffusionProblem(SC, C, D2, S2, u0, tspan)
        sol = solve(prob, Trapezoid(autodiff=false))
        @test isapprox(sol[end], [-0.7508504346928699, -0.7818932132970585, -0.6562225354742912, -0.7235558365810306, 0.3311916198615287, 0.8192511341454509, 0.9084243118142071, 1.0413090506857137, 1.0834513751978398, 1.1088748754483329, 1.1067791284797306])
    end
end