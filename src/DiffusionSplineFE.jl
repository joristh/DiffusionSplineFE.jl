module DiffusionSplineFE

using BSplineKit, LinearAlgebra, OrdinaryDiffEq, SparseArrays, MakieCore, RecipesBase, ColorSchemes, PlotUtils

include("splines.jl")
include("problem.jl")
#include("solve.jl")
include("plot.jl")

export SplineComplex, approximate_nonlinear, LinearDiffusionProblem, NonlinearDiffusionProblem, solve, plot

end
