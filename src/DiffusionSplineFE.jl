module DiffusionSplineFE

using BSplineKit, LinearAlgebra, OrdinaryDiffEq, SparseArrays, MakieCore, ColorSchemes, PlotUtils

include("splines.jl")
include("problem.jl")
include("plot.jl")

export SplineComplex, initial_coefficients, DiffusionProblem, GeneralDiffusionProblem

end
