function sol_to_matrix(SC::SplineComplex, sol::ODESolution, resolution::Int=100)
    x = range(minimum(SC.B.t), maximum(SC.B.t), length=resolution)
    data = zeros(length(sol.t), resolution)
    for (i, u) in enumerate(sol.u)
        S = Spline(SC.R, u)
        data[i, :] .= S.(x)
    end
    return data
end

function coeffs_to_fct(SC::SplineComplex, coeffs::Vector)
    if length(coeffs) != SC.R.M.M
        throw(DimensionMismatch("Coefficient vector (n=$(length(coeffs))) and spline complex (n=$(SC.R.M.M)) are not compatible."))
    end
    return Spline(SC.R, coeffs)
end

function eval_spline(SC::SplineComplex, coeffs::Vector, resolution::Int)
    return Spline(SC.R, coeffs).(range(minimum(SC.B.t), maximum(SC.B.t), length=resolution))
end

function get_grid(SC::SplineComplex, resolution::Int)
    return range(minimum(SC.B.t), maximum(SC.B.t), length=resolution)
end

MakieCore.convert_arguments(P::MakieCore.PointBased, SC::SplineComplex, coeffs::Vector, resolution::Int=100) = MakieCore.convert_arguments(P, get_grid(SC, resolution), eval_spline(SC, coeffs, resolution))

MakieCore.convert_arguments(P::Type{<:Heatmap}, SC::SplineComplex, sol::ODESolution, resolution::Int=100) = MakieCore.convert_arguments(P, sol_to_matrix(SC, sol, resolution))


MakieCore.@recipe(DiffusionPlot, SC, sol, resolution) do scene
    MakieCore.Theme(
        cmap=:viridis
    )
end

function MakieCore.plot!(diffusionplot::DiffusionPlot)
    sol = diffusionplot[:sol].val
    SC = diffusionplot[:SC].val
    color = cgrad(diffusionplot[:cmap].val)
    tmax = maximum(sol.t)
    for (i, t) in enumerate(sol.t)
        lines!(diffusionplot, SC, sol.u[i], diffusionplot[:resolution], color=get(color, t / tmax))
    end
    diffusionplot
end