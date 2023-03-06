function sol_to_matrix(sc::SplineComplex, sol::ODESolution, resolution::Int)
    x = range(minimum(sc.B.t), maximum(sc.B.t), length=resolution)
    data = zeros(length(sol.t), resolution)
    for (i, u) in enumerate(sol.u)
        S = Spline(sc.R, u)
        data[i, :] .= S.(x)
    end
    return data
end

function eval_spline(sc::SplineComplex, coeffs::Vector, resolution::Int)
    return Spline(sc.R, coeffs).(range(minimum(sc.B.t), maximum(sc.B.t), length=resolution))
end

function get_grid(sc::SplineComplex, resolution::Int)
    return range(minimum(sc.B.t), maximum(sc.B.t), length=resolution)
end

MakieCore.convert_arguments(P::MakieCore.PointBased, sc::SplineComplex, coeffs::Vector, resolution::Int=100) = MakieCore.convert_arguments(P, get_grid(sc, resolution), eval_spline(sc, coeffs, resolution))

MakieCore.convert_arguments(P::Type{<:Heatmap}, sc::SplineComplex, sol::ODESolution, resolution::Int=100) = MakieCore.convert_arguments(P, sol_to_matrix(sc, sol, resolution))


MakieCore.@recipe(DiffusionPlot, sc, sol, resolution) do scene
    MakieCore.Theme(
        cmap = :viridis
    )
end

function MakieCore.plot!(diffusionplot::DiffusionPlot)
    sol = diffusionplot[:sol].val
    sc  = diffusionplot[:sc].val
    color = cgrad(diffusionplot[:cmap].val)
    tmax = maximum(sol.t)
    for (i,t) in enumerate(sol.t)
        lines!(diffusionplot, sc, sol.u[i], diffusionplot[:resolution], color=get(color, t/tmax))
    end
    diffusionplot
end

# default value for resolution??