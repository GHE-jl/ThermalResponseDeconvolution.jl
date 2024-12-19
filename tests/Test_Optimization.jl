using Optimization, OptimizationOptimJL
using SpecialFunctions

function example_Optimization()
    x₀ = [1.0, 1.0]
    p = P(5.0, -4.0)
    opt_fun = OptimizationFunction(obj_fun, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(opt_fun, x₀, p)
    sol = solve(prob, Optim.NelderMead(), g_tol=1e-3)
    return sol
end

function obj_fun(x::Vector{T}, p) where T<:Real
    return sum(abs.((p.x1 .* -expinti.(-p.x2 ./ collect(1:10))).-(x[1] .* -expinti.(-x[2] ./ collect(1:10)))))
end

struct P
    x1
    x2
end

@time sol = example_Optimization()
println("Optimal solution : ", sol.u)
println("Objective function value: ", sol.objective)