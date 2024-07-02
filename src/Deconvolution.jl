using Optim
using SpecialFunctions
using SparseArrays
using PCHIPInterpolation
using Plots

function deconvolution(
    t::Vector{Float64},
    f::Vector{Float64},
    temperature::Vector{Float64},
    n::Int=35,
    c::Int=2,
    show_ini::Bool=false,
    show_opt::Bool=false,
)
    #=
    Optimization algorithm to perform deconvolution on the experimental data
    extracted from a TRT to recover a short-term g-function (STgF).The function
    computes the estimated STgF and the estimated convolved
    temperature variations.
    Inputs:
      - t: Time array, starting at 0, with constant time step (see [1]) (s)
      - f: Incremental perturbation function (see note [2] below) [degC, W, W/m]
      - TExp: Experimental temperature variation (T_exp = T_out-T_0) [degC]
      - [Optionnal Inputs]:
          - n: Number of nodes (default: 35)
          - c: Choice of constraints (default: 2):
              0: No constraint
              1: Positivity
              2: Positivity and negative derivative slope
         - show_ini: Display initial guess figure (default: false)
         - show_opt: Display final guess figure (default: false)
    Outputs:
        - g: Interpolated estimated STgF obtained by deconvolution [-]
       - T: Estimated convolved temperature variation (T_conv=f*ghat) [degC]
       - [Optionnal Outputs]:
           - Time: Computing time to converge [s]
           - [Id,gOpt]: Node positions and optimized transfer function values
           - Out: Output variables of the solver [-]
           - IterOut: Optimization variables at each iteration [-]
               - [Iteration,Time,fval,funccount]

     Notes:
     [1]: The inputs "f" and "TExp" must have constant time step for the
     deconvolution to work effectively. Otherwise, these vectors must be
     interpolated.
     [2]: The input "f" is the time derivative of the perturbation function, which
     can be $T_{in}-T_{out}$ [$^\circ$C], $Q$ [W] or $q$ [W/m]. The transfer function 
     will be defined for the units of input used.

     Reference: Dion, G., Pasquier, P., & Marcotte, D. (2022). Deconvolution of
     experimental thermal response test data to recover short-term g-function.
     Geothermics, 100, 102302. https://doi.org/10.1016/j.geothermics.2021.102302

     Author: Gabriel Dion
     Date: 04-2024
    =#

    # 0. Data validation and Initialization
    data_validation(t, f, temperature)
    option_validation(n, c, show_ini, show_opt)

    # 1. Initial parameters
    nt = length(temperature)
    idall = collect(1:nt)

    # 2. Compute initial solution
    g_0 = deconv_ini(idall, f, temperature)
    show_ini ? show_fig(t, f, g_0, temperature) : nothing

    # 3. Define nodes positions
    id = set_nodes(nt, n)

    # 4. Set weights for the multi-objective function
    w = set_weights(g_0, f, temperature)

    # 5. Set linear inequality constraints and bounds
    lb = zeros(n, 1)

    a, b = const_derivative(t, id, c)
    # Note: Usually, the 2 constraints give best results. Either using 1 or 0 constraints
    # may be usefull in some cases.

    # 6. Optimization
    x_opt = optimize(x -> obj_fun(x, id, idall, f, temperature, w, w_array),
        TwiceDifferentiableConstraints(a, b, lb, fill(Inf, length(lb))),
        g_0,
        LBFGS(),
        Optim.Options(
            iterations=1000,
            show_trace=true)
        )
    
    # 7. Final interpolation and convolution
    interp = Interpolator(id, x_opt.minimizer)
    g = interp(idall)
    T = convolution(f, g)
    
    # 8. Plot final result
    show_opt ? show_fig(t, f, g, T) : nothing
    return g, T, x_opt
end

function data_validation(t::Vector{Float64}, f::Vector{Float64}, T::Vector{Float64})
    #=

    =#

    # Check if vectors are all the same length
    if length(t) != length(f) || length(t) != length(T)
        error("inputs are not of the same lengths.")

        # Check for imaginary, NaN or infinite value
    elseif any(!isreal, t) || any(!isreal, f) || any(!isreal, T)
        error("inputs contain imaginary values.")

        # Check for NaN values
    elseif any(isnan, t) || any(isnan, f) || any(isnan, T)
        error("inputs contains NaN values.")
        # Check for NaN values
    elseif any(isinf, t) || any(isinf, f) || any(isinf, T)
        error("inputs contains infinite values.")
    end
end

function option_validation(n, c, show_ini, show_opt)
    #=

    =#
    if n < 0
        error("number of nodes must be a positive integer.")
    elseif n < 15
        @warn println(
            "the number of nodes should be at least larger than 20, ideally 30 to 60.",
        )
    elseif n > 100
        @warn println(
            "the number of nodes is large (>100) and should be arount 60 at the maximum.",
        )
    elseif c < 0 || c > 2
        error("optionnal input for constraints must be 0, 1 or 2.")
    elseif !isinteger(c)
        error("optionnal constraints must be an integer.")
        # elseif tol #TODO: when I have the optimization, add constraints for this
    elseif typeof(show_ini) != Bool
        error("input either `true` or `false` to show or not the initial fit figure.")
    elseif typeof(show_opt) != Bool
        error("input either `true` or `false` to show or not the optimized fit figure.")
    end
end

function deconv_ini(idall::Vector{Int}, f::Vector{Float64},
    temperature::Vector{Float64})
    #= 

    =#

    # Define objective function (RMSE) between model and experimental values.
    obj_fun(x) = sqrt(sum((convolution(f, x[1] .* expint.(x[2] ./ idall))
                           .-
                           (temperature .- temperature[1])) .^ 2) / length(idall))

    # Define the optimization on 2 variables
    x0 = [1.0, 1.0]
    x_opt = optimize(obj_fun, x0, NelderMead(), Optim.Options(show_trace=true))

    return x_opt.minimizer[1] .* expint.(x_opt.minimizer[2] ./ idall)
end

function set_nodes(nt::Int, n0::Int)
    #=
    Function that sets the position of the nodes on the transfer function based on
    the number of nodes asked by the user.
    Inputs:
        - nt: Total number of data in the input vectors [-]
        - n0: User defined number of nodes on the transfer function [-]
    Output:
        - id: A vector of length "n" of node positions on the transfer function [-]
    =#

    n_tmp = n0 - 1
    id = []

    while length(id) != n0
        id = unique(trunc.(Int,exp10.(range(0, stop=log10(nt), length=n_tmp))))
        n_tmp += 1
    end
    return id
end

function set_weights(g_0::Vector{Float64}, f::Vector{Float64}, temperature::Vector{Float64})
    #=

    =#
    # Define proportion for terms in the objective function
    w_0 = [0.7, 0.15, 0.15]

    # Compute initial transfer functionn derivatives
    dg_0 = diff([0; g_0])
    ddg_0 = diff(diff([0; g_0]))

    # Compute each terms of the objective function
    nt = length(g_0)
    e = zeros(3)
    e[1] = sqrt((sum((convolution(f, g_0) - temperature) .^ 2)) / nt)
    e[2] = sqrt((sum((dg_0) .^ 2)) / nt)
    e[3] = sqrt((sum((ddg_0) .^ 2)) / nt)

    # Define the objective functions weights
    prop = w_0 / sum(e)
    w = w_0 / prop

    # Set array weight
    w_array = [3 * ones(trunc(Int, nt * 0.05)); ones(trunc(Int, nt * 0.95))]

    return w, w_array
end

function const_derivative(t::Vector{Float64}, id::Vector{Int}, cnst::Int)
    #=

    =#

    # Input parameters
    n = length(id)
    e = ones(Float64, n)
    h = diff([0; id])

    # First constraint: positive first derivative
    a1 = spdiagm(0 => -e, 1 => e[1:n-1])
    a1 = -a1[1:end-1, :] ./ h[1:end-1]
    b1 = zeros(Float64, n - 1)

    # Second constraint: negative second derivative on a SpecialFunctions
    # Find time at around 3 hours of test
    if maximum(t) >= 3600 * 3
        id_nodes = argmin(abs.(3600 * 3 .- t[id]))
    else
        id_nodes = 0
    end

    @show id_nodes

    # Construct the second derivative a2 and b2
    c = zero(n-2-id_nodes)
    for kk in 1:n-2-id_nodes
        c[kk] = (id[id_nodes+kk+1] - id[id_nodes+kk])/(id[id_nodes+kk] - id[id_nodes+kk-1])
        # a2[kk, id_nodes+kk-1:id_nodes+kk+1] .= [c, -(c+1), 1]
    end
    # TODO: Try to do directly a sparse matrix like with a1
    a2 = spdiagm(0 => c, 1 => -(c + 1), 2 => ones(Float64, n - 2 - id_nodes))
    b2 = zeros(Float64, n - 2 - id_nodes)

    # Setting output constraints depending on the optional input `c`
    if cnst == 2
        A = [a1; a2]
        B = [b1, b2]
    elseif cnst == 1
        A = a1
        B = b1
    elseif cnst == 0
        A = []
        B = []
    else
        error("Optionnal constraint input must be 0, 1 or 2.")
    end

    # Verification
    if any(isnan, A) || any(isnan, B)
        error("Presence of NaN in the output.")
    end

    return A, B
end

function obj_fun(g_opt::Vector{Float64}, id::Vector{Int}, idall::Vector{Int},
    f::Vector{Float64}, T::Vector{Float64}, w, w_array)
    #=

    =#

    # Interpolate the transfer function
    interp = Interpolator(id, g_opt)
    g = interp(idall)

    # Compute initial transfer functionn derivatives
    dg = diff([0; g])
    ddg = diff(diff([0; g]))

    # Compute each terms of the objective function
    nt = length(g)
    e = zeros(3)
    e[1] = w[1].*sqrt((sum((w_array.*(convolution(f, g) - T)) .^ 2)) / nt)
    e[2] = w[2].*sqrt((sum((dg) .^ 2)) / nt)
    e[3] = w[3].*sqrt((sum((ddg) .^ 2)) / nt)

    return sum(e)
end

function show_fig(t::Vector{Float64}, f::Vector{Float64}, g::Vector{Float64},
    temperature::Vector{Float64})
    #= Function that prints results of a deconvolution process (either the initial results 
    or the optimized one).
    =#
    # Define the first plot with ĝ and ĝ'.
    p1 = plot(
        t / 3600 / 24,
        g,
        xaxis="Time (d)",
        yaxis="ĝ (-)",
        xscale=:log10,
        linewidth=1.5,
        linestyle=:solid,
        linecolor=:blue,
        label="ĝ",
        legend=:topleft,
    )
    p1 = plot!(
        twinx(),
        t / 3600 / 24,
        diff([0; g]),
        yaxis="ĝ' (-)",
        xscale=:log10,
        yscale=:log10,
        linewidth=1.5,
        linestyle=:dash,
        linecolor=:black,
        label="ĝ'",
        legend=:left,
    )
    p1 = plot!(
        framestyle=:box,
        grid=false,
        xlabelfontsize=8,
        ylabelfontsize=8,
        xtickfontsize=8,
        ytickfontsize=8,
        legendfontsize=8,
    )

    # Define the second plot with the experimental and reconstructed temperature signals
    p2 = plot(
        t / 3600 / 24,
        temperature,
        linewidth=1.5,
        linestyle=:solid,
        linecolor=:black,
        label="Reference",
    )
    p2 = plot!(
        t / 3600 / 24,
        convolution(f, g),
        linewidth=1.5,
        linestyle=:dash,
        linecolor=:cyan,
        label="Convolved",
    )
    p2 = plot!(
        framestyle=:box,
        grid=false,
        xlabel="Time (d)",
        ylabel="Temperature (°C)",
        xlabelfontsize=8,
        ylabelfontsize=8,
        xtickfontsize=8,
        ytickfontsize=8,
        legendfontsize=8,
    )

    # Join both plots and print them in a layout with adequate sizes
    display(plot(p1, p2, layout=(2, 1), size=(480, 340))) # ~[17cm,12cm]
    #savefig("ExampleDeconv_fig.pdf")
end
