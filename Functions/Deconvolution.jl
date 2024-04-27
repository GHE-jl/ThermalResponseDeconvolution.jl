using FFTW, Optim, SpecialFunctions, Plots
include("../Functions/Convolution.jl")


function deconvolution(t::Vector, f::Vector, temperature::Vector, n::Int=35, c::Int=2, tol::Float64=1e-6, show_ini::Bool=false, show_opt::Bool=false)
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
          - tol: Optimality tolerance for solver (default: 1e-6)
          - StepTol: Step tolerance for solver (default: 1e-15)
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
    option_validation(n, c, tol, show_ini, show_opt)

    # 1. Initial parameters
    nt = length(temperature)
    idall = collect(1:nt)

    # 2. Compute initial solution
    g_0 = deconv_ini(idall, f, temperature)
    show_ini ? show_fig(t, f, g_0, temperature) : nothing

    # 3. Define nodes positions
    id = set_nodes(nt, n)

    # 4. Set weights for the multi-objective function
    w = set_weights(g_0, f, temperature, n)

    # 5. Set linear inequality constraints and bounds
    lb = zeros(n, 1)

    if c == 2
        a = empty
        b = empty
    elseif c == 1
        #[a, b] = ConstDeriv(t, id, 1)
    elseif c == 2
        #[a, b] = ConstDeriv(t, id, 1)
    end

    # Show results

    # Optimization part
    return g_0
end

function data_validation(t, f, temperature)
    #=

    =#

    # Check if vectors are all the same length
    if length(t) != length(f) || length(t) != length(temperature)
        error("inputs are not of the same lengths.")

        # Check for imaginary, NaN or infinite value
    elseif any(!isreal, t) || any(!isreal, f) || any(!isreal, temperature)
        error("inputs contain imaginary values.")

        # Check for NaN values
    elseif any(isnan, t) || any(isnan, f) || any(isnan, temperature)
        error("inputs contains NaN values.")
        # Check for NaN values
    elseif any(isinf, t) || any(isinf, f) || any(isinf, temperature)
        error("inputs contains infinite values.")
    end
end

function option_validation(n, c, tol, show_ini, show_opt)
    #=

    =#
    if n < 0
        error("number of nodes must be a positive integer.")
    elseif n < 15
        @warn println("the number of nodes should be at least larger than 20, ideally 30 to 60.")
    elseif n > 100
        @warn println("the number of nodes is large (>100) and should be arount 60 at the maximum.")
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

function set_nodes(nt, n0)
    #=
    Function that sets the position of the nodes on the transfer function based on
    the number of nodes asked by the user.
    Inputs:
        - nt: Total number of data in the input vectors [-]
        - n0: User defined number of nodes on the transfer function [-]
    Output:
        - id: A vector of length "n" of node positions on the transfer function [-]
    =#

    n_tmp = n0-1
    id = []

    while length(id) != n0
        id = unique(round.(exp10.(range(0, stop=log10(nt), length=n_tmp))))
        n_tmp += 1
    end
    return id
end

function deconv_ini(idall, f, temperature)
    #= 

    =#

    # Define objective function (RMSE) between model and experimental values.
    # obj_fun(x) = @. sqrt(sum((convolution(f, x[1]*expinti(x[2]/idall)) - temperature)^2) / length(t))
    objfun(x) = sqrt(sum((convolution(f, x[1].*expinti.(x[2]./idall)).-(T_in.-T_in[1])).^2)/length(idall))

    # Define the optimization on 2 variables
    x0 = [1.0, 1.0]
    x_opt = optimize(obj_fun, x0, GradientDescent(), Optim.Options(show_trace=true))

    return x_opt[1]*expint(x_opt[2]/idall)
end

function show_fig(t, f, g, temperature)
    #= Function that prints results of a deconvolution process (either the initial results 
    or the optimized one).
    =#
    # Define the first plot with ĝ and ĝ'.
    p1 = plot(t/3600/24, g, xaxis="Time (d)", yaxis="ĝ (-)", xscale=:log10,
        linewidth=1.5, linestyle=:solid, linecolor=:blue, label="ĝ",legend=:topleft)
    p1 = plot!(twinx(), t[2:end]/3600/24, diff(g),
        yaxis="ĝ' (-)", xscale=:log10, yscale=:log10,
        linewidth=1.5, linestyle=:dash, linecolor=:black, label="ĝ'",legend=:left)
    p1 = plot!(framestyle=:box, grid=false,
        xlabelfontsize=8, ylabelfontsize=8, xtickfontsize=8, ytickfontsize=8,
        legendfontsize=8)

    # Define the second plot with the experimental and reconstructed temperature signals
    p2 = plot(t/3600/24, temperature, linewidth=1.5, linestyle=:solid, linecolor=:black,
        label="Reference")
    p2 = plot!(t/3600/24, convolution(f, g), linewidth=1.5, linestyle=:dash,
        linecolor=:cyan, label="Convolved")
    p2 = plot!(framestyle=:box, grid=false, xlabel="Time (d)", ylabel="Temperature (°C)",
        xlabelfontsize=8, ylabelfontsize=8, xtickfontsize=8, ytickfontsize=8,
        legendfontsize=8)

    # Join both plots and print them in a layout with adequate sizes
    display(plot(p1, p2, layout=(2, 1), size=(480, 340))) # ~[17cm,12cm]
    #savefig("ExampleDeconv_fig.pdf")
end
