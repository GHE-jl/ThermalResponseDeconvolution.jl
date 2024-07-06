using SpecialFunctions, Optim
include("../src/Convolution.jl")

using CSV, DataFrames, Plots

# Import synthetic data
DATA = CSV.read("Data/Data.csv", DataFrame)
G = CSV.read("Data/gRef.csv", DataFrame)

t = DATA[:, 1];
T_in = DATA[:, 2]
T_out = DATA[:, 3]
g = G[:, 2]

# Convolution algorithm
f = diff([0; T_in - T_out])

# test
idall = collect(1:length(t))

function obj_fun(x)
    return sqrt(
        sum((convolution(f, x[1] .* expint.(x[2] ./ t)) .- (T_out .- T_out[1])) .^ 2) /
        length(t),
    )
end

x0 = [1.0, 1.0]

opt = optimize(obj_fun, x0, NelderMead(), Optim.Options(; show_trace = true))
x_opt = opt.minimizer
g0 = x_opt[1] .* expint.(x_opt[2] ./ t)

rmse = sqrt(sum((convolution(f, g0) .- (T_out .- T_out[1])) .^ 2) / length(t))

println("RMSE: ", round(rmse; digits = 3), " °C")

function show_fig(t, f, g, temperature)
    #= Function that prints results of a deconvolution process (either the initial results 
    or the optimized one).
    =#
    # Define the first plot with ĝ and ĝ'.
    p1 = plot(
        t / 3600 / 24,
        g;
        xaxis = "Time (d)",
        yaxis = "ĝ (-)",
        xscale = :log10,
        linewidth = 1.5,
        linestyle = :solid,
        linecolor = :blue,
        label = "ĝ",
        legend = :topleft
    )
    p1 = plot!(
        twinx(),
        t[2:end] / 3600 / 24,
        diff(g);
        yaxis = "ĝ' (-)",
        xscale = :log10,
        yscale = :log10,
        linewidth = 1.5,
        linestyle = :dash,
        linecolor = :black,
        label = "ĝ'",
        legend = :left
    )
    p1 = plot!(;
        framestyle = :box,
        grid = false,
        xlabelfontsize = 8,
        ylabelfontsize = 8,
        xtickfontsize = 8,
        ytickfontsize = 8,
        legendfontsize = 8
    )

    # Define the second plot with the experimental and reconstructed temperature signals
    p2 = plot(
        t / 3600 / 24,
        temperature;
        linewidth = 1.5,
        linestyle = :solid,
        linecolor = :black,
        label = "Reference"
    )
    p2 = plot!(
        t / 3600 / 24,
        convolution(f, g);
        linewidth = 1.5,
        linestyle = :dash,
        linecolor = :cyan,
        label = "Convolved"
    )
    p2 = plot!(;
        framestyle = :box,
        grid = false,
        xlabel = "Time (d)",
        ylabel = "Temperature (°C)",
        xlabelfontsize = 8,
        ylabelfontsize = 8,
        xtickfontsize = 8,
        ytickfontsize = 8,
        legendfontsize = 8
    )
    p2 = annotate!((0.15, 0.9), text("RMSE: " * string(round(rmse; digits = 3)) * " °C", 8))

    # Join both plots and print them in a layout with adequate sizes
    return display(plot(p1, p2; layout = (2, 1), size = (480, 340))) # ~[17cm,12cm]
    #savefig("ExampleDeconv_fig.pdf")
end

show_fig(t, f, g0, T_out .- T_out[1])
