"""
Script that includes some general functions used throughout the pagkage.
"""

using Plots

function rms(x::Vector{T} where {T<:Real})
    rmse = sqrt(sum(y^2 for y in x) / length(x))
    return rmse
end

function show_fig(t::Vector{Float64},
    f::Vector{Float64},
    g::Vector{Float64},
    temperature::AbstractVector{Float64}
)
    """ Function that prints results of a deconvolution process (either the initial results 
    or the optimized one).
    """

    # Define the first plot with ĝ and ĝ'.
    p1 = plot(
        t / 3600 / 24,
        g;
        xaxis="Time (d)",
        yaxis="ĝ (-)",
        xscale=:log10,
        linewidth=1.5,
        linestyle=:solid,
        linecolor=:blue,
        label="ĝ",
        legend=:topleft
    )
    p1 = plot!(
        twinx(),
        t / 3600 / 24,
        diff([0; g]);
        yaxis="ĝ' (-)",
        xscale=:log10,
        yscale=:log10,
        linewidth=1.5,
        linestyle=:dash,
        linecolor=:black,
        label="ĝ'",
        legend=:left
    )
    p1 = plot!(;
        framestyle=:box,
        grid=false,
        xlabelfontsize=8,
        ylabelfontsize=8,
        xtickfontsize=8,
        ytickfontsize=8,
        legendfontsize=8
    )

    # Define the second plot with the experimental and reconstructed temperature signals
    p2 = plot(
        t / 3600 / 24,
        temperature;
        linewidth=1.5,
        linestyle=:solid,
        linecolor=:black,
        label="Reference"
    )
    p2 = plot!(
        t / 3600 / 24,
        convolution(f, g);
        linewidth=1.5,
        linestyle=:dash,
        linecolor=:cyan,
        label="Convolved"
    )
    p2 = plot!(;
        framestyle=:box,
        grid=false,
        xlabel="Time (d)",
        ylabel="Temperature (°C)",
        xlabelfontsize=8,
        ylabelfontsize=8,
        xtickfontsize=8,
        ytickfontsize=8,
        legendfontsize=8
    )

    # Join both plots and print them in a layout with adequate sizes
    return display(plot(p1, p2; layout=(2, 1), size=(480, 340))) # ~[17cm,12cm]
    #savefig("ExampleDeconv_fig.pdf")
end