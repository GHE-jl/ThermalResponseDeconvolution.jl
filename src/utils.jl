"""
Script that includes some general functions used throughout the pagkage.
"""

using Plots

function rms(x::AbstractVector{<:Real})
    return sqrt(sum(y^2 for y in x) / length(x))
end

"""
    set_nodes(nt, n₀)

Function that sets a logarithmic progression of node positions on a transfer function.
# Arguments
    - `nt`: Total number of data in the input vectors [-]
    - `n₀`: User defined number of nodes on the transfer function [-]
# Output
    - `id`: A vector of length "n₀" of node positions on the transfer function [-]
"""
function set_nodes(nt::Real, n₀::Integer)
    # Basic inputs
    n_tmp = n₀ - 1
    id = Vector{Integer}(undef, n_tmp)
    # Fill the vector with node positions
    while length(id) < n₀
        empty!(id)
        for x in range(0, stop=log10(nt), length=n_tmp)
            push!(id, round(Int, exp10(x)))
        end
        unique!(id)
        n_tmp += 1
    end
    return id
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
        t[2:end] / 3600 / 24,
        diff(g);
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