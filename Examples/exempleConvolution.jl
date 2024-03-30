"""
Exemple of convolution on thermal response test data.
In Julia, this script can be called from the project folder.
"""

include("../Functions/Convolution.jl")

using CSV, Plots

# Import synthetic data
DATA = CSV.read("Data/Data.csv", DataFrame)
G = CSV.read("Data/gRef.csv", DataFrame)

t = DATA[:, 1]
T_in = DATA[:, 2]
T_out = DATA[:, 3]
g = G[:, 2]

# Convolution algorithm
f = diff([0; T_in - T_out])

Time = @time T_tmp = Convolution(f, g) #todo: 
T_conv = T_tmp.+T_out[1]

# Print results
println("Temperature RMSE: $round(sqrt(mean((T_conv.-T_out).^2)), 2) [degC]")
# println("Convolution time: $(Time) [s]")

# Plot results
plot(t/3600/24, T_out, linewidth=1.5, linestyle=:solid, linecolor=:black, label="Reference")
plot!(t/3600/24, T_conv, linewidth=1.5, linestyle=:dash, linecolor=:cyan, label="Convolved")
plot!(framestyle=:box, grid=false, xlabel="Time (d)", ylabel="Temperature (°C)",
xlabelfontsize=11, ylabelfontsize=11, xtickfontsize=11, ytickfontsize=11,
legendfontsize=11)
