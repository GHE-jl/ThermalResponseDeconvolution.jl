"""
Exemple of convolution on thermal response test data.
In Julia, this script can be called from the project folder.
"""

using GHEDeconvolutions
using CSV, DataFrames
using Plots
using BenchmarkTools

# Import synthetic data
data = TRTData_import("data/DataCL_TRT.csv")
gG = CSV.read("data/gNumCL.csv", DataFrame, skipto=6)
DATA = CSV.read("data/data.csv", DataFrame)
G = CSV.read("data/gRef.csv", DataFrame)

t = DATA[:, 1]
T_in = DATA[:, 2]
T_out = DATA[:, 3]
g = G[:, 2]

# Convolution algorithm
f = diff([0; T_in - T_out])

@btime T_tmp = convolution(f, g)
T_conv = T_tmp .+ T_out[1]

# Print results
rmse = rms(T_conv - T_out)
println("Temperature RMSE: ", round(rmse; digits=2), " [degC]")

# Plot results
plot(t/3600/24, T_out; linewidth=1.5, linestyle=:solid, linecolor=:black, label="Reference")
plot!(t/3600/24, T_conv; linewidth=1.5, linestyle=:dash, linecolor=:deepskyblue, label="Convolved")
plot!(; framestyle=:box, grid=false, xlabel="Time (d)", ylabel="Temperature (°C)",
    xlabelfontsize=11, ylabelfontsize=11, xtickfontsize=11, ytickfontsize=11, legendfontsize=11)
annotate!((0.1, 0.1), text("RMSE = " * string(round(rmse; digits=2)) * "(°C)", :left, 11))
