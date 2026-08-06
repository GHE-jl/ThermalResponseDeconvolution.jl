# Script that demonstrates the deconvolution of a thermal response function from numerical data. The
# script reads in the data and reference thermal response function, performs the deconvolution,
# calculates the RMSE for both temperature and thermal response function, and plots the results.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate()

using DelimitedFiles
using CairoMakie
using ThermalResponseDeconvolution

# Load data
data_dir = joinpath(@__DIR__, "..", "data")
data, _ = readdlm(joinpath(data_dir, "data.csv"), ',', header=true)
gref, _ = readdlm(joinpath(data_dir, "gRef.csv"), ',', header=true)

# Prepare data
t = Float64.(data[:, 1])
Tin = Float64.(data[:, 2])
Tout = Float64.(data[:, 3])
gRef = Float64.(gref[:, 2])

# Perform deconvolution
f = diff([0.0; Tin .- Tout])                # Incremental heat load function (for impulse in °C)
Texp = Tout .- Tout[1]                      # Experimental temperature response

n = 50                                       # Number of nodes used in the optimization
ĝ, gOpt = deconvolution(t, f, Texp; n=n, c=2) # Deconvolved thermal response function
T̂ = collect(convolution(f, ĝ))              # Reconstructed temperature response

# Calculate RMSE for temperature and thermal response function
rmse_T = rms(T̂ .- Texp)
rmse_g = rms(ĝ .- gRef)
println("Temperature RMSE: $(round(rmse_T, digits=3)) [degC]")
println("Thermal Response Function RMSE: $(round(rmse_g, digits=3)) [-]")

# Calculate the first derivative of the reference and deconvolved thermal response functions
id = set_nodes(length(Texp), n)             # Node indices used in the optimization
dgRef = diff([0.0; gRef])
dĝ = diff([0.0; ĝ])

# Plot results
fig = Figure(size=(1000, 700))

ax1 = Axis(fig[1, 1], xlabel="Time (d)", ylabel="ĝ (-)", xscale=log10)
lines!(ax1, t ./ 86400, gRef, color=:black, label="Reference")
lines!(ax1, t ./ 86400, ĝ, color=:dodgerblue, linestyle=:dash, label="Deconvolved")
scatter!(ax1, t[id] ./ 86400, gOpt, color=:dodgerblue, markersize=8, label="Nodes")
axislegend(ax1, position=:lt)
text!(ax1, 0.02, 0.5, text="RMSE: $(round(rmse_g, digits=3)) (-)", space=:relative,
    align=(:left, :center))

ax2 = Axis(fig[1, 2], xlabel="Time (d)", ylabel="dĝ (-)", xscale=log10, yscale=log10)
lines!(ax2, t ./ 86400, dgRef, color=:black, label="Reference")
lines!(ax2, t ./ 86400, dĝ, color=:dodgerblue, linestyle=:dash, label="Deconvolved")
axislegend(ax2, position=:lt)

ax3 = Axis(fig[2, 1:2], xlabel="Time (d)", ylabel="Tout - T0 (°C)")
lines!(ax3, t ./ 86400, Texp, color=:black, label="Experimental")
lines!(ax3, t ./ 86400, T̂, color=:crimson, linestyle=:dash, label="Convolved")
axislegend(ax3, position=:lt)
text!(ax3, 0.98, 0.9, text="RMSE: $(round(rmse_T, digits=3)) (°C)", space=:relative,
    align=(:right, :top))

display(fig)
