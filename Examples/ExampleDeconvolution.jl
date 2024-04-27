"""
Exemple of convolution on thermal response test data.
In Julia, this script can be called from the project folder.
"""

include("../Functions/Deconvolution.jl")

using CSV, DataFrames, Plots

# Import synthetic data
DATA = CSV.read("Data/Data.csv", DataFrame)
G = CSV.read("Data/gRef.csv", DataFrame)

t = DATA[:, 1];
temperature_in = DATA[:, 2]
temperature_out = DATA[:, 3]
temperature_exp = temperature_out .- temperature_out[1]
g_ref = G[:, 2]

# Deconvolution algorithm
f = diff([0; temperature_in - temperature_out])
g = deconvolution(t, f, temperature_exp, 35, 2, 1e-6, true)

println("Done")