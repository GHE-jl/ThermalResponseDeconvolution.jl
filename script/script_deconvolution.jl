"""
Example of deconvolution on thermal response test data.
"""

using ThermalResponseDeconvolutions
using CSV
using DataFrames

data = TRTData_import(joinpath(@__DIR__, "..", "data", "DataCL_TRT.csv"))

# Deconvolution algorithm
g, T = deconvolution(data; n=35, c=2, show_opt=true)

# Reconstruction quality
rmse = rms(T .- data.Texp)
println("Temperature RMSE: ", round(rmse; digits=3), " [degC]")