"""
Exemple of convolution on thermal response test data.
In Julia, this script can be called from the project folder.
"""

include("../src/Deconvolution.jl")

using CSV
using DataFrames
using Plots
using Revise

function main()

    # Import synthetic data
    data = TRTData_import("data/DataCL_TRT.csv")
    gG = CSV.read("data/gNumCL.csv", DataFrame, skipto=6)
    DATA = CSV.read("data/data.csv", DataFrame)
    G = CSV.read("data/gRef.csv", DataFrame)
    g_ref = G[:, 2]

    # Deconvolution algorithm
    g = deconvolution(data, 35, 2)

    return println("Done")
end

main()
