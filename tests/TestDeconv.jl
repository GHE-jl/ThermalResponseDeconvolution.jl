"""
To debug the file, in the REPL:
cd("C:/Users/Gabriel/OneDrive - polymtlus/Documents/Coding/Julia/GHEDeconvolutions")
using Revise
using GHEDeconvolutions
includet("tests\\TestDeconv.jl")
@run main()
"""

using CSV
using DataFrames
using Revise

using GHEDeconvolutions

#function main()

    # Import synthetic data
    DATA = CSV.read("Data/Data.csv", DataFrame)
    G = CSV.read("Data/gRef.csv", DataFrame)

    t = tuple(DATA[:, 1])
    temperature_in = tuple(DATA[:, 2])
    temperature_out = tuple(DATA[:, 3])
    temperature_0 = temperature_out[1]
    temperature_exp = temperature_out - temperature_0
    g_ref = tuple(G[:, 2])

    # Deconvolution algorithm
    f = diff([0; temperature_in - temperature_out])
    @profview g = deconvolution(t, f, temperature_exp, 35, 2, true, true)

#    return println("Done")
#end

#main()
