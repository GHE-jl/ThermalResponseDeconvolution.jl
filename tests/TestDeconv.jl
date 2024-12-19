"""
To debug the file, in the REPL:
using GHEDeconvolutions
includet("tests\\TestDeconv.jl")
@run main()
"""

using GHEDeconvolutions

function main()

    # Import synthetic data
    data = TRTData_import("data/DataCL_TRT.csv")
    n, c = 35, 2
    #G = CSV.read("Data/gNumCL.csv", DataFrame)

    # Deconvolution algorithm    
    g, T = deconvolution(data, n=n, c=c)

    show_fig(data.t, data.Tin-data.Tout, g, data.Texp)

    return g, T, println("Done")
end

main()
