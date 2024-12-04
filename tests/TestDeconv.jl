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
    #G = CSV.read("Data/gNumCL.csv", DataFrame)

    # Deconvolution algorithm    
    g = deconvolution(data, 35, 2)

    return println("Done")
end

main()
