module ThermalResponseDeconvolutions

# Files to include in the package (order matters!)
include("structures.jl")
include("convolution.jl")
include("utils.jl")
include("trtdata.jl")
include("deconvolution.jl")

# Main API
export TRTData_import, convolution, deconvolution

# Data structures
export TRTData, f_FFT, DeconvParams, Deconv0Params

# Utilities
export rms, show_fig, convolution_g, set_nodes

end
