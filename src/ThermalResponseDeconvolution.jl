module ThermalResponseDeconvolution

# Files to include in the package
include("convolution.jl")
include("deconvolution.jl")

# Convolution
export convolution

# Deconvolution
export deconvolution, set_nodes, rms

end
