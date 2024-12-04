"""
GHEDeconvolutions.jl is a package used to apply a deconvolution algorithm to the data
obtained from ground heat exchanger (GHE). The goal to estimate the GHE transfer function
of the system (i.e., GHE) only from experimental data.
"""

module GHEDeconvolutions

# Other files to include in the package
include("TRTData.jl")
include("Convolution.jl")
include("Deconvolution.jl")
include("utils.jl")

# Basic export
export TRTData_import, convolution, deconvolution

# Utils export
export rms, show_fig

# Structures of data
export TRTData, f_FFT

# Optionnal export
export data_validation,
       option_validation,
       define_f,
       set_nodes,
       convolution_g,
       deconv_ini!,
       set_weights,
       const_derivative,
       obj_fun
end
