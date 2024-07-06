"""
GHEDeconvolutions.jl is a package used to apply a deconvolution algorithm to the data
obtained from ground heat exchanger (GHE). The goal to estimate the GHE transfer function
of the system (i.e., GHE) only from experimental data.
"""

module GHEDeconvolutions

# Other files to include in the package
include("Convolution.jl")
include("Deconvolution.jl")
include("utils.jl")

# Basic export
export deconvolution, convolution

# Optionnal export
export data_validation,
       option_validation,
       deconv_ini,
       set_nodes,
       set_weights,
       const_derivative,
       obj_fun,
       show_fig
end
