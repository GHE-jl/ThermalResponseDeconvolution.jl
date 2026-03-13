module ThermalResponseDeconvolution

# Files to include in the package (order is important!)
include("Structures.jl")
include("TRTData.jl")
include("Convolution.jl")
include("Deconvolution.jl")
include("utils.jl")

# Basic export
export TRTData_import, convolution, deconvolution

# Structures of data
export TRTData, f_FFT, deconv_optim₀, deconv_obj, deconv_optim

# Utils export
export rms, show_fig

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
