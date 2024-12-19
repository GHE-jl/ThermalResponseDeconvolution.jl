"""
File that defines structures used in GHEDeconvolutions.
"""

struct TRTData
    t::Vector{Float64}      # Time array [s]
    Q::Vector{Float64}      # Heating power [W]
    Tin::Vector{Float64}    # Inlet fluid temperature [degC]
    Tout::Vector{Float64}   # Outlet fluid temperature [degC]
    Texp::Vector{Float64}   # Experimental temperature to deconvolve [degC]
end

struct f_FFT
    f::Vector{Float64}      # Incremental load function
    f_pad::Vector{Float64}  # Padded incremental load function
    fft_plan                # Plan of FFT for optimized computation afterwards
    F_f                     # Frequency domain of the incremental load function
end

struct deconv_optim₀
    data::TRTData
    f_fft::f_FFT
end

mutable struct deconv_obj
    g::Vector{Float64}
    dg::Vector{Float64}
    ddg::Vector{Float64}
    e::Vector{Float64}
end

struct deconv_optim
    data::TRTData
    id::Vector{Integer}
    f_fft::f_FFT
    w::Vector{Float64}
    wₐ::Vector{Float64}
    cnst::Integer
    objfun::deconv_obj
end
