"""
File that defines the structures used in ThermalResponseDeconvolutions.
"""

"""
    TRTData

Immutable container for a thermal response test (TRT) dataset.
"""
struct TRTData
    t::Vector{Float64}      # Time array [s]
    Q::Vector{Float64}      # Heating power [W]
    Tin::Vector{Float64}    # Inlet fluid temperature [degC]
    Tout::Vector{Float64}   # Outlet fluid temperature [degC]
    Texp::Vector{Float64}   # Experimental temperature to deconvolve (Tout - Tout[1]) [degC]
end

"""
    f_FFT

Holds the incremental load function `f` together with its zero-padded version,
the cached real-FFT plan and the frequency-domain representation `F_f`. Caching
the plan/transform avoids recomputing them at every objective evaluation.
"""
struct f_FFT
    f::Vector{Float64}      # Incremental load function
    f_pad::Vector{Float64}  # Zero-padded incremental load function (length 2n-1)
    fft_plan                # Cached rFFT plan
    F_f                     # Frequency domain of the incremental load function
end

"""
    DeconvParams

Parameters passed to the objective and constraint functions during the main
constrained minimization. Everything here is constant w.r.t. the optimization
variables (the node values), so it is computed once before `solve`.
"""
struct DeconvParams
    Texp::Vector{Float64}       # Experimental temperature variation [degC]
    f_fft::f_FFT                # Cached FFT of the incremental load function
    idnodes::Vector{Float64}    # Node positions (indices) used for interpolation
    idall::Vector{Float64}      # Full index axis 1:nData
    w::Vector{Float64}          # 3 scalar weights of the multi-objective function
    wₐ::Vector{Float64}         # Array weight emphasizing the first samples (length nData)
    A::Matrix{Float64}          # Linear inequality matrix; constraints are A*g .<= 0
end

"""
    Deconv0Params

Parameters for the 2-parameter initial-guess fit (`deconv_ini`).
"""
struct Deconv0Params
    Texp::Vector{Float64}
    f_fft::f_FFT
    idall::Vector{Float64}
end
