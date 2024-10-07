using FFTW

function convolution(f::Vector{T}, g::Vector{T}) where T<:Real
    # Perform a stationary convolution with fft and ifft with 2xn -1 zero padding to avoid circular convolution.

    # Initialization
    n = length(f)
    pad = n - 1
    n_full = n+pad

    # Preallocate array
    f_pad = zeros(T, n_full)
    g_pad = zeros(T, n_full)

    # Filling the preallocated array
    f_pad[1:n] .= f
    g_pad[1:n] .= g

    # Preallocate FFT results
    F_f = similar(f_pad)
    F_g = similar(g_pad)
    F_results = similar(f_pad)

    # Compute FFT
    F_f = rfft(f_pad)
    F_g = rfft(g_pad)

    # Element-wise multiplication
    F_results = F_f.*F_g

    # Inverse FFT to get convolution result
    y = irfft(F_results, n+pad)

    # Verification of the presence of complex or NaN values in the output
    if any(isnan.(y))
        throw(ArgumentError("Presence of NaN"))
    end

    return y[1:n]
end

function convolution2(f::Vector{T}, g::Vector{T}) where T<:Real
    # Perform a stationary convolution with fft and ifft with 2xn -1 zero padding to avoid circular convolution.

    # Initialization
    n = length(f)
    Pad = n - 1

    # Perform non-circular convolution in the spectral domain using padded signals
    y = real(ifft(fft([f; zeros(Pad)]) .* fft([g; zeros(Pad)])))

    # Verification of the presence of complex or NaN values in the output
    if any(isnan.(y)) || any(imag.(y) .!= 0)
        throw(ArgumentError("Presence of NaN or Complex"))
    end

    return y[1:n]
end

function convolution3(f::Vector{T}, g::Vector{T}) where T<:Real
    # Perform a stationary convolution with fft and ifft with 2xn -1 zero padding to avoid circular convolution.

    # Initialization
    n = length(f)
    pad = n - 1
    n_full = n+pad

    # Preallocate array
    f_pad = zeros(T, n_full)
    g_pad = zeros(T, n_full)

    # Filling the preallocated array
    copyto!(f_pad, f)
    copyto!(g_pad, g)

    # Precompute FFT plans
    fft_plan = plan_rfft(f_pad)
    ifft_plan = plan_irfft(f_pad, n_full)

    # Compute FFTs in-place
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Element-wise multiplication
    F_f .*= F_g

    # Inverse FFT to get convolution result
    y = ifft_plan * F_f

    # Verification of NaN values
    any(isnan, result) && throw(ArgumentError("Presence of NaN"))

    return y[1:n]
end

function convolution4(f::Vector{T}, g::Vector{T}) where T<:Real
    n = length(f)
    Pad = n - 1
    total_length = n + Pad

    # Preallocate arrays
    f_pad = zeros(T, total_length)
    g_pad = zeros(T, total_length)

    # Fill preallocated arrays
    copyto!(f_pad, f)
    copyto!(g_pad, g)

    # Precompute FFT plans
    fft_plan = plan_rfft(f_pad)
    ifft_plan = plan_irfft(fft_plan * f_pad, total_length)

    # Compute FFTs
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Element-wise multiplication (in-place)
    F_f .*= F_g

    # Inverse FFT to get convolution result
    y = ifft_plan * F_f

    # Use view to avoid copying
    result = @view y[1:n]

    # Verification of NaN values
    any(isnan, result) && throw(ArgumentError("Presence of NaN"))

    return result
end