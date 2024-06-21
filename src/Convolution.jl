using FFTW

function convolution(f::Vector{Float64}, g::Vector{Float64})
    # Perform a stationary convolution with fft and ifft with 2xn -1zero padding to avoid circular convolution.

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
