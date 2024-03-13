using FFTW

function Convolution(f::Vector{T}, g::Vector{T}) where T
    # Perform a stationary convolution with fft and ifft with 2xn -1zero padding to avoid circular convolution.

    # Initialization
    n = length(f)
    Pad = 2*n-1

    # Perform non-circular convolution in the spectral domain using padded signals
    y = real(ifft(fft([f; zeros(T, Pad-n)]).*fft([g; zeros(T, Pad-n)])))

    # Verification of the presence of complex or NaN values in the output
    if any(isnan.(y)) || any(imag.(y) .!= 0)
        throw(ArgumentError("Presence of NaN or Complex"))
    end

    return y[1:n]
end