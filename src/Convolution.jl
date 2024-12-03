using FFTW

function convolution(f::Vector{T}, g::Vector{T}) where {T<:Real}
    n = length(f)
    #total_length = nextpow(2, n)
    total_length = 2*n-1

    # Preallocate arrays
    f_pad = zeros(T, total_length)
    g_pad = zeros(T, total_length)

    # Fill preallocated arrays
    copyto!(f_pad, f)
    copyto!(g_pad, g)

    # Precompute FFT plans
    fft_plan = plan_rfft(f_pad)

    # Compute FFTs
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Element-wise multiplication (in-place)
    F_f .*= F_g

    # Inverse FFT to get convolution resuls
    y = irfft(F_f, total_length)
    return @view y[1:n]
end