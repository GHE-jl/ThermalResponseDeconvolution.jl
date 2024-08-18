"""
Script that includes some general functions used throughout the pagkage.
"""

function rms(x::Vector{T} where T<:Real)
    rmse = sqrt(sum(y^2 for y in x) / length(x))
    return rmse
end
