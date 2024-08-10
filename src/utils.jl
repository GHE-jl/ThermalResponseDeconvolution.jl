"""
Script that includes some general functions used throughout the pagkage.
"""

function rms(x::AbstractVector{T} where T<:Number)
    rmse = sqrt(sum(y^2 for y in x) / length(x))
    return rmse
end
