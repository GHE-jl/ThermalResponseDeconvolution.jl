"""
Script that includes some general functions used throughout the pagkage.
"""

function rms(x::Vector{Float64})
    rmse = @. sqrt(sum(y^2)/length(x))
    return rmse
end