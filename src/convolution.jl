using DSP: conv

"""
    convolution(f, g)

Non-circular convolution of `f` and `g` (zero-padded internally to avoid circular aliasing),
truncated to `length(f)`. Thin wrapper around `DSP.conv`, matching `GroundHeatExchanger.jl`'s
`convolutionf`.
# Arguments
    - `f`: The first signal to convolve.
    - `g`: The second signal to convolve.
# Returns
    - The convolution of `f` and `g`, truncated to `length(f)`.
# Reference
    - Marcotte, D., & Pasquier, P. (2008). Fast fluid and ground temperature computation for
        geothermal ground-loop heat exchanger systems. Geothermics, 37(6), 651–665.
        https://doi.org/10.1016/j.geothermics.2008.08.003
    - Pasquier, P., & Marcotte, D. (2013). Efficient computation of heat flux signals to ensure the
        reproduction of prescribed temperatures at several interacting heat sources. Applied Thermal
        Engineering, 59(1–2), 515–526. https://doi.org/10.1016/j.applthermaleng.2013.06.018
"""
function convolution(f::AbstractVector{<:Real}, g::AbstractVector{<:Real})
    n = length(f)
    n == length(g) || throw(ArgumentError("f and g must have the same length"))
    return @view conv(f, g)[1:n]
end
