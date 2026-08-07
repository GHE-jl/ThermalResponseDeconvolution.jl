# API reference

This page is an alphabetical index of every exported symbol, followed by the full docstrings.

```@index
Modules = [ThermalResponseDeconvolution]
```

## By topic

### Deconvolution

- [`deconvolution`](@ref) — recover a thermal response function from paired temperature/heat-load data
- [`set_nodes`](@ref) — log-spaced node indices used by `deconvolution`
- [`rms`](@ref) — root-mean-square, e.g. for reporting a fit's residual

### Convolution

- [`convolution`](@ref) — non-circular convolution, e.g. to reconstruct a temperature response

!!! note "Internal functions"
    Everything else in the package (names starting with `_`) is an unexported implementation
    detail of `deconvolution` and is not part of the public API — see [Modeling theory](@ref) for
    what they do.

```@autodocs
Modules = [ThermalResponseDeconvolution]
```
