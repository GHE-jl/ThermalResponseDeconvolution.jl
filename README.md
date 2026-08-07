# ThermalResponseDeconvolution.jl

[![CI](https://github.com/GHE-jl/ThermalResponseDeconvolution.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GHE-jl/ThermalResponseDeconvolution.jl/actions/workflows/CI.yml)
[![Docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GHE-jl.github.io/ThermalResponseDeconvolution.jl/dev)
[![Docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GHE-jl.github.io/ThermalResponseDeconvolution.jl/stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A Julia package to recover a borehole outlet thermal response function of a ground heat exchanger by
deconvolution of paired fluid-temperature and heat-load data. The data can come from
either a thermal response test (TRT) or ground source heat pump (GSHP) operating data (the
method only assumes a constant time step, not a controlled experiment).

Deconvolution is posed as a constrained, multi-objective optimization problem (weighted temperature
misfit plus first- and second-derivative regularization on the response function), solved with
[Optimization.jl](https://github.com/SciML/Optimization.jl) using the NLopt SLSQP backend.

## Quick start

```julia
using ThermalResponseDeconvolution

# f: incremental perturbation function [degC], Texp: measured temperature variation [degC]
f = diff([0.0; Tin .- Tout])
Texp = Tout .- Tout[1]

ĝ, gOpt = deconvolution(t, f, Texp; n=50, c=2)  # ĝ: thermal response function, gOpt: node values
T̂ = collect(convolution(f, ĝ))                  # reconstructed temperature, for validation
rmse = rms(T̂ .- Texp)                           # Temperature root mean square error
```

## Convolution

| Function | Purpose |
|---|---|
| `convolution(f, g)` | Non-circular convolution of `f` and `g`, truncated to `length(f)`. Used to reconstruct the temperature response from a (measured or deconvolved) thermal response function. Same implementation as [`convolutionf`](https://github.com/GHE-jl/GroundHeatExchanger.jl/blob/main/src/temporal_superposition.jl) in GroundHeatExchanger.jl |

## Deconvolution

| Function | Purpose |
|---|---|
| `deconvolution(t, f, Texp; n=35, c=2)` | Recover a thermal response function `ĝ` (and the optimized node values `gOpt`, before interpolation) from a time array, an incremental perturbation function, and a measured temperature variation. |
| `set_nodes(nData, n0)` | Log-spaced node indices on `1:nData`, as used internally by `deconvolution` — call it with the same `nData`/`n` to recover where `gOpt` sits along `ĝ`. |
| `rms(x)` | Root-mean-square of a vector, e.g. for reporting a fit's temperature or response-function residual. |

The `c` keyword selects the inequality constraints on the response function: `0` (none), `1`
(positive first derivative), or `2` (adds a negative second derivative past roughly 3 h from the
start of the record).

## Scripts

Run from the package root with `julia --project=script/ script/<name>.jl`. First-time setup:
```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

| Script | What it shows |
|---|---|
| `script_deconvolution.jl` | End-to-end deconvolution on TRT data, with the reconstructed temperature and response-function derivative, validated against a reference solution. |

## Installation

The package is not yet registered. Install directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/ThermalResponseDeconvolution.jl")
```

Or in the Julia REPL package manager (`]`):

```
pkg> add https://github.com/GHE-jl/ThermalResponseDeconvolution.jl
```

## Dependencies

### Library

| Package | Used for |
|---|---|
| [DSP.jl](https://github.com/JuliaDSP/DSP.jl) | FFT-based convolution |
| [Optimization.jl](https://github.com/SciML/Optimization.jl) / [OptimizationNLopt.jl](https://github.com/SciML/Optimization.jl) | Constrained multi-objective optimization with the NLopt SLSQP backend |
| [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl) | Finite-difference gradients for the optimization |
| [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl) | Shape-preserving interpolation of the node values onto the full index axis |
| [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) | Exponential-integral fit used for the initial guess |

### Scripts only

| Package | Used in |
|---|---|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | Plotting in `script_deconvolution.jl` |
| [DelimitedFiles.jl](https://github.com/JuliaData/DelimitedFiles.jl) | Reading the sample CSV data in `script_deconvolution.jl` |

## References

- Dion, G., Pasquier, P., & Marcotte, D. (2022). Deconvolution of experimental thermal response
  test data to recover short-term g-function. *Geothermics*, 100, 102302.
  https://doi.org/10.1016/j.geothermics.2021.102302
- Dion, G., Pasquier, P., Marcotte, D., & Beaudry, G. (2023). Multi-deconvolution in non-stationary
  conditions applied to experimental thermal response test analysis to obtain short-term transfer
  functions. *Science and Technology for the Built Environment*, 30(3), 1–14.
  https://doi.org/10.1080/23744731.2023.2217729
- Dion, G., Pasquier, P., & Marcotte, D. (2024). Application of deconvolution to interpretation of
  distributed thermal response test (DTRT) and to determination of thermal conductivity profiles.
  *Applied Thermal Engineering*, 236, 121680. https://doi.org/10.1016/j.applthermaleng.2023.121680
- Marcotte, D., & Pasquier, P. (2008). Fast fluid and ground temperature computation for geothermal
  ground-loop heat exchanger systems. *Geothermics*, 37(6), 651–665.
  https://doi.org/10.1016/j.geothermics.2008.08.003
