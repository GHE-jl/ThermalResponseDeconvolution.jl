# ThermalResponseDeconvolution.jl

*Recover a borehole outlet thermal response function of a ground heat exchanger by deconvolution.*

`ThermalResponseDeconvolution.jl` estimates a borehole's thermal response function
directly from paired fluid-temperature and heat-load data, by solving a constrained multi-objective
optimization problem, rather than fitting a physical ground model (as
[ThermalResponseTest.jl](https://github.com/GHE-jl/ThermalResponseTest.jl) does). Because it makes no ground-model assumption, it applies equally to data
from a dedicated thermal response test (TRT) or from ordinary ground source heat pump (GSHP)
operating data logged over time (anywhere paired temperature and heat-load series with a constant
time step are available).

## Installation

The package is not yet registered. Install it directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/ThermalResponseDeconvolution.jl")
```

or, in the Pkg REPL mode (press `]`):

```
pkg> add https://github.com/GHE-jl/ThermalResponseDeconvolution.jl
```

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

## Manual outline

- **[Tutorial](@ref)** — a worked example on sample thermal response test data, step by step.
- **[Modeling theory](@ref)** — the multi-objective formulation, constraints, and initial guess.
- **[API reference](@ref)** — the complete docstring reference for every exported function.
- **[References](@ref)** — the bibliography underpinning the implementation.

## Conventions used throughout

| Symbol | Meaning | Unit |
|---|---|---|
| `t` | Time array, starting at 0, with constant time step | s |
| `f` | Incremental perturbation function, `diff([0; Tin .- Tout])` | degC |
| `Texp` | Measured temperature variation, `Tout .- Tout[1]` | degC |
| `ĝ` | Estimated thermal response function, interpolated at every index | - |
| `gOpt` | Optimized thermal response function values at the node indices, before interpolation | - |
| `T̂` | Reconstructed temperature variation, `T̂ = convolution(f, ĝ)` | degC |
| `n` | Number of nodes in the optimization problem | - |
| `c` | Choice of inequality constraints on `ĝ` (0, 1 or 2) | - |

## Ecosystem

`ThermalResponseDeconvolution.jl` is a model-free, standalone alternative to the model-based
inversion in [ThermalResponseTest.jl](https://github.com/GHE-jl/ThermalResponseTest.jl): the two
packages solve related but distinct problems and can be used independently or side by side:
deconvolution recovers the borehole thermal response function itself, while model inversion
recovers physical ground properties (e.g. thermal conductivity) from an assumed ground model.
Neither package depends on the other.
