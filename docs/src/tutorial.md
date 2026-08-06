# Tutorial

This tutorial walks through a complete deconvolution, from raw fluid-temperature and heat-load
data to a validated thermal response function, using the sample data shipped in `data/`. Every
function used here is documented in the [API reference](@ref).

## 1. Load the data

The only requirement on the input data is a **constant time step**. Here `t` is a time array
starting at 0 [s], `Tin`/`Tout` are the inlet/outlet fluid temperatures [degC], measured either
during a thermal response test or from normal GSHP operation:

```julia
using DelimitedFiles
using ThermalResponseDeconvolution

data, _ = readdlm("data/data.csv", ',', header=true)
t    = Float64.(data[:, 1])
Tin  = Float64.(data[:, 2])
Tout = Float64.(data[:, 3])
```

## 2. Build the deconvolution inputs

`deconvolution` takes an **incremental perturbation function** `f` (the discrete derivative of the
inlet/outlet temperature difference) and the **measured temperature variation** `Texp` (relative
to its initial value):

```julia
f    = diff([0.0; Tin .- Tout])
Texp = Tout .- Tout[1]
```

## 3. Run the deconvolution

```julia
ĝ, gOpt = deconvolution(t, f, Texp; n=50, c=2)
```

- `n` is the number of log-spaced nodes the optimization solves for (`ĝ` itself is returned at
  every index, interpolated from those nodes — see [`set_nodes`](@ref) to recover their
  positions and [Modeling theory](@ref) for why a reduced node set is used).
- `c=2` enforces both a positive first derivative and a negative second derivative past roughly
  3 h into the record (see [Modeling theory](@ref) for the full constraint set, `c=1` or `c=0`
  relax it).

## 4. Validate against the reconstructed temperature

Since `ĝ` is a thermal response function, convolving it back with `f` should reproduce the
measured temperature response — this is the standard sanity check, and does not require any
reference thermal response function to be known:

```julia
T̂ = collect(convolution(f, ĝ))

rmse_T = rms(T̂ .- Texp)
```

## 5. Locate the optimization nodes on the curve

`gOpt` holds the optimized values *before* interpolation, at the node indices returned by
`set_nodes` (called with the same `length(Texp)` and `n` used in step 3):

```julia
id = set_nodes(length(Texp), 50)

# ĝ[id] ≈ gOpt: the interpolation passes exactly through the optimized node values
```

`script/script_deconvolution.jl` runs this full example end to end, additionally plotting `ĝ`,
its derivative, and `T̂` against a reference solution.
