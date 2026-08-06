# Modeling theory

## Problem statement

The forward model relating the incremental perturbation function `f` to the temperature
variation `Texp` is a discrete, non-circular convolution with the thermal response function `g`:

```math
T_{exp} \approx f * g.
```

Deconvolution recovers `g` from measured `f` and `Texp` — the inverse of
[`convolution`](@ref). This is ill-posed as a direct inverse (small measurement noise in `Texp`
is amplified into large, unphysical oscillations in a naively recovered `g`), so
[`deconvolution`](@ref) instead poses it as a **constrained, regularized optimization problem**,
following Dion et al. (2022).

## Node parameterization

Rather than optimizing one value of `g` per data point (`nData`, typically thousands of samples),
the optimization solves for `g` at only `n` log-spaced node indices ([`set_nodes`](@ref) on
`1:nData`), then interpolates ([`PCHIPInterpolation.jl`](https://github.com/gerlero/PCHIPInterpolation.jl),
shape-preserving) onto the full index axis to obtain `ĝ`. Log-spacing matches the thermal
response function's own behavior: it varies quickly at early times and slowly at late times, so
fewer nodes are needed there.

## Initial guess

The optimization is initialized with a 2-parameter exponential-integral fit
(Dion et al. 2022, Eqs. 20–21):

```math
g_0(idx) = x_1 \, E_1(x_2 / idx),
```

with `E₁` the exponential integral, fit to `Texp` by a derivative-free Nelder-Mead search. No
physical ground model is required for this step.

## Multi-objective function

The optimization minimizes a weighted combination of three terms (following Dion et al. 2023):

```math
E = w_1 \, \mathrm{rms}\big(w_a \cdot (T_{conv} - T_{exp})\big) + w_2 \, \mathrm{rms}(g') + w_3 \, \mathrm{rms}(g''),
```

where `T_conv` is the temperature reconstructed from the current iterate of `g` by convolution,
`g'`/`g''` are its first and second discrete derivatives (regularization terms, penalizing
non-physical oscillation), and `wₐ` is an array weight emphasizing the leading 5% of samples,
where the response function changes fastest. The scalar weights `w₁, w₂, w₃` are set from the
initial guess so that the three terms contribute roughly 70%/15%/15% to `E`, following the
reference implementation.

## Constraints

The `c` keyword of [`deconvolution`](@ref) selects which physical constraints are imposed on `g`
as linear inequalities:

- `c = 0`: none.
- `c = 1`: a positive first derivative (`g` strictly growing).
- `c = 2`: adds a negative second derivative past roughly 3 h from the start of the record, where
  the thermal response function is expected to be concave.

## Solver

The constrained problem is solved with
[Optimization.jl](https://github.com/SciML/Optimization.jl), using NLopt's SLSQP algorithm
(sequential quadratic programming, via
[OptimizationNLopt.jl](https://github.com/SciML/Optimization.jl)) and finite-difference gradients
([FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl), since the FFT-based convolution in
the objective does not support automatic differentiation). As an active-set method, SLSQP handles
a starting point on or near a constraint boundary natively, so no interior-point-specific
safeguard on the initial guess is needed.
