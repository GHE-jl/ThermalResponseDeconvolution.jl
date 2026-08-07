# Modeling theory

## Problem statement

The forward model relating the incremental perturbation function `f` (also known as an incremental heat load or temperature function) to the temperature
variation `Texp` is a discrete, non-circular convolution with the thermal response function `g`:

```math
T_{exp} = (f * g)(t).
```

Deconvolution recovers `g` from measured `f` and `Texp` (the inverse of
[`convolution`](@ref)). This is ill-posed as a direct inverse (small measurement noise in `Texp`
is amplified into large, unphysical oscillations in a naively recovered `g`), so
[`deconvolution`](@ref) instead poses it as a **constrained, regularized optimization problem**,
following Dion et al. (2022, 2023, 2024).

The optimization below relies throughout on the root-mean-square, [`rms`](@ref), to turn a vector
into a scalar deviation. For a vector `x` of length `n`, it is:

```math
\mathrm{rms}(x) = \sqrt{\frac{1}{n}\sum_{i=1}^n x_i^2}
```

so the deviation between two vectors `x` and `y` is `rms(x .- y)`.

## Node parameterization

Rather than optimizing one value of `g` per data point (typically thousands of samples),
the optimization solves for `g` at only `n` log-spaced node indices `τ` ([`set_nodes`](@ref) on
`1:nData`), giving `ĝ(τ)` then interpolates ([`PCHIPInterpolation.jl`](https://github.com/gerlero/PCHIPInterpolation.jl),
shape-preserving) onto the full index axis to obtain `ĝ(t)`. Log-spacing matches the thermal
response function's own behavior: it varies quickly at early times and slowly at late times, so
fewer nodes are needed there.

## Initial guess

The optimization is initialized with a 2-parameter exponential-integral fit
(Dion et al. 2024, Eqs. 20–21):

```math
\hat x_{1,2} = \min_{x_1, x_2} \Big[\mathrm{rms}\Big((f*\hat g_0)(t) - T_{exp} \Big)\Big] \\

g_0(t) = x_1 \cdot E_1(x_2 / t),
```

with `E₁` the exponential integral, fit to `Texp` by a derivative-free Nelder-Mead search. This step is akin to fitting
an infinite line source model (through the exponential integral) to the experimental data to have a first fit.

## Multi-objective function

The main optimization, over the node values `ĝ(τ)`, minimizes a weighted combination of three
terms (following Dion et al. 2023):

```math
\hat g(\tau) = \arg\min_{\hat g(\tau)} E \quad \text{s.t. } C_1(\tau),\ C_2(\tau),
```

with the objective function being:

```math
E = w_1 \, \mathrm{rms}\big(w_a \cdot (T_{conv} - T_{exp})\big) + w_2 \, \mathrm{rms}(g') + w_3 \, \mathrm{rms}(g''),
```

where `T_conv` is the temperature reconstructed from the current iterate of `ĝ(t)` by convolution,
`g'`/`g''` are its first and second discrete derivatives (regularization terms, penalizing
non-physical oscillation), and `wₐ` is an array weight emphasizing the leading 5% of samples,
where the response function changes fastest. The scalar weights `w₁, w₂, w₃` are set from the
initial guess so that the three terms contribute roughly 70%/15%/15% to `E` (Dion et al., 2022, 2024).

## Constraints

The `c` keyword of [`deconvolution`](@ref) selects which physical constraints are imposed on `g`
as linear inequalities:

- `c = 0`: none.
- `c = 1`: a positive first derivative (`ĝ(τ)` strictly growing).
- `c = 2`: adds a negative second derivative past roughly 3 h from the start of the record, where
  the thermal response function is expected to be concave (`ĝ'(τ)` is strictly decreasing after a node `j`).

The constraints are set as linear inequality constraints in the optimization through a matrix formulation such that `A⋅x ≤ 0`, and are defined as:

```math
C_1 \rightarrow 0 \lt \hat g (\tau_i) \lt \hat g (\tau_{i+1})\ \forall i \in [1,\ n_{\tau} - 1] \\
C_2 \rightarrow 0 \lt \hat g' (\tau_{i+1}) \lt \hat g' (\tau_i)\ \forall i \in [1,\ n_{\tau} - 1]
```

## Solver

The constrained problem is solved with
[Optimization.jl](https://github.com/SciML/Optimization.jl), using NLopt's SLSQP algorithm
(sequential quadratic programming, via
[OptimizationNLopt.jl](https://github.com/SciML/Optimization.jl)) and finite-difference gradients
([FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl)). As an active-set method, SLSQP handles
a starting point on or near a constraint boundary natively, so no interior-point-specific
safeguard on the initial guess is needed.
