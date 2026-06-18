using Optimization, OptimizationOptimJL
using FiniteDiff
using SpecialFunctions
using PCHIPInterpolation
using LinearAlgebra
using FFTW

"""
    deconvolution(q, T; n=35, c=2, show_ini=false, show_opt=false)

Optimization algorithm to perform deconvolution on the experimental data extracted from a TRT to
recover a short-term g-function (STgF). The function computes the estimated STgF and the estimated
convolved temperature variations.

Arguments
- `q`: Incremental load function.
- `T`: Experimental temperature variations.
- `n`: Number of nodes (default: 35).
- `c`: Choice of constraints (default: 2):
    - 0: No constraint
    - 1: Positivity (lower bound) + positive first derivative
    - 2: Positivity + positive first derivative + negative second derivative
- `show_ini`: Display the initial-guess figure (default: false).
- `show_opt`: Display the final figure (default: false).

Outputs
- `ĝ`: Interpolated estimated STgF obtained by deconvolution [-]
- `Tconv`: Estimated convolved temperature variation (Tconv = f * ĝ) [degC]

Notes
The inputs must have a constant time step for the convolution to be valid.
Interpolation, derivatives and constraints are all expressed over the sample
*indices* `1:nData` (as in the reference MATLAB implementation), not over the
physical time axis.

Because the convolution relies on FFTW (which does not accept `ForwardDiff`
dual numbers), the gradients/Hessians required by the interior-point solver are
obtained by finite differences (`AutoFiniteDiff`), mirroring MATLAB's `fmincon`.

# Reference
Dion, G., Pasquier, P., & Marcotte, D. (2022). Deconvolution of experimental
thermal response test data to recover short-term g-function. Geothermics, 100,
102302. https://doi.org/10.1016/j.geothermics.2021.102302
"""
function deconvolution(q, T; n::Integer=35, c::Integer=2,
    show_ini::Bool=false, show_opt::Bool=false)

    # 0. Validation
    data_validation(q, T)
    option_validation(n, c)

    # 1. Initial parameters
    nData = length(T)
    idall = collect(1.0:nData)              # Index axis used everywhere

    # 2. Prepare the incremental load function and cache its FFT
    f_fft = define_f(q)

    # 3. Node positions (log-spaced integer indices on 1:nData)
    id = set_nodes(nData, n)
    idnodes = Float64.(id)

    # 4. Initial guess via the exponential-integral fit
    g₀ = deconv_ini(Deconv0Params(T, f_fft, idall))
    show_ini && show_fig(data.t, f_fft.f, g₀, convolution_g(f_fft, g₀))
    g₀i = max.(g₀[id], 1e-8)                # Strictly interior starting point

    # 5. Weights of the multi-objective function
    w, wₐ = set_weights(g₀, T, f_fft)

    # 6. Linear inequality constraints (A*g .<= 0) and bounds
    A = const_matrix(data.t, id, c)
    p = DeconvParams(T, f_fft, idnodes, idall, w, wₐ, A)

    lb = zeros(n)                           # Positivity lower bound
    ub = fill(Inf, n)                       # No upper bound

    # 7. Build and solve the constrained problem with Optim's interior point
    if c == 0
        optf = OptimizationFunction(obj_fun, Optimization.AutoFiniteDiff())
        prob = OptimizationProblem(optf, g₀i, p; lb=lb, ub=ub)
    else
        m = size(A, 1)
        optf = OptimizationFunction(obj_fun, Optimization.AutoFiniteDiff();
            cons=cons_derivative)
        prob = OptimizationProblem(optf, g₀i, p; lb=lb, ub=ub,
            lcons=fill(-Inf, m), ucons=zeros(m))
    end

    sol = solve(prob, IPNewton(); maxiters=1000)
    gOpt = sol.u

    # 8. Final interpolation and convolution
    interp = Interpolator(idnodes, gOpt)
    ĝ = interp.(idall)
    T̂ = collect(convolution_g(f_fft, ĝ))
    show_opt && show_fig(data.t, f_fft.f, ĝ, T̂)

    return ĝ, T̂
end

function data_validation(q, T)
    if length(q) != length(T)
        error("the incremental load and temperature vectors must have the same length.")
    elseif length(q) < 10
        @warn "the input data are very short; the deconvolution may be inaccurate."
    end
end

function option_validation(n, c)
    if n < 0
        error("number of nodes must be a positive integer.")
    elseif n <= 15
        @warn "the number of nodes should be larger than 20, ideally between 30 and 60"
    elseif n > 100
        @warn "the number of nodes is large (>100); it should be between 30 and 60"
    end

    if !(c in (0, 1, 2))
        error("optional input for constraints must be 0, 1 or 2.")
    end
end

"""
    define_f(data)

Build the incremental load function `f = diff([0; Tin - Tout])` and cache its
zero-padded real FFT for fast convolution.
"""
function define_f(data::TRTData)
    f = diff([0.0; data.Tin .- data.Tout])
    f_pad = zeros(Float64, 2 * length(data.t) - 1)
    copyto!(f_pad, f)
    fft_plan = plan_rfft(f_pad)
    F_f = fft_plan * f_pad
    return f_FFT(f, f_pad, fft_plan, F_f)
end

"""
    set_nodes(nData, n0)

Return `n0` log-spaced, unique integer node indices on `1:nData`.
"""
function set_nodes(nData::Integer, n0::Integer)
    n_tmp = n0 - 1
    id = Int[]
    while length(id) != n0
        id = unique(round.(Int, exp10.(range(0.0, stop=log10(nData), length=n_tmp))))
        n_tmp += 1
    end
    return id
end

"""
    convolution_g(f_fft, g)

Non-circular convolution of the cached `f` with `g`, using the zero-padded
spectral product and the cached rFFT plan. The result is truncated to the
original signal length.
"""
function convolution_g(f_fft::f_FFT, g::AbstractVector{Float64})
    n_full = length(f_fft.f_pad)
    g_pad = zeros(Float64, n_full)
    copyto!(g_pad, g)
    F_g = f_fft.fft_plan * g_pad
    F_g .*= f_fft.F_f
    y = irfft(F_g, n_full)
    return @view y[1:length(f_fft.f)]
end

# E1 exponential integral, expressed from Ei: E1(z) = -Ei(-z) for z > 0.
expint_E1(z) = -expinti(-z)

"""
    deconv_ini(p0)

Fit `g0(idx) = x1 * E1(x2 / idx)` to the experimental temperatures using a
derivative-free Nelder–Mead search (no gradients, so FFTW is safe here), and
return the resulting initial guess of the transfer function.
"""
function deconv_ini(p0::Deconv0Params)
    function obj_ini(x, p)
        g = x[1] .* expint_E1.(x[2] ./ p.idall)
        return rms(convolution_g(p.f_fft, g) .- p.Texp)
    end

    optf = OptimizationFunction(obj_ini)        # No AD: Nelder–Mead is derivative-free
    prob = OptimizationProblem(optf, [1.0, 1.0], p0)
    sol = solve(prob, Optim.NelderMead(); g_tol=1e-3)
    return sol.u[1] .* expint_E1.(sol.u[2] ./ p0.idall)
end

"""
    set_weights(g₀, Texp, f_fft)

Set the scalar weights and the array weight of the multi-objective function
from the initial guess, following the reference implementation.

Note: as in the MATLAB code, `w = w0 ./ (w0 ./ sum(e))` collapses to the scalar
`sum(e)` for all three terms; the relative weighting between the three objective
terms is therefore driven by their own magnitudes, while the leading-samples
emphasis comes from the array weight `wₐ`.
"""
function set_weights(g₀::Vector{Float64}, Texp::Vector{Float64}, f_fft::f_FFT)
    w₀ = [0.7, 0.15, 0.15]

    dg₀ = diff([0.0; g₀])
    ddg₀ = diff(dg₀)

    e = [
        rms(convolution_g(f_fft, g₀) .- Texp),
        rms(dg₀),
        rms(ddg₀),
    ]

    prop = w₀ ./ sum(e)
    w = w₀ ./ prop

    nt = length(g₀)
    n5 = round(Int, nt * 0.05)
    wₐ = [fill(3.0, n5); fill(1.0, nt - n5)]

    return w, wₐ
end

"""
    const_matrix(t, id, c)

Build the linear inequality matrix `A` such that the constraints are
`A * g .<= 0`:

  - First constraint (`c >= 1`): positive first derivative (strictly growing).
  - Second constraint (`c == 2`): negative second derivative on the section past
    roughly 3 h of test.

Returns a `0×n` matrix when `c == 0`.
"""
function const_matrix(t::Vector{Float64}, id::Vector{Int}, c::Integer)
    n = length(id)

    function const_1()
        e = ones(Float64, n)
        h = diff([0.0; id])
        a1 = diagm(0 => -e, 1 => e[1:(n-1)])
        a1 = -a1[1:(end-1), :] ./ h[1:(end-1)]
        return a1
    end

    function const_2()
        # Index of the node closest to 3 h of test (1-based), or skip if too short
        idnodes = maximum(t) >= 3600 * 3 ? argmin(abs.(3600 * 3 .- t[id])) : 0

        a2 = zeros(n - 2 - idnodes, n)
        cc = @. (id[(idnodes+2):(end-1)] - id[(idnodes+1):(end-2)]) /
                (id[(idnodes+1):(end-2)] - id[idnodes:(end-3)])
        for kk in 1:(n-2-idnodes)
            a2[kk, (idnodes+kk):(idnodes+kk+2)] = [cc[kk], -(1 + cc[kk]), 1.0]
        end
        return a2
    end

    A = c == 0 ? zeros(0, n) :
        c == 1 ? const_1() :
        [const_1(); const_2()]

    any(isnan, A) && error("NaN present in the constraint matrix.")
    return Matrix{Float64}(A)
end

"""
    cons_derivative(res, u, p)

In-place constraint function for Optimization.jl: `res .= A * u`, combined with
`-Inf .<= res .<= 0` to express `A * g .<= 0`.
"""
function cons_derivative(res, u, p::DeconvParams)
    mul!(res, p.A, u)
    return nothing
end

"""
    obj_fun(u, p)

Multi-objective function minimized by the deconvolution: weighted temperature
misfit plus first- and second-derivative regularization, with the node values
`u` interpolated over the full index axis.
"""
function obj_fun(u, p::DeconvParams)
    interp = Interpolator(p.idnodes, u)
    g = interp.(p.idall)

    dg = diff(vcat(zero(eltype(g)), g))
    ddg = diff(dg)

    Tconv = convolution_g(p.f_fft, g)

    e = p.w[1] * rms(p.wₐ .* (Tconv .- p.Texp))
    e += p.w[2] * rms(dg)
    e += p.w[3] * rms(ddg)
    return e
end
