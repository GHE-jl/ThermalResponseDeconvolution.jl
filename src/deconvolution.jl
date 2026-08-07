using Optimization, OptimizationNLopt, NLopt
using FiniteDiff
using SpecialFunctions: expinti
using PCHIPInterpolation: Interpolator
using LinearAlgebra: mul!, diagm

"""
    deconvolution(t, f, Texp; n=35, c=2)

Optimization algorithm to recover a thermal response function from paired fluid temperature and
heat load data, whether collected during a thermal response test (TRT) or from normal GSHP
operating data. The inputs must have a constant time step for the convolution to be valid. The
estimated convolved temperature variation can be recovered with `T̂ = convolution(f, ĝ)`.
# Arguments
    - `t`: Time array, starting at 0, with constant time step [s].
    - `f`: Incremental perturbation function (`diff([0; Tin .- Tout])`) [degC].
    - `Texp`: Measured temperature variation (`Tout .- Tout[1]`) [degC].
    - `n`: Number of nodes (default: 35).
    - `c`: Choice of constraints (default: 2):
        - 0: No constraint
        - 1: Positivity (positive first derivative)
        - 2: Positivity (positive first derivative) + negative second derivative
# Returns
    - `ĝ`: Interpolated estimated thermal response function obtained by deconvolution [-]
    - `gOpt`: Optimized thermal response function values at the node indices, before interpolation
        (see [`set_nodes`](@ref) to recover their positions) [-]
# Reference
    - Dion, G., Pasquier, P., & Marcotte, D. (2022). Deconvolution of experimental thermal response
        test data to recover short-term g-function. Geothermics, 100, 102302.
        https://doi.org/10.1016/j.geothermics.2021.102302
    - Dion, G., Pasquier, P., & Marcotte, D. (2024). Application of deconvolution to interpretation
        of distributed thermal response test (DTRT) and to determination of thermal conductivity
        profiles. Applied Thermal Engineering, 236, 121680.
        https://doi.org/10.1016/j.applthermaleng.2023.121680
"""
function deconvolution(t::AbstractVector{<:Real}, f::AbstractVector{<:Real},
    Texp::AbstractVector{<:Real}; n::Integer=35, c::Integer=2)

    # Initial setup
    nData = length(Texp)                # Number of data points
    idall = collect(1.0:nData)          # All indices of the data points
    id = set_nodes(nData, n)            # Log-spaced node indices for the optimization problem
    idnodes = Float64.(id)              # Node indices for the optimization problem

    # Initial guess of the thermal response function
    g₀ = _deconv_ini(f, idall, Texp)
    g₀i = g₀[id]

    # Set weights for the multi-objective function
    w, wₐ = _set_weights(g₀, Texp, f)

    # Build the linear inequality constraint matrix
    A = _const_matrix(t, id, c)
    p = (Texp=Texp, f=f, idnodes=idnodes, idall=idall, w=w, wₐ=wₐ, A=A)

    # Set lower and upper bounds for the optimization problem
    lb = zeros(n)
    ub = fill(Inf, n)

    # Set up the optimization problem and solve it as a function of the choice of constraints
    if c == 0
        optf = OptimizationFunction(_obj_fun, Optimization.AutoFiniteDiff())
        prob = OptimizationProblem(optf, g₀i, p; lb=lb, ub=ub)
    else
        m = size(A, 1)
        optf = OptimizationFunction(_obj_fun, Optimization.AutoFiniteDiff(); cons=_cons_derivative)
        prob = OptimizationProblem(optf, g₀i, p; lb=lb, ub=ub, lcons=fill(-Inf, m), ucons=zeros(m))
    end

    # Solve the optimization problem using NLopt's SLSQP algorithm
    sol = solve(prob, NLopt.LD_SLSQP(); maxiters=1000)
    gOpt = sol.u

    # Interpolate the optimized node values to obtain the estimated thermal response function at
    # all indices
    interp = Interpolator(idnodes, gOpt)
    ĝ = interp.(idall)

    return ĝ, gOpt
end

"""
    rms(x)

Root-mean-square of `x`.
# Arguments
    - `x`: Input vector [-].
# Returns
    - Root-mean-square of `x` [-].
"""
rms(x::AbstractVector{<:Real}) = sqrt(sum(y^2 for y in x) / length(x))

"""
    set_nodes(nData, n0)

Return `n0` log-spaced, unique integer node indices on `1:nData`.
# Arguments
    - `nData`: Total number of data points [-].
    - `n0`: Requested number of nodes [-].
# Returns
    - `id`: Log-spaced, unique integer node indices on `1:nData` [-].
"""
function set_nodes(nData::Integer, n0::Integer)
    # 1:nData holds only nData distinct integers, so n0 > nData can never be reached and would
    # otherwise loop forever below.
    n0 <= nData || throw(ArgumentError("n0 ($n0) cannot exceed nData ($nData)"))

    n_tmp = n0 - 1
    id = Int[]
    while length(id) != n0
        id = unique(round.(Int, exp10.(range(0.0, stop=log10(nData), length=n_tmp))))
        n_tmp += 1
    end
    return id
end

"""
    _deconv_ini(f, idall, Texp)

Fit `g0(idx) = x1 * E1(x2 / idx)` to the measured temperatures using a derivative-free
Nelder-Mead search, and return the resulting initial guess of the thermal response function. This
is the 2-parameter exponential-integral fit of Dion et al. (2022, Eq. 20-21). Here E₁ = -expinti(-z).
# Arguments
    - `f`: Incremental perturbation function [degC].
    - `idall`: All indices of the data points [-].
    - `Texp`: Measured temperature variation [degC].
# Returns
    - `g0`: Exponential-integral initial guess of the thermal response function, sampled at every
        index [-].
"""
function _deconv_ini(f::Vector{Float64}, idall::Vector{Float64}, Texp::Vector{Float64})
    function obj_ini(x, _p)
        g = x[1] .* -expinti.(-x[2] ./ idall)
        return rms(convolution(f, g) .- Texp)
    end

    optf = OptimizationFunction(obj_ini)
    prob = OptimizationProblem(optf, [1.0, 1.0], nothing)
    sol = solve(prob, NLopt.LN_NELDERMEAD(); maxiters=2000, xtol_rel=1e-6)
    return sol.u[1] .* -expinti.(-sol.u[2] ./ idall)
end

"""
    _set_weights(g₀, Texp, f)

Set the scalar weights and the array weight of the multi-objective function from the initial guess.
The scalar weights are set to balance the contributions of the temperature misfit, first- and
second-derivative terms.
# Arguments
    - `g₀`: Initial guess of the thermal response function, sampled at every index [-].
    - `Texp`: Measured temperature variation [degC].
    - `f`: Incremental perturbation function [degC].
# Returns
    - `w`: Scalar weights of the temperature-misfit, first- and second-derivative terms [-].
    - `wₐ`: Array weight emphasizing the leading 5% of samples [-].
# Reference
    - Dion, G., Pasquier, P., Marcotte, D., & Beaudry, G. (2023). Multi-deconvolution in
        non-stationary conditions applied to experimental thermal response test analysis to obtain
        short-term transfer functions. Science and Technology for the Built Environment, 30(3),
        1–14. https://doi.org/10.1080/23744731.2023.2217729
"""
function _set_weights(g₀::Vector{Float64}, Texp::Vector{Float64}, f::Vector{Float64})
    # Set the scalar weights of the multi-objective function.
    w₀ = [0.7, 0.15, 0.15]

    # First- and second-derivative terms
    dg₀ = diff([0.0; g₀])
    ddg₀ = diff(dg₀)

    # Calculate the root-mean-square of the three terms in the multi-objective function
    e = [rms(convolution(f, g₀) .- Texp), rms(dg₀), rms(ddg₀)]

    # Calculate the scalar weights to balance the contributions of the three terms
    prop = w₀ ./ sum(e)
    w = w₀ ./ prop

    # Set the array weight emphasizing the leading 5% of samples
    nt = length(g₀)
    n5 = round(Int, nt * 0.05)
    wₐ = [fill(3.0, n5); fill(1.0, nt - n5)]

    return w, wₐ
end

"""
    _const_matrix(t, id, c)

Build the linear inequality matrix `A` such that the constraints are
`A * g .<= 0`:
  - First constraint (`c >= 1`): positive first derivative (strictly growing).
  - Second constraint (`c == 2`): negative second derivative on the section past
    roughly 3 h from the start of the record.
# Arguments
    - `t`: Time array, starting at 0, with constant time step [s].
    - `id`: Node indices for the optimization problem [-].
    - `c`: Choice of constraints (0, 1 or 2, see [`deconvolution`](@ref)).
# Returns
    - `A`: Linear inequality constraint matrix, `0×n` when `c == 0` [-].
"""
function _const_matrix(t::Vector{Float64}, id::Vector{Int}, c::Integer)
    n = length(id)

    """
        const_1()

    Build the linear inequality constraint matrix for the first constraint (positive first
    derivative) following a finite-difference scheme in matrix form.
    # Returns
        - `a1`: Linear inequality constraint matrix for the first constraint [-].
    """
    function const_1()
        e = ones(Float64, n)
        h = diff([0.0; id])
        a1 = diagm(0 => -e, 1 => e[1:(n-1)])
        a1 = -a1[1:(end-1), :] ./ h[1:(end-1)]
        return a1
    end

    """
        const_2()

    Build the linear inequality constraint matrix for the second constraint (negative second
    derivative) following a finite-difference scheme in matrix form. The constraint is only applied
    to the section of the thermal response function past roughly 3 h from the start of the record,
    as the second derivative is expected to be positive (or varying) in the early part of the
    thermal response function.
    # Returns
        - `a2`: Linear inequality constraint matrix for the second constraint [-].
    """
    function const_2()
        idnodes = maximum(t) >= 3600 * 3 ? argmin(abs.(3600 * 3 .- t[id])) : 0

        a2 = zeros(n - 2 - idnodes, n)
        cc = @. (id[(idnodes+2):(end-1)] - id[(idnodes+1):(end-2)]) /
                (id[(idnodes+1):(end-2)] - id[idnodes:(end-3)])
        for kk in 1:(n-2-idnodes)
            a2[kk, (idnodes+kk):(idnodes+kk+2)] = [cc[kk], -(1 + cc[kk]), 1.0]
        end
        return a2
    end

    # Build the linear inequality constraint matrix `A` based on the choice of constraints
    A = c == 0 ? zeros(0, n) :
        c == 1 ? const_1() :
        [const_1(); const_2()]

    return Matrix{Float64}(A)
end

"""
    _cons_derivative(res, u, p)

In-place constraint function for Optimization.jl: `res .= A * u`, combined with `-Inf .<= res .<= 0`
to express `A * g .<= 0`. This is used to enforce the linear inequality constraints on the
optimization problem, where `A` is the linear inequality constraint matrix built by `_const_matrix`.
# Arguments
    - `res`: Output vector, overwritten in place with `A * u` [-].
    - `u`: Node values of the thermal response function at the current iterate [-].
    - `p`: Named tuple of problem parameters, including the constraint matrix `p.A`.
"""
function _cons_derivative(res, u, p)
    mul!(res, p.A, u)
    return nothing
end

"""
    _obj_fun(u, p)

Multi-objective function minimized by the deconvolution: weighted temperature misfit plus first- and
second-derivative regularization, with the node values `u` interpolated over the full index axis.
# Arguments
    - `u`: Node values of the thermal response function at the current iterate [-].
    - `p`: Named tuple of problem parameters (measured data, incremental perturbation
        function, weights, indices).
# Returns
    - `e`: Value of the multi-objective function [-].
"""
function _obj_fun(u, p)
    # Interpolate the node values to obtain the estimated thermal response function at all indices
    interp = Interpolator(p.idnodes, u)
    g = interp.(p.idall)

    # Calculate the first- and second-derivative terms
    dg = diff(vcat(zero(eltype(g)), g))
    ddg = diff(dg)

    # Calculate the reconstructed temperature response from the convolution of the incremental
    # perturbation function with the current thermal response function
    Tconv = convolution(p.f, g)

    # Calculate the multi-objective function value
    e = p.w[1] * rms(p.wₐ .* (Tconv .- p.Texp))
    e += p.w[2] * rms(dg)
    e += p.w[3] * rms(ddg)
    return e
end
