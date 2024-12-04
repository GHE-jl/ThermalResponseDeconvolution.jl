using Optim
using PCHIPInterpolation
using LinearAlgebra
using GHEDeconvolutions

function deconvolution(data::TRTData; n::Int=35, c::Int=2)
    """
        deconvolution(t, f, T, n=35, c=2, show_ini=false, show_opt=false)
    
    Optimization algorithm to perform deconvolution on the experimental data
    extracted from a TRT to recover a short-term g-function (STgF).The function
    computes the estimated STgF and the estimated convolved temperature variations.
    Inputs:
        - t: Time array, starting at 0, with constant time step (see [1]) (s)
        - f: Incremental perturbation function (see note [2] below) [degC, W, W/m]
        - TExp: Experimental temperature variation (T_exp = T_out-T_0) [degC]
        - [Optionnal Inputs]:
          - n: Number of nodes (default: 35)
          - c: Choice of constraints (default: 2):
              0: No constraint
              1: Positivity
              2: Positivity and negative derivative slope
         - show_ini: Display initial guess figure (default: false)
         - show_opt: Display final guess figure (default: false)
    Outputs:
        - g: Interpolated estimated STgF obtained by deconvolution [-]
        - T_conv: Estimated convolved temperature variation (T_conv=f*ghat) [degC]
        - [Optionnal Outputs]:
           - Time: Computing time to converge [s]
           - [Id,gOpt]: Node positions and optimized transfer function values
           - Out: Output variables of the solver [-]
           - IterOut: Optimization variables at each iteration [-]
               - [Iteration,Time,fval,funccount]

     Notes:
     [1]: The inputs "f" and "TExp" must have constant time step for the
     deconvolution to work effectively, because of how convolution works. Otherwise, these
     vectors must be interpolated.
     [2]: The input "f" is the time derivative of the perturbation function, which
     can be T_{in}-T_{out} [degC], Q [W] or q [W/m]. The transfer function 
     will be defined for the units of input used.

     Reference: Dion, G., Pasquier, P., & Marcotte, D. (2022). Deconvolution of
     experimental thermal response test data to recover short-term g-function.
     Geothermics, 100, 102302. https://doi.org/10.1016/j.geothermics.2021.102302

     Author: Gabriel Dion
     Date: 2024-12
    """

    # 0. Data validation
    data_validation(data)
    option_validation(n, c)

    # 1. Preallocation
    g_0 = similar(data.t)
    g = similar(data.t)
    T_conv = similar(data.t)

    # 2. Prepare the incremental impulse function (here as T_in-T_out)
    f_fft = define_f(data)

    # 3. Define nodes positions
    id = set_nodes(length(data.t), n)
    g_0_id = similar(id)

    # 4. Compute initial solution
    deconv_ini!(g_0, data, f_fft)
    copyto!(g_0_id, g_0[id])

    # 5. Set weights for the multi-objective function
    w, w_array = set_weights(g_0, f, temperature)

    # 6. Set linear inequality constraints and bounds
    lb = fill(0.0, n)

    a, b = const_derivative(t, id, c)
    # Note: Usually, the 2 constraints give best results. Either using 1 or 0 constraints
    # may be usefull in some cases.

    # 7. Optimization
    x_opt = optimize(
        (x) -> obj_fun(x, id, idall, f, temperature, w, w_array),
        #TwiceDifferentiableConstraints(a, b, lb, fill(Inf, length(lb))),
        g_0_id,
        #IPNewton(),
        Optim.Options(; iterations=1000, show_trace=true)
    )

    # 8. Final interpolation and convolution
    interp = Interpolator(id, x_opt.minimizer)
    ĝ = interp.(idall)
    T_conv = convolution(f, ĝ)

    # 8. Function outputs
    return ĝ, T_conv, x_opt
end

function data_validation(data::TRTData)
    """
        data_validation(data)
    
    Function that ensures that the input data respect correct form. The input is based on the
    structure constructed with the TRT dataset.
    """

    # Check for imaginary, NaN or infinite value
    if any(!isreal(getfield(data, field)) for field in fieldnames(typeof(data)))
        error("inputs contain imaginary values.")
        # Check for NaN values
    elseif any(isnan.(getfield(data, field)) for field in fieldnames(typeof(data)) if fieldtype(typeof(data), field) <: AbstractFloat)
        error("inputs contains NaN values.")
        # Check for NaN values
    elseif any(isinf(getfield(data, field)) for field in fieldnames(typeof(data)) if fieldtype(typeof(data), field) <: AbstractFloat)
        error("inputs contains infinite values.")
    end
end

function option_validation(n, c)
    """
        option_validation(n, c)
    
    Function that ensures that the optional inputs are correctly set for the deconvolution function.
    """
    # Checking on the number of nodes
    if n < 0
        error("number of nodes must be a positive integer.")
    elseif n <= 15
        @warn "the number of nodes should be larger than 20, ideally between 30 to 60"
    elseif n > 100
        @warn "the number of nodes is large (>100). It should be between 30 and 60"
    end

    # Checking on the constraint selection
    if c < 0 || c > 2
        error("optionnal input for constraints must be 0, 1 or 2.")
    elseif !isinteger(c)
        error("optionnal constraints must be an integer.")
    end
end

function define_f(data::TRTData)
    f = diff([0; data.Tin - data.Tout])
    f_pad = zeros(Float64, 2 * length(data.t) - 1)
    copyto!(f_pad, f)
    fft_plan = plan_rfft(f_pad)
    F_f = fft_plan * f_pad
    return f_FFT(f, f_pad, fft_plan, F_f)
end

function set_nodes(nt::Int, n0::Int)
    """
        set_nodes(nt, n0)
    Function that sets the position of the nodes on the transfer function based on
    the number of nodes asked by the user.
    Inputs:
        - nt: Total number of data in the input vectors [-]
        - n0: User defined number of nodes on the transfer function [-]
    Output:
        - id: A vector of length "n" of node positions on the transfer function [-]
    """

    n_tmp = n0 - 1
    id = Vector{Integer}(undef, n_tmp)

    while length(id) != n0
        id = unique(round.(Int, exp10.(range(0, stop=log10(nt), length=n_tmp))))
        n_tmp += 1
    end
    return id
end

function convolution_g(f_fft::f_FFT, g::Vector{Float64})
    """
        convolution_g(f_pad, fft_plan, F_f, g)

    Function that convolves an already prepared vector "f" to the transfer function "g" in the
    spectral domain.
    """
    n_full = length(f_fft.f_pad)
    # Prepare the transfer function
    g_pad = zeros(Float64, n_full)
    copyto!(g_pad, g)
    F_g = f_fft.fft_plan * g_pad

    # Element-wise multiplication (in-place)
    F_g .*= f_fft.F_f

    # Inverse FFT to get convolution resuls
    y = irfft(F_g, n_full)
    return @view y[1:length(f_fft.f)]
end

function deconv_ini!(g_0::Vector{Float64}, data::TRTData, f_fft::f_FFT)
    """
        deconv_ini!(g_0, data, f_fft)
    
    Function that ajust an exponential integral equation to obtain a first approximation of
    the ground heat exchanger transfer function. The algorithm uses a convolution product
    to compute the objective function.
    Inputs:
        - idall: Numerized time vector [-]
        - f: Incremental perturbation function [degC, W, W/m]
        - temperature: Experimental temperature variation at the borehole outlet [degC]
    Output:
        - g_0: Initial guess of the transfer function based on an exponential integral [-]
    """

    # Define objective function (RMSE) between model and experimental values.
    function obj_fun(x)
        """Objective function: RMSE between experimental and computed temperatures."""
        return rms(convolution_g(f_fft, x[1] .* expint.(x[2] ./ data.t)) .- data.Texp)
    end

    # Define the optimization on 2 variables
    x0 = [1.0, 1.0]
    x_opt = optimize(obj_fun, x0)
    g_0 .= x_opt.minimizer[1] .* expint.(x_opt.minimizer[2] ./ data.t)
end

function set_weights(g_0::Vector{Float64}, data::TRTData, f_fft::f_FFT)
    """
        set_weights(g_0, data, f_fft)
    
    Function that computes the weights of each term in the deconvolution objective function
    based on the weights obtained using the initial transfer function guess. The desired
    weights are so that e[1]~0.7, e[2]~0.15, e[3]~0.15.
    Inputs:
        - g_0: Initial transfer function guess based on an exponential integral [-]
        - f: Incremental perturbation function [degC, W, W/m]
        - temperature: Experimental temperature variation at the borehole outlet [degC]
    Outputs:
        - w: Weights of each term in the objective function (3x1 vector) [-]
        - w_array: Array weights to emphazise the reconstruction of the transfer function's
            initial impulse. 
    """

    # Define proportion for terms in the objective function
    w_0 = [0.7, 0.15, 0.15]

    # Compute initial transfer functionn derivatives
    dg_0 = diff([0; g_0])
    ddg_0 = diff(diff([0; g_0]))

    # Compute each terms of the objective function
    nt = length(g_0)
    e = zeros(3)
    e[1] = rms(convolution_g(f_fft, g_0) - data.Texp)
    e[2] = rms(dg_0 .^ 2)
    e[3] = rms(ddg_0 .^ 2)

    # Define the objective functions weights
    prop = w_0 ./ sum(e)
    w = w_0 ./ prop

    # Set array weight
    w_array = [3 * ones(round(Int, nt * 0.05)); ones(round(Int, nt * 0.95))]

    return w, w_array
end

function const_derivative(t::Vector{Float64}, id::Vector{Integer}, cnst::Int)
    """
        const_derivative(t, id, cnst)
    Set the linear inequality constraints for the problem Ax <= b used to constrain the
    deconvolution algorithm. There are two different constraints applied on the transfer 
    function:
        1. A positive first derivative constraint (strictly growing function)
        2. A negative second derivative constraint after an inflexion point
    Inputs:
        - t: Time vector starting at 0, with constant time step [s]
        - id: Nodes position on the transfer function [-]
        - cnst: Constraint choice (default to 2):
            0: No constraint applied to the deconvolution algorithm
            1: Only the first constraint is applied
            2: The two available constraints are applied.
    Outputs:
        - a: Linear inequality constraint matrix in Ax <= b of size (m x n), where "m" is
            the number of inequalities and "n" is the number of nodes (lenght(id)).
        - b: A constraint vector of (m x 1) in Ax <= b to respect in the optimization. 
    """

    # Input parameters
    n = length(id)

    function const_1(id::Vector{Integer}, n::Int)
        """
        First constraint: positive first derivative
        """
        # Construct parameters
        e = ones(Float64, n)
        h = diff([0; id])
        # Build matrix and constraint vector
        a1 = diagm(0 => -e, 1 => e[1:(n-1)])
        a1 = -a1[1:(end-1), :] ./ h[1:(end-1)]
        b1 = zeros(Float64, n - 1)
        return a1, b1
    end

    function const_2(t::Vector{Float64}, id::Vector{Integer}, n::Int)
        """
        Second constraint: negative second derivative on a SpecialFunctions
        """
        # Find time at around 3 hours of test
        if maximum(t) >= 3600 * 3
            id_nodes = argmin(abs.(3600 * 3 .- t[id]))
        else
            id_nodes = 0
        end

        # Construct the second derivative a2 and b2
        c = @. (id[(id_nodes+2):(end-1)] - id[(id_nodes+1):(end-2)]) /
               (id[(id_nodes+1):(end-2)] - id[id_nodes:(end-3)])

        a2 = zeros(n - 2 - id_nodes, n)
        for ii in 1:(n-2-id_nodes)
            a2[ii, (id_nodes+ii):(id_nodes+ii+2)] = [c[ii], -(c[ii] .+ 1), 1.0]
        end

        b2 = zeros(Float64, n - 2 - id_nodes)
        return a2, b2
    end

    # Setting output constraints depending on the optional input `c`
    if cnst == 2
        a1, b1 = const_1(id, n)
        a2, b2 = const_2(t, id, n)
        a = [a1; a2]
        b = [b1; b2]
    elseif cnst == 1
        a, b = const_1(id, n)
    elseif cnst == 0
        a = []
        b = []
    else
        error("Constraint number must be either 0, 1 or 2.")
    end

    # Verification
    if any(isnan, a) || any(isnan, b)
        error("Presence of NaN in the output.")
    end

    return a, b
end

function obj_fun(g_opt::Vector{Real}, id::Vector{Int}, idall::Vector{Int}, f::Vector{Real},
    T::Vector{Real},
    w,
    w_array
)
    """
        obj_fun(g_opt, id, idall, f, T, w, w_array)
    Computes the multi-objective function that allows to optimize the position of nodes on
    the transfer function, so that experimental temperature are recreated.
    Inputs:
        - g_opt: The noded transfer function that is being iterated on [-]
        - id: Nodes position on the transfer function [-]
        - idall: Numerized time vector [-]
        - f: Incremental perturbation function [degC, W, W/m]
        - T: Experimental temperature variation at the borehole outlet [degC]
        - w: Weights of each term in the objective function (3x1 vector) [-]
        - w_array: Array weights to emphazise the reconstruction of the transfer function's
            initial impulse.
    Output:
        - sum(e): The value of the objective function to minimize by the deconvolution
    """

    # Interpolate the transfer function
    interp = Interpolator(id, g_opt)
    g = interp.(idall)

    # Compute initial transfer functionn derivatives
    dg = diff([0; g])
    ddg = diff(diff([0; g]))

    # Compute each terms of the objective function
    nt = length(g)
    e = zeros(3)
    e[1] = w[1] .* sqrt((sum((w_array .* (convolution(f, g) - T)) .^ 2)) / nt)
    e[2] = w[2] .* sqrt((sum((dg) .^ 2)) / nt)
    e[3] = w[3] .* sqrt((sum((ddg) .^ 2)) / nt)

    return sum(e)
end