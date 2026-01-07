module MATLABDiffEq

using Reexport
@reexport using DiffEqBase
using MATLAB, ModelingToolkit
using PrecompileTools

# MATLAB only supports Float64 arrays. Check if a type is MATLAB-compatible.
# Note: We specifically accept standard Julia integer types that MATLAB can convert,
# but NOT BigInt since MATLAB doesn't support arbitrary precision integers.
_is_matlab_compatible_eltype(::Type{Float64}) = true
_is_matlab_compatible_eltype(::Type{<:Union{Int8, Int16, Int32, Int64, Int128}}) = true
_is_matlab_compatible_eltype(::Type{<:Union{UInt8, UInt16, UInt32, UInt64, UInt128}}) = true
_is_matlab_compatible_eltype(::Type{<:Complex{Float64}}) = true
_is_matlab_compatible_eltype(::Type) = false

function _check_matlab_compatible(u0, tspan)::Nothing
    T = eltype(u0)
    if !_is_matlab_compatible_eltype(T)
        throw(
            ArgumentError(
                "MATLABDiffEq.jl requires Float64-compatible element types. " *
                    "Got eltype(u0) = $T. MATLAB does not support arbitrary precision " *
                    "(BigFloat) or GPU arrays (JLArrays, CuArrays). Please convert your " *
                    "initial conditions to Float64: u0 = Float64.(u0)"
            )
        )
    end
    tT = eltype(tspan)
    if !_is_matlab_compatible_eltype(tT)
        throw(
            ArgumentError(
                "MATLABDiffEq.jl requires Float64-compatible time span types. " *
                    "Got eltype(tspan) = $tT. MATLAB does not support arbitrary precision " *
                    "(BigFloat). Please use Float64 for tspan: tspan = Float64.(tspan)"
            )
        )
    end
    # Check that the array type itself is a standard Julia array
    if !(u0 isa Array || u0 isa Number)
        @warn "MATLABDiffEq.jl works best with standard Julia Arrays. " *
            "Got $(typeof(u0)). The array will be converted to a standard Array " *
            "before being sent to MATLAB."
    end
    return nothing
end

# Handle ModelingToolkit API changes: states -> unknowns
if isdefined(ModelingToolkit, :unknowns)
    const mtk_states = ModelingToolkit.unknowns
else
    const mtk_states = ModelingToolkit.states
end

abstract type MATLABAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
struct ode23 <: MATLABAlgorithm end
struct ode45 <: MATLABAlgorithm end
struct ode113 <: MATLABAlgorithm end
struct ode23s <: MATLABAlgorithm end
struct ode23t <: MATLABAlgorithm end
struct ode23tb <: MATLABAlgorithm end
struct ode15s <: MATLABAlgorithm end
struct ode15i <: MATLABAlgorithm end

function DiffEqBase.__solve(
        prob::DiffEqBase.AbstractODEProblem{uType, tupType, isinplace},
        alg::AlgType,
        timeseries = [],
        ts = [],
        ks = [];
        saveat = eltype(tupType)[],
        timeseries_errors = true,
        reltol = 1.0e-3,
        abstol = 1.0e-6,
        callback = nothing,
        kwargs...
    ) where {uType, tupType, isinplace, AlgType <: MATLABAlgorithm}
    # Validate that input types are MATLAB-compatible
    _check_matlab_compatible(prob.u0, prob.tspan)

    tType = eltype(tupType)

    if prob.tspan[end] - prob.tspan[1] < tType(0)
        error("final time must be greater than starting time. Aborting.")
    end

    f = prob.f
    u0 = prob.u0

    if saveat isa Number
        tspan = Array(prob.tspan[1]:saveat:prob.tspan[2])
        tspan = sort(unique([prob.tspan[1]; tspan; prob.tspan[2]]))
    else
        tspan = sort(unique([prob.tspan[1]; saveat; prob.tspan[2]]))
    end

    sizeu = size(prob.u0)

    if uType <: AbstractArray
        u0 = vec(prob.u0)
    else
        u0 = prob.u0
    end

    callback !== nothing && error("Callbacks are not supported in MATLABDiffEq.jl")

    sys = modelingtoolkitize(prob)

    matstr = ModelingToolkit.build_function(
        map(x -> x.rhs, equations(sys)),
        mtk_states(sys),
        parameters(sys),
        independent_variables(sys)[1],
        target = ModelingToolkit.MATLABTarget()
    )

    # Send the variables
    put_variable(get_default_msession(), :tspan, tspan)
    put_variable(get_default_msession(), :u0, u0)
    put_variable(get_default_msession(), :internal_var___p, prob.p)
    put_variable(get_default_msession(), :reltol, reltol)
    put_variable(get_default_msession(), :abstol, abstol)

    # Define the ifelse helper function in MATLAB
    # MATLAB doesn't have a built-in ifelse, so we need to define one
    # This allows symbolic ifelse expressions to be evaluated properly
    eval_string("ifelse = @(cond, a, b) cond .* a + ~cond .* b;")

    # Send the function over
    eval_string(matstr)

    eval_string("options = odeset('RelTol',reltol,'AbsTol',abstol);")
    algstr = string(typeof(alg).name.name)
    eval_string("mxsol = $(algstr)(diffeqf,tspan,u0,options);")
    eval_string("mxsolstats = struct(mxsol.stats);")
    solstats = get_variable(:mxsolstats)
    eval_string("t = mxsol.x;")
    ts = jvector(get_mvariable(:t))
    eval_string("u = mxsol.y';")
    timeseries_tmp = jarray(get_mvariable(:u))

    # Reshape the result if needed
    if uType <: AbstractArray
        timeseries = Vector{uType}(undef, length(ts))
        for i in 1:length(ts)
            timeseries[i] = @view timeseries_tmp[i, :]
        end
    else
        timeseries = timeseries_tmp
    end

    stats = buildDEStats(solstats)

    return DiffEqBase.build_solution(
        prob,
        alg,
        ts,
        timeseries,
        timeseries_errors = timeseries_errors,
        stats = stats
    )
end

"""
    buildDEStats(solverstats::Dict{String, <:Any}) -> DiffEqBase.Stats

Convert MATLAB ODE solver statistics dictionary to DiffEqBase.Stats.

The function extracts statistics from the MATLAB solver output and maps them
to the corresponding fields in DiffEqBase.Stats. Missing keys default to 0.
"""
function buildDEStats(solverstats::Dict{String, <:Any})::DiffEqBase.Stats
    destats = DiffEqBase.Stats(0)
    destats.nf = Int(get(solverstats, "nfevals", 0))
    destats.nreject = Int(get(solverstats, "nfailed", 0))
    destats.naccept = Int(get(solverstats, "nsteps", 0))
    destats.nsolve = Int(get(solverstats, "nsolves", 0))
    destats.njacs = Int(get(solverstats, "npds", 0))
    destats.nw = Int(get(solverstats, "ndecomps", 0))
    return destats
end

@setup_workload begin
    # Precompile algorithm struct instantiations and buildDEStats
    @compile_workload begin
        # Instantiate algorithm structs - this precompiles their constructors
        _ = ode23()
        _ = ode45()
        _ = ode113()
        _ = ode23s()
        _ = ode23t()
        _ = ode23tb()
        _ = ode15s()
        _ = ode15i()

        # Precompile buildDEStats with typical MATLAB stats dictionaries
        test_stats = Dict{String, Any}(
            "nfevals" => 100,
            "nfailed" => 5,
            "nsteps" => 95,
            "nsolves" => 50,
            "npds" => 10,
            "ndecomps" => 8
        )
        _ = buildDEStats(test_stats)

        # Also precompile with missing keys (common case)
        _ = buildDEStats(Dict{String, Any}())

        # Precompile type compatibility checks
        _ = _is_matlab_compatible_eltype(Float64)
        _ = _is_matlab_compatible_eltype(Int64)
        _ = _is_matlab_compatible_eltype(UInt64)
        _ = _is_matlab_compatible_eltype(Complex{Float64})
        _ = _is_matlab_compatible_eltype(BigFloat)
        _ = _is_matlab_compatible_eltype(BigInt)
        _ = _check_matlab_compatible([1.0, 2.0], (0.0, 1.0))
        _ = _check_matlab_compatible(1.0, (0.0, 1.0))
    end
end

end # module
