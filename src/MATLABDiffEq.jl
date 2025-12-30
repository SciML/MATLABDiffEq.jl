module MATLABDiffEq

using Reexport
@reexport using DiffEqBase
using MATLAB, ModelingToolkit
using PrecompileTools

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
        reltol = 1e-3,
        abstol = 1e-6,
        callback = nothing,
        kwargs...
) where {uType, tupType, isinplace, AlgType <: MATLABAlgorithm}
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

    DiffEqBase.build_solution(
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
    destats
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
    end
end

end # module
