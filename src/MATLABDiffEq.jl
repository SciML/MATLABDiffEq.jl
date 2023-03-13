module MATLABDiffEq

using Reexport
@reexport using DiffEqBase
using MATLAB, ModelingToolkit

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
    prob::DiffEqBase.AbstractODEProblem{uType,tupType,isinplace},
    alg::AlgType,
    timeseries = [],
    ts = [],
    ks = [];
    saveat = eltype(tupType)[],
    timeseries_errors = true,
    reltol = 1e-3,
    abstol = 1e-6,
    callback = nothing,
    kwargs...,
) where {uType,tupType,isinplace,AlgType<:MATLABAlgorithm}

    tType = eltype(tupType)

    if prob.tspan[end] - prob.tspan[1] < tType(0)
        error("final time must be greater than starting time. Aborting.")
    end

    f = prob.f
    u0 = prob.u0

    if typeof(saveat) <: Number
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
        states(sys),
        parameters(sys),
        independent_variables(sys)[1],
        target = ModelingToolkit.MATLABTarget(),
    )

    # Send the variables
    put_variable(get_default_msession(), :tspan, tspan)
    put_variable(get_default_msession(), :u0, u0)
    put_variable(get_default_msession(), :internal_var___p, prob.p)
    put_variable(get_default_msession(), :reltol, reltol)
    put_variable(get_default_msession(), :abstol, abstol)

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
        for i = 1:length(ts)
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
        stats = stats,
    )
end

function buildDEStats(solverstats::Dict)

    destats = DiffEqBase.Stats(0)
    destats.nf = if (haskey(solverstats, "nfevals"))
        solverstats["nfevals"]
    else
        0
    end
    destats.nreject = if (haskey(solverstats, "nfailed"))
        solverstats["nfailed"]
    else
        0
    end
    destats.naccept = if (haskey(solverstats, "nsteps"))
        solverstats["nsteps"]
    else
        0
    end
    destats.nsolve = if (haskey(solverstats, "nsolves"))
        solverstats["nsolves"]
    else
        0
    end
    destats.njacs = if (haskey(solverstats, "npds"))
        solverstats["npds"]
    else
        0
    end
    destats.nw = if (haskey(solverstats, "ndecomps"))
        solverstats["ndecomps"]
    else
        0
    end
    destats
end

end # module
