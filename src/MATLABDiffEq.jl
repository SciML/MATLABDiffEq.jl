module MATLABDiffEq

using DiffEqBase, MATLAB
import DiffEqBase: solve

abstract MATLABAlgorithm <: AbstractODEAlgorithm
immutable ode23 <: MATLABAlgorithm end
immutable ode45 <: MATLABAlgorithm end
immutable ode113 <: MATLABAlgorithm end
immutable ode23s <: MATLABAlgorithm end
immutable ode23t <: MATLABAlgorithm end
immutable ode23tb <: MATLABAlgorithm end
immutable ode15s <: MATLABAlgorithm end
immutable ode15i <: MATLABAlgorithm end

function solve{uType,tType,isinplace,AlgType<:MATLABAlgorithm,F}(
    prob::AbstractODEProblem{uType,tType,isinplace,F},
    alg::AlgType,timeseries=[],ts=[],ks=[];
    saveat=tType[],timeseries_errors=true,reltol = 1e-3, abstol = 1e-6,
    kwargs...)

    if !(typeof(prob.f) <: AbstractParameterizedFunction)
      error("Functions must be defined via ParameterizedFunctions.jl to work with this package.")
    end

    if prob.tspan[end]-prob.tspan[1]<tType(0)
        error("final time must be greater than starting time. Aborting.")
    end
    
    f = prob.f
    u0 = prob.u0

    tspan = sort(unique([prob.tspan[1];saveat;prob.tspan[2]]))

    sizeu = size(prob.u0)

    if uType <: AbstractArray
        u0 = vec(prob.u0)
    else
        u0 = prob.u0
    end

    strs = [string(f.fex.args[2i].args[2]) for i in 1:length(f.syms)]
    matstr = ""
    for i in 1:length(strs)
      matstr *= strs[i]
      i < length(strs) && (matstr *= "; ")
    end
    matstr = replace(matstr,"[","(")
    matstr = replace(matstr,"]",")")
    matstr = "f = @(t,u) ["*matstr*"];"

    # Send the function over
    eval_string(matstr)
    # Send the variables
    put_variable(get_default_msession(),:tspan,tspan)
    put_variable(get_default_msession(),:u0,u0)
    put_variable(get_default_msession(),:reltol,reltol)
    put_variable(get_default_msession(),:abstol,abstol)

    eval_string("options = odeset('RelTol',reltol,'AbsTol',abstol);")
    algstr = string(typeof(alg).name.name)
    #algstr = replace(string(typeof(alg)),"MATLABDiffEq.","")
    eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
    ts = jvector(get_mvariable(:t))
    timeseries_tmp = jarray(get_mvariable(:u))

    # Reshape the result if needed
    if uType <: AbstractArray
        timeseries = Vector{uType}(0)
        for i=1:length(ts)
            push!(timeseries,timeseries_tmp[i,:])
        end
    else
        timeseries = timeseries_tmp
    end

    build_solution(prob,alg,ts,timeseries,
                 timeseries_errors = timeseries_errors)
end

end # module
