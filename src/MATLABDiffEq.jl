__precompile__()

module MATLABDiffEq

using Reexport
@reexport using DiffEqBase
using MATLAB

abstract type MATLABAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
struct ode23 <: MATLABAlgorithm end
struct ode45 <: MATLABAlgorithm end
struct ode113 <: MATLABAlgorithm end
struct ode23s <: MATLABAlgorithm end
struct ode23t <: MATLABAlgorithm end
struct ode23tb <: MATLABAlgorithm end
struct ode15s <: MATLABAlgorithm end
struct ode15i <: MATLABAlgorithm end

function DiffEqBase.__solve{uType,tupType,isinplace,AlgType<:MATLABAlgorithm}(
    prob::DiffEqBase.AbstractODEProblem{uType,tupType,isinplace},
    alg::AlgType,timeseries=[],ts=[],ks=[];
    saveat=eltype(tupType)[],timeseries_errors=true,reltol = 1e-3, abstol = 1e-6,
    kwargs...)

    tType = eltype(tupType)

    if !(typeof(prob.f) <: DiffEqBase.AbstractParameterizedFunction)
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
    matstr = "f = @(t,internal_var___u) ["*matstr*"];"

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
        timeseries = Vector{uType}(undef,length(ts))
        for i=1:length(ts)
            timeseries[i] = timeseries_tmp[i,:]
        end
    else
        timeseries = timeseries_tmp
    end

    DiffEqBase.build_solution(prob,alg,ts,timeseries,
                 timeseries_errors = timeseries_errors)
end

end # module
