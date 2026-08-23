module MATLABDiffEq

using MATLAB: eval_string, get_default_msession, get_mvariable, get_variable, jarray,
    jvector, put_variable
using ModelingToolkitBase: equations, independent_variables, modelingtoolkitize, parameters,
    unknowns
using PrecompileTools: @compile_workload, @setup_workload
using SciMLBase: AbstractODEAlgorithm, AbstractODEProblem, DEStats, build_solution
import SciMLBase: __solve
using SciMLPublic: @public
using Symbolics: MATLABTarget, build_function

# The SciML common interface that MATLABDiffEq reexports (see the `export` block below),
# so that `using MATLABDiffEq` on its own is enough to build an ODE problem, solve it,
# and inspect the result -- the workflow the README and docs/src/index.md document. Every
# name stays owned and documented upstream.
using SciMLBase: EnsembleAnalysis, EnsembleDistributed, EnsembleProblem, EnsembleSerial,
    EnsembleSolution, EnsembleSplitThreads, EnsembleSummary, EnsembleThreads,
    NullParameters, ODEFunction, ODEProblem, ODESolution, ReturnCode, remake, solve,
    successful_retcode

# Reexported SciML common interface; approved via `reexports_allow` in test/qa/qa.jl.
# `DEStats` is imported above.
export DEStats, EnsembleAnalysis, EnsembleDistributed, EnsembleProblem, EnsembleSerial,
    EnsembleSolution, EnsembleSplitThreads, EnsembleSummary, EnsembleThreads,
    NullParameters, ODEFunction, ODEProblem, ODESolution, ReturnCode, remake, solve,
    successful_retcode

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

"""
    MATLABAlgorithm

Abstract supertype for the ODE algorithm markers provided by MATLABDiffEq.

Concrete subtypes select the MATLAB routine used by `SciMLBase.solve`. Downstream code may
dispatch on `MATLABAlgorithm` to identify MATLAB-backed ODE algorithms. External subtyping
is not supported: the bridge derives the MATLAB routine name from the concrete Julia type,
so MATLABDiffEq must own and test every subtype.

# Examples

```jldoctest
julia> MATLABDiffEq.ode45() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
abstract type MATLABAlgorithm <: AbstractODEAlgorithm end

"""
    ode23()

Select MATLAB's low-order explicit Runge-Kutta ODE solver.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode23() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode23 <: MATLABAlgorithm end

"""
    ode45()

Select MATLAB's variable-step explicit Runge-Kutta `(4, 5)` ODE solver.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode45() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode45 <: MATLABAlgorithm end

"""
    ode113()

Select MATLAB's variable-order Adams-Bashforth-Moulton ODE solver.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode113() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode113 <: MATLABAlgorithm end

"""
    ode23s()

Select MATLAB's low-order Rosenbrock solver for stiff ODEs.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode23s() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode23s <: MATLABAlgorithm end

"""
    ode23t()

Select MATLAB's trapezoidal-rule solver for moderately stiff ODEs.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode23t() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode23t <: MATLABAlgorithm end

"""
    ode23tb()

Select MATLAB's TR-BDF2 solver for stiff ODEs.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode23tb() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode23tb <: MATLABAlgorithm end

"""
    ode15s()

Select MATLAB's variable-order BDF/NDF solver for stiff ODEs.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode15s() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode15s <: MATLABAlgorithm end

"""
    ode15i()

Select MATLAB's solver for fully implicit differential equations.

Pass the resulting marker to `SciMLBase.solve`; configure tolerances and saved output with
the usual SciML solve keywords.

# Examples

```jldoctest
julia> MATLABDiffEq.ode15i() isa MATLABDiffEq.MATLABAlgorithm
true
```
"""
struct ode15i <: MATLABAlgorithm end

@public MATLABAlgorithm, ode23, ode45, ode113, ode23s, ode23t, ode23tb, ode15s, ode15i

function __solve(
        prob::AbstractODEProblem{uType, tupType, isinplace},
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

    matstr = build_function(
        map(x -> x.rhs, equations(sys)),
        unknowns(sys),
        parameters(sys),
        independent_variables(sys)[1],
        target = MATLABTarget()
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

    return build_solution(
        prob,
        alg,
        ts,
        timeseries,
        timeseries_errors = timeseries_errors,
        stats = stats
    )
end

"""
    buildDEStats(solverstats::Dict{String, <:Any}) -> DEStats

Convert MATLAB ODE solver statistics dictionary to SciMLBase.DEStats.

This internal bridge helper maps the counters returned by a MATLAB ODE solver to the
corresponding fields in `DEStats`. Missing keys default to zero.

# Arguments

- `solverstats`: Solver statistics keyed by the MATLAB field names `"nfevals"`,
  `"nfailed"`, `"nsteps"`, `"nsolves"`, `"npds"`, and `"ndecomps"`.

# Returns

- A `DEStats` value populated from the available MATLAB counters.
"""
function buildDEStats(solverstats::Dict{String, <:Any})::DEStats
    destats = DEStats(0)
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
