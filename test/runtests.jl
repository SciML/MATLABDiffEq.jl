using Pkg
using Test
using DiffEqBase, MATLABDiffEq, ParameterizedFunctions

const GROUP = get(ENV, "GROUP", "Core")

# QA (Aqua + JET) runs in an isolated environment (test/qa) so its tooling deps
# never enter the main test target's resolve. On Julia < 1.11 the [sources] table
# is ignored, so develop the package by path to test the PR branch code.
function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

if GROUP == "Core" || GROUP == "All"
    # Interface tests - these test type validation without needing MATLAB runtime
    include("interface_tests.jl")

    # Static type-stability tests - these also run without MATLAB
    include("jet_tests.jl")

    # The following tests require MATLAB runtime to be available
    # They test the actual ODE solving functionality

    f = @ode_def_bare LotkaVolterra begin
        dx = a * x - b * x * y
        dy = -c * y + d * x * y
    end a b c d
    p = [1.5, 1, 3, 1]
    tspan = (0.0, 10.0)
    u0 = [1.0, 1.0]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, MATLABDiffEq.ode45())

    function lorenz(du, u, p, t)
        du[1] = 10.0(u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        return du[3] = u[1] * u[2] - (8 / 3) * u[3]
    end
    u0 = [1.0; 0.0; 0.0]
    tspan = (0.0, 100.0)
    prob = ODEProblem(lorenz, u0, tspan)
    sol = solve(prob, MATLABDiffEq.ode45())

    algs = [
        MATLABDiffEq.ode23
        MATLABDiffEq.ode45
        MATLABDiffEq.ode113
        MATLABDiffEq.ode23s
        MATLABDiffEq.ode23t
        MATLABDiffEq.ode23tb
        MATLABDiffEq.ode15s
    ]

    for alg in algs
        sol = solve(prob, alg())
    end
end

if GROUP == "QA"
    activate_qa_env()
    include("qa/qa.jl")
end
