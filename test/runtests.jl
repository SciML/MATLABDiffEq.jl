using DiffEqBase, MATLABDiffEq, ParameterizedFunctions, Test

# Interface tests - these test type validation without needing MATLAB runtime
include("interface_tests.jl")

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
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, MATLABDiffEq.ode45())

algs = [MATLABDiffEq.ode23
        MATLABDiffEq.ode45
        MATLABDiffEq.ode113
        MATLABDiffEq.ode23s
        MATLABDiffEq.ode23t
        MATLABDiffEq.ode23tb
        MATLABDiffEq.ode15s]

for alg in algs
    sol = solve(prob, alg())
end
