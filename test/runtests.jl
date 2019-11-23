using DiffEqBase, MATLABDiffEq, ParameterizedFunctions, Test

f = @ode_def_bare LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1,3,1]
tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan,p)

sol = solve(prob,MATLABDiffEq.ode45())

algs = [MATLABDiffEq.ode23
        MATLABDiffEq.ode45
        MATLABDiffEq.ode113
        MATLABDiffEq.ode23s
        MATLABDiffEq.ode23t
        MATLABDiffEq.ode23tb
        MATLABDiffEq.ode15s]

for alg in algs
  sol = solve(prob,alg())
end
