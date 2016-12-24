using DiffEqBase, MATLABDiffEq, ParameterizedFunctions
using Base.Test


f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=1.5 b=1 c=3 d=1

tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan)

algs = [MATLABDiffEq.ode23
MATLABDiffEq.ode45
MATLABDiffEq.ode113
MATLABDiffEq.ode23s
MATLABDiffEq.ode23t
MATLABDiffEq.ode23tb
MATLABDiffEq.ode15s
MATLABDiffEq.ode15i]

for alg in algs
  sol = solve(prob,alg())
end
