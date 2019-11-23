# MATLABDiffEq.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

MATLABDiffEq.jl is a common interface binding for the [MATLAB](https://www.mathworks.com/products/matlab.html) ordinary differential equation solvers. It uses the [MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) interop in order to
send the differential equation over to MATLAB and solve it. Note that this
package requires the differential equation function to be defined using
[ParameterizedFunctions.jl](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl).

Note that this package isn't for production use and is mostly just for benchmarking. For well-developed differential equation package, see
[DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## Installation

To install MATLABDiffEq.jl, use the following:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/MATLABDiffEq.jl")
```

## Using MATLABDiffEq.jl

MATLABDiffEq.jl is simply a solver on the DiffEq common interface, so for details see the [DifferentialEquations.jl documentation](https://juliadiffeq.github.io/DiffEqDocs.jl/dev/). However, there are three things to know:

1. The only options implemented are those for error calculations (`timeseries_error`), `saveat` and tolerances.
2. The input function must be defined by a `ParameterizedFunction`
3. The input function must not use parameters

Note that the algorithms are defined to have the same name as the MATLAB algorithms, but are not exported. Thus to use `ode45`, you would specify the algorithm as `MATLABDiffEq.ode45()`.

## Example

```julia
using DiffEqBase, MATLABDiffEq, ParameterizedFunctions


f = @ode_def LotkaVolterra begin
  dx = 1.5x - x*y
  dy = -3y + x*y
end

tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,MATLABDiffEq.ode45())
```

## Measuring Overhead

To measure the overhead of over the wrapper, note that the variables
from the session will be still stored in MATLAB after the computation
is done. Thus you can simply call the same ODE function and time it
directly. This is done by:

```julia
@time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
```

To be even more pedantic, you can play around in the actual MATLAB
session by using

```
MATLABDiffEq.show_msession()
```

## Overhead Amount

Generally, for long enough problems the overhead is minimal. Example:

```julia
using DiffEqBase, ParameterizedFunctions, MATLABDiffEq
f = @ode_def_bare RigidBodyBench begin
  dy1  = -2*y2*y3
  dy2  = 1.25*y1*y3
  dy3  = -0.5*y1*y2 + 0.25*sin(t)^2
end
prob = ODEProblem(f,[1.0;0.0;0.9],(0.0,100.0))
alg = MATLABDiffEq.ode45()
algstr = string(typeof(alg).name.name)
```

For this, we get the following:

```julia
julia> @time sol = solve(prob,alg);
  0.063918 seconds (38.84 k allocations: 1.556 MB)

julia> @time sol = solve(prob,alg);
  0.062600 seconds (38.84 k allocations: 1.556 MB)

julia> @time sol = solve(prob,alg);
  0.061438 seconds (38.84 k allocations: 1.556 MB)

julia> @time sol = solve(prob,alg);
  0.065460 seconds (38.84 k allocations: 1.556 MB)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.058249 seconds (11 allocations: 528 bytes)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.060367 seconds (11 allocations: 528 bytes)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.060171 seconds (11 allocations: 528 bytes)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.058928 seconds (11 allocations: 528 bytes)
```

## Benchmark

The following benchmarks demonstrate a **100x performance advantage for the
pure-Julia methods over the MATLAB ODE solvers** across a range of stiff and
non-stiff ODEs. These were ran with Julia 1.2 and MATLAB 2019B after verifying
negligible overhead on interop.

```julia
using ParameterizedFunctions, MATLABDiffEq, OrdinaryDiffEq, ODEInterface,
      ODEInterfaceDiffEq, Plots, Sundials
using DiffEqDevTools
using LinearAlgebra

## Non-Stiff Problem 1: Lotka-Volterra

f = @ode_def_bare LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1,3,1]
tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

setups = [Dict(:alg=>DP5())
          Dict(:alg=>dopri5())
          Dict(:alg=>BS5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>Vern6())
          Dict(:alg=>Vern7())
          Dict(:alg=>MATLABDiffEq.ode45())
          Dict(:alg=>MATLABDiffEq.ode113())
          Dict(:alg=>CVODE_Adams())
  ]
abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)
wp = WorkPrecisionSet(prob,abstols,reltols,setups;appxsol=test_sol,dense=false,save_everystep=false,numruns=100,maxiters=10000000,ttimeseries_errors=false,verbose=false)
plot(wp,title="Non-stiff 1: Lotka-Volterra")
savefig("benchmark1.png")
```

![benchmark1](https://user-images.githubusercontent.com/1814174/69478405-f9dac480-0dbf-11ea-8f8e-5572b86d2a97.png)

```julia
## Non-Stiff Problem 2: Rigid Body

f = @ode_def_bare RigidBodyBench begin
  dy1  = -2*y2*y3
  dy2  = 1.25*y1*y3
  dy3  = -0.5*y1*y2 + 0.25*sin(t)^2
end
prob = ODEProblem(f,[1.0;0.0;0.9],(0.0,100.0))
sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

setups = [Dict(:alg=>DP5())
          Dict(:alg=>dopri5())
          Dict(:alg=>BS5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>Vern6())
          Dict(:alg=>Vern7())
          Dict(:alg=>MATLABDiffEq.ode45())
          Dict(:alg=>MATLABDiffEq.ode113())
          Dict(:alg=>CVODE_Adams())
  ]
abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)
wp = WorkPrecisionSet(prob,abstols,reltols,setups;appxsol=test_sol,dense=false,save_everystep=false,numruns=100,maxiters=10000000,ttimeseries_errors=false,verbose=false)
plot(wp,title="Non-stiff 2: Rigid-Body")
savefig("benchmark2.png")
```

![benchmark2](https://user-images.githubusercontent.com/1814174/69478406-fe9f7880-0dbf-11ea-8f15-a0a0f3e8717f.png)


```julia
## Stiff Problem 1: ROBER

rober = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),[0.04,3e7,1e4])
sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);
setups = [Dict(:alg=>Rosenbrock23()),
          Dict(:alg=>TRBDF2()),
          Dict(:alg=>rodas()),
          Dict(:alg=>radau()),
          Dict(:alg=>RadauIIA5()),
          Dict(:alg=>MATLABDiffEq.ode23s()),
          Dict(:alg=>MATLABDiffEq.ode15s())
          ]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=100)
plot(wp,title="Stiff 1: ROBER")
savefig("benchmark3.png")
```

![benchmark3](https://user-images.githubusercontent.com/1814174/69478407-ffd0a580-0dbf-11ea-9311-909bd3938b42.png)


```julia
## Stiff Problem 2: HIRES

f = @ode_def Hires begin
  dy1 = -1.71*y1 + 0.43*y2 + 8.32*y3 + 0.0007
  dy2 = 1.71*y1 - 8.75*y2
  dy3 = -10.03*y3 + 0.43*y4 + 0.035*y5
  dy4 = 8.32*y2 + 1.71*y3 - 1.12*y4
  dy5 = -1.745*y5 + 0.43*y6 + 0.43*y7
  dy6 = -280.0*y6*y8 + 0.69*y4 + 1.71*y5 -
           0.43*y6 + 0.69*y7
  dy7 = 280.0*y6*y8 - 1.81*y7
  dy8 = -280.0*y6*y8 + 1.81*y7
end

u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057
prob = ODEProblem(f,u0,(0.0,321.8122))

sol = solve(prob,Rodas5(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);
setups = [Dict(:alg=>Rosenbrock23()),
          Dict(:alg=>TRBDF2()),
          Dict(:alg=>rodas()),
          Dict(:alg=>radau()),
          Dict(:alg=>RadauIIA5()),
          Dict(:alg=>MATLABDiffEq.ode23s()),
          Dict(:alg=>MATLABDiffEq.ode15s())
          ]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=100)
plot(wp,title="Stiff 2: Hires")
savefig("benchmark4.png")
```

![benchmark4](https://user-images.githubusercontent.com/1814174/69478409-019a6900-0dc0-11ea-9fce-c11a5e4de9a4.png)


Together this demonstrates that the algorithms like `ode45` and `ode15s` are
not competitive performance-wise against the Fortran and Julia methods.
This shows that being able to run MATLAB ODE algorithms with MATLAB functions
is cute, but does not really have a practical use due to MATLAB's lack of
performance (and its [pass by copy](https://www.mathworks.com/matlabcentral/answers/152-can-matlab-pass-by-reference) for functions).
