# MATLABDiffEq.jl

MATLABDiffEq.jl is a common interface binding for the [MATLAB](https://www.mathworks.com/products/matlab.html)
ordinary differential equation solvers. It uses the [MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) interop in order to
send the differential equation over to MATLAB and solve it. Note that this
package requires the differential equation function to be defined using
[ParameterizedFunctions.jl](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl).

Note that this package isn't for production use and is mostly just for benchmarking. For well-developed differential equation package, see
[DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## Installation

To install MATLABDiffEq.jl, use the following:

```julia
Pkg.clone("")
```

## Using MATLABDiffEq.jl

MATLABDiffEq.jl is simply a solver on the DiffEq common interface, so for details see the [DifferentialEquations.jl documentation](https://juliadiffeq.github.io/DiffEqDocs.jl/latest/). However, there are two things to know:

1. The only options implemented are those for error calculations (`timeseries_error`), `saveat` and tolerances.
2. The input function must be defined by a `ParameterizedFunction`

Note that the algorithms are defined to have the same name as the MATLAB algorithms, but are not exported. Thus to use `ode45`, you would specify the algorithm as `MATLABDiffEq.ode45()`.

## Example

```julia
using DiffEqBase, MATLABDiffEq, ParameterizedFunctions


f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=1.5 b=1 c=3 d=1

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
f = @ode_def_bare RigidBodyBench begin
  dy1  = I1*y2*y3
  dy2  = I2*y1*y3
  dy3  = I3*y1*y2 + 0.25*sin(t)^2
end I1=-2 I2=1.25 I3=-.5
prob = ODEProblem(f,[1.0;0.0;0.9],(0.0,300.0))
```

For this, we get the following:

```julia
julia> @time sol = solve(prob,MATLABDiffEq.ode45());
  0.467800 seconds (329.98 k allocations: 12.903 MB)

julia> @time sol = solve(prob,MATLABDiffEq.ode45());
  0.469933 seconds (329.98 k allocations: 12.903 MB, 1.06% gc time)

julia> @time sol = solve(prob,MATLABDiffEq.ode45());
  0.468583 seconds (329.98 k allocations: 12.903 MB)

julia> @time sol = solve(prob,MATLABDiffEq.ode45());
  0.475089 seconds (329.98 k allocations: 12.903 MB, 0.99% gc time)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.451407 seconds (11 allocations: 528 bytes)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.452138 seconds (11 allocations: 528 bytes)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.445815 seconds (11 allocations: 528 bytes)

julia> @time MATLABDiffEq.eval_string("[t,u] = $(algstr)(f,tspan,u0,options);")
  0.456558 seconds (11 allocations: 528 bytes)
```
