```@meta
CurrentModule = MATLABDiffEq
```

# MATLABDiffEq.jl

MATLABDiffEq provides SciML algorithm markers that solve ODE problems with MATLAB's ODE
routines. It is primarily intended for comparing implementations and migrating MATLAB
models; native Julia solvers are the recommended choice for production workloads.

## Basic Usage

`using MATLABDiffEq` brings the SciML common interface -- `ODEProblem`, `solve` and the
solution types -- into scope along with the algorithms, so pass a qualified MATLABDiffEq
algorithm to `solve`:

```julia
using MATLABDiffEq

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

prob = ODEProblem(lorenz!, [1.0, 0.0, 0.0], (0.0, 100.0))
sol = solve(prob, MATLABDiffEq.ode45(); reltol = 1.0e-6, abstol = 1.0e-8)
```

The MATLAB engine must be installed and available to MATLAB.jl. MATLABDiffEq accepts
`Float64`, standard integer, and `Complex{Float64}` state values. Callbacks are not
supported.

## Public API

Algorithm types are public but not exported, which avoids collisions with similarly named
algorithms from other solver packages. Use them through the `MATLABDiffEq` namespace. The
SciML common interface that `using MATLABDiffEq` reexports alongside them is listed on the
[API page](api.md).

```@docs
MATLABAlgorithm
ode23
ode45
ode113
ode23s
ode23t
ode23tb
ode15s
ode15i
```
