# These tests exercise the actual ODE solving functionality, which requires a
# usable MATLAB engine. They are gated on the engine being startable: on CI
# runners (and any machine without a licensed/working MATLAB) starting an
# MSession is impossible, so the engine-dependent smoke tests are skipped with a
# notice rather than failing. Where a working MATLAB engine is present the full
# solve path runs and must succeed.
using DiffEqBase, MATLABDiffEq, ParameterizedFunctions, MATLAB, Test

# `libmx_size == 0` means MATLAB.jl was built without a MATLAB installation
# (the build-time CI fallback), so the engine libraries were never loaded.
# Even when libraries are present, starting an MSession can still fail (no
# license, no display), so probe an actual session before running the solves.
function matlab_engine_available()
    MATLAB.libmx_size > 0 || return false
    return try
        MATLAB.MSession(0)
        true
    catch
        false
    end
end

@testset "MATLAB runtime solve" begin
    if !matlab_engine_available()
        @info "MATLAB engine not available on this machine; skipping engine-dependent solve tests."
    else
        f = @ode_def_bare LotkaVolterra begin
            dx = a * x - b * x * y
            dy = -c * y + d * x * y
        end a b c d
        p = [1.5, 1, 3, 1]
        tspan = (0.0, 10.0)
        u0 = [1.0, 1.0]
        prob = ODEProblem(f, u0, tspan, p)
        sol = solve(prob, MATLABDiffEq.ode45())
        @test length(sol.t) > 0

        function lorenz(du, u, p, t)
            du[1] = 10.0(u[2] - u[1])
            du[2] = u[1] * (28.0 - u[3]) - u[2]
            return du[3] = u[1] * u[2] - (8 / 3) * u[3]
        end
        u0 = [1.0; 0.0; 0.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(lorenz, u0, tspan)
        sol = solve(prob, MATLABDiffEq.ode45())
        @test length(sol.t) > 0

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
            @test length(sol.t) > 0
        end
    end
end
