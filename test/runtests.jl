using DiffEqBase, MATLABDiffEq, ParameterizedFunctions, Test
using MATLAB

# Interface tests - these test type validation without needing MATLAB runtime
include("interface_tests.jl")

# JET static analysis tests - these also run without MATLAB
include("jet_tests.jl")

# The following tests require a MATLAB installation to be available at build
# time. MATLAB.jl exposes `libmx_size`, which is 0 when no MATLAB was found
# (e.g. CI without a MATLAB license). In that case the function pointers into
# the MATLAB engine are never bound, so any solve would throw an UndefRefError
# when it tries to open an engine session. Gate the runtime tests on a real
# MATLAB being present rather than letting the suite error out.
const HAS_MATLAB = MATLAB.libmx_size > 0

if HAS_MATLAB
    @testset "MATLAB runtime solving" begin
        f = @ode_def_bare LotkaVolterra begin
            dx = a * x - b * x * y
            dy = -c * y + d * x * y
        end a b c d
        p = [1.5, 1, 3, 1]
        tspan = (0.0, 10.0)
        u0 = [1.0, 1.0]
        prob = ODEProblem(f, u0, tspan, p)
        sol = solve(prob, MATLABDiffEq.ode45())
        @test length(sol.t) > 1
        @test length(sol.u) == length(sol.t)
        @test sol.t[1] == tspan[1]
        @test sol.t[end] == tspan[2]

        function lorenz(du, u, p, t)
            du[1] = 10.0(u[2] - u[1])
            du[2] = u[1] * (28.0 - u[3]) - u[2]
            return du[3] = u[1] * u[2] - (8 / 3) * u[3]
        end
        u0 = [1.0; 0.0; 0.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(lorenz, u0, tspan)
        sol = solve(prob, MATLABDiffEq.ode45())
        @test length(sol.t) > 1
        @test length(sol.u[1]) == 3

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
            @test length(sol.t) > 1
        end
    end
else
    @info "No MATLAB installation detected (MATLAB.libmx_size == 0); " *
        "skipping the MATLAB runtime solving tests. The type-validation and " *
        "static-analysis tests above still run and must pass."
end
