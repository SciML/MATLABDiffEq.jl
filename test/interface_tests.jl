# Interface compatibility tests for MATLABDiffEq.jl
# These tests verify type checking and interface compliance

using Test

@testset "Interface Compatibility" begin
    @testset "MATLAB type compatibility checks" begin
        # Test _is_matlab_compatible_eltype function
        @test MATLABDiffEq._is_matlab_compatible_eltype(Float64) == true
        @test MATLABDiffEq._is_matlab_compatible_eltype(Float32) == false
        @test MATLABDiffEq._is_matlab_compatible_eltype(Int64) == true
        @test MATLABDiffEq._is_matlab_compatible_eltype(Int32) == true
        @test MATLABDiffEq._is_matlab_compatible_eltype(Complex{Float64}) == true
        @test MATLABDiffEq._is_matlab_compatible_eltype(Complex{Float32}) == false
        @test MATLABDiffEq._is_matlab_compatible_eltype(BigFloat) == false
        @test MATLABDiffEq._is_matlab_compatible_eltype(BigInt) == false
        @test MATLABDiffEq._is_matlab_compatible_eltype(Rational{Int}) == false
    end

    @testset "_check_matlab_compatible validation" begin
        # Valid Float64 inputs should pass
        @test MATLABDiffEq._check_matlab_compatible([1.0, 2.0], (0.0, 1.0)) === nothing
        @test MATLABDiffEq._check_matlab_compatible(1.0, (0.0, 1.0)) === nothing
        @test MATLABDiffEq._check_matlab_compatible([1, 2, 3], (0, 10)) === nothing  # Integers OK

        # Complex Float64 should pass
        @test MATLABDiffEq._check_matlab_compatible([1.0 + 2.0im], (0.0, 1.0)) === nothing

        # BigFloat u0 should throw ArgumentError
        @test_throws ArgumentError MATLABDiffEq._check_matlab_compatible(
            BigFloat[1.0, 2.0], (0.0, 1.0)
        )

        # BigFloat tspan should throw ArgumentError
        @test_throws ArgumentError MATLABDiffEq._check_matlab_compatible(
            [1.0, 2.0], (BigFloat(0.0), BigFloat(1.0))
        )

        # Float32 should throw ArgumentError (MATLAB expects Float64)
        @test_throws ArgumentError MATLABDiffEq._check_matlab_compatible(
            Float32[1.0, 2.0], (0.0, 1.0)
        )
    end

    @testset "Error messages are helpful" begin
        # Test that error messages contain useful information
        try
            MATLABDiffEq._check_matlab_compatible(BigFloat[1.0], (0.0, 1.0))
            @test false  # Should not reach here
        catch e
            @test e isa ArgumentError
            @test occursin("BigFloat", e.msg)
            @test occursin("Float64", e.msg)
            @test occursin("MATLABDiffEq", e.msg)
        end

        try
            MATLABDiffEq._check_matlab_compatible([1.0], (BigFloat(0.0), BigFloat(1.0)))
            @test false  # Should not reach here
        catch e
            @test e isa ArgumentError
            @test occursin("tspan", lowercase(e.msg))
        end
    end

    @testset "buildDEStats is type-generic" begin
        # Test that buildDEStats works with different Dict types
        stats1 = Dict{String, Any}("nfevals" => 100, "nsteps" => 50)
        result1 = MATLABDiffEq.buildDEStats(stats1)
        @test result1.nf == 100
        @test result1.naccept == 50

        # Test with empty dict
        stats2 = Dict{String, Any}()
        result2 = MATLABDiffEq.buildDEStats(stats2)
        @test result2.nf == 0
        @test result2.naccept == 0

        # Test with all fields
        stats3 = Dict{String, Any}(
            "nfevals" => 200,
            "nfailed" => 10,
            "nsteps" => 190,
            "nsolves" => 100,
            "npds" => 20,
            "ndecomps" => 15
        )
        result3 = MATLABDiffEq.buildDEStats(stats3)
        @test result3.nf == 200
        @test result3.nreject == 10
        @test result3.naccept == 190
        @test result3.nsolve == 100
        @test result3.njacs == 20
        @test result3.nw == 15
    end

    @testset "Algorithm structs instantiation" begin
        # Test that all algorithm structs can be instantiated
        @test MATLABDiffEq.ode23() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode45() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode113() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode23s() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode23t() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode23tb() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode15s() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode15i() isa MATLABDiffEq.MATLABAlgorithm
    end
end
