using MATLABDiffEq, Aqua, JET
using Test

@testset "Aqua" begin
    Aqua.test_all(MATLABDiffEq)
end

@testset "JET static analysis" begin
    rep = JET.report_package(MATLABDiffEq; target_modules = (MATLABDiffEq,))
    @test length(JET.get_reports(rep)) == 0
end
