using SafeTestsets

@safetestset "Aqua" begin
    using MATLABDiffEq, Aqua
    Aqua.test_all(MATLABDiffEq)
end

@safetestset "JET static analysis" begin
    using MATLABDiffEq, JET, Test
    rep = JET.report_package(MATLABDiffEq; target_modules = (MATLABDiffEq,))
    @test length(JET.get_reports(rep)) == 0
end
