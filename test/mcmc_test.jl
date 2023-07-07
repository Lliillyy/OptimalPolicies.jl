# [test/foo_test.jl]

@testset "U test" begin
    gam = 4
    x = -2.0
    U4 = OptimalPolicies.curried(gam)
    @test U4(x) == OptimalPolicies.U(gam, x)
    @test OptimalPolicies.U(gam, x) == 36.0
end