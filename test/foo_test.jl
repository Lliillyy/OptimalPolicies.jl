# [test/foo_test.jl]
@testset "foo test" begin
    v = OptimalPolicies.foo(10,5)
    @test v[1] == 10
    @test v[2] == 5
    @test eltype(v) == Int
    v = OptimalPolicies.foo(10.0, 5)
    @test v[1] == 10
    @test v[2] == 5
    @test eltype(v) == Float64
end

@testset "U test" begin
    gam = 4
    x = -2.0
    U4 = OptimalPolicies.curried(gam)
    @test U4(x) == OptimalPolicies.U(gam, x)
    @test OptimalPolicies.U(gam, x) == 36.0
end