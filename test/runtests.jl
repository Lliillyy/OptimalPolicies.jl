using OptimalPolicies
using Test

# Test scripts
@testset "OptimalPolicies test" begin
	@testset "foo_test" begin
		include("foo_test.jl")
	end
end