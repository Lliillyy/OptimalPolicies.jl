using OptimalPolicies
using Test

# Test scripts
@testset "OptimalPolicies test" begin
    @testset "mcmc_test" begin
        include("mcmc_test.jl")
    end
end
