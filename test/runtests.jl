using OptimalPolicies
using Test

# Test scripts
@testset "OptimalPolicies test" begin
    @testset "pt_test" begin
        include("pt_test.jl")
    end

    # @testset "mcmc_test" begin
    #     include("mcmc_test.jl")
    # end
end
