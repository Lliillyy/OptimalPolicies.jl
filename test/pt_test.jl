# [test/pt_test.jl]
using Distributions
using StatsPlots
# [test/pt_test.jl]

@testset "pigeons" begin
    include("../src/pt_v2.jl")
end
