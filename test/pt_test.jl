# [test/pt_test.jl]
using Distributions
using StatsPlots
# [test/pt_test.jl]

@testset "pigeons" begin
    include("../src/run_pt.jl")
end
