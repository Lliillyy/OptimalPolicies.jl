# [test/foo_test.jl]
using Distributions, StatsPlots

@testset "U test" begin
    gam = 4
    x = -2.0
    U4 = OptimalPolicies.curried(gam)

    @test U4(x) == OptimalPolicies.U(gam, x)
    @test OptimalPolicies.U(gam, x) == 36.0

    # Plot potential function
    p1 = plot(U4, -2, 2, title = "Potential function, U(x)", legend = false, fmt = :png)

    # Plot density function
    p2 = plot(
        x -> exp(-U4(x)),
        -2,
        2,
        title = "Unnormalized density function, exp(-U(x))",
        legend = false,
        fmt = :png,
    )

    # Combine the plots
    plot(p1, p2, layout = (2, 1))
    savefig("plots.png")  # Save the combined plots as plots.png
end

@testset "chain test" begin
    # Global settings
    temps = 2 .^ (0:3)
    iters = Int(1e5)

    # Run chains and store results
    mat = hcat([OptimalPolicies.chain(OptimalPolicies.curried(gam)) for gam in temps]...)
    colnames(mat) = ["gamma=$gam" for gam in temps]
    OptimalPolicies.print_summary(mat)

    # 5 chains
    mat1 = OptimalPolicies.chains(temps = temps, iters = iters)
    OptimalPolicies.print_summary(mat1)
end
