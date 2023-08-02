# [test/mcmc_test.jl]
using Distributions
using StatsPlots
using Pigeons


@testset "U" begin
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

@testset "chain" begin
    # Global settings
    temps = 2 .^ (0:3)
    iters = Int(1e5)

    # Run chains and store results
    mat = hcat(
        [
            OptimalPolicies.chain(OptimalPolicies.curried(gam), iters = iters) for
            gam in temps
        ]...,
    )
    OptimalPolicies.print_summary(mat, temps)

    for i in eachindex(temps)
        OptimalPolicies.plot_chain(mat[:, i], temps[i], dir = "chain_result/")
    end

    # 5 chains
    mat1 = OptimalPolicies.chains(temps = temps, iters = iters)
    OptimalPolicies.print_summary(mat1, temps)

    for i in eachindex(temps)
        OptimalPolicies.plot_chain(mat1[:, i], temps[i], dir = "chains_result/")
    end
end

@testset "chain_bb" begin
    # Global settings
    temps = 2 .^ (0:3)
    iters = Int(1e5)

    # Run chains and store results
    mat = hcat(
        [
            OptimalPolicies.chain_bb(OptimalPolicies.curried(gam), iters = iters) for
            gam in temps
        ]...,
    )
    OptimalPolicies.print_summary(mat, temps)

    for i in eachindex(temps)
        OptimalPolicies.plot_chain(mat[:, i], temps[i], dir = "chain_bb_result/")
    end
end
