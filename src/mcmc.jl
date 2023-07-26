# [src/mcmc.jl]

using Distributions, StatsPlots, StatsBase
gr(fmt = :png)

"""
    U(gam::Number, x::Number)

The potential function for some given potential barrier height `gam`
"""
function U(gam::Number, x::Number)
    gam * (x^2 - 1)^2
end

"""
    curried(gam::Number)

Returns a potential function for the given potential barrier height `gam`
"""
function curried(gam::Number)
    println("Returning a function for gamma =", gam)
    x -> U(gam, x)
end

"""
    density(gam::Number, x::Number)

The density function for some given potential barrier height `gam`
"""
function density(gam::Number, x::Number)
    exp(-U(gam, x))
end

"""
    chain(target; tune = 0.1, init = 1.0, iters = Int(1e3))

The Metropolis algorithm targeting a particular potential function `target`
"""
function chain(target; tune = 0.1, init = 1.0, iters = Int(1e3))
    x = init
    xvec = Vector{Float64}(undef, iters)
    for i = 1:iters
        can = x + randn() * tune
        logA = target(x) - target(can)
        if log(rand()) < logA
            x = can
        end
        xvec[i] = x
    end
    xvec
end

function chain_bb(target; init = 1, iters = Int(1e3))
    alpha = 0.2    # Shape parameter of the Beta distribution
    beta = 0.25     # Shape parameter of the Beta distribution
    n = 10         # Number of trials
    bb_dist = BetaBinomial(n, alpha, beta)

    # # Define the probability mass function (PMF)
    # pmf = [0.2, 0.4, 0.3, 0.1]

    # # Define the discrete distribution
    # dist = Categorical(pmf)

    x = init
    xvec = Vector{Int}(undef, iters)
    for i = 1:iters
        can = rand(bb_dist)
        logA = target(can) - target(x)
        if log(rand()) < logA
            x = can
        end
        xvec[i] = x
    end
    xvec
end

"""
    print_summary(mat, temps)

Computes and displays the statistics of the given matrix `mat`
"""
function print_summary(mat, temps)
    summary_stats = Dict(
        "gamma" => temps,
        "mean" => mean(mat, dims = 1),
        "std" => std(mat, dims = 1),
        "minimum" => minimum(mat, dims = 1),
        "maximum" => maximum(mat, dims = 1),
    )

    println(summary_stats)
end


"""
    plot_chain(xvec, gam; dir = "result/")

Plot the value, autocorrelation, and density of a chain `xvec` given the value
of `gam`, store the plot in the directory `dir`
"""
function plot_chain(xvec, gam; dir = "result/")
    # Compute autocorrelation
    acf_values = autocor(xvec)

    # Plot ACF
    p1 = plot(
        xvec,
        xlabel = "Iteration",
        ylabel = "Value",
        label = false,
        title = "gamma=$gam",
    )

    p2 = plot(
        acf_values,
        marker = :circle,
        markersize = 4,
        ylims = (0, 1),
        xlabel = "Lag",
        ylabel = "ACF",
        label = false,
        title = "gamma=$gam",
    )

    # Plot density histogram
    p3 = histogram(
        xvec,
        xlabel = "Value",
        ylabel = "Density",
        label = false,
        title = "gamma=$gam",
        density = true,
        bins = 100,
    )

    # Combine the plots
    plot(p1, p2, p3, layout = (1, 3), size = (1500, 300))

    if !isdir(dir)
        mkdir(dir)
    end
    savefig(dir * "gamma=$gam.png")
end

"""
    chains(; pot = U, tune = 0.1, init = 1.0, iters = Int(1e3), temps = [1, 2])

Generates chains based on the `temps`
"""
function chains(; pot = U, tune = 0.1, init = 1.0, iters = Int(1e3), temps = [1, 2])
    x = fill(init, length(temps))
    xmat = zeros(iters, length(temps))
    for i = 1:iters
        can = x + randn(length(temps)) * tune
        logA = [pot(gam, x[j]) - pot(gam, can[j]) for (j, gam) in enumerate(temps)]
        accept = log.(rand(length(temps))) .< logA
        x[accept] .= can[accept]
        # swap states
        swap = sample(1:length(temps), 2)
        logA =
            pot(temps[swap[1]], x[swap[1]]) + pot(temps[swap[2]], x[swap[2]]) -
            pot(temps[swap[1]], x[swap[2]]) - pot(temps[swap[2]], x[swap[1]])
        if log(rand()) < logA
            x[swap] = reverse(x[swap])
        end
        # end swapping
        xmat[i, :] = x
    end
    xmat
end
