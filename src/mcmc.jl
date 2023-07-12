# [src/mcmc.jl]

using Distributions, StatsPlots
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
    chain(target, tune=0.1, init=1, iters=1e5)

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

"""
    print_summary(mat)

Computes and displays the statistics of the given matrix `mat`
"""
function print_summary(mat)
    summary_stats = Dict(
        "mean" => mean(mat, dims = 1),
        "std" => std(mat, dims = 1),
        "minimum" => minimum(mat, dims = 1),
        "maximum" => maximum(mat, dims = 1),
    )

    println(summary_stats)
end

"""
    chains(pot = U, tune = 0.1, init = 1.)

Generates 5 chains at once
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
    colnames(xmat) = ["gamma=$gam" for gam in temps]
    xmat
end
