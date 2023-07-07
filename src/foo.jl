# [src/foo.jl]
# """
#     foo(x, y)

# Creates a 2-element static array from the scalars `x` and `y`.
# """
# function foo(x::Number, y::Number)
#     SA[x, y]
# end

# using Pkg
# Pkg.add("Distributions")
# Pkg.add("Plots")
# Pkg.add("StatsPlots")
# Pkg.add("StatsBase")
using Distributions, StatsPlots
gr(fmt=:png)

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

# Plot potential and density functions
U4 = curried(4)

p1 = plot(
    U4, -2, 2,
    title = "Potential function, U(x)",
    legend = false,
    fmt = :png
)

# Plot density function
p2 = plot(
    x -> exp(-U4(x)), -2, 2,
    title = "Unnormalized density function, exp(-U(x))",
    legend = false,
    fmt = :png
)

# Combine the plots
plot(p1, p2, layout = (2, 1))
savefig("plots.png")  # Save the combined plots as plots.png

"""
    chain(target, tune=0.1, init=1, iters=1e5)

The Metropolis algorithm targeting a particular potential function `target`
"""
function chain(target, tune=0.1, init=1)
    x = init
    xvec = Vector{Float64}(undef, iters)
    for i in 1:iters
        can = x + randn() * tune
        logA = target(x) - target(can)
        if log(rand()) < logA
            x = can
        end
        xvec[i] = x
    end
    xvec
end

# Global settings
temps = 2 .^(0:3)
iters = Int(1e5)

# Run chains and store results
mat = hcat([chain(curried(gam)) for gam in temps]...)
colnames(mat) = ["gamma=$gam" for gam in temps]

# Compute summary statistics
summary_stats = Dict(
    "mean" => mean(mat, dims = 1),
    "std" => std(mat, dims = 1),
    "minimum" => minimum(mat, dims = 1),
    "maximum" => maximum(mat, dims = 1)
)

println(summary_stats)


# Generate 5 chains at once
function chains(pot = U, tune = 0.1, init = 1)
    x = fill(init, length(temps))
    xmat = zeros(iters, length(temps))
    for i in 1:iters
        can = x + randn(length(temps)) * tune
        logA = [pot(gam, x[j]) - pot(gam, can[j]) 
                for (j, gam) in enumerate(temps)]
        accept = log.(rand(length(temps))) .< logA
        x[accept] .= can[accept]
        xmat[i, :] = x
    end
    colnames(xmat) = ["gamma=$gam" for gam in temps]
    xmat
end

mat = chains()
summary_stats = Dict(
    "mean" => mean(mat, dims = 1),
    "std" => std(mat, dims = 1),
    "minimum" => minimum(mat, dims = 1),
    "maximum" => maximum(mat, dims = 1)
)

println(summary_stats)