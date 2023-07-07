# [src/OptimalPolicies.jl]
module OptimalPolicies

# using StaticArrays, Pigeons
using Distributions, StatsPlots, Statistics
include("mcmc.jl")

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

# Global settings
temps = 2 .^(0:3)
iters = Int(1e5)

# Run chains and store results
mat = hcat([chain(curried(gam)) for gam in temps]...)
colnames(mat) = ["gamma=$gam" for gam in temps]
print_summary(mat)

# 5 chains
mat1 = chains()
print_summary(mat1)

end # module OptimalPolicies
