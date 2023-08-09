
using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("ising.jl")

Base.copy(state::IsingState) =
    IsingState(state.matrix, state.sum_pair_products, state.base_length)
Base.length(state::IsingState) = 1
Base.iterate(state::IsingState) = iterate([state], 1)

# perform sampling - sanity check: log(Z) ≈ true value of around 33.3 for this example
pt = pigeons(target = IsingLogPotential(1.0, 5))

# # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# using Plots
# plot(pt.shared.tempering.communication_barriers.localbarrier)



# perform sampling - sanity check: log(Z) ≈ true value of around 33.3 for this example
pt = pigeons(
    target = IsingLogPotential(1.0, 5),
    record = [traces; index_process; record_default()],
)

samples = Chains(sample_array(pt), variable_names(pt))

# https://stackoverflow.com/questions/25028539/how-to-implement-an-iterator-in-julia
