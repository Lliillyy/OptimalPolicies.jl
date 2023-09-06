
using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("ising.jl")


# perform sampling - sanity check: log(Z) ≈ true value of around 33.3 for this example
# pt = pigeons(target = IsingLogPotential(1.0, 5))

# # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# using Plots
# plot(pt.shared.tempering.communication_barriers.localbarrier)

Pigeons.extract_sample(state::IsingState, log_potential) = copy(state.matrix)

# perform sampling - sanity check: log(Z) ≈ true value of around 33.3 for this example
pt = pigeons(
    target = IsingLogPotential(1.0, 5),
    record = [traces; index_process; record_default()],
    # on = ChildProcess(
    #         dependencies = ["examples/ising.jl"],
    #         n_local_mpi_processes = 4,
    #         n_threads = 1),
)

vector = get_sample(pt)
print(length(vector), eltype(vector))
# https://stackoverflow.com/questions/25028539/how-to-implement-an-iterator-in-julia
