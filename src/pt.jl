using Pigeons
using Random
using Distributions
import Base.@kwdef

@kwdef struct ULogPotential
    beta::Float64
end

mutable struct UState
    x::Float64
    U_value::Float64
end

function UState(x::Float64)
    result = UState(x, -1)
    result.U_value = recompute_U(result)
    return result
end

recompute_U(state::UState) = (state.x^2 - 1)^2

Base.copy(state::UState) = UState(state.x)
# Base.length(state::UState) = 1

# Make ULogPotential conform the log_potential informal interface
(log_potential::ULogPotential)(state::UState) = log_potential.beta * state.U_value

# Reference distribution uses beta = 0
Pigeons.default_reference(log_potential::ULogPotential) = ULogPotential(0.0)

# Initialization
Pigeons.initialization(log_potential::ULogPotential, ::AbstractRNG, ::Int) = UState(0.0)

# MCMC explorer 
# This struct should not contain state that is replica-specific 
#     and/or changed in the inner sampling loop (use an Augmentation if 
#     the explorer needs replica-specific auxiliary state information)
@kwdef struct UMetropolis
    # n_steps::Int = 1
end

Pigeons.default_explorer(lp::ULogPotential) = UMetropolis()

# Perform explorer MCMC step
function Pigeons.step!(explorer::UMetropolis, replica, shared)
    # sampler = Normal(0.0, 0.1) # rand(replica.rng, sampler)
    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)

    # propose a change and record the change in log probability
    log_pr_before = log_potential(replica.state)
    prev_x = replica.state.x
    replica.state.x = prev_x + randn(replica.rng) * 0.2
    replica.state.U_value = recompute_U(replica.state)
    log_pr_after = log_potential(replica.state)

    # accept-reject step 
    accept_ratio = exp(log_pr_after - log_pr_before)
    if accept_ratio < 1 && rand(replica.rng) > accept_ratio
        # reject: revert the move we just proposed
        replica.state = UState(prev_x)
    end # (nothing to do if accept, we work in-place)
end

# perform sampling - sanity check: log(Z) ≈ true value of around 33.3 for this example
pt = pigeons(target = ULogPotential(8.0), record = [traces; index_process])

# pt = pigeons(
#     target = ULogPotential(8.0),
#     reference = ULogPotential(0.0001),
#     record = [traces; index_process],
#     # checked_round = 3,
#     # checkpoint = true,
#     # on = ChildProcess(
#     #         n_local_mpi_processes = 4,
#     #         n_threads = 1),
# )

using StatsPlots
# using MCMCChains
# samples = Chains(get_sample(pt), variable_names(pt))
vector = get_sample(pt)
print(length(vector))
x_vector = [element.x for element in vector]
u_vector = [element.U_value for element in vector]
p1 = plot(x_vector, xlabel = "Iteration", ylabel = "x")
p2 = plot(u_vector, xlabel = "Iteration", ylabel = "U")
chain_plot = plot(p1, p2, layout = (2, 1), size = (1500, 600))
# chain_plot = plot(samples)
savefig(chain_plot, "chain_plot.png");

# # # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# using Plots

# plot2 = plot(pt.reduced_recorders.index_process);
# savefig(plot2, "U_index_process_plot1.png");

# # plotlyjs() this line creates error
# plot1 = plot(pt.shared.tempering.communication_barriers.localbarrier)
# savefig(plot1, "U_localbarrier1.png")
