using Pigeons
using Random
using StatsPlots
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
Base.length(state::UState) = 1

# Make IsingLogPotential conform the log_potential informal interface
(log_potential::ULogPotential)(state::UState) = log_potential.beta * state.U_value


# sample
# function iid_bb!(state::UState)
#     @assert recompute_U(state) == state.U_value
#     alpha = 2   # Shape parameter of the Beta distribution
#     beta = 2     # Shape parameter of the Beta distribution
#     n = 500         # Number of trials
#     bb_dist = BetaBinomial(n, alpha, beta)
#     state.x = rand(bb_dist)
#     state.U_value = recompute_U(state)
#     return nothing
# end

# Reference distribution uses beta = 0...
Pigeons.default_reference(log_potential::ULogPotential) = ULogPotential(0.0)
# ... so that we can do i.i.d. sampling of Bernoullis at the reference:
# function Pigeons.sample_iid!(reference_log_potential::ULogPotential, replica, shared)
#     # @assert reference_log_potential.beta == 0.0
#     iid_bb!(replica.state)
# end

# Initialization: all entries to zeros (falses)
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
    # alpha = 2    # Shape parameter of the Beta distribution
    # beta = 2     # Shape parameter of the Beta distribution
    # n = 500         # Number of trials
    # bb_dist = BetaBinomial(n, alpha, beta) # rand(bb_dist)
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
pt = pigeons(
    target = ULogPotential(8.0),
    explorer = UMetropolis(),
    record = [traces; index_process],
)

# pt = pigeons(
#     target = ULogPotential(8.0),
#     # reference = ULogPotential(0.0001),
#     record = [traces; index_process],
#     # checked_round = 3,
#     # checkpoint = true,
#     # on = ChildProcess(
#     #         n_local_mpi_processes = 4,
#     #         n_threads = 1),
# )

# samples = Chains(sample_array(pt), variable_names(pt))
# chain_plot = plot(samples)
# savefig(chain_plot, "julia_posterior_densities_and_traces.png");

# # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# using Plots

plot2 = plot(pt.reduced_recorders.index_process);
savefig(plot2, "U_index_process_plot.png");

# plotlyjs() this line creates error
plot1 = plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(plot1, "U_localbarrier.png")
