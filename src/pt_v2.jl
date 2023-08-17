using Pigeons
using Random
using Distributions
using StatsPlots
import Base.@kwdef

@kwdef struct ULogPotential
    beta::Float64
end

mutable struct UState
    x::Float64
end

Base.copy(state::UState) = UState(state.x)

# Make ULogPotential conform the log_potential informal interface
(log_potential::ULogPotential)(state::UState) = log_potential.beta * (state.x^2 - 1)^2

# Reference distribution uses beta = 1.0
Pigeons.default_reference(log_potential::ULogPotential) = ULogPotential(0.0)

# Initialization
Pigeons.initialization(log_potential::ULogPotential, ::AbstractRNG, ::Int) = UState(0.0)

# MCMC explorer 
# This struct should not contain state that is replica-specific 
#     and/or changed in the inner sampling loop (use an Augmentation if 
#     the explorer needs replica-specific auxiliary state information)
@kwdef struct UMetropolis end

Pigeons.default_explorer(lp::ULogPotential) = UMetropolis()

# Perform explorer MCMC step
function Pigeons.step!(explorer::UMetropolis, replica, shared)
    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)

    # propose a change and record the change in log probability
    log_pr_before = log_potential(replica.state)
    prev_x = replica.state.x
    replica.state.x = prev_x + randn(replica.rng) * 0.1
    log_pr_after = log_potential(replica.state)

    # accept-reject step 
    accept_ratio = exp(log_pr_after - log_pr_before)
    if accept_ratio < 1 && rand(replica.rng) > accept_ratio
        # reject: revert the move we just proposed
        replica.state.x = prev_x
    end # (nothing to do if accept, we work in-place)
end

# perform sampling
pt = pigeons(
    target = ULogPotential(8.0),
    reference = ULogPotential(1.0),
    n_chains = 5,
    record = [traces; index_process],
)


# function target_chains(pt::PT) 
#     n = n_chains(pt.inputs)
#     return filter(i -> is_target(pt.shared.tempering.swap_graphs, i), 1:n)
# end

vector = get_sample(pt)
# print(length(vector))
x_vector = [element.x for element in vector]
p1 = StatsPlots.plot(x_vector, xlabel = "Iteration", ylabel = "x")
savefig(p1, "chain_plot1.png");

# chain_plot = plot(p1, p2, layout = (2, 1), size = (1500, 600))
# chain_plot = plot(samples)
# savefig(chain_plot, "chain_plot1.png");

# # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# using Plots

plot2 = StatsPlots.plot(pt.reduced_recorders.index_process);
savefig(plot2, "U_index_process_plot1.png");

# plotlyjs() this line creates error
plot1 = StatsPlots.plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(plot1, "U_localbarrier1.png")
