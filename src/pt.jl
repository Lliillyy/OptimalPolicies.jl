using Pigeons
using SplittableRandoms
import Base.@kwdef
using Distributions

struct ULogPotential 
    beta::Float64
end 

mutable struct UState 
    x::Float64 
    U_value::Int
end

function UState(x::Float64)
    result = UState(x, -1)
    result.U_value = recompute_U(result) 
    return result
end

function recompute_U(state::UState) 
    return (state.x^2 - 1)^2
end

# sample all binary variable i.i.d. Bern(1/2)
function iid_bb!(state::UState)
    @assert recompute_U(state) == state.U_value
    alpha = 0.2    # Shape parameter of the Beta distribution
    beta = 0.25     # Shape parameter of the Beta distribution
    n = 10         # Number of trials
    bb_dist = BetaBinomial(n, alpha, beta)
    state.x = rand(bb_dist)
    state.U_value = recompute_U(state)
    return nothing
end


# Make IsingLogPotential conform the log_potential informal interface
(log_potential::ULogPotential)(state::UState) = log_potential.beta * state.U_value

# Reference distribution uses beta = 0...
Pigeons.create_reference_log_potential(log_potential::ULogPotential, ::Inputs) = ULogPotential(0.0)
# ... so that we can do i.i.d. sampling of Bernoullis at the reference:
function Pigeons.sample_iid!(reference_log_potential::ULogPotential, replica, shared)
    @assert reference_log_potential.beta == 0.0
    iid_bb!(replica.state)
end

# Initialization: all entries to zeros (falses)
Pigeons.create_state_initializer(my_potential::ULogPotential, ::Inputs) = my_potential
Pigeons.initialization(log_potential::ULogPotential, ::SplittableRandom, ::Int) = UState(0.0)

# MCMC explorer 
# This struct should not contain state that is replica-specific 
#     and/or changed in the inner sampling loop (use an Augmentation if 
#     the explorer needs replica-specific auxiliary state information)
@kwdef struct UMetropolis
    n_steps::Int = 3
end 

Pigeons.default_explorer(lp::ULogPotential) = UMetropolis()

# Perform explorer MCMC step
function Pigeons.step!(explorer::UMetropolis, replica, shared)
    alpha = 0.2    # Shape parameter of the Beta distribution
    beta = 0.25     # Shape parameter of the Beta distribution
    n = 10         # Number of trials
    bb_dist = BetaBinomial(n, alpha, beta)
    # Note: the log_potential is an InterpolatedLogPotential between two IsingLogPotential's
    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)
    for k in 1:explorer.n_steps
        # propose a change and record the change in log probability
        log_pr_before = log_potential(replica.state)
        prev_x = replica.state.x
        replica.state.x = rand(bb_dist)
        replica.state.U_value = recompute_U(replica.state)
        log_pr_after = log_potential(replica.state) 
        # accept-reject step 
        accept_ratio = exp(log_pr_after - log_pr_before) 
        if accept_ratio < 1 && rand(replica.rng) > accept_ratio 
            # reject: revert the move we just proposed
            replica.state.x = prev_x
            replica.state.U_value = recompute_U(replica.state)
        end # (nothing to do if accept, we work in-place)
end
end

# perform sampling - sanity check: log(Z) ≈ true value of around 33.3 for this example
pt = pigeons(target = ULogPotential(1.0), recorder_builders = [index_process])

# # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# import Pkg; Pkg.add("Plots")
using Plots

plot2 = plot(pt.reduced_recorders.index_process);
savefig(plot2, "U_index_process_plot.png");

# plotlyjs()
plot1 = plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(plot1, "U_localbarrier.png")