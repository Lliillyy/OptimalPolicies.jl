using Pkg
Pkg.activate(".")

using Pigeons
using Random
using Distributions
using StatsPlots
import Base.@kwdef
include("mcmc.jl")

# not application specific, just like this
@kwdef struct ULogPotential
    beta::Float64
end

# in this struct we hold the state of the system, which is going to be very application-specific
# here, it's just a number, but this might be itself a large object, e.g. a Matrix, or a lists of lists, etc.
# ? where do we store static data objects that enter computation?
mutable struct UState
    x::Float64
    # ? what about other state variables?
end

# we need to define a copy function for the state, so that we can copy it when we need to
Base.copy(state::UState) = UState(state.x)

# ! This is crucial, here we are computing the objective function as a function of the state (x)
# Make ULogPotential conform the log_potential informal interface
(log_potential::ULogPotential)(state::UState) = -log_potential.beta * (state.x^2 - 1)^2

# Reference distribution 
# this is the lowest beta used in the tempering schedule, i.e. the reference distribution
# if we don't specify `reference=ULogPotential(beta)` as an argument for pt() we will be default use this
Pigeons.default_reference(log_potential::ULogPotential) = ULogPotential(1.0)

# Initialization
# initial conditions for the state of the system (x)
Pigeons.initialization(log_potential::ULogPotential, ::AbstractRNG, ::Int) = UState(1.0)

# MCMC explorer 
# This struct should not contain state that is replica-specific 
#     and/or changed in the inner sampling loop (use an Augmentation if 
#     the explorer needs replica-specific auxiliary state information)
# ? should we store anything here?
@kwdef struct UMetropolis end

# Make UMetropolis conform the explorer informal interface
Pigeons.default_explorer(lp::ULogPotential) = UMetropolis()


# Perform explorer MCMC step
# update within each chain
# glossary:
# explorer = object that does this updating
# replica = roughly corresponds to a chain. 
#    It's an object that holds the state of the system + 
#    random number generator state + 
#    some global parameters
# shared = even higher level, knows about other replicas
function Pigeons.step!(explorer::UMetropolis, replica, shared)

    # find the log potential for this replica
    # this is a function
    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)

    # ! key part to customize
    ### propose a change and record the change in log probability

    # current objective value 
    # * if in a given application this is an expensive operation, might also store it as a field somewhere to avoid recomputing it
    log_pr_before = log_potential(replica.state)

    # copy/store current state
    # * if in a given application this is a large object, will instead change it in place, and then change it back if we reject the change (see Ising example)
    prev_x = replica.state.x

    # update the state with a new candidate
    replica.state.x = prev_x + randn(replica.rng) * 0.1

    # new objective value
    log_pr_after = log_potential(replica.state)

    # accept-reject step 
    accept_ratio = log_pr_after - log_pr_before
    if log(rand(replica.rng)) >= accept_ratio
        # reject: revert the move we just proposed
        # * can do this in place if x is large using the inverse function of the candidate proposal function
        replica.state.x = prev_x
    end # (nothing to do if accept, we work in-place)
end

# perform pt
pt = pigeons(
    target = ULogPotential(8.0),
    # reference = ULogPotential(1.0),
    n_chains = 5,
    n_rounds = 17, # number of iteration = 2^n_rounds
    record = [traces; index_process],
    # on = ChildProcess(
    #         n_local_mpi_processes = 4,
    #         n_threads = 1)
)

vector = get_sample(pt)
print(length(vector), eltype(vector))
x_vector = [element.x for element in vector]
p1 = StatsPlots.plot(x_vector, xlabel = "Iteration", ylabel = "x")
plot_chain(x_vector, 8, dir = "pt_v3_result/")

for i = 1:5
    pt.reduced_recorders.index_process[i] =
        pt.reduced_recorders.index_process[i][10000:10050]
end

plot2 = StatsPlots.plot(pt.reduced_recorders.index_process) #|> display
savefig(plot2, "pt_v3_result/U_index_process_plot.png");

plot1 = StatsPlots.plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(plot1, "pt_v3_result/U_localbarrier.png")
