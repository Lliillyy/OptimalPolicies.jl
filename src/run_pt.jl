using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("pt_v2.jl")

Pigeons.extract_sample(state::UState, log_potential) = copy(state.x) #or copy(state)

pt = pigeons(
    target = ULogPotential(8.0),
    reference = ULogPotential(1.0),
    n_chains = 5,
    n_rounds = 14, # number of iteration = 2^n_rounds
    record = [traces; index_process],
    on = ChildProcess(
            dependencies = ["src/pt_v2.jl"],
            n_local_mpi_processes = 4,
            n_threads = 1),
)

vector = get_sample(pt)
print(length(vector), eltype(vector))
# x_vector = [element.x for element in vector]
# p1 = StatsPlots.plot(x_vector, xlabel = "Iteration", ylabel = "x")
# plot_chain(x_vector, 8, dir = "pt_result/")

# # sanity check: the local communication barrier has a peak near the predicted phase transition log(1+sqrt(2))/2
# # using Plots

# plot2 = StatsPlots.plot(pt.reduced_recorders.index_process);
# savefig(plot2, "pt_result/U_index_process_plot.png");

# # plotlyjs() this line creates error
# plot1 = StatsPlots.plot(pt.shared.tempering.communication_barriers.localbarrier)
# savefig(plot1, "pt_result/U_localbarrier.png")
