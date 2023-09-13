using Pkg
Pkg.activate(".")
Pkg.instantiate()

# import the classes and functions from pt_v2
include("../src/pt_v2.jl")

# perform pt
pt = load(
    pigeons(
        target = ULogPotential(8.0),
        reference = ULogPotential(1.0),
        n_chains = 5,
        n_rounds = 17, # number of iteration = 2^n_rounds
        record = [traces; index_process],
        on = ChildProcess(
            dependencies = ["../src/pt_v2.jl"],
            n_local_mpi_processes = 4,
            n_threads = 1,
        ),
        checkpoint = true,
    ),
)

# extract the target distribution and plot results
vector = get_sample(pt)
print(length(vector), eltype(vector))
x_vector = [element.x for element in vector]
p1 = StatsPlots.plot(x_vector, xlabel = "Iteration", ylabel = "x")
plot_chain(x_vector, 8, dir = "run_pt_result/")

plot2 = StatsPlots.plot(pt.reduced_recorders.index_process);
savefig(plot2, "run_pt_result/U_index_process_plot.png");

plot1 = StatsPlots.plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(plot1, "run_pt_result/U_localbarrier.png")
