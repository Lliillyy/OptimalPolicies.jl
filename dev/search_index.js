var documenterSearchIndex = {"docs":
[{"location":"api.html","page":"API","title":"API","text":"CurrentModule = OptimalPolicies ","category":"page"},{"location":"api.html","page":"API","title":"API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api.html#API","page":"API","title":"API","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"This page is a dump of all the docstrings found in the code. ","category":"page"},{"location":"api.html","page":"API","title":"API","text":"Modules = [OptimalPolicies]\nOrder = [:module, :type, :function, :macro]","category":"page"},{"location":"api.html#OptimalPolicies.UMetropolis","page":"API","title":"OptimalPolicies.UMetropolis","text":"Metropolis\n\nMCMC explorer, should not contain state that is replica-specific  and/or changed in the inner sampling loop (use an Augmentation if  the explorer needs replica-specific auxiliary state information)\n\n\n\n\n\n","category":"type"},{"location":"api.html#OptimalPolicies.U-Tuple{Number, Number}","page":"API","title":"OptimalPolicies.U","text":"U(gam::Number, x::Number)\n\nThe potential function for some given potential barrier height gam\n\n\n\n\n\n","category":"method"},{"location":"api.html#OptimalPolicies.chain-Tuple{Any}","page":"API","title":"OptimalPolicies.chain","text":"chain(target; tune = 0.1, init = 1.0, iters = Int(1e3))\n\nThe Metropolis algorithm targeting a particular potential function target\n\n\n\n\n\n","category":"method"},{"location":"api.html#OptimalPolicies.chains-Tuple{}","page":"API","title":"OptimalPolicies.chains","text":"chains(; pot = U, tune = 0.1, init = 1.0, iters = Int(1e3), temps = [1, 2])\n\nGenerates chains based on the temps\n\n\n\n\n\n","category":"method"},{"location":"api.html#OptimalPolicies.curried-Tuple{Number}","page":"API","title":"OptimalPolicies.curried","text":"curried(gam::Number)\n\nReturns a potential function for the given potential barrier height gam\n\n\n\n\n\n","category":"method"},{"location":"api.html#OptimalPolicies.density-Tuple{Number, Number}","page":"API","title":"OptimalPolicies.density","text":"density(gam::Number, x::Number)\n\nThe density function for some given potential barrier height gam\n\n\n\n\n\n","category":"method"},{"location":"api.html#OptimalPolicies.plot_chain-Tuple{Any, Any}","page":"API","title":"OptimalPolicies.plot_chain","text":"plot_chain(xvec, gam; dir = \"result/\")\n\nPlot the value, autocorrelation, and density of a chain xvec given the value of gam, store the plot in the directory dir\n\n\n\n\n\n","category":"method"},{"location":"api.html#OptimalPolicies.print_summary-Tuple{Any, Any}","page":"API","title":"OptimalPolicies.print_summary","text":"print_summary(mat, temps)\n\nComputes and displays the statistics of the given matrix mat\n\n\n\n\n\n","category":"method"},{"location":"index.html#OptimalPolicies.jl","page":"Introduction","title":"OptimalPolicies.jl","text":"","category":"section"},{"location":"index.html#Overview","page":"Introduction","title":"Overview","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Parallel Tempering is a computational method to sample complex probability  distributions, which is common in simulating physical and chemical processes.  The Parallel Tempering algorithm runs multiple Markov chain Monte Carlo (MCMC)  simulations at different temperatures in parallel to explore high-dimensional  spaces, where the temperature refers to a parameter that controls the  acceptance rate of the MCMC simulation. ","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Optimal policy search in economics refers to finding policies that maximize or  minimize certain objectives, such as maximizing profit and minimizing costs.  The economic models in this context can be highly complex,  involving numerous variables and parameters. Traditional optimization methods  could be struggled by those models due to their nature of high-dimensionality.  However, Parallel Tempering can be adapted to solve those problems by running  multiple simulations of the economic model with different parameter settings  defined as different temperatures. Therefore, Parallel Tempering has the  potential to be a systematic approach of economic policy search for policymakers. This Julia package aims to utilize the Parallel Tempering algorithm to find  optimal transportation networks.","category":"page"}]
}