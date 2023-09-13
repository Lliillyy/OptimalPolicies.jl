# OptimalPolicies.jl

## Overview

Parallel Tempering is a computational method to sample complex probability 
distributions, which is common in simulating physical and chemical processes. 
The Parallel Tempering algorithm runs multiple Markov chain Monte Carlo (MCMC) 
simulations at different temperatures in parallel to explore high-dimensional 
spaces, where the temperature refers to a parameter that controls the 
acceptance rate of the MCMC simulation. 

Optimal policy search in economics refers to finding policies that maximize or 
minimize certain objectives, such as maximizing profit and minimizing costs. 
The economic models in this context can be highly complex, 
involving numerous variables and parameters. Traditional optimization methods 
could be struggled by those models due to their nature of high-dimensionality. 
However, Parallel Tempering can be adapted to solve those problems by running 
multiple simulations of the economic model with different parameter settings 
defined as different temperatures. Therefore, Parallel Tempering has the 
potential to be a systematic approach of economic policy search for policymakers.
This Julia package aims to utilize the Parallel Tempering algorithm to find 
optimal transportation networks.