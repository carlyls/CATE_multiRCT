# CATE_multiRCT

This repository includes code for performing machine learning methods that estimate the CATE by combining data from multiple RCTs.

R files uploaded include those that apply the methods, and those that run simulation iterations to compare the methods in terms of performance.

Comparing_methods_functions.R performs the ensembling approaches based on those by Tan et al., 2021 (<https://github.com/ellenxtan/ifedtree>). These methods are adapted from Tan et al.'s ifedtree package to be applied to all trial datasets, rather than one coordinating dataset as required in Tan et al.'s data setting. The methods include X-Learner, S-Learner, and Causal Forest base methods, and ensemble regression tree, random forest, and lasso.

Simulation_MLOptions.R generates data for simulation based on a few different data generation setups, and it then performs each method on the simulated data and saves the MSEs for each method.

Finally, sim_array.R sets up the simulation for parallel computing in a high-performance computing cluster.
