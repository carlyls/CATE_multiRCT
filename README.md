# CATE_multiRCT

This repository includes code for performing machine learning methods that estimate the conditional average treatment effect (CATE) by combining data from multiple RCTs.

R files uploaded include those that apply the methods, and those that run simulation iterations to compare the methods in terms of performance.

In the R folder:

-   Comparing_methods_functions.R performs the ensembling approaches based on those by Tan et al., 2021 (<https://github.com/ellenxtan/ifedtree>). These methods are adapted from Tan et al.'s ifedtree package to be applied to all trial datasets, rather than one coordinating dataset as required in Tan et al.'s data setting. The methods include X-Learner, S-Learner (Kunzel et al., 2019 - causalToolbox package), and Causal Forest (Athey et al., 2019 - grf package) base methods, and ensemble regression tree, random forest, and lasso.

-   Simulation_MLOptions.R generates data for simulation based on a few different data generation setups, and it then performs each method on the simulated data and saves the MSEs for each method.

-   Sim_Array.R sets up the simulation for parallel computing in a high-performance computing cluster.

-   Sim_Combine_Results.Rmd reads in the simulation results and creates boxplot summaries of MSE.

-   Modifier_Simulation.R generates a single simulated datasets and assesses the CATE estimates and potential moderation.
