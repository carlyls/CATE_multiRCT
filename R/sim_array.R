### Array Method for Parallel Computing ###

library(tidyverse)
library(causalToolbox)
library(rsample)
library(rpart)
library(ranger)
library(glmnet)
library(grf)
library(fastDummies)
library(lme4)


#settings matrix with all the parameters that vary in each parallel run
settings_ind <- expand.grid(K=c(10),
                      n_mean=c(500),
                      n_sd=c(0),
                      study_mean=c(0),
                      study_inter_mean=c(0),
                      study_sds=c("0.5,0", "1,0", "1,0.5", 
                                  "1,1", "3,1"),
                      scenario=c("1a","1b")) %>%
  separate(study_sds, into=c("study_sd", "study_inter_sd"), sep=",") %>%
  mutate(study_sd = as.numeric(study_sd),
         study_inter_sd = as.numeric(study_inter_sd),
         scenario = factor(scenario, levels=c("1a","1b","2"))) %>%
  rbind(c(K=10, n_mean=500, n_sd=0, study_mean=0, study_inter_mean=0, study_sd=NA, study_inter_sd=NA, scenario=2))

settings <- do.call("rbind", replicate(1000, settings_ind, simplify = FALSE)) %>%
  cbind(iteration = rep(1:1000, each=11))


#sets the row of the settings that you will use
i=as.numeric(Sys.getenv('SGE_TASK_ID'))

iteration <- settings$iteration[i]
K <- settings$K[i]
n_mean <- settings$n_mean[i]
n_sd <- settings$n_sd[i]
study_mean <- settings$study_mean[i]
study_inter_mean <- settings$study_inter_mean[i]
study_sd <- settings$study_sd[i]
study_inter_sd <- settings$study_inter_sd[i]
scenario <- settings$scenario[i]
seed <- i + 100000
honesty <- F


#now code
source("Comparing_methods_functions.R", local=T)
source("Simulation_MLOptions.R", local=T)

set.seed(seed)
sim_mses <- compare_mse(K, n_mean, n_sd, study_mean, study_inter_mean,
                        study_sd, study_inter_sd, scenario, honesty)
sim_output <- data.frame(sim_mses, seed=seed, K=K, n_mean=n_mean, n_sd=n_sd, study_sd=study_sd,
                         study_inter_sd=study_inter_sd, scenario=scenario, honesty=honesty, iteration=iteration)

save(sim_output, file=paste(paste("sim_output",seed,"nmean",n_mean,"nsd",n_sd,"studysd",study_sd,
                                  "studyintersd",study_inter_sd,"scenario",scenario,"honesty",honesty,
                                  "iter",iteration,sep="_"),".Rdata", sep=""))

