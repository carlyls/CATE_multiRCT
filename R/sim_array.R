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


## Settings matrix with all the parameters that vary in each parallel run
settings_og <- expand.grid(K=c(10),
                           ns=c("same"),
                           cov_shift=c("no"),
                           study_sds=c("0.5,0", "1,0", "1,0.5", 
                                       "1,1", "3,1"),
                           scenario=c("1a","1b")) %>%
  separate(study_sds, into=c("study_sd", "study_inter_sd"), sep=",") %>%
  mutate(ns = factor(ns, levels=c("same","one large","half and half")),
         cov_shift = factor(cov_shift, levels=c("no","yes")),
         scenario = factor(scenario, levels=c("1a","1b","2"))) %>%
  rbind(c(K=10, ns="same", cov_shift="no", study_sd=NA, study_inter_sd=NA, scenario="2"))

settings_new <- expand.grid(K=c(10), 
                            ns=c("one large", "half and half"), 
                            cov_shift=c("no"), 
                            study_sd=1, 
                            study_inter_sd=0.5, 
                            scenario=c("1a","1b")) %>%
  mutate(ns = factor(ns, levels=c("same","one large","half and half")),
         cov_shift = factor(cov_shift, levels=c("no","yes")),
         scenario = factor(scenario, levels=c("1a","1b","2"))) %>%
  rbind(c(K=10, ns="one large", cov_shift="no", study_sd=NA, study_inter_sd=NA, scenario="2"),
        c(K=10, ns="half and half", cov_shift="no", study_sd=NA, study_inter_sd=NA, scenario="2"),
        c(K=10, ns="same", cov_shift="yes", study_sd=1, study_inter_sd=0.5, scenario="1a"),
        c(K=10, ns="one large", cov_shift="yes", study_sd=1, study_inter_sd=0.5, scenario="1a"),
        c(K=10, ns="same", cov_shift="yes", study_sd=1, study_inter_sd=0.5, scenario="1b"),
        c(K=10, ns="one large", cov_shift="yes", study_sd=1, study_inter_sd=0.5, scenario="1b"),
        c(K=30, ns="same", cov_shift="no", study_sd=1, study_inter_sd=0.5, scenario="1a"),
        c(K=30, ns="same", cov_shift="no", study_sd=1, study_inter_sd=0.5, scenario="1b"))

settings_combo <- settings_og %>%
  rbind(settings_new) %>%
  mutate(K = as.numeric(K),
         study_sd = as.numeric(study_sd),
         study_inter_sd = as.numeric(study_inter_sd))
  
settings <- do.call("rbind", replicate(1000, settings_combo, simplify = FALSE)) %>%
  cbind(iteration = rep(1:1000, each=nrow(settings_combo)))


## Set parameters for given iteration
i=as.numeric(Sys.getenv('SGE_TASK_ID'))

iteration <- settings$iteration[i]
study_mean <- 0
study_inter_mean <- 0
K <- settings$K[i]
ns <- settings$ns[i]
cov_shift <- settings$cov_shift[i]
study_sd <- settings$study_sd[i]
study_inter_sd <- settings$study_inter_sd[i]
scenario <- settings$scenario[i]
seed <- i
honesty <- F

# NOTE: some iterations did not run in cluster because of memory spikes in the causal forest procedure;
# therefore, some seeds had to be rerun or added on after original iterations

## Run code
source("R/Comparing_methods_functions.R", local=T)
source("R/Simulation_MLOptions.R", local=T)

set.seed(seed)
sim_mses <- compare_mse(K, ns, cov_shift, study_mean, study_inter_mean,
                        study_sd, study_inter_sd, scenario, honesty)
sim_output <- data.frame(sim_mses, seed=seed, K=K, ns=ns, cov_shift=cov_shift,
                         study_sd=study_sd, study_inter_sd=study_inter_sd, 
                         scenario=scenario, honesty=honesty, iteration=iteration)

save(sim_output, file=paste(paste("sim_output",seed,"K",K,"ns",ns,"covshift",cov_shift,
                                  "studysd",study_sd,"studyintersd",study_inter_sd,
                                  "scenario",scenario,"honesty",honesty,
                                  "iter",iteration,sep="_"),".Rdata", sep=""))

