## Assessing moderators in simulated data ##

library(tidyverse)
library(causalToolbox)
library(rsample)
library(rpart)
library(ranger)
library(glmnet)
library(grf)
library(fastDummies)
library(lme4)

source("R/Comparing_methods_functions.R", local=T)
source("R/Simulation_MLOptions.R", local=T)

#create a single simulated dataset according to scenario 1a
sim_dat <- gen_data(K=10, n_mean=500, n_sd=0, study_mean=0, study_inter_mean=0,
                    study_sd=1, study_inter_sd=0.5, scenario="1a")

covars <- grep("^X", names(sim_dat), value=TRUE)
feat <- select(sim_dat, c(S,all_of(covars))) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

feat_nostudy <- select(sim_dat, all_of(covars))
tr <- sim_dat$W
y <- sim_dat$Y
tau_true <- sim_dat$tau


#fit causal forest with pooling with trial indicator
tau_forest <- causal_forest(X=feat, Y=y, W=tr, num.threads=3, honesty=F, num.trees=1000)
tau_hat <- c(tau_forest$predictions)
# causal_studyind <- mean((tau_hat - tau_true)^2)
sim_preds <- sim_dat %>%
  mutate(tau_hat=tau_hat)


#assess moderation
#plot of all CATE estimates
sim_preds %>%
  ggplot(aes(x=tau_hat)) +
  geom_histogram()

#plot of CATE by X1
sim_preds %>%
  ggplot(aes(x=X1, y=tau_hat, color=S)) +
  geom_point() + ylab("CATE Estimate") +
  labs(color="Trial")

#test calibration
test_calibration(tau_forest) #shows evidence of heterogeneity

#best linear projection
blp <- best_linear_projection(tau_forest, feat)
blp
jtools::plot_summs(blp) +
  ggtitle("Best Linear Projection of CATE")
