### Simulation Approach for Comparing ML Options with Multiple RCTs ###

library(dplyr)
library(causalToolbox)
library(rsample)
library(rpart)
library(ranger)
library(glmnet)
library(grf)
library(fastDummies)
library(lme4)

source("Comparing_methods_functions.R")

## Data Generation Function
## Tau based on Tan et al., 2021; Wager and Athey, 2018; Kunzel et al., 2019

gen_data <- function (K, n_mean, n_sd, study_mean, study_inter_mean,
                      study_sd, study_inter_sd, scenario,
                      ncovar=5, sd=sqrt(0.01)) {
  
  all_dat <- data.frame()
  
  n_study <- floor(rnorm(K, mean=n_mean, sd=n_sd))
  
  if (scenario %in% c("1a","1b")) {
    study_main <- rnorm(K, mean=study_mean, sd=study_sd)
    study_inter <- rnorm(K, mean=study_inter_mean, sd=study_inter_sd)
  }
  
  for (k in 1:K) {
    
    n <- n_study[k]
    
    #sample covariates
    dat <- data.frame(matrix(rnorm(n*ncovar), nrow=n, ncol=ncovar))
    colnames(dat) <- paste0("X", seq(1,ncovar))
    
    #treatment
    dat$W <- rbinom(n, size=1, prob=0.5)
    
    #study and id
    dat$S <- rep(k, n)
    dat$id <- seq(1, n)
    
    #noise
    dat$eps <- rnorm(n, mean=0, sd=sd)
    
    all_dat <- bind_rows(all_dat, dat)
  }
  
  #tau and Y
  if (scenario == "1a") {
    all_dat$m <- all_dat$X1/2 + all_dat$X2 + all_dat$X3 + all_dat$X4 +
      study_main[all_dat$S] + study_inter[all_dat$S]*all_dat$X1
    all_dat$tau <- all_dat$X1*(all_dat$X1>0) + study_main[all_dat$S] + 
      study_inter[all_dat$S]*all_dat$X1
  }
  
  if (scenario == "1b") {
    all_dat$m <- 0
    all_dat$tau <- (2/(1+exp(-12*(all_dat$X1-1/2))))*(2/(1+exp(-12*(all_dat$X2-1/2)))) +
      study_main[all_dat$S] + study_inter[all_dat$S]*all_dat$X1
  }
  
  if (scenario == "2") {
    all_dat$m <- all_dat$X1/2 + all_dat$X2 + all_dat$X3 + all_dat$X4
    all_dat$tau <- ifelse(all_dat$S %in% c(1:4), (2/(1+exp(-12*(all_dat$X1-1/2))))*(2/(1+exp(-12*(all_dat$X2-1/2)))),
                          ifelse(all_dat$S %in% c(5:8), all_dat$X1*(all_dat$X1>0), 0))
  }
  
  all_dat$Y <- all_dat$m + (2*all_dat$W-1)/2*all_dat$tau + all_dat$eps
  
  all_dat <- all_dat %>%
    select(-eps,-m) %>%
    mutate(S = factor(S)) %>%
    relocate(S, id, W, X1, X2, X3, X4, X5, Y, tau)
  return(all_dat)
}


## Comparing Methods
compare_mse <- function (K, n_mean, n_sd, study_mean, study_inter_mean,
                         study_sd, study_inter_sd, scenario, honesty) {
  
  #generate data
  sim_dat <- gen_data(K, n_mean, n_sd, study_mean, study_inter_mean,
                      study_sd, study_inter_sd, scenario)
  
  covars <- grep("^X", names(sim_dat), value=TRUE)
  feat <- select(sim_dat, c(S,all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  
  feat_nostudy <- select(sim_dat, all_of(covars))
  tr <- sim_dat$W
  y <- sim_dat$Y
  tau_true <- sim_dat$tau
  
  
  ### causal forest options ###
  
  #no study
  tau_forest_nostudy <- causal_forest(X=feat_nostudy, Y=y, W=tr, num.threads=3, honesty=honesty, num.trees=1000)
  tau_hat_nostudy <- predict(tau_forest_nostudy)$predictions
  causal_nostudy <- mean((tau_hat_nostudy - tau_true)^2)
  
  rm(list = c("tau_forest_nostudy", "tau_hat_nostudy"))
  
  #study indicator
  tau_forest <- causal_forest(X=feat, Y=y, W=tr, num.threads=3, honesty=honesty, num.trees=1000)
  tau_hat <- predict(tau_forest)$predictions
  causal_studyind <- mean((tau_hat - tau_true)^2)
  
  rm(list = c("tau_forest", "tau_hat"))
  
  #tan
  causal_tan <- tan_preds(K, sim_dat, covars, "causalforest", honesty=honesty)
  causal_tree <- causal_tan["mse_tree"]
  causal_forest <- causal_tan["mse_forest"]
  causal_lasso <- causal_tan["mse_lasso"]
  
  rm(list = c("causal_tan"))
  
  
  ### x-learner options ###
  
  #no study
  x_rf_nostudy <- X_RF(feat=feat_nostudy, tr=tr, yobs=y, nthread=3)
  cate_x_rf_nostudy <- EstimateCate(x_rf_nostudy, feat_nostudy)
  x_nostudy <- mean((cate_x_rf_nostudy - tau_true)^2)
  
  rm(list = c("x_rf_nostudy", "cate_x_rf_nostudy"))
  
  #study indicator
  x_rf_studyind <- X_RF(feat = feat, tr = tr, yobs = y, nthread=3)
  cate_x_rf_studyind <- EstimateCate(x_rf_studyind, feat)
  x_studyind <- mean((cate_x_rf_studyind - tau_true)^2)
  
  rm(list = c("x_rf_studyind", "cate_x_rf_studyind"))
  
  #tan
  x_tan <- tan_preds(K, sim_dat, covars, "xlearner")
  x_tree <- x_tan["mse_tree"]
  x_forest <- x_tan["mse_forest"]
  x_lasso <- x_tan["mse_lasso"]
  
  rm(list = c("x_tan"))
  
  
  ### s-learner options ###
  
  #no study
  s_rf_nostudy <- S_RF(feat=feat_nostudy, tr=tr, yobs=y, nthread=3)
  cate_s_rf_nostudy <- EstimateCate(s_rf_nostudy, feat_nostudy)
  s_nostudy <- mean((cate_s_rf_nostudy - tau_true)^2)
  
  rm(list = c("s_rf_nostudy", "cate_s_rf_nostudy"))
  
  #study indicator
  s_rf_studyind <- S_RF(feat = feat, tr = tr, yobs = y, nthread=3)
  cate_s_rf_studyind <- EstimateCate(s_rf_studyind, feat)
  s_studyind <- mean((cate_s_rf_studyind - tau_true)^2)
  
  rm(list = c("s_rf_studyind", "cate_s_rf_studyind"))
  
  #tan
  s_tan <- tan_preds(K, sim_dat, covars, "slearner")
  s_tree <- s_tan["mse_tree"]
  s_forest <- s_tan["mse_forest"]
  s_lasso <- s_tan["mse_lasso"]
  
  rm(list = c("s_tan"))
  
  
  ### meta-analysis ###
  
  mixed_mod <- lmer(Y ~ W*X1 + W*X2 + W*X3 + W*X4 + W*X5 + (1|S), data=sim_dat) #naive
  sim_dat_trt <- mutate(sim_dat, W=1)
  sim_dat_cntrl <- mutate(sim_dat, W=0)
  cate_ma <- predict(mixed_mod, sim_dat_trt) - predict(mixed_mod, sim_dat_cntrl)
  ma <- mean((cate_ma - tau_true)^2)
  
  
  mses <- data.frame(x_nostudy, x_studyind, x_tree, x_forest, x_lasso,
                     causal_nostudy, causal_studyind, causal_tree, causal_forest, causal_lasso,
                     s_nostudy, s_studyind, s_tree, s_forest, s_lasso, ma)
  
  return(mses)
}


