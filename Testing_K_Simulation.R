## Testing Differences in K (Number of Trials) ##

library(tidyverse)

source("R/Comparing_methods_functions.R")
source("R/Simulation_MLOptions.R", local=T)

# for four different values of K:
# save indicators of top coefficients (most different studies)
# fit causal forest with pooling with trial indicator
# fit causal forest with ensemble forest
# figure out variable importance of the top coefficients
# calculate percent of those that get picked up by the two methods


test_k <- function(K) {
  
  # generate data 
  all_dat <- data.frame()
  n_study <- rep(500, K)
  study_main <- rnorm(K, mean=0, sd=0.5)
  study_inter <- rnorm(K, mean=0, sd=0)
  
  for (k in 1:K) {
    n <- n_study[k]
    #sample covariates
    Sigma <- matrix(.2, nrow=5, ncol=5)
    diag(Sigma) <- 1
    dat <- MASS::mvrnorm(n=n, mu=rep(0,5), Sigma=Sigma) %>% data.frame()
    #treatment
    dat$W <- rbinom(n, size=1, prob=0.5)
    #study and id
    dat$S <- rep(k, n)
    dat$id <- seq(1, n)
    #noise
    dat$eps <- rnorm(n, mean=0, sd=sqrt(0.01))
    all_dat <- bind_rows(all_dat, dat)
  }
  
  #tau and Y
  all_dat$m <- all_dat$X1/2 + all_dat$X2 + all_dat$X3 + all_dat$X4 +
      study_main[all_dat$S] + study_inter[all_dat$S]*all_dat$X1
  all_dat$tau <- all_dat$X1*(all_dat$X1>0) + study_main[all_dat$S] + 
      study_inter[all_dat$S]*all_dat$X1
  
  all_dat$Y <- all_dat$m + (2*all_dat$W-1)/2*all_dat$tau + all_dat$eps
  
  sim_dat <- all_dat %>%
    select(-eps,-m) %>%
    mutate(S = factor(S)) %>%
    relocate(S, id, W, X1, X2, X3, X4, X5, Y, tau)
  
  #save top coefficients
  #find farthest 20% from mean
  dist <- abs(study_main - mean(study_main))
  top_coefs <- which(dist > quantile(dist, .8))
  top_coef_vals <- data.frame(top_coef_vals=dist[top_coefs]) %>%
    arrange(desc(top_coef_vals)) %>%
    t()
  colnames(top_coef_vals) <- paste0("top_coef_var", 1:ncol(top_coef_vals))
  
  #fit models
  covars <- grep("^X", names(sim_dat), value=TRUE)
  feat <- select(sim_dat, c(S,all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  tr <- sim_dat$W
  y <- sim_dat$Y
  tau_true <- sim_dat$tau
  
  #study indicator
  tau_forest <- causal_forest(X=feat, Y=y, W=tr, num.threads=3, honesty=honesty, num.trees=1000)
  varimp <- data.frame(varimp = variable_importance(tau_forest))
  rownames(varimp) <- paste0(colnames(feat), "_imp")
  
  #report average importances
  x_imp <- varimp %>%
    slice(1:5) %>% 
    t()
  study_imp <- varimp %>%
    slice(-c(1:5)) %>%
    slice(top_coefs)
  avg_study_imp <- mean(study_imp$varimp)
  
  return(data.frame(x_imp, avg_study_imp, top_coef_vals))
}

setwd("~/Dropbox/Moderation/MLSimulations_CLS/Papers/SIM_Revision/comparingk_results")
set.seed(222)
res10 <- map_dfr(.x = c(rep(10, 50)),
                 .f = test_k)
saveRDS(res10, "k10_imp.rds")

set.seed(222)
res15 <- map_dfr(.x = c(rep(15, 50)),
                 .f = test_k)
saveRDS(res15, "k15_imp.rds")

set.seed(222)
res20 <- map_dfr(.x = c(rep(20, 50)),
                 .f = test_k)
saveRDS(res20, "k20_imp.rds")

set.seed(222)
res25 <- map_dfr(.x = c(rep(25, 50)),
                 .f = test_k)
saveRDS(res25, "k25_imp.rds")

set.seed(222)
res30 <- map_dfr(.x = c(rep(30, 50)),
                 .f = test_k)
saveRDS(res30, "k30_imp.rds")


#compare
all_res <- data.frame(colMeans(res10),
           colMeans(res20)[1:8],
           colMeans(res25)[1:8],
           colMeans(res30)[1:8]) %>%
  t()
View(all_res)

#get avg across all non-moderating covariates
mean(all_res[1,2:5])

#get sds
sapply(res10, sd)
sd(c(res10[,2],res10[,3],res10[,4],res10[,5]))
