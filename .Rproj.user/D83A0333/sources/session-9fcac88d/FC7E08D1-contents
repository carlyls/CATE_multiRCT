## Functions for Simulation_MLOptions_28Jul2022.R

## Tan method functions
EnsemTreeAll <- function(aug_data, site, covars) {
  
  fml <- as.formula(paste("tau_pred ~ ", site, " + ", paste(covars, collapse="+")))
  
  #adaptive ET
  myfit <- rpart(fml, data=aug_data)
  myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])
  
  return(myfit)
}

EnsemForestAll <- function(aug_data, site, covars) {
  
  fml <- as.formula(paste("tau_pred ~ ", site, " + ", paste(covars, collapse="+")))
  
  myfit <- ranger(fml, data=aug_data, respect.unordered.factors='order', 
                  importance="impurity", num.threads=3)
  
  return(myfit)
}

EnsemLassoAll <- function(aug_data, sim_dat_mod, site, covars) {
  
  fml <- as.formula(paste(paste(" ~ ", site, " + ", paste(covars, collapse="+")), 
                          paste(site, covars, sep="*", collapse=" + "), sep=" + "))
  grid <- 10^seq(10,-2,length=100)
  
  X <- model.matrix(fml, data=aug_data)
  myfit <- glmnet(X, aug_data$tau_pred, alpha=1, lambda=grid)
  cv.out <- cv.glmnet(X, aug_data$tau_pred, alpha=1)
  bestlam <- cv.out$lambda.min
  
  sim_dat_mod_matrix <- model.matrix(fml, data=sim_dat_mod)
  pred_tau_lasso <- predict(myfit, s=bestlam,
                            newx=sim_dat_mod_matrix)
  
  return(pred_tau_lasso)
}


## Augmented dataset
XLearner_aug <- function(K, sim_dat, covars) {
  
  aug_data <- c()
  
  for (k in 1:K) {
    #fit model to study
    df <- filter(sim_dat, S==k)
    fit <- X_RF(feat=select(df, all_of(covars)), tr=df$W, yobs=df$Y, nthread=3)
    
    #apply model to all data
    df_mod <- sim_dat %>%
      mutate(tau_pred = EstimateCate(fit, select(sim_dat, all_of(covars))),
             Model_Site=k)
    aug_data <- bind_rows(aug_data, df_mod)
    
    rm(fit)
  }
  
  aug_data$Model_Site <- factor(aug_data$Model_Site)
  return(aug_data)
}

SLearner_aug <- function(K, sim_dat, covars) {
  
  aug_data <- c()
  
  for (k in 1:K) {
    #fit model to study
    df <- filter(sim_dat, S==k)
    fit <- S_RF(feat=select(df, all_of(covars)), tr=df$W, yobs=df$Y, nthread=3)
    
    #apply model to all data
    df_mod <- sim_dat %>%
      mutate(tau_pred = EstimateCate(fit, select(sim_dat, all_of(covars))),
             Model_Site=k)
    aug_data <- bind_rows(aug_data, df_mod)
    
    rm(fit)
  }
  
  aug_data$Model_Site <- factor(aug_data$Model_Site)
  return(aug_data)
}

CausalForest_aug <- function(K, sim_dat, covars, honesty) {
  
  aug_data <- c()
  
  for (k in 1:K) {
    #fit model to study
    df <- filter(sim_dat, S==k)
    fit <- grf::causal_forest(X=select(df,all_of(covars)), Y=df$Y, W=df$W, num.threads=3, honesty=honesty, num.trees=1000)
    
    #apply model to all data
    #first get OOB predictions from same site (training data)
    df_mod1 <- df %>%
      mutate(tau_pred = predict(fit)$predictions,
             Model_Site = k)
    #then predict on rest of the sites (testing data)
    other_dat <- filter(sim_dat, S!=k)
    df_mod2 <- other_dat %>%
      mutate(tau_pred = predict(fit, select(other_dat, all_of(covars)))$predictions,
             Model_Site = k)
    
    rm(fit)
    df_mod <- rbind(df_mod1, df_mod2)
    aug_data <- bind_rows(aug_data, df_mod)
  }
  
  aug_data$Model_Site <- factor(aug_data$Model_Site)
  return(aug_data)
}


## Main function
tan_preds <- function(K, sim_dat, covars, method, honesty=F) {
  
  sim_dat_mod <- mutate(sim_dat, 
                        Model_Site = S,
                        Model_Site = factor(Model_Site))
  
  if (method=="xlearner") {
    aug_data <- XLearner_aug(K, sim_dat, covars)
  }
  
  if (method=="causalforest") {
    aug_data <- CausalForest_aug(K, sim_dat, covars, honesty)
  }
  
  if (method=="slearner") {
    aug_data <- SLearner_aug(K, sim_dat, covars)
  }
  
  #ensemble tree
  et_fit <- EnsemTreeAll(aug_data, "Model_Site", covars)
  tan_pred_tau_et <- predict(et_fit,sim_dat_mod)
  mse_tree <- mean((tan_pred_tau_et - sim_dat$tau)^2)
  
  rm(list = c("et_fit", "tan_pred_tau_et"))
  
  #ensemble forest
  ef_fit <- EnsemForestAll(aug_data, "Model_Site", covars)
  tan_pred_tau_ef <- predict(ef_fit,sim_dat_mod,num.threads=3)$predictions
  mse_forest <- mean((tan_pred_tau_ef - sim_dat$tau)^2)
  
  rm(list = c("ef_fit", "tan_pred_tau_ef"))
  
  #ensemble lasso
  pred_tau_lasso <- EnsemLassoAll(aug_data, sim_dat_mod, "Model_Site", covars)
  mse_lasso <- mean((pred_tau_lasso - sim_dat$tau)^2)
  
  rm(pred_tau_lasso)
  
  return(c(mse_tree=mse_tree, mse_forest=mse_forest, mse_lasso=mse_lasso))
}
