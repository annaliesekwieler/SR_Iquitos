# This file creates the multivariate prior for the baseline phase based on the marginal priors
# of the parameters and historical force of infection

#### load libraries ####

library(dplyr)
library(mvtnorm)
library(fitdistrplus)

source("functions.R")

new_prior_df = data.frame(qu=NA, au=NA, tau=NA, d=NA, gu=NA, lambda = NA, 
                          n=NA, X=NA, c=NA,  b=NA)


while(nrow(new_prior_df) < 100){
  print(nrow(new_prior_df))
  N=10000
  qu_true_candidates = rgamma(N, 2, 4)
  overall_time_between_bites = rgamma(N, 9, 2)
  time_to_bite = runif(N, overall_time_between_bites * .25, overall_time_between_bites * .75)
  au_true_candidates = 1/time_to_bite
  time_to_digest = overall_time_between_bites - time_to_bite
  d_true_candidates = 1 / time_to_digest
  tau_true_candidates = rgamma(N, 2*3.3, 4)
  gu_true_candidates = rbeta(N, .1, 18)
  lambda_true_candidates = rgamma(N, 1, 10)
  n_true_candidates = rgamma(N, 40, 40/14)
  X_true_candidates = rbeta(N, 1, 7)
  c_true_candidates = rbeta(N, 1, 1)
  b_true_candidates = rbeta(N, 1, 1)
  
  prior_df = data.frame(qu=qu_true_candidates, au=au_true_candidates, tau=tau_true_candidates,
                        d=d_true_candidates, gu=gu_true_candidates, lambda = lambda_true_candidates, 
                        n=n_true_candidates, X=X_true_candidates, c=c_true_candidates,
                        b=b_true_candidates)
  
  prior_df$foi_control = rep(NA, nrow(prior_df))
  for(i in 1:nrow(prior_df)){
    prior_df$foi_control[i] = force_of_infection_untreated(prior_df$au[i], prior_df$qu[i], prior_df$tau[i],
                                              prior_df$d[i], rho = 1, q_mult = 1, a_mult = 1, C = 0,
                                              prior_df$gu[i], g_mult = 0, prior_df$b[i], prior_df$lambda[i],
                                              prior_df$c[i],prior_df$X[i],prior_df$n[i])
    
  }
  
  to_keep = prior_df$foi > 0 & prior_df$foi < 0.33/365
  if( sum(to_keep== T) != 0){
    new_prior_df = rbind(new_prior_df, prior_df[to_keep, 1:10])
  }
  if(any(is.na(new_prior_df))){
    new_prior_df = new_prior_df[-which(is.na(new_prior_df)),]
  }
  
}

new_prior_df_og = new_prior_df


new_prior_df = new_prior_df_og

# get mean and var-covar matrix
new_prior_df = new_prior_df[-1,]


#### to be added later
# add in the other parameters that weren't involved in the calculation
new_prior_df$phi = rexp(nrow(new_prior_df), 1)
new_prior_df$catch_prop = rbeta(nrow(new_prior_df), 1, 1)



# convert to log and inverse sigmoid scale to do inference
new_prior_df_log = log(new_prior_df)
new_prior_df_inv_sig = inverse_sigmoid(new_prior_df)

# select the correct transformation for each parameter
new_prior_df_trans = cbind(new_prior_df_log[,c("qu", "au", "tau", "d", "lambda", "n", "phi")],
                           new_prior_df_inv_sig[c("gu", "X", "c", "b", "catch_prop")])

# re-order the parameters
new_prior_df_trans = new_prior_df_trans %>%
  dplyr::select(qu, au, tau, d, gu,
                lambda, catch_prop, phi,
                n, X, c, b)

# fit a multivariate normal distribution
mu_mvn = apply(new_prior_df_trans, 2, mean)
sigma_mvn = var(new_prior_df_trans)

# save
save(mu_mvn, sigma_mvn, file="../results/MultivariatePrior.RData")


#### For running altogether or only on intervention data ####

new_prior_df$rho =rbeta(nrow(new_prior_df), 1, 1)
new_prior_df$a_mult = runif(nrow(new_prior_df), 0.5, 1.5)
new_prior_df$g_mult = runif(nrow(new_prior_df), 0.5, 1.5)
new_prior_df$q_mult = runif(nrow(new_prior_df), 0.5, 1.5)

new_prior_df = new_prior_df %>%
  dplyr::select(qu, au, tau, d, gu,
                lambda, catch_prop, phi,
                n, X, c, b, rho, a_mult, g_mult, q_mult)


new_prior_df_log = log(new_prior_df)
new_prior_df_inv_sig = inverse_sigmoid(new_prior_df)
new_prior_df_transformed= cbind(new_prior_df_log[,1:4], new_prior_df_inv_sig[,5],
                          new_prior_df_log[,6], new_prior_df_inv_sig[,7],
                          new_prior_df_log[,8:9], new_prior_df_inv_sig[,c(10, 11, 12, 13)],
                          new_prior_df_log[,14:16])



mu_mvn = apply(new_prior_df_transformed, 2, mean)
sigma_mvn = var(new_prior_df_transformed)

save(mu_mvn, sigma_mvn, file="../results/MultivariatePrior_All.RData")
