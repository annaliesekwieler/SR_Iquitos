library(dplyr)
library(scales)
library(mvtnorm)
library(fitdistrplus)

# load the data sets and parameter values used to make the simulation
load("../../results/Simulation/Simulation_InterventionOnly_NoParity_1.RData")

# thin
result_thinned = list()
for(i in 1:length(result)){
  result_temp = result[[i]][[1]]
  
  result_thinned[[i]] = rbind(result_temp[seq(1, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(2, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(3, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),])
}

result_all = result_thinned



load("../../results/Simulation/Simulation_InterventionOnly_NoParity_2.RData")
result_thinned = list()
for(i in 1:length(result)){
  result_temp = result[[i]][[1]]
  
  result_thinned[[i]] = rbind(result_temp[seq(1, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(2, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(3, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),])
}
result_all = c(result_all, result_thinned)


load("../../results/Simulation/Simulation_InterventionOnly_NoParity_3.RData")
result_thinned = list()
for(i in 1:length(result)){
  result_temp = result[[i]][[1]]
  
  result_thinned[[i]] = rbind(result_temp[seq(1, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(2, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(3, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),])
}
result_all = c(result_all, result_thinned)


load("../../results/Simulation/Simulation_InterventionOnly_NoParity_4.RData")
result_thinned = list()
for(i in 1:length(result)){
  result_temp = result[[i]][[1]]
  
  result_thinned[[i]] = rbind(result_temp[seq(1, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(2, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(3, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),])
}
result_all = c(result_all, result_thinned)


load("../../results/Simulation/Simulation_InterventionOnly_NoParity_5.RData")
result_thinned = list()
for(i in 1:length(result)){
  result_temp = result[[i]][[1]]
  
  result_thinned[[i]] = rbind(result_temp[seq(1, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(2, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),],
                              result_temp[seq(3, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 300),])
}
result_all = c(result_all, result_thinned)

setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/Simulation")
load("100SimulatedDatasets.RData")
result = result_all

result_length = length(result)

combos_to_try = combos_to_try %>%
  dplyr::select(qu_true, au_true, tau_true, d_true, gu_true, lambda_true,
                catch_prop_true, phi_true, n_true, X_true, c_true, b_true, rho_true,
                a.mult_true, g.mult_true, q.mult_true)
# combos_to_try = combos_to_try[c(61:100),]
#### Comparison of posterior to true value ####

# non-rounded
# determine what proportion of the posterior contained the true value for each parameter

prop_true_value = rep(0,16)
for(par in 1:16){
  for(i in 1:result_length){
    if (min(result[[i]][,par]) <= combos_to_try[i,par] & 
        max(result[[i]][,par]) >= combos_to_try[i,par]){
      prop_true_value[par] = prop_true_value[par]+1
    }
  }
}
prop_true_value = prop_true_value / result_length

# print the proportion correct for each data set
prop_true_by_ds = rep(NA, 100)
for(ds in 1:length(result)){
  num_correct = 0
  for(par in 1:16){
    true_val = combos_to_try[ds,par]
    post = result[[ds]][,par]
    if(min(post) < true_val & max(post) > true_val){
      num_correct = num_correct + 1
    }
  }
  prop_true_by_ds[ds] = num_correct
}




#### correlation between true and median inferred value for each parameter ####
corrs = corrs_corr = rep(NA, 16)
for(par in 1:16){
  par_vals_inferred = par_vals_true = rep(NA, 40)
  for(i in 1:result_length){
    par_vals_inferred[i] = median(result[[i]][,par])
    par_vals_true[i] = combos_to_try[i,par]
  }
  corrs[par] = cor(par_vals_inferred, par_vals_true)
  mu_x = mean(par_vals_inferred)
  mu_y = mean(par_vals_true)
  var_x = var(par_vals_inferred)
  var_y = var(par_vals_true)
  
  corrs_corr[par] = (2 * corrs[par] * sqrt(var_x) * sqrt(var_y)) / (var_x + var_y + (mu_x - mu_y)^2)
}


# directionality


a.mult_correct_direction = g.mult_correct_direction = q.mult_correct_direction = 0


for(i in 1:length(result)){
  a.mult_true = combos_to_try[i,"a.mult_true"]
  g.mult_true = combos_to_try[i,"g.mult_true"]
  q.mult_true = combos_to_try[i, "q.mult_true"]
  if((min(result[[i]]$a.mult) < a.mult_true) & (max(result[[i]]$a.mult) >
                                                a.mult_true)){
    a.mult_correct_direction = a.mult_correct_direction + 1
  }
  if((min(result[[i]]$g.mult) < g.mult_true) & (max(result[[i]]$g.mult) >
                                                g.mult_true)){
    g.mult_correct_direction = g.mult_correct_direction + 1
  }
  
  if((min(result[[i]]$q.mult) < q.mult_true) & (max(result[[i]]$q.mult) >
                                                q.mult_true)){
    q.mult_correct_direction = q.mult_correct_direction + 1
  }

}



a.mult_pe_correct_direction = g.mult_pe_correct_direction =
  q.mult_pe_correct_direction = 0


for(i in 1:length(result)){
  a.mult_true = median(combos_to_try[i,"a.mult_true"])
  g.mult_true = median(combos_to_try[i,"g.mult_true"])
  q.mult_true = median(combos_to_try[i,"q.mult_true"])
  if(sign(a.mult_true - 1) == sign(median(result[[i]]$a.mult) - 1)){
    a.mult_pe_correct_direction = a.mult_pe_correct_direction + 1
  }
  if(sign(g.mult_true - 1) == sign(median(result[[i]]$g.mult) - 1)){
    g.mult_pe_correct_direction = g.mult_pe_correct_direction + 1
  }
  
  if(sign(q.mult_true - 1) == sign(median(result[[i]]$q.mult) - 1)){
    q.mult_pe_correct_direction = q.mult_pe_correct_direction + 1
  }

}

a.mult_pe_correct_direction / 100
g.mult_pe_correct_direction / 100
q.mult_pe_correct_direction / 100


#### visualize ####

# draw from priors
sigmoid = function(x){
  1/(1+exp(-x))
}


inverse_sigmoid = function(x){
  log(x / (1-x))
}

setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/")
load("MultivariatePrior_All.RData")

prior_samples = data.frame(rmvnorm(10000, mu_mvn, sigma_mvn))
df_names = names(prior_samples)
prior_samples_exp = exp(prior_samples)
prior_samples_sig = sigmoid(prior_samples)

prior_samples = cbind(prior_samples_exp[,1:4], prior_samples_sig[,5],
                      prior_samples_exp[,6],
                      prior_samples_sig[,7], prior_samples_exp[,8:9],
                      prior_samples_sig[,10:13], prior_samples_exp[,14:16])

names(prior_samples) = names(result[[1]])

# switch y and x axes

png("../../results/ManuscriptFigures/Simulation_SR_Effects.png", height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)

par(mfrow=c(2,2), mar=c(3,4,4,5), cex=1.5)

lows = c(0, 0, 0, 0,
         0, 0, 0, 0,
         8, 0, 0, 0,
         0, 0, 0, 0)
tops = c(3, 3, 6, 2.5,
         .8, 2.5, 1, 5
         ,20,1, 1, 1,
         1, 2, 3, 3)
dens_mult = c(1, 1, 6, 1,
              .05, .1, .5, 4,
              100, .05, .4, .4,
              .5, 1, 1.8, 1.8)

param_names = c("Baseline \nexit rate", "Baseline \nbiting rate", "Entrance rate", "Digestion rate",
                "Mortality Rate",  "Emergence", "Catch Proportion", "Dispersion",
                 "EIP", "Prevalence", "Transmission \nprob 1",  "Transmission \nprob 2",
                "Repellency", "SR effect biting", "SR effect mortality", "SR effect exiting")

for(par in 13:16){
  
 
  
  plot(0, type="n", xlim=c(lows[par], tops[par]), ylim=c(lows[par], tops[par]),
       xlab=xlab,
       ylab=ylab,
       main=param_names[par])
  
  if(par == 13){
    mtext("True", side = 2, line= 2, cex = 1.5)
  }
  
  if(par == 13){
    mtext("Inferred", side = 1, line= 3, cex = 1.5)
  }
 
  for(i in 1:result_length){
    points(median(result[[i]][,par]), combos_to_try[i,par], col="red",
           pch=19)
  }
  abline(0, 1)
  # add bars for posterior distribution
  for(i in 1:result_length){
    segments(y0=combos_to_try[i,par],y1=combos_to_try[i,par], 
             x0 =min(result[[i]][,par]),
             x1 = max(result[[i]][,par]), col=alpha("black", .5))
  }
  dens_ob = density(prior_samples[,par])
  lines(dens_ob$x, dens_ob$y * dens_mult[par])
}

dev.off()



png("../../results/ManuscriptFigures/Supplement/Simulation_NonSR.png", height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)

par(mfrow=c(3, 4), cex=1.2, mar=c(2.5, 3, 3, 1))

lows = c(0, 0, 0, 0,
         0, 0, 0, 0,
         8, 0, 0, 0,
         0, 0, 0, 0)
tops = c(3, 3, 6, 2.5,
         .8, 2.5, 1, 5
         ,20,1, 1, 1,
         1, 2, 3, 3)
dens_mult = c(1, 1, 6, 1,
              .05, .1, .5, 4,
              100, .05, .4, .4,
              .5, 1, 1.8, 1.8)

param_names = c("Baseline \nexit rate", "Baseline \nbiting rate", "Entrance rate", "Digestion rate",
                "Mortality Rate",  "Emergence", "Catch Proportion", "Dispersion",
                "EIP", "Prevalence", "Transmission \nProb 1",  "Transmission \nProb 2",
                "Repellency", "SR effect biting", "SR effect mortality", "SR effect exiting")

for(par in 1:12){
  
  
  
  plot(0, type="n", xlim=c(lows[par], tops[par]), ylim=c(lows[par], tops[par]),
       xlab=xlab,
       ylab=ylab,
       main=param_names[par])
  
  if(par == 1){
    mtext("True", side = 2, line= 2, cex = 1.5)
  }
  
  if(par == 1){
    mtext("Inferred", side = 1, line= 2, cex = 1.5)
  }
  
  for(i in 1:result_length){
    points(median(result[[i]][,par]), combos_to_try[i,par], col="red",
           pch=19)
  }
  abline(0, 1)
  # add bars for posterior distribution
  for(i in 1:result_length){
    segments(y0=combos_to_try[i,par],y1=combos_to_try[i,par], 
             x0 =min(result[[i]][,par]),
             x1 = max(result[[i]][,par]), col=alpha("black", .5))
  }
  dens_ob = density(prior_samples[,par])
  lines(dens_ob$x, dens_ob$y * dens_mult[par])
}

dev.off()
#### load functions and trial data ####

source("../functions.R")

# load trial data 
bloodmeal_intervention = read.csv("../../data/bloodmeal_intervention.csv")
bloodmeal_baseline = read.csv("../../data/bloodmeal_baseline.csv")

# take total full to be the sum of hald and full
bloodmeal_intervention = bloodmeal_intervention %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

bloodmeal_baseline = bloodmeal_baseline %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

abundance_intervention_control = read.csv("../../data/abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("../../data/abundance_intervention_treatment.csv")


parity_intervention = read.csv("../../data/parity_intervention.csv")

epi = read.csv("../../data/epi.csv")

load("../../data/treatedclusters.RData")

C=.44
p=5.2

#### check to ensure that all the quantities matched ####

prop_bloodfed_control = prop_bloodfed_control_lower = prop_bloodfed_control_upper =
  prop_bloodfed_treatment =prop_bloodfed_treatment_lower =prop_bloodfed_treatment_upper =
  true_prop_bloodfed_control = true_prop_bloodfed_treatment =
  rep(NA, length(result))


for(i in 1:length(result)){
  print(i)
  res = result[[i]]
  blood_control = blood_trt = rep(NA, nrow(res))
  for(j in 1:nrow(res)){
    res_row = res[j,]
    blood_control[j] = with(res_row, prop_bloodfed_untreated(au, qu, tau, d, rho,
                                                         q.mult, a.mult, C=0))
    
    blood_trt[j] = with(res_row, prop_bloodfed_treated(au, qu, tau, d, rho,
                                                             q.mult, a.mult, C))
  }
  
  prop_bloodfed_control[i] = median(blood_control)
  prop_bloodfed_control_lower[i] = quantile(blood_control, .025)
  prop_bloodfed_control_upper[i] = quantile(blood_control, .975)
  
  prop_bloodfed_treatment[i] = median(blood_trt)
  prop_bloodfed_treatment_lower[i] = quantile(blood_trt, .025)
  prop_bloodfed_treatment_upper[i] = quantile(blood_trt, .975)
  
  true_prop_bloodfed_control[i] = with(combos_to_try[i,], prop_bloodfed_untreated(au_true, qu_true, tau_true, d_true, rho_true,
                                                                        q.mult_true, a.mult_true, C=0))
  true_prop_bloodfed_treatment[i] = with(combos_to_try[i,], prop_bloodfed_treated(au_true, qu_true, tau_true, d_true, rho_true,
                                                                                  q.mult_true, a.mult_true, C))
  
}

plot(prop_bloodfed_control ~ true_prop_bloodfed_control, pch=19,
     xlab="True", ylab="Predicted", main= "Proportion Bloodfed, Control")


#### Sensitivity Analysis ####

source("../functions.R")

true_total_effect = rep(NA, 100)
inferred_total_effect = list()

C = .44
p = 5.2


i= sample(1:length(result), size = 1)

  res = result[[i]]
  
  true_params = combos_to_try[i,]
  
  true_foi_c = force_of_infection_untreated(true_params$au_true, true_params$qu_true, true_params$tau_true,
                                            true_params$d_true, 0, 1,
                                            1, C=0, true_params$gu_true,1, true_params$b_true,
                                            true_params$lambda_true,
                                            true_params$c_true, true_params$X_true, 
                                            true_params$n_true)
  
  true_foi_u = force_of_infection_untreated(true_params$au_true, true_params$qu_true,
                                            true_params$tau_true, true_params$d_true,
                                            true_params$rho_true, true_params$q.mult_true,
                                            true_params$a.mult_true, C = C, true_params$gu_true,
                                            true_params$g.mult_true, true_params$b_true,
                                            true_params$lambda_true,
                                            true_params$c_true, true_params$X_true,
                                            true_params$n_true)
  
  true_foi_t = force_of_infection_treated(true_params$au_true, true_params$qu_true,
                                            true_params$tau_true, true_params$d_true,
                                            true_params$rho_true, true_params$q.mult_true,
                                            true_params$a.mult_true, C = C, true_params$gu_true,
                                            true_params$g.mult_true, true_params$b_true,
                                            true_params$lambda_true,
                                            true_params$c_true, true_params$X_true,
                                            true_params$n_true)
  
  inferred_foi_c = inferred_foi_u = inferred_foi_t = rep(NA, nrow(res))
  
  for(j in 1:nrow(res)){
    inferred_foi_c[j] =  force_of_infection_untreated(res[j,]$au, res[j,]$qu, res[j,]$tau,
                                                      res[j,]$d, 0, 1,
                                                      1, C=0, res[j,]$gu,1, res[j,]$b,
                                                      res[j,]$lambda,
                                                      res[j,]$c, res[j,]$X, 
                                                      res[j,]$n)
    
    inferred_foi_u[j] =  force_of_infection_untreated(res[i,]$au, res[i,]$qu,
                                                      res[i,]$tau, res[i,]$d,
                                                      res[i,]$rho, res[i,]$q.mult,
                                                      res[i,]$a.mult, C = C, res[i,]$gu,
                                                      res[i,]$g.mult, res[i,]$b,
                                                      res[i,]$lambda,
                                                      res[i,]$c, res[i,]$X,
                                                      res[i,]$n)
    
    inferred_foi_t[j] =  force_of_infection_treated(res[i,]$au, res[i,]$qu,
                                                      res[i,]$tau, res[i,]$d,
                                                      res[i,]$rho, res[i,]$q.mult,
                                                      res[i,]$a.mult, C = C, res[i,]$gu,
                                                      res[i,]$g.mult, res[i,]$b,
                                                      res[i,]$lambda,
                                                      res[i,]$c, res[i,]$X,
                                                      res[i,]$n)
  }
  
  inferred_total_effect =  1 - (inferred_foi_t / inferred_foi_c)
  
  true_total_effect[i] =  1 - (true_foi_t / true_foi_c)
  
  
  plot(density(inferred_total_effect[[i]]), main = "Inferred Total Effect for Simulated Data",
       xlab = "Total Effect")
  polygon(density(inferred_total_effect[[i]]), col="grey")
  abline(v = true_total_effect[i], lty="dashed", col="red")
  
  
plot(inferred_total_effect[[i]] ~ res$au, xlab = "Mortality",
     ylab="Total Effect", pch=19)
abline(h = true_total_effect, lty="dashed", col = "red", lwd=3)


summary(lm(inferred_total_effect ~ res$q.mult))

plot(inferred_total_effect ~ res$gu, xlab = "Mortality",
     ylab="Total Effect", pch=19)
abline(h = true_total_effect, lty="dashed", col = "red", lwd=3)

