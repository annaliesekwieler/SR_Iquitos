# This script analyzes the 5 data sets for leaving out each data 
# types

#### load packages ####

library(mvtnorm)
library(tidyverse)

#### load results from leaving out each data type ####

load("../../results/TrialData/Intervention_Only.RData")

result_thinned = rbind(result[seq(1, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                                 by = 300),],
                              result[seq(2, nrow(result), by=3), ][seq(1,
                                                                       nrow(result)/3,
                                                                                 by = 300),],
                              result[seq(3, nrow(result), by=3), ][seq(1,
                                                                       nrow(result)/3,
                                                                              by = 300),])

result_all = result_thinned
head(result_thinned)


load("../../results/TrialData/InterventionOnly_NoParity.RData")
result_thinned = rbind(result[seq(1, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(2, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(3, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),])
result_no_parity = result_thinned
head(result_thinned)

load("../../results/TrialData/InterventionOnly_NoAspiration.RData")
result_thinned = rbind(result[seq(1, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(2, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(3, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),])
result_no_aspiration = result_thinned


load("../../results/TrialData/InterventionOnly_NoBloodmeal.RData")
result_thinned = rbind(result[seq(1, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(2, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(3, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),])
result_no_blood = result_thinned

load("../../results/TrialData/InterventionOnly_NoSC.RData")
result_thinned = rbind(result[seq(1, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(2, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),],
                       result[seq(3, nrow(result), by=3), ][seq(1,
                                                                nrow(result)/3,
                                                                by = 200),])
result_no_epi = result_thinned

#### load various priors ####

source("../../code/functions.R")

load("../../results/MultivariatePrior_All.RData")

# samples from the prior

prior_samples = data.frame(rmvnorm(10000, mu_mvn, sigma_mvn))
df_names = names(prior_samples)
prior_samples_exp = exp(prior_samples)
prior_samples_sig = sigmoid(prior_samples)

prior_samples = cbind(prior_samples_exp[,1:4], prior_samples_sig[,5],
                      prior_samples_exp[,6], prior_samples_sig[,7],
                      prior_samples_exp[,8:9], prior_samples_sig[,10:13],
                      prior_samples_exp[,14:16])
names(prior_samples) =names(result)


#### Parameter Estimates ####

png("../../results/ManuscriptFigures/Supplement/LeavingOutParamEstimates.png",height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)

par(mfrow=c(3, 4), mar=c(2, 2, 4, 1), cex = 1.2)

plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1.5), 
     ylab="Parameter Estimate",
     main="Exit Rate, Untreated",
     xaxt="n")
boxplot(result_all$qu, at = 1, add=T,col="red")

boxplot(result_no_parity$qu, at = 2, add=T,col="red")

boxplot(result_no_aspiration$qu, at = 3, add=T,col="red")

boxplot(result_no_blood$qu, at = 4, add=T,col="red")

boxplot(result_no_epi$qu, at = 5, add=T,col="red")

boxplot(prior_samples$qu, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)

plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1.5),
     ylab="Parameter Estimate",
     main="Biting Rate, \nUntreated",
     xaxt="n")
boxplot(result_all$au, at = 1, add=T,col="red")

boxplot(result_no_parity$au, at = 2, add=T,col="red")

boxplot(result_no_aspiration$au, at = 3, add=T,col="red")

boxplot(result_no_blood$au, at = 4, add=T,col="red")

boxplot(result_no_epi$au, at = 5, add=T,col="red")

boxplot(prior_samples$au, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 5),
     ylab="Parameter Estimate",
     main="Entrance Rate, \nUntreated",
     xaxt="n")
boxplot(result_all$tau, at = 1, add=T,col="red")

boxplot(result_no_parity$tau, at = 2, add=T,col="red")

boxplot(result_no_aspiration$tau, at = 3, add=T,col="red")

boxplot(result_no_blood$tau, at = 4, add=T,col="red")

boxplot(result_no_epi$tau, at = 5, add=T,col="red")

boxplot(prior_samples$tau, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1.5), 
     ylab="Parameter Estimate",
     main="Digestion Rate",
     xaxt="n")
boxplot(result_all$d, at = 1, add=T,col="red")

boxplot(result_no_parity$d, at = 2, add=T,col="red")

boxplot(result_no_aspiration$d, at = 3, add=T,col="red")

boxplot(result_no_blood$d, at = 4, add=T,col="red")

boxplot(result_no_epi$d, at = 5, add=T,col="red")

boxplot(prior_samples$d, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, .2),
     ylab="Parameter Estimate",
     main="Mortality, Untreated",
     xaxt="n")
boxplot(result_all$gu, at = 1, add=T,col="red")

boxplot(result_no_parity$gu, at = 2, add=T,col="red")

boxplot(result_no_aspiration$gu, at = 3, add=T,col="red")

boxplot(result_no_blood$gu, at = 4, add=T,col="red")

boxplot(result_no_epi$gu, at = 5, add=T,col="red")

boxplot(prior_samples$gu, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, .25),
     ylab="Parameter Estimate",
     main="Emergence Rate",
     xaxt="n")
boxplot(result_all$lambda, at = 1, add=T,col="red")

boxplot(result_no_parity$lambda, at = 2, add=T,col="red")

boxplot(result_no_aspiration$lambda, at = 3, add=T,col="red")

boxplot(result_no_blood$lambda, at = 4, add=T,col="red")

boxplot(result_no_epi$lambda, at = 5, add=T,col="red")

boxplot(prior_samples$lambda, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1), 
     ylab="Parameter Estimate",
     main="Catch Proportion",
     xaxt="n")
boxplot(result_all$catch_prop, at = 1, add=T,col="red")

boxplot(result_no_parity$catch_prop, at = 2, add=T,col="red")

boxplot(result_no_aspiration$catch_prop, at = 3, add=T,col="red")

boxplot(result_no_blood$catch_prop, at = 4, add=T,col="red")

boxplot(result_no_epi$catch_prop, at = 5, add=T,col="red")

boxplot(prior_samples$catch_prop, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 3), 
     ylab="Parameter Estimate",
     main="Dispersion, Aspiration",
     xaxt="n")
boxplot(result_all$phi, at = 1, add=T,col="red")

boxplot(result_no_parity$phi, at = 2, add=T,col="red")

boxplot(result_no_aspiration$phi, at = 3, add=T,col="red")

boxplot(result_no_blood$phi, at = 4, add=T,col="red")

boxplot(result_no_epi$phi, at = 5, add=T,col="red")

boxplot(prior_samples$phi, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(5, 20),
     ylab="Parameter Estimate",
     main="EIP",
     xaxt="n")
boxplot(result_all$n, at = 1, add=T,col="red")

boxplot(result_no_parity$n, at = 2, add=T,col="red")

boxplot(result_no_aspiration$n, at = 3, add=T,col="red")

boxplot(result_no_blood$n, at = 4, add=T,col="red")

boxplot(result_no_epi$n, at = 5, add=T,col="red")

boxplot(prior_samples$n, at = 6, add=T)
axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 0.6),
     ylab="Parameter Estimate",
     main="Prevalence",
     xaxt="n")
boxplot(result_all$X, at = 1, add=T,col="red")

boxplot(result_no_parity$X, at = 2, add=T,col="red")

boxplot(result_no_aspiration$X, at = 3, add=T,col="red")

boxplot(result_no_blood$X, at = 4, add=T,col="red")

boxplot(result_no_epi$X, at = 5, add=T,col="red")

boxplot(prior_samples$X, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)

plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1), 
     ylab="Parameter Estimate",
     main="Transmission \nProbability (b)",
     xaxt="n")
boxplot(result_all$b, at = 1, add=T,col="red")

boxplot(result_no_parity$b, at = 2, add=T,col="red")

boxplot(result_no_aspiration$b, at = 3, add=T,col="red")

boxplot(result_no_blood$b, at = 4, add=T,col="red")

boxplot(result_no_epi$b, at = 5, add=T,col="red")

boxplot(prior_samples$b, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)

plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1),
     ylab="Parameter Estimate",
     main="Transmission \nProbability (c)",
     xaxt="n")
boxplot(result_all$c, at = 1, add=T,col="red")

boxplot(result_no_parity$c, at = 2, add=T,col="red")

boxplot(result_no_aspiration$c, at = 3, add=T,col="red")

boxplot(result_no_blood$c, at = 4, add=T,col="red")

boxplot(result_no_epi$c, at = 5, add=T,col="red")

boxplot(prior_samples$c, at = 6, add=T)

axis("N", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("A", side = 1, at = 3, las = 1)
axis("B", side = 1, at = 4, las = 1)
axis("E", side = 1, at = 5, las = 1)
axis("Pr", side = 1, at = 6, las = 1)

dev.off()

png("../../results/ManuscriptFigures/LeavingOut_SREffects.png",
    height=8.5, width=11, units = "in", pointsize = 12,
    res = 72)

par(mfrow=c(2,2), cex = 1.2, mar= c(4,4,4,2))

plot(-1, xlim=c(0.5, 6.5), ylim=c(0, 1), xlab="Left Out",
     ylab="Parameter Estimate",
     main="SR Effect, Repellency",
     xaxt="n")
boxplot(result_all$rho, at = 1, add=T,col="red")

boxplot(result_no_parity$rho, at = 2, add=T,col="red")

boxplot(result_no_aspiration$rho, at = 3, add=T,col="red")

boxplot(result_no_blood$rho, at = 4, add=T,col="red")

boxplot(result_no_epi$rho, at = 5, add=T,col="red")

boxplot(prior_samples$rho, at = 6, add=T)

axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
axis("Prior", side = 1, at = 6, las = 1)


plot(-1, xlim=c(0.5, 6.5), ylim=c(0.2, 2), xlab="Left Out",
     ylab="Parameter Estimate",
     main="SR Effect, Biting",
     xaxt="n")
boxplot(result_all$a_mult, at = 1, add=T,col="red")

boxplot(result_no_parity$a_mult, at = 2, add=T,col="red")

boxplot(result_no_aspiration$a_mult, at = 3, add=T,col="red")

boxplot(result_no_blood$a_mult, at = 4, add=T,col="red")

boxplot(result_no_epi$a_mult, at = 5, add=T,col="red")

boxplot(prior_samples$a_mult, at = 6, add=T)

axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
axis("Prior", side = 1, at = 6, las = 1)

abline(h = 1, lty="dashed")



plot(-1, xlim=c(0.5, 6.5), ylim=c(0.2, 2), xlab="Left Out",
     ylab="Parameter Estimate",
     main="SR Effect, Mortality",
     xaxt="n")
boxplot(result_all$g_mult, at = 1, add=T,col="red")

boxplot(result_no_parity$g_mult, at = 2, add=T,col="red")

boxplot(result_no_aspiration$g_mult, at = 3, add=T,col="red")

boxplot(result_no_blood$g_mult, at = 4, add=T,col="red")

boxplot(result_no_epi$g_mult, at = 5, add=T,col="red")

boxplot(prior_samples$g_mult, at = 6, add=T)

axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
axis("Prior", side = 1, at = 6, las = 1)
abline(h = 1, lty="dashed")



plot(-1, xlim=c(0.5, 6.5), ylim=c(0.2, 2), xlab="Left Out",
     ylab="Parameter Estimate",
     main="SR Effect, Exiting",
     xaxt="n")
boxplot(result_all$q_mult, at = 1, add=T,col="red")

boxplot(result_no_parity$q_mult, at = 2, add=T,col="red")

boxplot(result_no_aspiration$q_mult, at = 3, add=T,col="red")

boxplot(result_no_blood$q_mult, at = 4, add=T,col="red")

boxplot(result_no_epi$q_mult, at = 5, add=T,col="red")

boxplot(prior_samples$q_mult, at = 6, add=T)

axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
axis("Prior", side = 1, at = 6, las = 1)
abline(h = 1, lty="dashed")


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

abundance_baseline = read.csv("../../data/abundance_baseline.csv")
abundance_intervention_control = read.csv("../../data/abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("../../data/abundance_intervention_treatment.csv")

parity_baseline = read.csv("../../data/parity_baseline.csv")
parity_intervention = read.csv("../../data/parity_intervention.csv")

epi = read.csv("../../data/epi.csv")

load("../../data/treatedclusters.RData")

C=.44
p=5.2


#### PPD ####


#### All data included ####


result = result_all

result$a_mult

prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau,
                                                                   d, rho, q_mult,
                                                                   a_mult, C = 0))
}

prop_fed_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_treated[i] =  with(result[i,],
                              prop_bloodfed_treated(au, qu, tau,
                                                    d, rho, q_mult,
                                                    a_mult, C))
}

f_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,], expected_catch_number_control(au, qu, tau, d, 
                                                                gu, lambda,
                                                                catch_prop))
}

f_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_treatment[i] = with(result[i,], expected_catch_number_treated(au, qu, tau, d, rho, q_mult, a_mult,
                                                                  C, gu, g_mult, lambda, catch_prop))
}

parity_rate_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,], parity(au, qu,
                                                   tau, d,
                                                   rho, q_mult,
                                                   a_mult, C = 0,
                                                   gu, g_mult))
}

parity_rate_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_treatment[i] = with(result[i,], parity(au, qu,
                                                     tau, d,
                                                     rho, q_mult,
                                                     a_mult, C,
                                                     gu, g_mult))
}

foi_c = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_c[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C= 0,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

foi_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_t[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

prob_not_sc_control = prob_sc_control = list()
for(i in 1:nrow(result)){
  prob_not_sc_control [[i]] = exp(-foi_c[i] *
                                    epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  prob_sc_control[[i]] = 1 - prob_not_sc_control[[i]]
}

prob_not_sc_treatment = prob_sc_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_treatment [[i]] = exp(-foi_c[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_treatment[[i]] = 1 - prob_not_sc_treatment[[i]]
}


parity_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_control_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                              filter(cluster_trt == "C") %>%
                                              pull(total_caught)),
                                     prob = parity_rate_control[i])
}

parity_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_treatment_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                                filter(cluster_trt == "T") %>%
                                                pull(total_caught)),
                                       prob = parity_rate_treatment[i])
}

bloodmeal_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_control_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                 filter(cluster_trt == "C") %>%
                                                 pull(total_evaluated)),
                                        prob = prop_fed_untreated[i])
}

bloodmeal_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_treatment_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                   filter(cluster_trt == "T") %>%
                                                   pull(total_evaluated)),
                                          prob = prop_fed_treated[i])
}

aspiration_int_control_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_control_ppd[[i]] = rnbinom(nrow(abundance_intervention_control),
                                            mu=f_control[i],
                                            size = result$phi[i])
}

aspiration_int_treatment_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_treatment_ppd[[i]] = rnbinom(nrow(abundance_intervention_treatment),
                                              mu=f_treatment[i],
                                              size = result$phi[i])
}

sc_control_ppd = sc_treatment_ppd = list()
for(i in 1:nrow(result)){
  sc_control_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "C",]),
                               size = 1, prob=prob_sc_control[[i]])
  
  sc_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_treatment[[i]])
}

prop_parous_int_control_ppd_all = parity_int_control_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_caught))

prop_parous_int_treatment_ppd_all = parity_int_treatment_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_caught))

prop_bloodfed_int_control_ppd_all = bloodmeal_int_control_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_evaluated))


prop_bloodfed_int_treated_ppd_all = bloodmeal_int_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

mean_aspiration_int_control_ppd_all = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_control_ppd_all[i] = mean(aspiration_int_control_ppd[[i]])
}


mean_aspiration_int_treatment_ppd_all = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_treatment_ppd_all[i] = mean(aspiration_int_treatment_ppd[[i]])
}

prop_sc_control_all = prop_sc_treatment_all = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_sc_control_all[i] = mean(sc_control_ppd[[i]])
  prop_sc_treatment_all[i] = mean(sc_treatment_ppd[[i]])
}

#### Parity ####


result = result_no_parity

result$a_mult = result$alpha

prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau,
                                                  d, rho, q_mult,
                                                  a_mult, C = 0))
}

prop_fed_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_treated[i] =  with(result[i,],
                              prop_bloodfed_treated(au, qu, tau,
                                              d, rho, q_mult,
                                                    a_mult, C))
}

f_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,], expected_catch_number_control(au, qu, tau, d, 
                                                                gu, lambda,
                                                                catch_prop))
}

f_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_treatment[i] = with(result[i,], expected_catch_number_treated(au, qu, tau, d, rho, q_mult, a_mult,
                                                                  C, gu, g_mult, lambda, catch_prop))
}

parity_rate_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,], parity(au, qu,
                                                   tau, d,
                                                   rho, q_mult,
                                                   a_mult, C = 0,
                                                   gu, g_mult))
}

parity_rate_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_treatment[i] = with(result[i,], parity(au, qu,
                                                   tau, d,
                                                   rho, q_mult,
                                                   a_mult, C,
                                                   gu, g_mult))
}

foi_c = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_c[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C= 0,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

foi_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_t[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

prob_not_sc_control = prob_sc_control = list()
for(i in 1:nrow(result)){
  prob_not_sc_control [[i]] = exp(-foi_c[i] *
                                    epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  prob_sc_control[[i]] = 1 - prob_not_sc_control[[i]]
}

prob_not_sc_treatment = prob_sc_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_treatment [[i]] = exp(-foi_c[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_treatment[[i]] = 1 - prob_not_sc_treatment[[i]]
}


parity_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_control_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                              filter(cluster_trt == "C") %>%
                                              pull(total_caught)),
                                     prob = parity_rate_control[i])
}

parity_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_treatment_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                                filter(cluster_trt == "T") %>%
                                                pull(total_caught)),
                                       prob = parity_rate_treatment[i])
}

bloodmeal_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_control_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                 filter(cluster_trt == "C") %>%
                                                 pull(total_evaluated)),
                                        prob = prop_fed_untreated[i])
}

bloodmeal_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_treatment_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                   filter(cluster_trt == "T") %>%
                                                   pull(total_evaluated)),
                                          prob = prop_fed_treated[i])
}

aspiration_int_control_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_control_ppd[[i]] = rnbinom(nrow(abundance_intervention_control),
                                            mu=f_control[i],
                                            size = result$phi[i])
}

aspiration_int_treatment_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_treatment_ppd[[i]] = rnbinom(nrow(abundance_intervention_treatment),
                                              mu=f_treatment[i],
                                              size = result$phi[i])
}

sc_control_ppd = sc_treatment_ppd = list()
for(i in 1:nrow(result)){
  sc_control_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "C",]),
                               size = 1, prob=prob_sc_control[[i]])
  
  sc_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_treatment[[i]])
}

prop_parous_int_control_ppd_no_parity = parity_int_control_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_caught))

prop_parous_int_treatment_ppd_no_parity = parity_int_treatment_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_caught))

prop_bloodfed_int_control_ppd_no_parity = bloodmeal_int_control_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_evaluated))


prop_bloodfed_int_treated_ppd_no_parity = bloodmeal_int_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

mean_aspiration_int_control_ppd_no_parity = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_control_ppd_no_parity[i] = mean(aspiration_int_control_ppd[[i]])
}


mean_aspiration_int_treatment_ppd_no_parity = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_treatment_ppd_no_parity[i] = mean(aspiration_int_treatment_ppd[[i]])
}

prop_sc_control_no_parity = prop_sc_treatment_no_parity = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_sc_control_no_parity[i] = mean(sc_control_ppd[[i]])
  prop_sc_treatment_no_parity[i] = mean(sc_treatment_ppd[[i]])
}

#### Aspiration ####

result = result_no_aspiration

result$a_mult = result$alpha

prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau,
                                                                   d, rho, q_mult,
                                                                   a_mult, C = 0))
}

prop_fed_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_treated[i] =  with(result[i,],
                              prop_bloodfed_treated(au, qu, tau,
                                                    d, rho, q_mult,
                                                    a_mult, C))
}

f_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,], expected_catch_number_control(au, qu, tau, d, 
                                                                gu, lambda,
                                                                catch_prop))
}

f_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_treatment[i] = with(result[i,], expected_catch_number_treated(au, qu, tau, d, rho, q_mult, a_mult,
                                                                  C, gu, g_mult, lambda, catch_prop))
}

parity_rate_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,], parity(au, qu,
                                                   tau, d,
                                                   rho, q_mult,
                                                   a_mult, C = 0,
                                                   gu, g_mult))
}

parity_rate_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_treatment[i] = with(result[i,], parity(au, qu,
                                                     tau, d,
                                                     rho, q_mult,
                                                     a_mult, C,
                                                     gu, g_mult))
}

foi_c = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_c[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C= 0,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

foi_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_t[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

prob_not_sc_control = prob_sc_control = list()
for(i in 1:nrow(result)){
  prob_not_sc_control [[i]] = exp(-foi_c[i] *
                                    epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  prob_sc_control[[i]] = 1 - prob_not_sc_control[[i]]
}

prob_not_sc_treatment = prob_sc_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_treatment [[i]] = exp(-foi_c[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_treatment[[i]] = 1 - prob_not_sc_treatment[[i]]
}


parity_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_control_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                              filter(cluster_trt == "C") %>%
                                              pull(total_caught)),
                                     prob = parity_rate_control[i])
}

parity_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_treatment_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                                filter(cluster_trt == "T") %>%
                                                pull(total_caught)),
                                       prob = parity_rate_treatment[i])
}

bloodmeal_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_control_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                 filter(cluster_trt == "C") %>%
                                                 pull(total_evaluated)),
                                        prob = prop_fed_untreated[i])
}

bloodmeal_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_treatment_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                   filter(cluster_trt == "T") %>%
                                                   pull(total_evaluated)),
                                          prob = prop_fed_treated[i])
}

aspiration_int_control_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_control_ppd[[i]] = rnbinom(nrow(abundance_intervention_control),
                                            mu=f_control[i],
                                            size = result$phi[i])
}

aspiration_int_treatment_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_treatment_ppd[[i]] = rnbinom(nrow(abundance_intervention_treatment),
                                              mu=f_treatment[i],
                                              size = result$phi[i])
}

sc_control_ppd = sc_treatment_ppd = list()
for(i in 1:nrow(result)){
  sc_control_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "C",]),
                               size = 1, prob=prob_sc_control[[i]])
  
  sc_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_treatment[[i]])
}

prop_parous_int_control_ppd_no_aspiration = parity_int_control_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_caught))

prop_parous_int_treatment_ppd_no_aspiration = parity_int_treatment_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_caught))

prop_bloodfed_int_control_ppd_no_aspiration = bloodmeal_int_control_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_evaluated))


prop_bloodfed_int_treated_ppd_no_aspiration = bloodmeal_int_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

mean_aspiration_int_control_ppd_no_aspiration = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_control_ppd_no_aspiration[i] = mean(aspiration_int_control_ppd[[i]])
}


mean_aspiration_int_treatment_ppd_no_aspiration = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_treatment_ppd_no_aspiration[i] = mean(aspiration_int_treatment_ppd[[i]])
}

prop_sc_control_no_aspiration = prop_sc_treatment_no_aspiration = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_sc_control_no_aspiration[i] = mean(sc_control_ppd[[i]])
  prop_sc_treatment_no_aspiration[i] = mean(sc_treatment_ppd[[i]])
}

#### Bloodmeal ####

result = result_no_blood

result$a_mult = result$alpha

prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau,
                                                                   d, rho, q_mult,
                                                                   a_mult, C = 0))
}

prop_fed_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_treated[i] =  with(result[i,],
                              prop_bloodfed_treated(au, qu, tau,
                                                    d, rho, q_mult,
                                                    a_mult, C))
}

f_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,], expected_catch_number_control(au, qu, tau, d, 
                                                                gu, lambda,
                                                                catch_prop))
}

f_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_treatment[i] = with(result[i,], expected_catch_number_treated(au, qu, tau, d, rho, q_mult, a_mult,
                                                                  C, gu, g_mult, lambda, catch_prop))
}

parity_rate_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,], parity(au, qu,
                                                   tau, d,
                                                   rho, q_mult,
                                                   a_mult, C = 0,
                                                   gu, g_mult))
}

parity_rate_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_treatment[i] = with(result[i,], parity(au, qu,
                                                     tau, d,
                                                     rho, q_mult,
                                                     a_mult, C,
                                                     gu, g_mult))
}

foi_c = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_c[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C= 0,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

foi_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_t[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

prob_not_sc_control = prob_sc_control = list()
for(i in 1:nrow(result)){
  prob_not_sc_control [[i]] = exp(-foi_c[i] *
                                    epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  prob_sc_control[[i]] = 1 - prob_not_sc_control[[i]]
}

prob_not_sc_treatment = prob_sc_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_treatment [[i]] = exp(-foi_c[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_treatment[[i]] = 1 - prob_not_sc_treatment[[i]]
}


parity_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_control_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                              filter(cluster_trt == "C") %>%
                                              pull(total_caught)),
                                     prob = parity_rate_control[i])
}

parity_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_treatment_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                                filter(cluster_trt == "T") %>%
                                                pull(total_caught)),
                                       prob = parity_rate_treatment[i])
}

bloodmeal_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_control_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                 filter(cluster_trt == "C") %>%
                                                 pull(total_evaluated)),
                                        prob = prop_fed_untreated[i])
}

bloodmeal_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_treatment_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                   filter(cluster_trt == "T") %>%
                                                   pull(total_evaluated)),
                                          prob = prop_fed_treated[i])
}

aspiration_int_control_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_control_ppd[[i]] = rnbinom(nrow(abundance_intervention_control),
                                            mu=f_control[i],
                                            size = result$phi[i])
}

aspiration_int_treatment_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_treatment_ppd[[i]] = rnbinom(nrow(abundance_intervention_treatment),
                                              mu=f_treatment[i],
                                              size = result$phi[i])
}

sc_control_ppd = sc_treatment_ppd = list()
for(i in 1:nrow(result)){
  sc_control_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "C",]),
                               size = 1, prob=prob_sc_control[[i]])
  
  sc_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_treatment[[i]])
}

prop_parous_int_control_ppd_no_blood = parity_int_control_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_caught))

prop_parous_int_treatment_ppd_no_blood = parity_int_treatment_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_caught))

prop_bloodfed_int_control_ppd_no_blood = bloodmeal_int_control_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_evaluated))


prop_bloodfed_int_treated_ppd_no_blood = bloodmeal_int_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

mean_aspiration_int_control_ppd_no_blood = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_control_ppd_no_blood[i] = mean(aspiration_int_control_ppd[[i]])
}


mean_aspiration_int_treatment_ppd_no_blood = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_treatment_ppd_no_blood[i] = mean(aspiration_int_treatment_ppd[[i]])
}

prop_sc_control_no_blood = prop_sc_treatment_no_blood = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_sc_control_no_blood[i] = mean(sc_control_ppd[[i]])
  prop_sc_treatment_no_blood[i] = mean(sc_treatment_ppd[[i]])
}


#### Epi ####

result = result_no_epi

result$a_mult = result$alpha

prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau,
                                                                   d, rho, q_mult,
                                                                   a_mult, C = 0))
}

prop_fed_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_treated[i] =  with(result[i,],
                              prop_bloodfed_treated(au, qu, tau,
                                                    d, rho, q_mult,
                                                    a_mult, C))
}

f_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,], expected_catch_number_control(au, qu, tau, d, 
                                                                gu, lambda,
                                                                catch_prop))
}

f_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_treatment[i] = with(result[i,], expected_catch_number_treated(au, qu, tau, d, rho, q_mult, a_mult,
                                                                  C, gu, g_mult, lambda, catch_prop))
}

parity_rate_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,], parity(au, qu,
                                                   tau, d,
                                                   rho, q_mult,
                                                   a_mult, C = 0,
                                                   gu, g_mult))
}

parity_rate_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_treatment[i] = with(result[i,], parity(au, qu,
                                                     tau, d,
                                                     rho, q_mult,
                                                     a_mult, C,
                                                     gu, g_mult))
}

foi_c = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_c[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C= 0,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

foi_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_t[i] = with(result[i,],force_of_infection(au, qu, tau,
                                                d, rho, q_mult,
                                                a_mult, C,
                                                gu, g_mult,
                                                b, lambda,
                                                c, X, n))
}

prob_not_sc_control = prob_sc_control = list()
for(i in 1:nrow(result)){
  prob_not_sc_control [[i]] = exp(-foi_c[i] *
                                    epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  prob_sc_control[[i]] = 1 - prob_not_sc_control[[i]]
}

prob_not_sc_treatment = prob_sc_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_treatment [[i]] = exp(-foi_c[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_treatment[[i]] = 1 - prob_not_sc_treatment[[i]]
}


parity_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_control_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                              filter(cluster_trt == "C") %>%
                                              pull(total_caught)),
                                     prob = parity_rate_control[i])
}

parity_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_int_treatment_ppd[i] = rbinom(1, sum(parity_intervention %>%
                                                filter(cluster_trt == "T") %>%
                                                pull(total_caught)),
                                       prob = parity_rate_treatment[i])
}

bloodmeal_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_control_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                 filter(cluster_trt == "C") %>%
                                                 pull(total_evaluated)),
                                        prob = prop_fed_untreated[i])
}

bloodmeal_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_treatment_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                   filter(cluster_trt == "T") %>%
                                                   pull(total_evaluated)),
                                          prob = prop_fed_treated[i])
}

aspiration_int_control_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_control_ppd[[i]] = rnbinom(nrow(abundance_intervention_control),
                                            mu=f_control[i],
                                            size = result$phi[i])
}

aspiration_int_treatment_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_treatment_ppd[[i]] = rnbinom(nrow(abundance_intervention_treatment),
                                              mu=f_treatment[i],
                                              size = result$phi[i])
}

sc_control_ppd = sc_treatment_ppd = list()
for(i in 1:nrow(result)){
  sc_control_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "C",]),
                               size = 1, prob=prob_sc_control[[i]])
  
  sc_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_treatment[[i]])
}

prop_parous_int_control_ppd_no_epi = parity_int_control_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_caught))

prop_parous_int_treatment_ppd_no_epi = parity_int_treatment_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_caught))

prop_bloodfed_int_control_ppd_no_epi = bloodmeal_int_control_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_evaluated))


prop_bloodfed_int_treated_ppd_no_epi = bloodmeal_int_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

mean_aspiration_int_control_ppd_no_epi = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_control_ppd_no_epi[i] = mean(aspiration_int_control_ppd[[i]])
}


mean_aspiration_int_treatment_ppd_no_epi = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_treatment_ppd_no_epi[i] = mean(aspiration_int_treatment_ppd[[i]])
}

prop_sc_control_no_epi = prop_sc_treatment_no_epi = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_sc_control_no_epi[i] = mean(sc_control_ppd[[i]])
  prop_sc_treatment_no_epi[i] = mean(sc_treatment_ppd[[i]])
}





#### Plot ####
par(mfrow=c(1,1))
plot(-1, xlim = c(0, 6), ylim=c(0, .85), main="Parity", xlab="Left Out",
     ylab="PPD", xaxt="n")
axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)

abline(h = sum(parity_intervention$parous[parity_intervention$cluster_trt == "C"]) / 
         sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "C"]),
       col="black", lty = "dashed")

abline(h = sum(parity_intervention$parous[parity_intervention$cluster_trt == "T"]) / 
         sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "T"]),
       col="red", lty = "dashed")
boxplot(prop_parous_int_control_ppd_all, add = T, at = .75)
boxplot(prop_parous_int_treatment_ppd_all, add = T, at = 1.25, col="red")
boxplot(prop_parous_int_control_ppd_no_parity, add = T, at = 1.75)
boxplot(prop_parous_int_treatment_ppd_no_parity, add = T, at = 2.25, col="red")
boxplot(prop_parous_int_control_ppd_no_aspiration, add = T, at = 2.75)
boxplot(prop_parous_int_treatment_ppd_no_aspiration, add = T, at = 3.25, col="red")
boxplot(prop_parous_int_control_ppd_no_blood, add = T, at = 3.75)
boxplot(prop_parous_int_treatment_ppd_no_blood, add = T, at = 4.25, col="red")
boxplot(prop_parous_int_control_ppd_no_epi, add = T, at = 4.75)
boxplot(prop_parous_int_treatment_ppd_no_epi, add = T, at = 5.25, col="red")
legend("bottomright", c("Control", "Treatment"), col=c("black", "red"),
       lwd= 5)




plot(-1, xlim = c(0, 6), ylim=c(0, .6), main="Aspirations", xlab="Left Out",
     ylab="PPD", xaxt="n")
axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
abline(h = mean(abundance_intervention_control$total_caught),
       col="black", lty = "dashed")

abline(h = mean(abundance_intervention_treatment$total_caught),
       col="red", lty = "dashed")
boxplot(mean_aspiration_int_control_ppd_all, add = T, at = .75)
boxplot(mean_aspiration_int_treatment_ppd_all, add = T, at = 1.25, col="red")
boxplot(mean_aspiration_int_control_ppd_no_parity, add = T, at = 1.75)
boxplot(mean_aspiration_int_treatment_ppd_no_parity, add = T, at = 2.25, col="red")
boxplot(mean_aspiration_int_control_ppd_no_aspiration, add = T, at = 2.75)
boxplot(mean_aspiration_int_treatment_ppd_no_aspiration, add = T, at = 3.25, col="red")
boxplot(mean_aspiration_int_control_ppd_no_blood, add = T, at = 3.75)
boxplot(mean_aspiration_int_treatment_ppd_no_blood, add = T, at = 4.25, col="red")
boxplot(mean_aspiration_int_control_ppd_no_epi, add = T, at = 4.75)
boxplot(mean_aspiration_int_treatment_ppd_no_epi, add = T, at = 5.25, col="red")


plot(-1, xlim = c(0, 6), ylim=c(0.2, .7), main="Bloodmeal", xlab="Left Out",
     ylab="PPD", xaxt="n")
axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
abline(h = sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "C"]) /
         sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "C"]),
       col="black", lty = "dashed")

abline(h = sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "T"]) /
         sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "T"]),
       col="red", lty = "dashed")
boxplot(prop_bloodfed_int_control_ppd_all, add = T, at = .75)
boxplot(prop_bloodfed_int_treated_ppd_all, add = T, at = 1.25, col="red")
boxplot(prop_bloodfed_int_control_ppd_no_parity, add = T, at = 1.75)
boxplot(prop_bloodfed_int_treated_ppd_no_parity, add = T, at = 2.25, col="red")
boxplot(prop_bloodfed_int_control_ppd_no_aspiration, add = T, at = 2.75)
boxplot(prop_bloodfed_int_treated_ppd_no_aspiration, add = T, at = 3.25, col="red")
boxplot(prop_bloodfed_int_control_ppd_no_blood, add = T, at = 3.75)
boxplot(prop_bloodfed_int_treated_ppd_no_blood, add = T, at = 4.25, col="red")
boxplot(prop_bloodfed_int_control_ppd_no_epi, add = T, at = 4.75)
boxplot(prop_bloodfed_int_treated_ppd_no_epi, add = T, at = 5.25, col="red")



plot(-1, xlim = c(0, 6), ylim=c(0, 1), main="Epi", xlab="Left Out",
     ylab="PPD", xaxt="n")
axis("None", side = 1, at = 1, las = 1)
axis("Pa", side = 1, at = 2, las = 1)
axis("Asp", side = 1, at = 3, las = 1)
axis("Blood", side = 1, at = 4, las = 1)
axis("Epi", side = 1, at = 5, las = 1)
abline(h = mean(epi$outcome[epi$Cluster_Allocation == "C"]),
       col="black", lty = "dashed")
abline(h = mean(epi$outcome[epi$Cluster_Allocation == "T"]),
       col="red", lty = "dashed")
boxplot(prop_sc_control_all, add = T, at = .75)
boxplot(prop_sc_treatment_all, add = T, at = 1.25, col="red")
boxplot(prop_sc_control_no_parity, add = T, at = 1.75)
boxplot(prop_sc_treatment_no_parity, add = T, at = 2.25, col="red")
boxplot(prop_sc_control_no_aspiration, add = T, at = 2.75)
boxplot(prop_sc_treatment_no_aspiration, add = T, at = 3.25, col="red")
boxplot(prop_sc_control_no_blood, add = T, at = 3.75)
boxplot(prop_sc_treatment_no_blood, add = T, at = 4.25, col="red")
boxplot(prop_sc_control_no_epi, add = T, at = 4.75)
boxplot(prop_sc_treatment_no_epi, add = T, at = 5.25, col="red")