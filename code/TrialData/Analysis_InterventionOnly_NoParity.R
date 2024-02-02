
# This script analyzes CRC results for only intervention data

#### load packages ####

library(dplyr)
library(scales)
library(mvtnorm)
library(fitdistrplus)
library(PerformanceAnalytics)
library(corrplot)
library(GGally)



#### load the result ####

load("../../results/TrialData/InterventionOnly_NoParity.RData")
head(result)

# plot the chains

chain_1 = result[seq(1, nrow(result), by=3), ]
chain_2 = result[seq(2, nrow(result), by=3), ]
chain_3 = result[seq(3, nrow(result), by=3), ]


plot(chain_1$q_mult, type = "l", xlab = "Iteration", ylab = "Baseline Biting Rate",
     main = "Exit Effect")
lines(chain_2$q_mult, col = "red")
lines(chain_3$q_mult, type = "l", col = "green")


result = rbind(result[seq(1, nrow(result), by=3), ][seq(1, nrow(result)/3,
                                                        by = 350),],
               result[seq(2, nrow(result), by=3), ][seq(1, nrow(result)/3,
                                                        by = 350),],
               result[seq(3, nrow(result), by=3), ][seq(1, nrow(result)/3,
                                                        by = 350),])


chain_1_thin = result[1:38, ]
chain_2_thin = result[39:77,]
chain_3_thin = result[78:nrow(result),]


plot(chain_1_thin$au, type = "l", xlab = "Iteration", ylab = "Baseline Biting Rate",
     main = "Before thinning", ylim= c(0, 2))
lines(chain_2_thin$au, col = "red")
lines(chain_3_thin$au, type = "l", col = "green")

head(result)

#### plot densities of the parameters ####

source("../functions.R")

load("../../results/MultivariatePrior_All.RData")



prior_samples = data.frame(rmvnorm(10000, mu_mvn, sigma_mvn))
df_names = names(result)
prior_samples_exp = exp(prior_samples)
prior_samples_sig = sigmoid(prior_samples)

prior_samples = cbind(prior_samples_exp[,1:4], prior_samples_sig[,5],
                      prior_samples_exp[,6],
                      prior_samples_sig[,7], prior_samples_exp[,8:9],
                      prior_samples_sig[,10:13], prior_samples_exp[,14:16])

names(prior_samples) = names(result)

png("../../results/ManuscriptFigures/Non_SR_Effects.png",
    height=8.5, width=11, units = "in", pointsize = 12, res = 72)


par(mfrow=c(3, 4), cex=1.2, mar=c(2, 3, 3, 1))


plot(density(result$qu), main="Exit Rate, \nUntreated", xlab="", xlim=c(0, 2.5), ylim=c(0, 2.5))
lines(density(prior_samples$qu))
polygon(density(prior_samples$qu),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$qu), col=rgb(red=1, 
                                    green=0,
                                    blue=0,
                                    alpha=.5))

legend("topright", c("Prior", "Posterior"), col = c(rgb(red=0,  green=0,  blue=1,alpha=.5),
                                                    rgb(red=1, 
                                                        green=0,
                                                        blue=0,
                                                        alpha=.5)), pch=19)


plot(density(result$au), main="Biting Rate, \nUntreated",
     xlab="", xlim=c(0, 2), ylim= c(0, 2.5))
lines(density(prior_samples$au))
polygon(density(prior_samples$au),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$au), col=rgb(red=1, 
                                    green=0,
                                    blue=0,
                                    alpha=.5))


plot(density(result$tau), main="Entrance Rate",
     xlab="", ylim=c(0, .8), xlim=c(0, 6))
lines(density(prior_samples$tau))
polygon(density(prior_samples$tau),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$tau), col=rgb(red=1, 
                                     green=0,
                                     blue=0,
                                     alpha=.5))


plot(density(result$d), main="Digestion",
     xlab="", xlim=c(0, 1.5), ylim=c(0, 5))
lines(density(prior_samples$d))
polygon(density(prior_samples$d),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$d), col=rgb(red=1, 
                                   green=0,
                                   blue=0,
                                   alpha=.5))


plot(density(result$gu), main="Mortality Rate, \nUntreated",
     xlab="", xlim=c(0, .3))
lines(density(prior_samples$gu))
polygon(density(prior_samples$gu),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$gu), col=rgb(red=1, 
                                    green=0,
                                    blue=0,
                                    alpha=.5))


plot(density(result$lambda), main="Mosquito \nEmergence",
     xlab="", xlim=c(0, .15), ylim=c(0,80))
lines(density(prior_samples$lambda[prior_samples$lambda< 1]))
polygon(density(prior_samples$lambda[prior_samples$lambda< 1]),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$lambda), col=rgb(red=1, 
                                        green=0,
                                        blue=0,
                                        alpha=.5))


plot(density(result$catch_prop), main="Aspiration \nCatch Proportion",
     xlab="", ylim=c(0, 2), xlim=c(0,1))
lines(density(prior_samples$catch_prop))
polygon(density(prior_samples$catch_prop),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$catch_prop), col=rgb(red=1, 
                                            green=0,
                                            blue=0,
                                            alpha=.5))


plot(density(result$phi), main="Dispersion,\n Aspiration",
     xlab="", xlim=c(.15, .18),)
lines(density(prior_samples$phi))
polygon(density(prior_samples$phi),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$phi), col=rgb(red=1, 
                                     green=0,
                                     blue=0,
                                     alpha=.5))


plot(density(result$n), main="EIP",
     xlab="", xlim=c(8, 22), ylim=c(0, .25))
lines(density(prior_samples$n))
polygon(density(prior_samples$n),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$n), col=rgb(red=1, 
                                   green=0,
                                   blue=0,
                                   alpha=.5))


plot(density(result$X), main="Prevalence",
     xlab="", xlim=c(0, .4), ylim=c(0, 15))
lines(density(prior_samples$X))
polygon(density(prior_samples$X),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$X), col=rgb(red=1, 
                                   green=0,
                                   blue=0,
                                   alpha=.5))


plot(density(result$c), main="Transmission Prob 1",
     xlab="", xlim=c(0, 1), ylim=c(0, 2.2))
lines(density(prior_samples$c))
polygon(density(prior_samples$c),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$c), col=rgb(red=1, 
                                   green=0,
                                   blue=0,
                                   alpha=.5))


plot(density(result$b), main="Transmission Prob 2",
     xlab="", xlim=c(0, 1), ylim=c(0, 2.3))
lines(density(prior_samples$b))
polygon(density(prior_samples$b),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
polygon(density(result$b), col=rgb(red=1, 
                                   green=0,
                                   blue=0,
                                   alpha=.5))


dev.off()


png("../../results/ManuscriptFigures/SR_Effects.png", height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)

par(mfrow=c(2,2), mar=c(3,4,4,4), cex=1.5)
plot(density(result$rho, bw=.05), main="Repellency",
     xlab="", xlim=c(0, 1))
lines(density(prior_samples$rho))
polygon(density(result$rho, bw=.05), col=rgb(red=1, 
                                             green=0,
                                             blue=0,
                                             alpha=.5))
polygon(density(prior_samples$rho),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
abline(v=0, lty="dashed", lwd=2)


plot(density(result$a_mult), main="SR Biting Effect",
     xlab="", xlim=c(0.5, 2), ylim=c(0, 10))
lines(density(prior_samples$a_mult))
polygon(density(result$a_mult), col=rgb(red=1, 
                                        green=0,
                                        blue=0,
                                        alpha=.5))
polygon(density(prior_samples$a_mult),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
abline(v = 1, lty="dashed", lwd = 2)

legend("topright", c("Prior", "Posterior"), col = c(rgb(red=0,  green=0,  blue=1,alpha=.5),
                                                    rgb(red=1, 
                                                        green=0,
                                                        blue=0,
                                                        alpha=.5)), pch=19)

plot(density(result$g_mult), main="SR Effect, Mortality",
     xlab="", xlim=c(0, 2.5), ylim=c(0, 2))
lines(density(prior_samples$g_mult))
polygon(density(result$g_mult), col=rgb(red=1, 
                                        green=0,
                                        blue=0,
                                        alpha=.5))
polygon(density(prior_samples$g_mult),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
abline(v = 1, lty="dashed", lwd = 2)

plot(density(result$q_mult, bw=.03), main="SR Effect, Exit Rate",
     xlab="", xlim=c(0, 2))
lines(density(prior_samples$q_mult))
polygon(density(result$q_mult, bw=.03), col=rgb(red=1, 
                                                green=0,
                                                blue=0,
                                                alpha=.5))
polygon(density(prior_samples$q_mult),
        col=rgb(red=0,  green=0,  blue=1,alpha=.5))
abline(v = 1, lty="dashed", lwd = 2)


dev.off()


# 95% credible intervals
quantile(result$rho, c(0.025, .5, .975))
1 - quantile(result$a_mult, c(0.025, .5, .975))
(1 - quantile(result$g_mult, c(0.025, .5, .975))) * 100
(1 - quantile(result$q_mult, c(0.025, .5, .975))) * 100

# effect on overall biting and mortality rates

ac_calc_control =ac_calc_treatment = rep(NA, nrow(result))

for(i in 1:nrow(result)){
  ac_calc_control[i] = with(result[i,], ac_calc(au, qu, tau, d, rho, q_mult, a_mult, C=0))
  ac_calc_treatment[i] = with(result[i,], ac_calc(au, qu, tau, d, rho, q_mult, a_mult, C))
  
}

ac_calc_treatment / ac_calc_control

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


#### PPD: Calculate quantities ####


prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau, d,
                                                                   rho, q_mult = 1,
                                                                   a_mult = 1, C = 0))
}

prop_fed_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_treated[i] = with(result[i,],
                             prop_bloodfed_treated(au, qu, tau, d, rho, q_mult, a_mult, C))
}

prop_fed_untreated_in_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated_in_treated[i] = with(result[i,], prop_bloodfed_untreated(au, qu, tau, d,
                                                                   rho, q_mult,
                                                                   a_mult, C))
}


f_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,],
                      expected_catch_number_untreated(au, qu, tau, d, rho, q_mult,
                                                      a_mult, C=0, gu, g_mult, lambda, catch_prop))
}


f_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_treatment[i] = with(result[i,],
                        expected_catch_number_treated(au, qu, tau, d, rho, q_mult,
                                                      a_mult, C, gu, g_mult, lambda, catch_prop))
}

f_untreated_in_treated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_untreated_in_treated[i] = with(result[i,],
                      expected_catch_number_untreated(au, qu, tau, d, rho, q_mult,
                                                      a_mult, C, gu, g_mult, lambda, catch_prop))
}


parity_rate_control = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,],
                                parity(au, qu, tau, d, rho = 0, q_mult = 1, a_mult = 1, C = 0, gu, g_mult = 1))
}

parity_rate_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_treatment[i] = with(result[i,],
                                  parity(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult))
}

foi_c = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_c[i] = with(result[i,], force_of_infection_untreated(au, qu, tau, d, rho = 0, q_mult = 1,
                                                           a_mult = 1, C = 0, gu, g_mult = 1,
                                                           b, lambda, c, X, n))
}






foi_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_t[i] = with(result[i,], force_of_infection_treated(au, qu, tau, d, rho, q_mult,
                                                         a_mult, C, gu, g_mult,
                                                         b, lambda, c, X, n))
}


foi_u_in_t = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  foi_u_in_t[i] = with(result[i,], force_of_infection_untreated(au, qu, tau, d, rho, q_mult,
                                                           a_mult, C, gu, g_mult,
                                                           b, lambda, c, X, n))
}


prob_not_sc_control = prob_sc_control = list()
for(i in 1:nrow(result)){
  prob_not_sc_control [[i]] = exp(-foi_c[i] *
                                    epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  prob_sc_control[[i]] = 1 - prob_not_sc_control[[i]]
}

prob_not_sc_treatment = prob_sc_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_treatment [[i]] = exp(-foi_t[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_treatment[[i]] = 1 - prob_not_sc_treatment[[i]]
}

prob_not_sc_untreated_in_treatment = prob_sc_untreated_in_treatment = list()
for(i in 1:nrow(result)){
  prob_not_sc_untreated_in_treatment [[i]] = exp(-foi_t[i] *
                                      epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
  prob_sc_untreated_in_treatment[[i]] = 1 - prob_not_sc_untreated_in_treatment[[i]]
}



#### PPD: calculate data sets ####

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

bloodmeal_int_untreated_in_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_int_untreated_in_treatment_ppd[i] = rbinom(1, sum(bloodmeal_intervention %>%
                                                   filter(cluster_trt == "T") %>%
                                                   pull(total_evaluated)),
                                          prob = prop_fed_untreated_in_treated[i])
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

aspiration_int_untreated_in_treatment_ppd = list()
for(i in 1:nrow(result)){
  aspiration_int_untreated_in_treatment_ppd[[i]] = rnbinom(nrow(abundance_intervention_treatment),
                                              mu=f_untreated_in_treated[i],
                                              size = result$phi[i])
}

sc_control_ppd = sc_treatment_ppd = sc_untreated_in_treatment_ppd =  list()
for(i in 1:nrow(result)){
  sc_control_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "C",]),
                               size = 1, prob=prob_sc_control[[i]])

  sc_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_treatment[[i]])

  sc_untreated_in_treatment_ppd[[i]] = rbinom(nrow(epi[epi$Cluster_Allocation == "T",]),
                                 size = 1, prob=prob_sc_untreated_in_treatment[[i]])


}


#### calculate summary statistics ####

prop_parous_int_control_ppd = parity_int_control_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_caught))

prop_parous_int_treatment_ppd = parity_int_treatment_ppd / 
  sum(parity_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_caught))


prop_bloodfed_int_control_ppd = bloodmeal_int_control_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "C") %>%
        pull(total_evaluated))


prop_bloodfed_int_treated_ppd = bloodmeal_int_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

prop_bloodfed_int_untreated_in_treated_ppd = bloodmeal_int_untreated_in_treatment_ppd / 
  sum(bloodmeal_intervention %>%
        filter(cluster_trt == "T") %>%
        pull(total_evaluated))

mean_aspiration_int_control_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_control_ppd[i] = mean(aspiration_int_control_ppd[[i]])
}


mean_aspiration_int_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_treatment_ppd[i] = mean(aspiration_int_treatment_ppd[[i]])
}

mean_aspiration_int_untreated_in_treatment_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_int_untreated_in_treatment_ppd[i] = 
    mean(aspiration_int_untreated_in_treatment_ppd[[i]])
}

prop_sc_control = prop_sc_treatment = prop_sc_untreated_in_treatment = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_sc_control[i] = mean(sc_control_ppd[[i]])
  prop_sc_treatment[i] = mean(sc_treatment_ppd[[i]])
  prop_sc_untreated_in_treatment[i] = mean(sc_untreated_in_treatment_ppd[[i]])
}


#### plot ####

confidence_interval_proportion = function(successes, sample_size){
  
  p_hat = successes / sample_size
  
  z = 1.96
  
  lower_limit = p_hat - z * sqrt(p_hat * (1 - p_hat) / sample_size)
  upper_limit = p_hat + z * sqrt(p_hat * (1 - p_hat) / sample_size)
  
  return(c(lower_limit, upper_limit))
}

confidence_interval_count = function(events){
  return(c(mean(events) - 1.96 * sqrt(mean(events) / length(events)),
           mean(events) + 1.96 * sqrt(mean(events) / length(events))))
}

png("../../results/ManuscriptFigures/PPD.png", height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)

par(mfrow=c(5, 3), mar=c(1, 5, 1, 2), cex=1.2)

boxplot(prop_parous_int_control_ppd, main="Control", ylim=c(.5, 1), ylab="Parity")
points(1, sum(parity_intervention$parous[parity_intervention$cluster_trt == "C"]) / 
         sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "C"]),
       col="red", pch=19)
segments(x0 = 1, x1 = 1, y0 = confidence_interval_proportion(parity_intervention %>%
                                                               filter(cluster_trt == "C") %>%
                                                               pull(parous) %>%
                                                               sum(), parity_intervention %>%
                                                               filter(cluster_trt == "C") %>%
                                                               pull(total_caught) %>%
                                                               sum())[1],
         y1 = confidence_interval_proportion(parity_intervention %>%
                                               filter(cluster_trt == "C") %>%
                                               pull(parous) %>%
                                               sum(), parity_intervention %>%
                                               filter(cluster_trt == "C") %>%
                                               pull(total_caught) %>%
                                               sum())[2], col = "red")



boxplot(prop_parous_int_treatment_ppd, main="Treatment", ylim=c(.5, 1))
points(1, sum(parity_intervention$parous[parity_intervention$cluster_trt == "T"]) / 
         sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "T"]),
       col="red", pch=19)

segments(x0 = 1, x1 = 1, y0 = confidence_interval_proportion(parity_intervention %>%
                                                               filter(cluster_trt == "T") %>%
                                                               pull(parous) %>%
                                                               sum(), parity_intervention %>%
                                                               filter(cluster_trt == "T") %>%
                                                               pull(total_caught) %>%
                                                               sum())[1],
         y1 = confidence_interval_proportion(parity_intervention %>%
                                               filter(cluster_trt == "T") %>%
                                               pull(parous) %>%
                                               sum(), parity_intervention %>%
                                               filter(cluster_trt == "T") %>%
                                               pull(total_caught) %>%
                                               sum())[2], col = "red")





boxplot(prop_parous_int_treatment_ppd, main="Untreated", ylim=c(.5, 1))




boxplot(prop_bloodfed_int_control_ppd, ylab="Bloodmeal", ylim=c(.6, .66))
points(1, sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "C"]) /
         sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "C"]),
       col="red", pch=19)

segments(x0 = 1, x1 = 1, y0 = confidence_interval_proportion(bloodmeal_intervention %>%
           filter(cluster_trt == "C") %>%
           pull(total_fed) %>%
           sum(), bloodmeal_intervention %>%
                                            filter(cluster_trt == "C") %>%
                                            pull(total_evaluated) %>%
                                            sum())[1],
         y1 = confidence_interval_proportion(bloodmeal_intervention %>%
                                               filter(cluster_trt == "C") %>%
                                               pull(total_fed) %>%
                                               sum(), bloodmeal_intervention %>%
                                               filter(cluster_trt == "C") %>%
                                               pull(total_evaluated) %>%
                                               sum())[2], col = "red")


boxplot(prop_bloodfed_int_treated_ppd, ylim=c(.6, .66))
points(1, sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "T"]) /
         sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "T"]),
       col="red", pch=19)
segments(x0 = 1, x1 = 1, y0 = confidence_interval_proportion(bloodmeal_intervention %>%
                                                               filter(cluster_trt == "T") %>%
                                                               pull(total_fed) %>%
                                                               sum(), bloodmeal_intervention %>%
                                                               filter(cluster_trt == "T") %>%
                                                               pull(total_evaluated) %>%
                                                               sum())[1],
         y1 = confidence_interval_proportion(bloodmeal_intervention %>%
                                               filter(cluster_trt == "T") %>%
                                               pull(total_fed) %>%
                                               sum(), bloodmeal_intervention %>%
                                               filter(cluster_trt == "T") %>%
                                               pull(total_evaluated) %>%
                                               sum())[2], col = "red")



boxplot(prop_bloodfed_int_untreated_in_treated_ppd, ylim=c(.6, .66))


boxplot(mean_aspiration_int_control_ppd, ylab="Aspiration", ylim=c(.2, .6))
points(1, mean(abundance_intervention_control$total_caught),
       col="red", pch=19)
segments(x0=1, x1=1, y0 = confidence_interval_count(abundance_intervention_control$total_caught)[1],
         y1 = confidence_interval_count(abundance_intervention_control$total_caught)[2] )



boxplot(mean_aspiration_int_treatment_ppd,  ylim=c(.2, .6))
points(1, mean(abundance_intervention_treatment$total_caught),
       col="red", pch=19)
  segments(x0=1, x1=1, y0 = confidence_interval_count(abundance_intervention_treatment$total_caught)[1],
         y1 = confidence_interval_count(abundance_intervention_treatment$total_caught)[2] )

boxplot(mean_aspiration_int_untreated_in_treatment_ppd,  ylim=c(.2, .6))



boxplot(foi_c, ylab="Force of \nInfection", ylim=c(0, .0015))


boxplot(foi_t, ylim=c(0, .0015))


boxplot(foi_u_in_t, ylim=c(0, .0015))



boxplot(prop_sc_control, ylab="Proportion \nSeroconverting", ylim=c(0.1, .4))
points(1, mean(epi$outcome[epi$Cluster_Allocation == "C"]), pch=19, col="red")

segments(x0 = 1, x1 = 1, y0 = confidence_interval_count(epi %>%
                                                               filter(Cluster_Allocation == "C") %>%
                                                          pull(outcome))[1],
         y1 = confidence_interval_count(epi %>%
                                          filter(Cluster_Allocation == "C") %>%
                                          pull(outcome))[2], col="red")




boxplot(prop_sc_treatment, ylim=c(0.1, .4))
points(1, mean(epi$outcome[epi$Cluster_Allocation == "T"]), pch=19, col="red")


segments(x0 = 1, x1 = 1, y0 = confidence_interval_count(epi %>%
                                                          filter(Cluster_Allocation == "T") %>%
                                                          pull(outcome))[1],
         y1 = confidence_interval_count(epi %>%
                                          filter(Cluster_Allocation == "T") %>%
                                          pull(outcome))[2], col="red")

dev.off()


#### scatter plots of variables ####
png("../../results/ManuscriptFigures/Correlations.png",height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)


par(cex = 1.5)
pairs(result %>%
        rename(Mortality = gu, Biting = au, Digestion = d, "Prop'n Caught" = catch_prop, Repellency = rho,
               "Exit Effect" = q_mult, "Biting Effect" = a_mult) %>%
        dplyr::select(Mortality, Biting, Digestion, "Prop'n Caught", Repellency, "Exit Effect", 
                      "Biting Effect"))
dev.off()

#### varying coverage and its effect on quantities ####

coverages = seq(0, 1, by=.01)

direct_effect = indirect_effect = total_effect = overall_effect =
  rep(NA, length(coverages))

direct_effect_lower = indirect_effect_lower = total_effect_lower =
  overall_effect_lower = direct_effect_upper = indirect_effect_upper = 
  total_effect_upper = overall_effect_upper = rep(NA, length(coverages))

result = result[sample(1:nrow(result), size = 200),]

for(i in 1:length(coverages)){
  
  print(i)
  
  coverage = coverages[i]
  
  direct_effects = indirect_effects = total_effects = overall_effects = rep(NA, nrow(result))
  
  for(j in 1:nrow(result)){
    parms = result[j,]
    
    
    foi_c = force_of_infection_untreated(parms$au, parms$qu, parms$tau, parms$d, 0, 1,
                                         1, C=0, parms$gu,1, parms$b, parms$lambda,
                                         parms$c, parms$X, parms$n)
    
    foi_u = force_of_infection_untreated(au = parms$au, qu = parms$qu, tau = parms$tau, d = parms$d,
                                         rho = parms$rho, q_mult = parms$q_mult,
                                         a_mult = parms$a_mult, C = coverage, gu = parms$gu, g_mult = parms$g_mult,
                                         b = parms$b, lambda = parms$lambda,
                                         c = parms$c, X = parms$X, n = parms$n)
    foi_t = force_of_infection_treated(parms$au, parms$qu, parms$tau, parms$d, parms$rho, parms$q_mult,
                                       parms$a_mult, C = coverage, parms$gu, parms$g_mult, parms$b, parms$lambda,
                                       parms$c, parms$X, parms$n)
    
    direct_effects[j] = 1 - (foi_t / foi_u)
    indirect_effects[j] = 1 - (foi_u / foi_c)
    total_effects[j] = 1 - (foi_t / foi_c)
    overall_effects[j] = 1 - ((coverage * foi_t + (1-coverage)*foi_u) / foi_c)
    
  }
  
  direct_effect[i] = median(direct_effects)
  indirect_effect[i] = median(indirect_effects)
  total_effect[i] = median(total_effects)
  overall_effect[i] = median(overall_effects)
  
  direct_effect_lower[i] = quantile(direct_effects, .025, na.rm=T)
  indirect_effect_lower[i] = quantile(indirect_effects, .025, na.rm=T)
  total_effect_lower[i] = quantile(total_effects, .025, na.rm=T)
  overall_effect_lower[i] = quantile(overall_effects, .025, na.rm=T)
  
  direct_effect_upper[i] = quantile(direct_effects, .975, na.rm=T)
  indirect_effect_upper[i] = quantile(indirect_effects, .975, na.rm=T)
  total_effect_upper[i] = quantile(total_effects, .975, na.rm=T)
  overall_effect_upper[i] = quantile(overall_effects, .975, na.rm=T)
}

png("../../results/ManuscriptFigures/Coverage.png",height=8.5,
    width=11, units = "in", pointsize = 12, res = 72)

par(mfrow=c(1,1), mar=c(4,5,3,4), cex=1.5)
plot(direct_effect ~ coverages, type="l", ylab="Effect Size",
     xlab="Proportion Covered", lwd=3, ylim=c(-2.5, .6))
lines(indirect_effect[-101] ~ coverages[-101], col="red", lwd=3)
lines(total_effect ~ coverages, col="blue", lwd=3)
lines(overall_effect ~ coverages, col="green", lwd=3)
legend("bottomleft", c("Direct T:U", "Indirect U:C", "Total T:C", "Overall"),
       col=c("black", "red", "blue", "green"), lty=1, lwd=3)

lines(direct_effect_lower ~ coverages, lty="dashed", col="black")
lines(direct_effect_upper ~ coverages, lty="dashed", col="black")


lines(indirect_effect_lower[-101] ~ coverages[-101], lty="dashed", col="red")
lines(indirect_effect_upper[-101] ~ coverages[-101], lty="dashed", col="red")

lines(total_effect_lower ~ coverages, lty="dashed", col="blue")
lines(total_effect_upper ~ coverages, lty="dashed", col="blue")

lines(overall_effect_lower ~ coverages, lty="dashed", col="green")
lines(overall_effect_upper ~ coverages, lty="dashed", col="green")


abline(v=.44,  col = "grey")

dev.off()
