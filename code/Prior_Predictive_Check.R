# for baseline ento data

source("functions.R")

# load trial data 
bloodmeal_intervention = read.csv("../data/bloodmeal_intervention.csv")
bloodmeal_baseline = read.csv("../data/bloodmeal_baseline.csv")

# take total full to be the sum of hald and full
bloodmeal_intervention = bloodmeal_intervention %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

bloodmeal_baseline = bloodmeal_baseline %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

abundance_baseline = read.csv("../data/abundance_baseline.csv")
abundance_intervention_control = read.csv("../data/abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("../data/abundance_intervention_treatment.csv")

parity_baseline = read.csv("../data/parity_baseline.csv")
parity_intervention = read.csv("../data/parity_intervention.csv")

epi = read.csv("../data/epi.csv")

load("../data/treatedclusters.RData")
p=5.2


#### PPD: Calculate quantities ####

load("../results/MultivariatePrior.RData")

prior_samples = data.frame(rmvnorm(1000, mu_mvn, sigma_mvn))
df_names = names(prior_samples)
prior_samples_exp = exp(prior_samples)
prior_samples_sig = sigmoid(prior_samples)

prior_samples = cbind(prior_samples_exp[,1:4], prior_samples_sig[,5],
                      prior_samples_exp[,6],
                      prior_samples_sig[,7], prior_samples_exp[,8:9],
                      prior_samples_sig[,10:12])
names(prior_samples) = df_names

result = prior_samples

prop_fed_untreated = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  prop_fed_untreated[i] = with(result[i,],
                               prop_bloodfed_untreated(au, qu, tau, d, rho = 0,
                                                       q_mult = 1, a_mult = 1, C = 0))
}


f_control= rep(NA, nrow(result))
for(i in 1:nrow(result)){
  f_control[i] = with(result[i,],
                      expected_catch_number_untreated(au, qu, tau, d, rho=0, q_mult=1,
                                                      a_mult=1, C=0, gu, g_mult=1, lambda, catch_prop))
}

parity_rate_control= rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_rate_control[i] = with(result[i,],
                                parity(au, qu, tau, d, rho = 0, q_mult = 1,
                                       a_mult = 1, C = 0, gu, g_mult = 1))
}


#### PPD: calculate data sets ####

parity_baseline_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  parity_baseline_ppd[i] = rbinom(1, sum(parity_baseline$total_caught),
                                  prob = parity_rate_control[i])
}

bloodmeal_baseline_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  bloodmeal_baseline_ppd[i] = rbinom(1, sum(bloodmeal_baseline$total_evaluated),
                                     prob = prop_fed_untreated[i])
}

aspiration_baseline_ppd = list()
for(i in 1:nrow(result)){
  aspiration_baseline_ppd[[i]] = rnbinom(nrow(abundance_baseline),
                                         mu=f_control[i],
                                         size = result$phi[i])
}

#### calculate summary statistics ####

prop_parous_baseline_ppd = parity_baseline_ppd / 
  sum(parity_baseline$total_caught)

prop_bloodfed_baseline_ppd = bloodmeal_baseline_ppd / 
  sum(bloodmeal_baseline$total_evaluated)

mean_aspiration_baseline_ppd = rep(NA, nrow(result))
for(i in 1:nrow(result)){
  mean_aspiration_baseline_ppd[i] = mean(aspiration_baseline_ppd[[i]])
}

#### plot ####

par(mfrow=c(1,3), mar=c(1, 3, 4, 1), cex=1.2)

boxplot(prop_parous_baseline_ppd, main="Parity")
points(1, sum(parity_baseline$parous) / sum(parity_baseline$total_caught),
       col="red", pch=19)

boxplot(prop_bloodfed_baseline_ppd, main="Bloodmeal", ylim=c(0.48, .56))
points(1, sum(bloodmeal_baseline$total_fed) /
         sum(bloodmeal_baseline$total_evaluated),
       col="red", pch=19)

boxplot(mean_aspiration_baseline_ppd, main="Aspiration", ylim=c(0,20))
points(1, mean(abundance_baseline$total_caught),
       col="red", pch=19)




##### Intervention phase data- from Multivariate Prior All#####

source("functions.R")

# load trial data 
bloodmeal_intervention = read.csv("../data/bloodmeal_intervention.csv")
bloodmeal_baseline = read.csv("../data/bloodmeal_baseline.csv")

# take total full to be the sum of hald and full
bloodmeal_intervention = bloodmeal_intervention %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

bloodmeal_baseline = bloodmeal_baseline %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

abundance_baseline = read.csv("../data/abundance_baseline.csv")
abundance_intervention_control = read.csv("../data/abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("../data/abundance_intervention_treatment.csv")

parity_baseline = read.csv("../data/parity_baseline.csv")
parity_intervention = read.csv("../data/parity_intervention.csv")

epi = read.csv("../data/epi.csv")

load("../data/treatedclusters.RData")

C=.44
p=5.2

load("../results/MultivariatePrior_All.RData")
prior_samples = data.frame(rmvnorm(1000, mu_mvn, sigma_mvn))


prior_samples_exp = exp(prior_samples)
prior_samples_sig = sigmoid(prior_samples)

prior_samples = cbind(prior_samples_exp[,1:4], prior_samples_sig[,5],
                      prior_samples_exp[,6],
                      prior_samples_sig[,7], prior_samples_exp[,8:9],
                      prior_samples_sig[,10:13], prior_samples_exp[,14:16])


names(prior_samples) = c("qu", "au", "tau", "d", "gu", "lambda", "catch_prop",
                         "phi", "n", "X", "c", "b", "rho", "a_mult", "g_mult",
                         "q_mult")

result = prior_samples


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


par(mfrow=c(5, 3), mar=c(1, 4, 1, 1), cex=1.2)

boxplot(prop_parous_int_control_ppd, main="Control", ylim=c(.5, 1), ylab="Parity")
points(1, sum(parity_intervention$parous[parity_intervention$cluster_trt == "C"]) / 
         sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "C"]),
       col="red", pch=19)



boxplot(prop_parous_int_treatment_ppd, main="Treatment", ylim=c(.5, 1))
points(1, sum(parity_intervention$parous[parity_intervention$cluster_trt == "T"]) / 
         sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "T"]),
       col="red", pch=19)



boxplot(prop_parous_int_treatment_ppd, main="Untreated", ylim=c(.5, 1))




boxplot(prop_bloodfed_int_control_ppd, ylab="Bloodmeal", ylim=c(0, 1))
points(1, sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "C"]) /
         sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "C"]),
       col="red", pch=19)



boxplot(prop_bloodfed_int_treated_ppd, ylim=c(0,1))
points(1, sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "T"]) /
         sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "T"]),
       col="red", pch=19)



boxplot(prop_bloodfed_int_untreated_in_treated_ppd, ylim=c(0, 1))


boxplot(mean_aspiration_int_control_ppd, ylab="Aspiration", ylim=c(0, 100))
points(1, mean(abundance_intervention_control$total_caught),
       col="red", pch=19)







boxplot(mean_aspiration_int_treatment_ppd,  ylim=c(0, 100))
points(1, mean(abundance_intervention_treatment$total_caught),
       col="red", pch=19)



boxplot(mean_aspiration_int_untreated_in_treatment_ppd,  ylim=c(.2, .6))



boxplot(foi_c, ylab="Force of Infection", ylim=c(0, .0015))
# points(1, mean(epi$outcome[epi$Cluster_Allocation == "C"]), pch=19, col="red")

boxplot(foi_t, ylab="Force of Infection", ylim=c(0, .0015))
# boxplot(prop_sc_treatment, ylim=c(0.1, .4))
# points(1, mean(epi$outcome[epi$Cluster_Allocation == "T"]), pch=19, col="red")

boxplot(foi_u_in_t, ylab="Force of Infection", ylim=c(0, .0015))
# boxplot(prop_sc_untreated_in_treatment, ylim=c(0.1, .4))



boxplot(prop_sc_control, ylab="Seroconversions", ylim=c(0, 1))
points(1, mean(epi$outcome[epi$Cluster_Allocation == "C"]), pch=19, col="red")






boxplot(prop_sc_treatment, ylim=c(0, 1))
points(1, mean(epi$outcome[epi$Cluster_Allocation == "T"]), pch=19, col="red")

