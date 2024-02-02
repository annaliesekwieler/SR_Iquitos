# This script creates 100 data sets from parameter combinations
# that make epidemiological sense given transmission patterns in 
# the past


#### load libraries ####

library(dplyr)
library(fitdistrplus)

#### import trial data ####
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


#### function to generate a single data set ####
generate_data = function(au_true, d_true, gu_true, qu_true, tau_true, n_true,
                         b_true, c_true,
                         X_true, a_mult_true, g_mult_true, rho_true,
                         q_mult_true, lambda_true,
                         phi_true, catch_prop_true){

  
  
  # find the proportion fed 
  prop_fed_control_true = prop_bloodfed_untreated(au_true, qu_true, tau_true,
                                                  d_true, rho = 0,
                                                  q_mult = 1, a_mult = 1,
                                                   C = 0)
  prop_fed_treated_true = prop_bloodfed_treated(au_true, qu_true, tau_true,
                                                d_true, rho_true,
                                                q_mult_true, a_mult_true, C)
  
  # find the force of infection (daily)
  foi_control_true = force_of_infection_untreated(au_true, qu_true, tau_true,
                                        d_true, rho = 0, q_mult = 1,
                                        a_mult = 1, C = 0, gu_true,
                                        g_mult = 1, b_true, lambda_true, 
                                        c_true, X_true, n_true)
  foi_treatment_true = force_of_infection_treated(au_true, qu_true, tau_true,
                                          d_true, rho_true, q_mult_true,
                                          a_mult_true, C, gu_true,
                                          g_mult_true, b_true, lambda_true, 
                                          c_true, X_true, n_true)
  
  # find the parity rates
  parity_control_true = parity(au_true, qu_true, tau_true, d_true, rho = 0,
                               q_mult = 1, a_mult = 1, C = 0, gu_true,
                               g_mult = 1)
  parity_treatment_true = parity(au_true, qu_true, tau_true, d_true, rho_true,
                                 q_mult_true, a_mult_true, C, gu_true,
                                 g_mult_true)
  
  # find the expected number of mosquitoes located inside a single house
  f_control_true = expected_catch_number_untreated(au_true, qu_true, tau_true, d_true, rho=0, q_mult=1,
                                                   a_mult=1, C=0, gu_true, g_mult=1, lambda_true,
                                                   catch_prop_true)
  f_treatment_true = expected_catch_number_treated(au_true, qu_true, tau_true,
                                                   d_true, 
                                                   rho_true, q_mult_true,
                                                   a_mult_true, 
                                                   C, gu_true,
                                                   g_mult_true, lambda_true,
                                                   catch_prop_true)
  
  # sample size for epidemiology- vector of followup times
  followup_days_control = epi %>%
    filter(Cluster_Allocation == "C") %>%
    pull(followup_time) *365
  
  followup_days_treatment = epi %>%
    filter(Cluster_Allocation == "T") %>%
    pull(followup_time) *365
  
  # get the probability of seroconverting for each person during this time
  prob_seroconverting_control = 1 - exp(-foi_control_true*followup_days_control)
  prob_seroconverting_treatment = 1 - exp(-foi_treatment_true*followup_days_treatment)
  
  # for each person, get a binomial yes/no outcome
  sc_data_control = data.frame(outcome = rbinom(length(followup_days_control), size=1,
                                                prob=prob_seroconverting_control),
                               followup_days = followup_days_control)
  sc_data_treatment = data.frame(outcome= rbinom(length(followup_days_treatment), size=1,
                                                 prob=prob_seroconverting_treatment),
                                 followup_days = followup_days_treatment)
  
  #### generate parity data
  parity_data_baseline = rbinom(n=1, size=sum(parity_baseline$total_caught), prob=parity_control_true)
  parity_data_intervention_control = rbinom(n=1,
                                            size=sum(parity_intervention[parity_intervention$cluster_trt == "C",]$total_caught),
                                            prob=parity_control_true)
  parity_data_intervention_treatment = rbinom(n=1,
                                              size=sum(parity_intervention[parity_intervention$cluster_trt == "T",]$total_caught),
                                              prob=parity_treatment_true)
  
  #### generate Prokopack data
  catch_data_baseline = rnbinom(nrow(abundance_baseline), mu = f_control_true, 
                                size=phi_true)
  catch_data_intervention_control= rnbinom(nrow(abundance_intervention_control),
                                           mu = f_control_true, 
                                           size=phi_true)
  catch_data_intervention_treatment= rnbinom(nrow(abundance_intervention_treatment),
                                             mu = f_treatment_true, 
                                             size=phi_true)
  
  #### generate blood feeding data
  
  feeding_data_baseline = rbinom(1, size=sum(bloodmeal_baseline$total_evaluated),
                                 prob=prop_fed_control_true)
  feeding_data_intervention_control = rbinom(1,
                                             size=sum(bloodmeal_intervention[bloodmeal_intervention$cluster_trt ==
                                                                               "C",]$total_evaluated),
                                             prob=prop_fed_control_true)
  feeding_data_intervention_treatment = rbinom(1,
                                               size=sum(bloodmeal_intervention[bloodmeal_intervention$cluster_trt ==
                                                                                 "T",]$total_evaluated),
                                               prob=prop_fed_treated_true)
  
  # return a list of all data
  return(list(
    # epidemiology 
    list(sc_data_control, sc_data_treatment),
    list(parity_data_baseline, parity_data_intervention_control, parity_data_intervention_treatment),
    list(catch_data_baseline, catch_data_intervention_control, catch_data_intervention_treatment),
    list(feeding_data_baseline, feeding_data_intervention_control, 
         feeding_data_intervention_treatment)
  ))
  
}

#### generate parameter values to generate from ####

N=40000
# These come from ten bosch et al. She only gave default values, so I used these
# as the mean with a pretty narrow 95% CI

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

a.mult_true_candidates = runif(N, 0.5, 1.5)
g.mult_true_candidates = runif(N, .5, 1.5)
rho_true_candidates = rbeta(N, 1, 1)
q.mult_true_candidates = runif(N, .5, 1.5)

# This is informed from "a new cost effective" article on Prokopack
catch_prop_true_candidates = rbeta(N, 1, 1)

# These are pretty flat because i don't know
lambda_true_candidates= rexp(N, 10)
phi_true_candidates = rexp(N,1)



candidate_df = data.frame(au_true=au_true_candidates, d_true=d_true_candidates,
                          gu_true=gu_true_candidates,
                          qu_true=qu_true_candidates, tau_true=tau_true_candidates,
                          n_true=n_true_candidates,
                          b_true=b_true_candidates, 
                          X_true=X_true_candidates, 
                          a.mult_true=a.mult_true_candidates, 
                          g.mult_true=g.mult_true_candidates,
                          rho_true=rho_true_candidates, q.mult_true=q.mult_true_candidates,
                          lambda_true = lambda_true_candidates, c_true = c_true_candidates,
                          phi_true = phi_true_candidates,
                          catch_prop_true = catch_prop_true_candidates)

# for each parameter set, calculate the mosquito locations and overall death and biting rates
candidate_df$foi = 0
for(i in 1:nrow(candidate_df)){
  print(i)
 
  candidate_df$foi[i] = with(candidate_df[i,], force_of_infection_untreated(au_true, qu_true,
                                                                  tau_true, d_true,
                                                                  rho = 0, q_mult = 1,
                                                                  a_mult = 1, C = 0,
                                                                  gu_true, g_mult = 1,
                                                                  b_true, lambda_true,
                                                                  c_true, X_true, n_true))
}


# take only those between 0.2 and 2
combos_to_try = candidate_df[candidate_df$foi > 0 & candidate_df$foi < 0.33/365,]

combos_to_try_ppp = combos_to_try %>%
  filter(a.mult_true > 1, g.mult_true > 1, q.mult_true > 1) %>%
  head(12)
combos_to_try_npp = combos_to_try %>%
  filter(a.mult_true < 1, g.mult_true > 1, q.mult_true > 1) %>%
  head(12)
combos_to_try_pnp = combos_to_try %>%
  filter(a.mult_true > 1, g.mult_true < 1, q.mult_true > 1) %>%
  head(14)
combos_to_try_nnp = combos_to_try %>%
  filter(a.mult_true < 1, g.mult_true < 1, q.mult_true > 1) %>%
  head(12)
combos_to_try_ppn = combos_to_try %>%
  filter(a.mult_true > 1, g.mult_true > 1, q.mult_true < 1) %>%
  head(12)


combos_to_try_npn = combos_to_try %>%
  filter(a.mult_true < 1, g.mult_true > 1, q.mult_true < 1) %>%
  head(12)
combos_to_try_pnn = combos_to_try %>%
  filter(a.mult_true > 1, g.mult_true < 1, q.mult_true < 1) %>%
  head(14)
combos_to_try_nnn = combos_to_try %>%
  filter(a.mult_true < 1, g.mult_true < 1, q.mult_true < 1) %>%
  head(12)

# take balanced combinations of SR effects
combos_to_try = rbind(combos_to_try_ppp, combos_to_try_npp, combos_to_try_pnp, 
                      combos_to_try_nnp, combos_to_try_ppn, combos_to_try_npn, 
                      combos_to_try_pnn, combos_to_try_nnn)

nrow(combos_to_try)
#### generate 100 datasets ####

simulated_datasets = list()
for(i in 1:100){
  print(i)
  simulated_datasets[[i]] = generate_data(au_true=combos_to_try$au_true[i],
                                          d_true=combos_to_try$d_true[i],
                                          gu_true=combos_to_try$gu_true[i],
                                          qu_true=combos_to_try$qu_true[i],
                                          tau_true=combos_to_try$tau_true[i],
                                          n_true=combos_to_try$n_true[i],
                                          b_true=combos_to_try$b_true[i],
                                          X_true=combos_to_try$X_true[i], 
                                          a_mult_true=combos_to_try$a.mult_true[i],
                                          g_mult_true=combos_to_try$g.mult_true[i],
                                          rho_true=combos_to_try$rho_true[i],
                                          q_mult_true=combos_to_try$q.mult_true[i],
                                          lambda_true=combos_to_try$lambda_true[i],
                                          c_true = combos_to_try$c_true[i],
                                          phi_true=combos_to_try$phi_true[i],
                                          catch_prop_true = combos_to_try$catch_prop_true[i])
}
save(simulated_datasets, combos_to_try, file="../../results/Simulation/100SimulatedDatasets.RData")


min(combos_to_try$gu_true)

