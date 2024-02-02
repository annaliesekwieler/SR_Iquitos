library(tidyverse)
library(fitdistrplus)
library(BayesianTools)
library(postpack)
library(coda)
library(doParallel)
library(mvtnorm)

#### load functions and trial data ####

source("functions.R")

# load trial data 
bloodmeal_intervention = read.csv("bloodmeal_intervention.csv")


# take total full to be the sum of hald and full
bloodmeal_intervention = bloodmeal_intervention %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))


abundance_intervention_control = read.csv("abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("abundance_intervention_treatment.csv")

parity_intervention = read.csv("parity_intervention.csv")

epi = read.csv("epi.csv")

load("treatedclusters.RData")

C=.44
p=5.2



load("MultivariatePrior_All.RData")
prior_density = function(params){
  
  return(dmvnorm(params, mu_mvn, sigma_mvn, log=T))
}

prior_sampler = function(n=1){
  return(rmvnorm(n, mu_mvn, sigma_mvn))
}

#### function to perform inference on a data set ####

infer_parameters = function(simulated_data,  index){
  
  # unpack the simulated data
  sc_control_sim = simulated_data[[1]][[1]]
  sc_treatment_sim = simulated_data[[1]][[2]]
  
  
  parity_intervention_control_sim = simulated_data[[2]][[2]]
  parity_intervention_treatment_sim = simulated_data[[2]][[3]]
  # 
  
  catch_intervention_control_sim = simulated_data[[3]][[2]]
  catch_intervention_treatment_sim = simulated_data[[3]][[3]]
  # 
  
  feeding_intervention_control_sim = simulated_data[[4]][[2]]
  feeding_intervention_treatment_sim = simulated_data[[4]][[3]]
  
  likelihood = function(params){
    
    
    qu=exp(params[1])
    au=exp(params[2])
    tau=exp(params[3])
    d = exp(params[4])
    gu = sigmoid(params[5])
    lambda = exp(params[6])
    catch_prop = sigmoid(params[7])
    phi=exp(params[8])
    
    n=exp(params[9])
    X= sigmoid(params[10])
    c=sigmoid(params[11])
    b=sigmoid(params[12])
    
    rho = sigmoid(params[13])
    a_mult = exp(params[14])
    g_mult = exp(params[15])
    q_mult = exp(params[16])
    
    
    prop_fed_untreated = prop_bloodfed_untreated(au, qu, tau, d, rho, q_mult, a_mult, C =0)
    prop_fed_treated = prop_bloodfed_treated(au, qu, tau, d, rho,
                                             q_mult, a_mult, C)
    
    # expected catch number 
    f_control = expected_catch_number_untreated(au, qu, tau, d, 
                                                rho, q_mult, a_mult, C=0, gu, g_mult, lambda, catch_prop)
    f_treated = expected_catch_number_treated(au, qu, tau, d, rho, q_mult,
                                              a_mult, C, gu, g_mult, lambda,
                                              catch_prop)
    
    # parity rates
    parity_rate_control = parity(au, qu, tau, d, rho=0, q_mult = 1,
                                 a_mult = 1, C = 0, gu, g_mult= 1)
    parity_rate_treatment = parity(au, qu, tau, d,
                                   rho, q_mult, a_mult,
                                   C, gu, g_mult)
    
    # force of infection
    foi_c = force_of_infection_untreated(au, qu, tau, d, rho = 0,
                                         q_mult = 1, a_mult = 1, C = 0, gu,
                                         g_mult = 1, b, lambda, c, X,
                                         n)
    foi_t = force_of_infection_treated(au, qu, tau, d, rho, q_mult,
                                       a_mult, C, gu, g_mult, b,
                                       lambda, c, X, n)
    # 
    prob_not_sc_control = exp(-foi_c * 
                                epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
    prob_sc_control =1 - prob_not_sc_control
    # 
    prob_not_sc_treatment = exp(-foi_t * epi$followup_time[epi$Cluster_Allocation == "T"] * 365)
    prob_sc_treatment =1 - prob_not_sc_treatment
    # 
    # 
    # # take care of some zeroes
    # 
    
    if(is.nan(prop_fed_treated)){
      prop_fed_treated = .0000000001
    }
    
    
    if(prop_fed_treated == 0){
      prop_fed_treated = .0000000001
    }
    
    if(prop_fed_untreated == 0){
      prop_fed_untreated = .0000000001
    }
    # 
    # 
    if(f_treated == 0){
      f_treated = .0000000001
    }
    # 
    if(f_control == 0){
      f_untreated = .0000000001
    }
    # 
    # 
    # 
    prob_sc_control[prob_sc_control == 0] = 0.00000001
    
    prob_sc_treatment[prob_sc_treatment == 0] = 0.00000001
    # 
    prob_sc_control[prob_sc_control == 1] = .999999999
    prob_sc_treatment[prob_sc_treatment == 1] = .999999999
    # 
    bloodmeal_likelihood_intervention_treated = dbinom(feeding_intervention_treatment_sim,
                                                       size=sum(bloodmeal_intervention %>%
                                                                  filter(cluster_trt == "T") %>%
                                                                  pull(total_evaluated)),
                                                       prob=prop_fed_treated, log=T)
    
    bloodmeal_likelihood_intervention_control = dbinom(feeding_intervention_control_sim,
                                                       size=sum(bloodmeal_intervention %>%
                                                                  filter(cluster_trt == "C") %>%
                                                                  pull(total_evaluated)),
                                                       prob=prop_fed_untreated, log=T)
    
    
    abundance_likelihood_intervention_treated = sum(dnbinom(catch_intervention_treatment_sim,
                                                            mu=f_treated,
                                                            size=phi, log=T))
    abundance_likelihood_intervention_control = sum(dnbinom(catch_intervention_control_sim,
                                                            mu=f_control,
                                                            size=phi, log=T))
    
    parity_likelihood_intervention_control = dbinom(parity_intervention_control_sim,
                                                    size=sum(parity_intervention %>%
                                                               filter(cluster_trt == "C") %>%
                                                               pull(total_caught)),
                                                    prob=parity_rate_control, log=T)
    
    parity_likelihood_intervention_treatment = dbinom(parity_intervention_treatment_sim,
                                                      size=sum(parity_intervention %>%
                                                                 filter(cluster_trt == "T") %>%
                                                                 pull(total_caught)),
                                                      prob=parity_rate_treatment, log=T)
    
    epi_likelihood_control = sum(dbinom(sc_control_sim$outcome, size=1,
                                        prob=prob_sc_control, log=T))
    
    epi_likelihood_treatment = sum(dbinom(sc_treatment_sim$outcome, size=1,
                                          prob=prob_sc_treatment, log=T))
    
    return(
      bloodmeal_likelihood_intervention_treated +
        bloodmeal_likelihood_intervention_control +
        abundance_likelihood_intervention_treated + 
        abundance_likelihood_intervention_control +
        # parity_likelihood_intervention_control +
        # parity_likelihood_intervention_treatment +
        epi_likelihood_control + 
        epi_likelihood_treatment
    )
    
  }
  
  bayesianSetup = createBayesianSetup(likelihood, prior_density, prior_sampler)
  
  convergence = F
  
  output= runMCMC(bayesianSetup, settings=list(iterations = 100000),
                  sampler="DEzs")
  total_iterations_so_far = 100000
  
  while(!convergence){
    
    output= runMCMC(output, settings=list(iterations=1000000))
    gelman_output = gelmanDiagnostics(output)
    
    
    if (max(gelman_output$psrf[,2])< 1.1){
      convergence = T
      print("We have converged")
      
      result = getSample(output, start = 1)
      result = tail(result, 1000000)
      result = as.data.frame(result)
      
      
      # convert to natural units
      result_exp = exp(result)
      result_sigmoid = sigmoid(result)
      result = cbind(result_exp[,1:4], result_sigmoid[,5], result_exp[,6],
                     result_sigmoid[,7],result_exp[,8:9],
                     result_sigmoid[,10:13],result_exp[,14:16])
      names(result) =c("qu","au","tau","d", "gu", "lambda", "catch_prop", "phi",
                       "n", "X", "c", "b", "rho", "a.mult", "g.mult", "q.mult")
      # 
      
      
      total_iterations_so_far = total_iterations_so_far + 100000
      
    } else{
      total_iterations_so_far = total_iterations_so_far + 100000
      print(paste("Did not converge on", total_iterations_so_far, "iterations"))
      print(max(gelman_output$psrf[,2]))
      
      
    }
  }
  
  return(list(result, total_iterations_so_far))
}

#### set up and run ####


# import simulated datasets (big file)
load("100SimulatedDatasets.RData")

simulated_datasets = simulated_datasets[61:80]


#### run the algorithm on each simulated data set ####
number_cores = 20
registerDoParallel(cores=number_cores)
cl=makeCluster(20)
result = foreach(i=1:20,.errorhandling = "pass") %dopar% infer_parameters(simulated_data= simulated_datasets[[i]],
                                                                          index=i)
stopCluster(cl)

save(result, file=paste0("Simulation_InterventionOnly_NoParity_4.RData"))