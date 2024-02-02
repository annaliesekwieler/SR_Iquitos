
##### Dengue
#### load trial data ####


setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/code")
source("functions.R")

setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/data")
# load trial data 
bloodmeal_intervention = read.csv("bloodmeal_intervention.csv")
bloodmeal_baseline = read.csv("bloodmeal_baseline.csv")

# take total full to be the sum of hald and full
bloodmeal_intervention = bloodmeal_intervention %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

bloodmeal_baseline = bloodmeal_baseline %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

abundance_baseline = read.csv("abundance_baseline.csv")
abundance_intervention_control = read.csv("abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("abundance_intervention_treatment.csv")

parity_baseline = read.csv("parity_baseline.csv")
parity_intervention = read.csv("parity_intervention.csv")

epi = read.csv("epi.csv")

load("treatedclusters.RData")

C=.44
p=5.2


#### load posterior samples ####


setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/TrialData")
load("Intervention_Only.RData")
head(result)


result_thinned = list()
for(i in 1:length(result)){
  result_temp = result
  
  result_thinned[[i]] = rbind(result_temp[seq(1, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 200),],
                              result_temp[seq(2, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 200),],
                              result_temp[seq(3, nrow(result_temp), by=3), ][seq(1, nrow(result_temp)/3,
                                                                                 by = 200),])
}

result = result_thinned[[1]]


#### estimate foi from data ####

# total followup time
control_fut = epi %>%
  filter(Cluster_Allocation == "U")

##### code to plot the quantities vs the trial data for a single ####
# parameter set

ppd = function(qu, au, tau, d, gu, lambda, catch_prop, phi, 
               n, X, c, b, rho, a_mult, g_mult, q_mult){
  
  
  # calculate the quantities for control and treatment 
  parity_control_post = parity(au, qu, tau, d, rho, q_mult, a_mult, C=0,
                               gu, g_mult)
  
  parity_treatment_post = parity(au, qu, tau, d, rho, q_mult, a_mult, C,
                               gu, g_mult)
  
  expected_catch_control_post = expected_catch_number_untreated(au, qu, tau, d, rho,
                                                                q_mult,
                                                                a_mult, C=0,
                                                                gu, g_mult, lambda,
                                                                catch_prop)
  expected_catch_treatment_post = expected_catch_number_treated(au, qu, tau,
                                                                d, rho, q_mult,
                                                                a_mult, C, gu,
                                                                g_mult, lambda,
                                                                catch_prop)
  
 bloodfed_control_post = prop_bloodfed_untreated(au, qu, tau, d, rho,
                                                 q_mult, a_mult, C=0)
 
 bloodfed_treatment_post = prop_bloodfed_treated(au, qu, tau, d, rho, q_mult,
                                                 a_mult, C)
  
  foi_control_post = force_of_infection_untreated(au, qu, tau, d, rho, q_mult,
                                        a_mult, C = 0, gu, g_mult, b,
                                        lambda, c, X, n)
  
  foi_treatment_post = force_of_infection_treated(au, qu, tau, d, rho, q_mult,
                                        a_mult, C, gu, g_mult, b,
                                        lambda, c, X, n)
  
  # convert this to probability of seroconverting during the trial
  prop_sc_control_post = 1 - exp(-foi_control_post * epi$followup_time[epi$Cluster_Allocation == "C"] * 365)
  
  
  par(mfrow=c(2,2))
  
  plot(.5, parity_control_post, col = "red", pch=19,
       xlim=c(0, 3), xaxt="n", main="Parity", xlab="",
       ylim=c(0,1), ylab = "Parity")
  points(x=1, y=sum(parity_intervention$parous[parity_intervention$cluster_trt == "C"]) / 
           sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "C"]),
         col="black", pch=19)
  
  points(2, parity_treatment_post, col = "red", pch=19,
         xlim=c(0, 3), xaxt="n")
  points(x=2.5, y=sum(parity_intervention$parous[parity_intervention$cluster_trt == "T"]) / 
           sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "T"]),
         col="black", pch=19)
  abline(v = 1.5, lty="dashed")
  mtext("Control",1, at = 1)
  mtext("Treatment",1, at =2)
  
  plot(.5, bloodfed_control_post, col = "red", pch=19,
       xlim=c(0, 3), xaxt="n", main="Proportion Bloodfed", xlab="",
       ylim=c(0,1), ylab="Prop Bloodfed")
  points(x=1, y=sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "C"]) /
           sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "C"]),
         col="black", pch=19)
  
  points(2, bloodfed_treatment_post, col = "red", pch=19,
         xlim=c(0, 3), xaxt="n")
  points(x=2.5, y=sum(bloodmeal_intervention$total_fed[bloodmeal_intervention$cluster_trt == "T"]) /
           sum(bloodmeal_intervention$total_evaluated[bloodmeal_intervention$cluster_trt == "T"]),
         col="black", pch=19)
  abline(v = 1.5, lty="dashed")
  mtext("Control",1, at = 1)
  mtext("Treatment",1, at =2)
  
  plot(.5, expected_catch_control_post, col = "red", pch=19,
       xlim=c(0, 3), xaxt="n", main="Aspiration",
       ylim=c(0, 3), ylab="Expected Catch")
  
  points(x=1, y=mean(abundance_intervention_control$total_caught),
         col="black", pch=19)
  
  points(1.75, expected_catch_treatment_post, col = "red", pch=19)
  points(x=2.25, y=mean(abundance_intervention_treatment$total_caught),
         col="black", pch=19)
  abline(v = 1.5, lty="dashed")
  mtext("Control",1, at = 1)
  mtext("Treatment",1, at =2)
  
  
  plot(.5, foi_control_post, col = "red", pch=19,
       xlim=c(0, 3), xaxt="n", main="Force of Infection",
       ylim=c(0, .002), ylab="Force of Infection")
  points(x=1, y=foi_control_data,
         col="black", pch=19)
  
  points(1.75,foi_treatment_post, col = "red", pch=19)
  points(x=2.25, y=foi_treatment_data,
         col="black", pch=19)
  abline(v = 1.5, lty="dashed")
  mtext("Control",1, at = 1)
  mtext("Treatment",1, at =2)
  
  legend("topright", c("Posterior", "Data"),
         col=c("red", "black"), pch=19)
  
  
  
  
  
}

#### likelihood function ####



likelihood = function(qu, au, tau, d, gu, lambda, catch_prop, phi, n, X, 
                      c, b, rho, a_mult, g_mult, q_mult){

  
  
  prop_fed_untreated = prop_bloodfed_untreated(au, qu, tau, d, rho = 0,
                                               q_mult = 1, a_mult = 1,
                                               C = 0)
  prop_fed_treated = prop_bloodfed_treated(au, qu, tau, d, rho,
                                           q_mult, a_mult, C)
  
  # expected catch number 
  f_control = expected_catch_number_untreated(au, qu, tau, d, rho, q_mult,
                                              a_mult, C=0, gu, g_mult, lambda, catch_prop)
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
  bloodmeal_likelihood_intervention_treated = dbinom(sum(bloodmeal_intervention %>%
                                                           filter(cluster_trt == "T") %>%
                                                           pull(total_fed)),
                                                     size=sum(bloodmeal_intervention %>%
                                                                filter(cluster_trt == "T") %>%
                                                                pull(total_evaluated)),
                                                     prob=prop_fed_treated, log=T)
  
  bloodmeal_likelihood_intervention_control = dbinom(sum(bloodmeal_intervention %>%
                                                           filter(cluster_trt == "C") %>%
                                                           pull(total_fed)),
                                                     size=sum(bloodmeal_intervention %>%
                                                                filter(cluster_trt == "C") %>%
                                                                pull(total_evaluated)),
                                                     prob=prop_fed_untreated, log=T)
  
  
  abundance_likelihood_intervention_treated = sum(dnbinom(abundance_intervention_treatment$total_caught,
                                                          mu=f_treated,
                                                          size=phi, log=T))
  abundance_likelihood_intervention_control = sum(dnbinom(abundance_intervention_control$total_caught,
                                                          mu=f_control,
                                                          size=phi, log=T))
  
  parity_likelihood_intervention_control = dbinom(sum(parity_intervention %>%
                                                        filter(cluster_trt == "C") %>%
                                                        pull(parous)),
                                                  size=sum(parity_intervention %>%
                                                             filter(cluster_trt == "C") %>%
                                                             pull(total_caught)),
                                                  prob=parity_rate_control, log=T)
  
  parity_likelihood_intervention_treatment = dbinom(sum(parity_intervention %>%
                                                          filter(cluster_trt == "T") %>%
                                                          pull(parous)),
                                                    size=sum(parity_intervention %>%
                                                               filter(cluster_trt == "T") %>%
                                                               pull(total_caught)),
                                                    prob=parity_rate_treatment, log=T)
  
  epi_likelihood_control = sum(dbinom(epi$outcome[epi$Cluster_Allocation == "C"], size=1,
                                      prob=prob_sc_control, log=T))
  
  epi_likelihood_treatment = sum(dbinom(epi$outcome[epi$Cluster_Allocation == "T"], size=1,
                                        prob=prob_sc_treatment, log=T))
  
  return(
    bloodmeal_likelihood_intervention_treated +
      bloodmeal_likelihood_intervention_control +
      abundance_likelihood_intervention_treated + 
      abundance_likelihood_intervention_control +
      parity_likelihood_intervention_control +
      parity_likelihood_intervention_treatment +
      epi_likelihood_control + 
      epi_likelihood_treatment
  )
  
}



#### find the highest likelihood posterior sample ####

# likelihoods = rep(NA, nrow(result))
# 
# for(i in 1:nrow(result)){
#     res = result[i,]
#     likelihoods[i] = likelihood(res$qu, res$au, res$tau, res$d, res$gu, res$lambda, res$catch_prop, res$phi, res$n, res$X, 
#                                   res$c, res$b, res$rho, res$a_mult, res$g_mult, res$q_mult)
# }
# 
# which.max(likelihoods)

random_row = 312
result_random = result[random_row,]
qu = as.numeric(result_random["qu"])
au = as.numeric(result_random["au"])
tau = as.numeric(result_random["tau"])
d = as.numeric(result_random["d"])
gu = as.numeric(result_random["gu"])
lambda = as.numeric(result_random["lambda"])
catch_prop = as.numeric(result_random["catch_prop"])
phi = as.numeric(result_random["phi"])
n = as.numeric(result_random["n"])
X = as.numeric(result_random["X"])
c = as.numeric(result_random["c"])
b = as.numeric(result_random["b"])
rho = as.numeric(result_random["rho"])
a_mult = as.numeric(result_random["a_mult"])
g_mult = as.numeric(result_random["g_mult"])
q_mult = as.numeric(result_random["q_mult"])


#### Build shiny app ####

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Posteriors of Parity, HBR, Prevalence and Incidence"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("qu_slider", "Exit Rate, Untreated",
                  min=0, max=3, value=qu, step = .0001),
      sliderInput("au_slider", "Biting Rate, Untreated",
                  min=0, max=5, value=au, step = .0001),
      sliderInput("tau_slider", "Entrance Rate",
                  min=0, max=5, value=tau, step = .0001),
      sliderInput("d_slider", "Digestion Rate",
                  min=0, max=3, value=d, step = .0001),
      sliderInput("gu_slider", "Mortality Rate, Untreated",
                  min=0, max=1, value=gu, step = .0001),
      sliderInput("lambda_slider", "Mosquito Emergence",
                  min=0.01, max=150, value=lambda, step = .0001),
      sliderInput("catch_prop_slider", "Aspiration Catch Proportion",
                  min=0, max=.99, value=catch_prop, step = .0001),
      sliderInput("phi_slider", "Dispersion, Aspiration",
                  min=0, max=3, value=phi, step = .0001),
      sliderInput("n_slider", "EIP",
                  min=10, max=25, value=n, step = .0001),
      sliderInput("X_slider", "Prevalence",
                  min=0, max=1, value=X, step = .0001),
      sliderInput("c_slider", "Transmission Probability, Human to Mosquito",
                  min=0, max=1, value=c, step = .0001),
      sliderInput("b_slider", "Transmission Probability, Mosquito to Human",
                  min=0, max=1, value=b, step = .0001),
      sliderInput("rho_slider", "Repellency",
                  min=0, max=1, value=rho, step = .0001),
      sliderInput("a_mult_slider", "SR Effect, Biting",
                  min=0.1, max=2, value=a_mult, step = .0001),
      sliderInput("g_mult_slider", "SR Effect, Mortality",
                  min=0, max=2, value=g_mult, step = .0001),
      sliderInput("q_mult_slider", "SR Effect, Exiting",
                  min=.1, max=2, value=q_mult, step = .0001),
      actionButton("reset", "Reset")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput("parameter_set"),
      plotOutput("ppdplot", height="800px"),
      textOutput("overall_likelihood")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  output$ppdplot = renderPlot({
    ppd(input$qu_slider, input$au_slider, input$tau_slider, input$d_slider, input$gu_slider,
        input$lambda_slider, input$catch_prop_slider, input$phi_slider, input$n_slider,
        input$X_slider, input$c_slider, input$b_slider, input$rho_slider, input$a_mult_slider,
        input$g_mult_slider, input$q_mult_slider)
  })
  
  output$overall_likelihood = renderText({
    paste("Overall Likelihood: ", likelihood(input$qu_slider, input$au_slider, input$tau_slider, input$d_slider, input$gu_slider,
                                               input$lambda_slider, input$catch_prop_slider, input$phi_slider, input$n_slider,
                                               input$X_slider, input$c_slider, input$b_slider, input$rho_slider, input$a_mult_slider,
                                               input$g_mult_slider, input$q_mult_slider))
  })
  
  observeEvent(input$reset,{
    updateSliderInput(session,'qu_slider',value = qu)
    updateSliderInput(session,'au_slider',value = au)
    updateSliderInput(session,'tau_slider',value = tau)
    updateSliderInput(session,'d_slider',value = d)
    updateSliderInput(session,'gu_slider',value = gu)
    updateSliderInput(session,'lambda_slider',value = lambda)
    updateSliderInput(session,'catch_prop_slider',value = catch_prop)
    updateSliderInput(session,'phi_slider',value = phi)
    updateSliderInput(session,'n_slider',value = n)
    updateSliderInput(session,'X_slider',value = X)
    updateSliderInput(session,'c_slider',value = c)
    updateSliderInput(session,'b_slider',value = b)
    updateSliderInput(session,'rho_slider',value =rho)
    updateSliderInput(session,'a_mult_slider',value = a_mult)
    updateSliderInput(session,'g_mult_slider',value = g_mult)
    updateSliderInput(session,'q_mult_slider',value = q_mult)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)