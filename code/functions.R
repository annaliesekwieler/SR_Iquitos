#### Intro ####

# This script contains the functions for the model


#### Load packages ####
library(msm)
library(ohenery)

#### Main mosquito location model ####

location_probability = function(au, qu, tau, d, rho, q_mult, a_mult, C) {
  inf.matrix = matrix(c(
    
    # untreated, blood fed
    -(d+qu), d, qu, 0, 0, 0,
    
    # untreated, not blood fed
    au, -(au + qu), 0, qu, 0, 0,
    
    # Transition, Blood Fed
    tau * (1-C), 0, -(tau*(1-C) + tau * (1-rho) * C + d), d, tau * C * (1-rho), 0,
    
    # Transition, not blood fed
    0, tau * (1-C), 0, -(tau*(1-C) + tau*C*(1-rho)), 0, tau * C * (1-rho),
    
    # Treated, blood fed
    0, 0, q_mult*qu, 0, -(qu*q_mult + d), d,
    
    # Treated, not blood fed
    0, 0, 0, q_mult * qu, a_mult * au, -(q_mult * qu + a_mult*au)
    
  ), nrow=6, ncol=6,byrow=T)
  prob_transition_matrix = MatrixExp(inf.matrix)
  eigen_vector = eigen(t(prob_transition_matrix))$vector[,1]
  normalized_eigenvector = as.numeric(normalize(eigen_vector))
  
  normalized_eigenvector[normalized_eigenvector < 10e-5] = 0
  
  return(normalized_eigenvector)
}

pi_U = function(au, qu, tau, d, rho, q_mult, a_mult, C){
  return(sum(location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)[1:2]))
}

pi_T = function(au, qu, tau, d, rho, q_mult, a_mult, C){
  return(sum(location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)[5:6]))
}



#### mosquito biting and mortality rates ####

ac_calc = function(au, qu, tau, d, rho, q_mult, a_mult, C){
  
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  pi.U = sum(locations[1:2])
  pi.T = sum(locations[5:6])
  
  au = au * d / (au + d)
  at = a_mult * au * d / ((a_mult*au) + d)
  
  return(pi.U*au + pi.T*at)
}


gc_calc = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult){
  
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  pi.U = sum(locations[1:2])
  pi.tau = sum(locations[3:4])
  pi.T = sum(locations[5:6])
  
  gt = gu * g_mult
  return(pi.T * gt + pi.U * gu + pi.tau*gu)
}


#### ratio of mosquitoes to humans ####

m = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult, lambda){
  return(lambda / gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult))
}

# expected catch number for aspirations

expected_catch_number_untreated = function(au, qu, tau, d, rho, q_mult,
                                           a_mult, C, gu, g_mult, lambda, catch_prop){
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  m = lambda / g
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  return(((m * 5.2 * sum(locations[1:2])) / (1-C)) * catch_prop)
}

expected_catch_number_treated = function(au, qu, tau, d, rho, q_mult,
                                         a_mult, C, gu, g_mult, lambda, catch_prop){
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  m = lambda / g
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  return(((m * 5.2 * sum(locations[5:6])) / C) * catch_prop)
}


m_u = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult, lambda){
  
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  m = lambda / g
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  if(C == 1){
    return(0)
  } else{
    return(((m * sum(locations[1:2])) / (1-C)))}
  
}

m_t = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult, lambda){
  
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  m = lambda / g
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  return(((m * sum(locations[5:6])) / C))
}

#### epidemiological function ####

force_of_infection_untreated = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu,
                                        g_mult, b, lambda, c, X, n){
  
  a = au*d / (au + d)
  
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  
  
  return((b * m_u(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult, lambda) *
            (a ^ 2) * c * X * exp(-g*n)) / ((g + (a*c*X))))
}



force_of_infection_treated = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu,
                                      g_mult, b, lambda, c, X, n){
  
  at = a_mult*au*d / (a_mult*au + d)
  
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  
  return((b * m_t(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult, lambda) *
            (at ^ 2) * c * X * exp(-g*n)) / ((g + (at*c*X))))
}

#### mosquito parity ####

parity = function(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult){
  
  a = ac_calc(au, qu, tau, d, rho, q_mult, a_mult, C)
  
  g = gc_calc(au, qu, tau, d, rho, q_mult, a_mult, C, gu, g_mult)
  
  return(a / (a+g))
}

#### proportion bloodfed ####

prop_bloodfed_untreated = function(au, qu, tau, d, rho, q_mult, a_mult, C){
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  return(locations[1]/sum(locations[1:2]))
}

prop_bloodfed_treated = function(au, qu, tau, d, rho, q_mult, a_mult, C){
  locations = location_probability(au, qu, tau, d, rho, q_mult, a_mult, C)
  return(locations[5]/sum(locations[5:6]))
}


# sigmoid and inverse sigmoid for transforming (0-1) and [0,1] parameters

sigmoid = function(x){
  1/(1+exp(-x))
}

inverse_sigmoid = function(x){
  log(x / (1-x))
}