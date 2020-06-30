# Fitting the Activation Suppression Race Model (Miller & Schwarz, 2020)
# 
# This file contains the top-level method asr_fit and associated utility methods
# to derive the parameters of the ASR model for a set of input data. Data are 
# presumed to be reaction times (in msec) for one congruent and one incongruent
# condition (see Miller & Schwarz, 2020).



#################################################################################
# asr_fit(congruent_rts, incongruent_rts, soa, resolution = 0.2, estimate_excitation = FALSE)

# Given vectors of reaction times from congruent and incongruent trials,
# this method returns best estimates for the six parameters of the Activation
# Suppression Race model (Miller and Schwarz, 2020)

# Addition parameters:
# resolution: step size for iteration through values for proportion of variance
# due to the exponential component in estimating the starting parameters for 
# fminsearch and selects the set of parameters that produce the lowest fmin.
# Default values are {0.1, 0.3, 0.5, 0.7, 0.9} i.e. seq(0.1, 0.9, 0.2)
# estimate_excitation: If TRUE, fminsearch optimises lambda_exc; if FALSE, lambda_exc
# is fixed at 0. (See Miller and Schwarz, 2020)
#################################################################################
asr_fit <- function(congruent_rts, incongruent_rts, soa, 
                    resolution = 0.2, estimate_excitation = FALSE,
                    display_params = TRUE)
{
  
  library(gamlss) # for dexGAUS
  library(pracma) # for fminsearch
  
  # Prepare to track best (minimum) value of fmin (error) returned by fminsearch
  min_fmin <- Inf
  best_opt <- NA
  
  # Prepare to gather estimates and error for each attempt at fminsearh
  xmin_list <- list()
  fmin_list <- list()
  
  # For loop runs through a series of fminsearches, differing by their starting
  # "best guess" parameter estimates. Search accuracy can, for some subjects, be
  # very sensitive to these values.
  
  count <- 0
  
  # Starting parameters vary by the proportion of data variance attributed to the
  # exponential terms in the model.
  # Grain size of iterations determined by resolution parameter
  prop_exp_variance_values <- seq(0.1, 0.9, resolution)
  
  # For each putative proportion of variance due to the exponential components
  for (prop_exp_variance in prop_exp_variance_values)
  {
    count <- count+1
    
    # Returns five or six params (wo or w lambda_exc) depending on
    # estimate_excitation The resulting vector is passed to fminsearch, which
    # returns an estimate for each provided parameter. Thus fminsearch returns
    # five or six values depending on estimate_excitation.
    starting_params <- asr_starting_params(congruent_rts, 
                                           incongruent_rts,
                                           prop_exp_variance,
                                           estimate_excitation)
    
    # fminsearch over the solution space of asr_error (see function asr_error below)
    # asr_error checks the number of params it gets and if there is no value supplied
    # for lambda_exc (i.e. it gets only five), it fixes lambda_exc to 0
    opt <- fminsearch(asr_error, 
                      starting_params, 
                      congruent_data_to_fit = congruent_rts, 
                      incongruent_data_to_fit = incongruent_rts,
                      soa = soa)
    
   
     # Store estimates and final error of this fminsearch for return.
    xmin_list[[count]] <- opt$xmin
    fmin_list[[count]] <- opt$fmin
    
    # Track minimum error
    current_fmin <- opt$fmin
    if (current_fmin < min_fmin)
    {
      min_fmin <- current_fmin
      best_opt <- opt
    }
  } # for each value of exponential proportion of variance
  
  # Grab best estimates
  best_estimates <- best_opt$xmin
  
  
  # Check for that sixth parameter
  if (length(best_estimates) == 6) {
    lambda_exc <- best_estimates[6]
  } else {
    lambda_exc = 0
  }
  
  # fminsearch may have used negative values for tau_A, tau_B and sigma_C, which are
  # inverted by asr_error. To insure meaningful output values, we take the abs of those
  # parameters here.
  parameter_estimates <- list(tau_A = abs(best_estimates[1]), 
                              tau_B = abs(best_estimates[2]), 
                              mu_C = abs(best_estimates[3]), 
                              sigma_C = abs(best_estimates[4]), 
                              lambda_inh = best_estimates[5],
                              lambda_exc = lambda_exc)
  
  # Gather other elements for return
  output = list(BestParameters = parameter_estimates, 
                MinError = min_fmin, 
                AllParameterEstimates = xmin_list,
                AllErrors = fmin_list)
  
  # Display !!! Remove in prod
  if (display_params) {
    print(as.data.frame(parameter_estimates))
  }
  
  return(output)
}

#################################################################################
# Utility functions
#################################################################################

##############################################################################
# f_rt_b_less_a: Function describing the density of t when B is less than A
# Formula A.3 in the manuscript
###############################################################################
f_rt_b_less_a <- function(t,alpha, beta, mu_C, sigma_C, lambda)
{
  dexGAUS(t,
          mu = mu_C + lambda,
          sigma = sigma_C,
          nu = 1/(alpha + beta))
}

##############################################################################
# f_rt_b_greater_a: Function describing the density of t when B is greater
# than A. Formula A.8 in the manuscript
###############################################################################
f_rt_b_greater_a <- function(t, alpha, beta, mu_C, sigma_C, soa)
{
  exg_01 <- dexGAUS(t,
                    mu = mu_C,
                    sigma = sigma_C,
                    nu = 1/beta)
  exg_02 <- dexGAUS(t,
                    mu = mu_C,
                    sigma = sigma_C,
                    nu = 1/(alpha + beta))
  
 # ((1 + (beta/alpha)) * exg_01) - ((beta/alpha) * exg_02)
  
 mult1 <- (alpha + beta)/(alpha + (beta * (1-exp(-1 * alpha * soa))))
 mult2 <- (beta * exp(-1 * alpha * soa))/(alpha + (beta * (1-exp(-1 * alpha * soa))))
 
 (mult1 * exg_01) - (mult2 * exg_02)
}

##############################################################################
# f_rt_all: Function describing the overall density of t 
# Formula A.6 in the manuscript
###############################################################################
# lambda_inh for incongruent trials and lambda_exc for congruent trials
f_rt_all <- function(t, alpha, beta, mu_C, sigma_C, lambda, soa)
{
  q_s <- (beta/(alpha + beta)) * exp(-alpha * soa)
  f_plus <- f_rt_b_less_a(t, alpha, beta, mu_C, sigma_C, lambda)
  f_minus <- f_rt_b_greater_a(t, alpha, beta, mu_C, sigma_C, soa)
  (q_s * f_plus) + ((1-q_s) * f_minus)
}

#############################################################################
# asr_error(asr_params, congruent_data_to_fit, incongruent_data_to_fit)

# Given a set of estimates for the six parameters of an ASR model, and a data
# set, this function computes a measure of goodness of fit. This is the function
# optimised by fminsearch in asr_fit.

# Error is the negative log of the standard likelihood function.
#############################################################################
asr_error <- function(asr_params, congruent_data_to_fit, incongruent_data_to_fit, soa)
{
  # We take the absolutevalue of the incoming estimates of tau_A, tau_B and sigma_C 
  # to prevent negative values, for which the density functions are undefined.
  
  # Pull out params to increase readability. Transform where required.
  tau_A <- asr_params[1]
  tau_A <- abs(tau_A)
  
  tau_B <- asr_params[2]
  tau_B <- abs(tau_B)
  
  mu_C <- asr_params[3]
  
  sigma_C <- asr_params[4]
  sigma_C <- abs(sigma_C)
  
  lambda_inh <- asr_params[5] 
  
  # As determined by estimate_excitation flag in asr_fit
  if(length(asr_params) == 6){
    lambda_exc <- asr_params[6] # User should provide a negative lambda_exc for facilitation
  } else {
    lambda_exc <- 0
  }


  # The definitional functions use the rate of the exponential components, rather than
  # the means (tau). Prepare here -- rate = 1/mean
  alpha <- 1/tau_A
  beta <- 1/tau_B
  
  # Get predicted densities, using the appropriate lambda term
  # Passing in the vector to f_rt_all is faster than a for loop or sapply
  predicted_densities_congruent <- f_rt_all(congruent_data_to_fit,
                                         alpha = alpha,
                                         beta = beta,
                                         mu_C = mu_C,
                                         sigma_C = sigma_C,
                                         lambda = lambda_exc,
                                         soa = soa)
  
  
  predicted_densities_incongruent <- f_rt_all(incongruent_data_to_fit,
                                            alpha = alpha,
                                            beta = beta,
                                            mu_C = mu_C,
                                            sigma_C = sigma_C,
                                            lambda = lambda_inh,
                                            soa = soa)
 
  
  predicted_density_vector <- c(predicted_densities_congruent, predicted_densities_incongruent)
  
  # Log and sum to find likelihood
  log_predicted_density_vector <- log(predicted_density_vector)
  total_log_predicted_density <- sum(log_predicted_density_vector)
  
  # Flip sign, as we will search for minimum error value via fminsearch
  error <- -total_log_predicted_density
  
  return(error)
} # end asr_error


#############################################################################
# asr_starting_params(congruent_rts, incongruent_rts, prop_exp_variance, estimate_excitation)
#
# Heuristic for estimating good starting values for each model parameter to be 
# passed to fminsearch. Values derived from observed means and standard
# deviations of the input data. A range of such sets is tested during the
# fitting process (see asr_fit, above)
#############################################################################
asr_starting_params <- function(congruent_rts, incongruent_rts, prop_exp_variance, estimate_excitation)
{

  # Compute descriptives
  mean_con <- mean(congruent_rts)
  mean_inc <- mean(incongruent_rts)
  sd_con <- sd(congruent_rts)
  sd_inc <- sd(incongruent_rts)
  
  # We estimate tau_A as a proportion of the total variability in congruent trials
  est_tau_A <- prop_exp_variance * sd_con
  
  # We estimate tau_B as equal to tau_A
  est_tau_B <- est_tau_A
  
  # We estimate mu_C as the mean of the congruent trials - the mean of the exponential
  # component B
  est_mu_C <- mean_con - est_tau_B
  
  # We estimate sigma_C as the sqrt of the total variance minus that accounted for by the
  # exponential component.
  est_sigma_C <- sqrt((sd_con^2) - (est_tau_B^2))
  
  # We estimate lambda_inh as the mean inhibitory effect
  est_lambda_inh <- 2 * (mean_inc - mean_con) 
  
  # We estimate no excitatory effect
  est_lambda_exc <- 0
  

  # Gather the starting extimates into a vector
  params <- c(est_tau_A,
              est_tau_B,
              est_mu_C,
              est_sigma_C,
              est_lambda_inh)
  
  # Provide either five or six parameters to fminsearh depending on how many are to
  # be optimised (ess asr_fit, above)
  if (estimate_excitation) {
    params <- c(params, est_lambda_exc)
  }
  
  return(params)
}




