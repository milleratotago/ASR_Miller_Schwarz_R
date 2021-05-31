# Fitting the Activation Suppression Race Model (Miller & Schwarz, 2021)
# 
# This file contains the top-level method asr_fit and associated utility methods
# to derive the parameters of the ASR model for a set of input data. Data are 
# presumed to be reaction times (in msec) for one congruent and one incongruent
# condition (see Miller & Schwarz, 2021).



#################################################################################
# asr_fit(congruent_rts, incongruent_rts, soa, resolution = 0.2,
#         estimate_excitation = FALSE, TwoSigmas = FALSE, display_params = TRUE)

# Given vectors of reaction times from congruent and incongruent trials,
# this method returns best estimates for the 5-8 parameters of the indicated
# version of the Activation Suppression Race model (Miller and Schwarz, 2021)
# The 4 possible versions are distinguished according to:
#    estimate_excitation  TwoSigmas  Parameters
#         FALSE             FALSE      5:  tau_A, tau_B, mu_C, sigma_C, lambda_inh                                      (assume lambda_exc=0 and sigma_C=sigma_Cexc=sigma_Cinh=sigma_C)
#         TRUE              FALSE      6:  tau_A, tau_B, mu_C, sigma_C, lambda_exc, lambda_inh                          (assume                  sigma_C=sigma_Cexc=sigma_Cinh=sigma_C)
#         FALSE             TRUE       6:  tau_A, tau_B, mu_C, sigma_C, lambda_inh, sigma_Cinh                          (assume lambda_exc=0 and sigma_Cexc=sigma_C)
#         TRUE              TRUE       8:  tau_A, tau_B, mu_C, sigma_C, lambda_exc, sigma_Cexc, lambda_inh, sigma_Cinh

# Note:
#  sigma_C    is the std deviation of the normal when A wins the race and irrelevant info is suppressed
#  sigma_Cexc is the std deviation of the normal when B wins the race and irrelevant info is congruent
#  sigma_Cinh is the std deviation of the normal when B wins the race and irrelevant info is incongruent

# Additional parameters:
# resolution: step size for iteration through values for proportion of variance
# due to the exponential component in estimating the starting parameters for 
# fminsearch and selects the set of parameters that produce the lowest fmin.
# Default values are {0.1, 0.3, 0.5, 0.7, 0.9} i.e. seq(0.1, 0.9, 0.2)
# estimate_excitation: If TRUE, fminsearch optimises lambda_exc; if FALSE, lambda_exc
# is fixed at 0. (See Miller and Schwarz, 2021)
#################################################################################
asr_fit <- function(congruent_rts, incongruent_rts, soa, 
                    resolution = 0.2, estimate_excitation = FALSE,
                    TwoSigmas = FALSE, display_params = TRUE)
{
  
  #  library(gamlss) # originally included for dexGAUS, but this was replaced by exgauss for better control over exp overflow problems.
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
    
    # Returns five starting parameter estimates based on con/inc RTs and the proportion of
    # RT variance allocated to the exponential component.
    ests <- asr_starting_params(congruent_rts, 
                                incongruent_rts,
                                prop_exp_variance)
    
    if (estimate_excitation && TwoSigmas) {
      starting_params = c(ests[1], ests[2], ests[3], ests[4], 0, ests[4], ests[5], ests[4]);  # Start lambda_exc at 0
    }
    else if (!estimate_excitation && TwoSigmas) {
      starting_params = c(ests[1], ests[2], ests[3], ests[4], ests[5], ests[4]); # (assume lambda_exc=0 and sigma_Cexc=sigma_C)
    }
    else if ( estimate_excitation && !TwoSigmas) {
      starting_params = c(ests[1], ests[2], ests[3], ests[4], 0, ests[5]);  # Start lambda_exc at 0
    }
    else if (!estimate_excitation && !TwoSigmas) {
      starting_params = c(ests[1], ests[2], ests[3], ests[4], ests[5]);
    }
    
    # fminsearch over the solution space of asr_error1or2sigmas (see function below)
    # asr_error checks TwoSigmas and the number of params to see what version
    # of the model is being fit.
    opt <- fminsearch(asr_error1or2sigmas, 
                      starting_params, 
                      congruent_data_to_fit = congruent_rts, 
                      incongruent_data_to_fit = incongruent_rts,
                      soa = soa,
                      TwoSigmas)
    
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
  
  # fminsearch may have used negative values for tau_A, tau_B and sigma_C's, which are
  # inverted by asr_error. To insure meaningful output values, we take the abs of those
  # parameters here.
  best_estimates[1:4] <- abs(best_estimates[1:4])
  if (estimate_excitation && TwoSigmas) {
    best_estimates[6] <- abs(best_estimates[6]);  # sigma_Cexc
    best_estimates[8] <- abs(best_estimates[8]);  # sigma_Cinh
  }
  else if (!estimate_excitation && TwoSigmas) {
    best_estimates[6] <- abs(best_estimates[6]);  # sigma_Cinh
  }
  
  if (estimate_excitation && TwoSigmas) {
    parameter_estimates <- list(tau_A = best_estimates[1],
                                tau_B = best_estimates[2], 
                                mu_C = best_estimates[3], 
                                sigma_C = best_estimates[4], 
                                lambda_exc = best_estimates[5],
                                sigma_Cexc = best_estimates[6],
                                lambda_inh = best_estimates[7],
                                sigma_Cinh = best_estimates[8])
  }
  else if (!estimate_excitation && TwoSigmas) {
    parameter_estimates <- list(tau_A = best_estimates[1], 
                                tau_B = best_estimates[2], 
                                mu_C = best_estimates[3], 
                                sigma_C = best_estimates[4], 
                                lambda_inh = best_estimates[5],
                                sigma_Cinh = best_estimates[6])
  }
  else if (estimate_excitation && !TwoSigmas) {
    parameter_estimates <- list(tau_A = best_estimates[1], 
                                tau_B = best_estimates[2], 
                                mu_C = best_estimates[3], 
                                sigma_C = best_estimates[4], 
                                lambda_exc = best_estimates[5],
                                lambda_inh = best_estimates[6])
  }
  else if (!estimate_excitation && !TwoSigmas) {
    parameter_estimates <- list(tau_A = best_estimates[1], 
                                tau_B = best_estimates[2], 
                                mu_C = best_estimates[3], 
                                sigma_C = best_estimates[4], 
                                lambda_inh = best_estimates[5])
  }
  
  
  # Gather other elements for return
  output = list(BestParameters = parameter_estimates, 
                MinError = min_fmin, 
                AllParameterEstimates = xmin_list,
                AllErrors = fmin_list)
  
  # Display parameters for this subject
  if (display_params) {
    print(as.data.frame(parameter_estimates))
  }
  
  return(output)
}

#################################################################################
# Utility functions
#################################################################################

################################################################################
# pdf of ex-Gaussian distribution
################################################################################
dexgauss <- function(X,mu,sigma,tau)
{
  rate <- 1/tau
  t1 <- -X*rate + mu*rate + 0.5*(sigma*rate)^2
  t2 <- (X - mu - sigma^2*rate) / sigma
  thispdf <- rate*exp( t1 + log(pnorm(t2)) )
  # The above is the theoretical definition of the ex-Gaussian pdf,
  # but there are numerical problems if t1 > 708 or so, because then exp(t1) overflows.
  # That happens when tau is small relative to sigma (so (sigma*rate)^2 is large).
  # In that case, though, the exG density is very close to a normal with mean mu+tau & variance sigma^2+tau^2,
  # so we can just use that corresponding normal pdf in those cases to avoid numerical problems:
  t1bad <- t1 > 708;  # vector indicating too-large t1's for which we should use the normal approximation.
  thispdf[t1bad] <- dnorm(X[t1bad],mu+tau,sqrt(sigma^2+tau^2))
  return(thispdf)
}

##############################################################################
# f_rt_b_less_a: Function describing the density of t when B is less than A
# Formula A.3 in the manuscript
###############################################################################
f_rt_b_less_a <- function(t,alpha, beta, mu_C, sigma_C, lambda)
{
  dexgauss(t,
           mu = mu_C + lambda,
           sigma = sigma_C,
           tau = 1/(alpha + beta))
}

##############################################################################
# f_rt_b_greater_a: Function describing the density of t when B is greater
# than A. Formula A.8 in the manuscript
###############################################################################
f_rt_b_greater_a <- function(t, alpha, beta, mu_C, sigma_C, soa)
{
  exg_01 <- dexgauss(t,
                     mu = mu_C,
                     sigma = sigma_C,
                     tau = 1/beta)
  exg_02 <- dexgauss(t,
                     mu = mu_C,
                     sigma = sigma_C,
                     tau = 1/(alpha + beta))
  
  # ((1 + (beta/alpha)) * exg_01) - ((beta/alpha) * exg_02)
  
  mult1 <- (alpha + beta)/(alpha + (beta * (1-exp(-1 * alpha * soa))))
  mult2 <- (beta * exp(-1 * alpha * soa))/(alpha + (beta * (1-exp(-1 * alpha * soa))))
  
  (mult1 * exg_01) - (mult2 * exg_02)
}

##############################################################################
# f_rt: Function describing the overall density of t 
# Formula A.6 in the manuscript
###############################################################################
# lambda_inh for incongruent trials and lambda_exc for congruent trials
# sigma_Cab applies when A wins and irrelevant info is suppressed.
# sigma_Cba applies when B wins and irrelevant info has an effect of lambda.
f_rt <- function(t, alpha, beta, mu_C, sigma_Cab, lambda, sigma_Cba, soa)
{
  # q_s = Pr(A > B + SOA) = Pr(B finishes before A)
  q_s <- (beta/(alpha + beta)) * exp(-alpha * soa)
  # density when B finishes before A
  f_plus <- f_rt_b_less_a(t, alpha, beta, mu_C, sigma_Cba, lambda)
  # density when A finishes before B
  f_minus <- f_rt_b_greater_a(t, alpha, beta, mu_C, sigma_Cab, soa)
  thispdf <- (q_s * f_plus) + ((1-q_s) * f_minus)
  # Due to numerical problems with some extreme parameter combinations,
  # values of thispdf can sometimes be less than zero (which is impossible).
  # To circumvent that problem, we set pdf values <= 0 to very small positive values.
  badpdf <- thispdf <= 0
  thispdf[badpdf] = .Machine$double.xmin
  return(thispdf)
}



#############################################################################
# asr_error1or2sigmas(asr_params, congruent_data_to_fit, incongruent_data_to_fit, soa, TwoSigmas)

# Compute likelihood of congruent & incongruent RTs as f(parameters in asr_params)
#  for 4 different ASR models depending on:
#    estimate_excitation or not (indicated by the number of parameters in asr_params)
#    whether or not SigmaC changes when there is excitation/inhibition associated with lambda (indicated by TwoSigmas)
# TwoSigmas is a Boolean indicating whether to use 1 or 2 sigmas:
# TwoSigmas FALSE version uses a single SigmaC.
#   asr_params can have 5 or 6 parameters:
#   5: TauA, TauB, MuC, SigmaC, LambdaInh,    (assume lambdaExc=0)
#   6: TauA, TauB, MuC, SigmaC, LambdaExc, LambdaInh
# TwoSigmas TRUE version uses two SigmaC's for each condition.
#   asr_params can have 6 or 8 parameters:
#   6: TauA, TauB, MuC, SigmaC, LambdaInh, SigmaCinh (assume lambdaExc=0 and SigmaCexc=SigmaC)
#   8: TauA, TauB, MuC, SigmaC, LambdaExc, SigmaCexc, LambdaInh, SigmaCinh
#
# If LambdaExc is used, negative values mean that there is facilitation (speeding) of RT
# Note using abs() of TauA, TauB, MuC, SigmaC's but not Lambda's

#############################################################################
asr_error1or2sigmas <- function(asr_params, congruent_data_to_fit, incongruent_data_to_fit, soa, TwoSigmas)
{
  
  # Pull out params to increase readability. Transform where required.
  tau_A <- abs(asr_params[1])
  tau_B <- abs(asr_params[2])
  mu_C <- abs(asr_params[3])
  sigma_C <- abs(asr_params[4])
  
  if (TwoSigmas) {  
    if(length(asr_params) == 6){
      lambda_exc <- 0
      sigma_Cexc <- sigma_C
      lambda_inh <- asr_params[5] 
      sigma_Cinh <- abs(asr_params[6])
    } else if (length(asr_params) == 8){
      lambda_exc <- asr_params[5] # User should provide a negative lambda_exc for facilitation
      sigma_Cexc <- abs(asr_params[6])
      lambda_inh <- asr_params[7] 
      sigma_Cinh <- abs(asr_params[8])
    } else {
      stop("Wrong number of parameters: asr_error1or2sigmas requires exactly 6 or 8 parameters when TwoSigmas is true.") 
    } # end TwoSigmas
  } else {
    # one sigma
    sigma_Cexc <- sigma_C
    sigma_Cinh <- sigma_C
    if(length(asr_params) == 5){
      lambda_exc <- 0
      lambda_inh <- asr_params[5] 
    } else if (length(asr_params) == 6){
      lambda_exc <- asr_params[5] # User should provide a negative lambda_exc for facilitation
      lambda_inh <- asr_params[6] 
    } else {
      stop("Wrong number of parameters: asr_error1or2sigmas requires exactly 5 or 6 parameters when TwoSigmas is false.") 
    }
  } # end if (TwoSigmas)
  
  # The definitional functions use the rate of the exponential components, rather than
  # the means (tau). Prepare here -- rate = 1/mean
  alpha <- 1/tau_A
  beta <- 1/tau_B
  
  # Get predicted densities, using the appropriate lambda term
  # Passing in the vector to f_rt2sigmas is faster than a for loop or sapply
  predicted_densities_congruent <- f_rt(congruent_data_to_fit,
                                        alpha = alpha,
                                        beta = beta,
                                        mu_C = mu_C,
                                        sigma_Cab = sigma_C,
                                        lambda = lambda_exc,
                                        sigma_Cba = sigma_Cexc,
                                        soa = soa)
  
  
  predicted_densities_incongruent <- f_rt(incongruent_data_to_fit,
                                          alpha = alpha,
                                          beta = beta,
                                          mu_C = mu_C,
                                          sigma_Cab = sigma_C,
                                          lambda = lambda_inh,
                                          sigma_Cba = sigma_Cinh,
                                          soa = soa)
  
  
  predicted_density_vector <- c(predicted_densities_congruent, predicted_densities_incongruent)
  
  # Log and sum to find likelihood
  log_predicted_density_vector <- log(predicted_density_vector)
  total_log_predicted_density <- sum(log_predicted_density_vector)
  
  # Flip sign, as we will search for minimum error value via fminsearch
  error <- -total_log_predicted_density
  
  return(error)
} # end asr_error1or2sigmas


#############################################################################
# asr_starting_params(congruent_rts, incongruent_rts, prop_exp_variance)
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
  tau_A <- prop_exp_variance * sd_con
  
  # We estimate tau_B as equal to tau_A
  tau_B <- tau_A
  
  # We estimate mu_C as the mean of the congruent trials - the mean of the exponential
  # component B
  mu_C <- mean_con - tau_B
  
  # We estimate sigma_C as the sqrt of the total variance minus that accounted for by the
  # exponential component.
  sigma_C <- sqrt((sd_con^2) - (tau_B^2))
  
  # We estimate lambda_inh as the mean inhibitory effect
  lambda_inh <- 2 * (mean_inc - mean_con) 
  
  # Gather the starting extimates into a vector
  params <- c(tau_A,
              tau_B,
              mu_C,
              sigma_C,
              lambda_inh)
  
  return(params)
}

