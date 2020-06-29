##########################################################
# Corresponding to R distribution family structure:
#
# dasr(t, ...) returns the densities of t_con and t_incon
# pasr(t, ...) returns the cumulative densities at t_con and t_incon
# qasr(p, ...) returns the values of t_con and t_incon at cumulative
#              density p
#
# Additional params for each method are tau_A, tau_B, mu_C,
# sigma_C, lambda_inh and lambda_exc
###########################################################
source("asr_fit.R")

# lambda_exc should be supplied by the user as a negative value for
# facilitiation (i.e. speeding reaction times)
dasr <- function(t, tau_A, tau_B, mu_C, sigma_C, lambda_inh, lambda_exc, soa)
{
  alpha <- 1/tau_A
  beta <- 1/tau_B
  d_con <- f_rt_all(t, alpha, beta, mu_C, sigma_C, lambda_exc, soa)
  d_incon <- f_rt_all(t, alpha, beta, mu_C, sigma_C, lambda_inh, soa)
  result <- c(d_con = d_con, d_incon = d_incon)
  return(result)
}

# lambda_exc should be supplied by the user as a negative value for
# facilitiation (i.e. speeding reaction times)
pasr <- function(t, tau_A, tau_B, mu_C, sigma_C, lambda_inh, lambda_exc, soa)
{
  alpha <- 1/tau_A
  beta <- 1/tau_B
  cdf_con <- integrate(function(t) {f_rt_all(t, alpha, beta, mu_C, sigma_C, lambda_exc, soa)},
                   lower = 1, upper = t)
  cdf_incon <- integrate(function(t) {f_rt_all(t, alpha, beta, mu_C, sigma_C, lambda_inh, soa)},
                       lower = 1, upper = t)
  result <- c(p_con = cdf_con$value, p_incon = cdf_incon$value)
  return(result)
}



f_to_zero <- function(t, target_p, tau_A, tau_B, mu_C, sigma_C, lambda, soa)
{
  alpha <- 1/tau_A
  beta <- 1/tau_B
  cdf_t <- integrate(function(t) {f_rt_all(t, alpha, beta, mu_C, sigma_C, lambda, soa)},
                                  lower = 1, upper = t)
  return(cdf_t$value - target_p)
}

# lambda_exc should be supplied by the user as a negative value for
# facilitiation (i.e. speeding reaction times)
qasr <- function(p, tau_A, tau_B, mu_C, sigma_C, lambda_inh, lambda_exc, soa)
{
  estimate <- tau_B + mu_C
  inv_cdf_con <- fzero(function(t) {f_to_zero(t, p, tau_A, tau_B, mu_C, sigma_C, lambda_exc, soa)},
                 x = estimate)
  inv_cdf_incon <- fzero(function(t) {f_to_zero(t, p, tau_A, tau_B, mu_C, sigma_C, lambda_inh, soa)},
                 x = estimate)
  result <- c(q_con <- inv_cdf_con$x, q_incon <- inv_cdf_incon$x)
  return(result)
}


