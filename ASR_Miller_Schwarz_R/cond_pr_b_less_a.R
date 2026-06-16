# cond_pr_b_less_a.R

# Conditional probability of B < A given RT = t.
# lambda differs for the congruent & incongruent conditions.
cond_pr_b_less_a <- function(t, tau_a, tau_b, mu_C, sigma_C, lambda, soa) {
  alpha <- 1 / tau_a
  beta <- 1 / tau_b
  alpha_plus_beta <- alpha + beta
  pr_b_less_a <- beta / alpha_plus_beta  # NEWJEFF: This only works for SOA = 0
  # result <- probability that B wins * density when B wins / total density regardless of winner
  result <- pr_b_less_a *
    f_rt_b_less_a(t, alpha = alpha, beta = beta, mu_C, sigma_C = sigma_C, lambda = lambda) /
    f_rt(t, alpha = alpha, beta = beta, mu_C, sigma_Cab = sigma_C, lambda, sigma_Cba = sigma_C, soa)
}
