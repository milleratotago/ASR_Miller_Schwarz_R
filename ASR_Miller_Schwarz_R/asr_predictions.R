# !!! Move plotting into its own method 23-06-20



source("asr_fit.R")
source("asr_distribution.R")

# lambda_exc should be supplied by the user as a negative value for
# facilitiation (i.e. speeding reaction times)

asr_predictions <- function(tau_A = 100, tau_B = 100, mu_C = 200, sigma_Cab = 40, sigma_Cba = 40,
                            lambda_inh = 50, lambda_exc = -10, soa = 0,
                            percentiles = seq(0.005,0.995,0.01),
                            show_plots = TRUE)
{
  
  # Compute t values
  
  # sapply returns a 1 x n_percentile matrix. 1st row is congruent, 2nd row is incongruent
  quantiles <- sapply(percentiles, FUN = qasr, tau_A, tau_B, mu_C, sigma_Cab, sigma_Cba, lambda_inh, lambda_exc, soa)
  quantiles_con <- quantiles[1,]
  quantiles_incon <- quantiles[2,]
  asr_df <- data.frame(p = percentiles, 
                       t_con = quantiles_con, 
                       t_incon = quantiles_incon)
  
  # We have different t values to inspect for con and incon
  
  # Get densities for the congruent t values (quantiles)
  pdf_t_con <- sapply(quantiles_con, FUN = dasr, tau_A, tau_B, mu_C, sigma_Cab, sigma_Cba, lambda_inh, lambda_exc, soa)
  pdf_con <- pdf_t_con["d_con",]
  
  # Get densities for the incongruent t values (quantiles)
  pdf_t_incon <- sapply(quantiles_incon, FUN = dasr, tau_A, tau_B, mu_C, sigma_Cab, sigma_Cba, lambda_inh, lambda_exc, soa)
  pdf_incon <- pdf_t_incon["d_incon", ] 
  
  # Add the new columns
  asr_df$pdf_con <- pdf_con
  asr_df$pdf_incon <- pdf_incon
  
  # Compute the means and deltas
  asr_df$mean_t <- (asr_df$t_con + asr_df$t_incon)/2
  asr_df$delta_t <- asr_df$t_incon - asr_df$t_con

  if (show_plots) {
    
    par(mfrow = c(3,1))
    
    # PDF plot
    # !!! Adjust lims
    plot(asr_df$t_con, asr_df$pdf_con, type = "l", xlab = "RT", ylab = "PDF")
    lines(asr_df$t_incon, asr_df$pdf_incon, col = 2)
    legend("topright", legend=c("Congruent", "Incongruent"),
           col=c(1, 2), lty=c(1,1))
    
    # CDF plot
    # !!! Adjust lims
    plot(asr_df$t_con, asr_df$p, type = "l", xlab = "RT", ylab = "CDF")
    lines(asr_df$t_incon, asr_df$p, col = 2)
    legend("bottomright", legend=c("Congruent", "Incongruent"),
           col=c(1, 2), lty=c(1,1))
    
    # Delta plot
    plot(asr_df$mean_t, asr_df$delta_t, type = "l",
         xlab = "RT", ylab = "Delta")
  }
  
  return(asr_df)
} # end asr_predictions


