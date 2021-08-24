# To run these examples:

# 1. Place file asr_sample_data.csv into the project's root directory
# 2. Source this file by clicking the "Source" button in the upper right of this pane
# 3. Call function asr_fit_demo()
# 4. Call function asr_predictions_demo()

###### Note re lambda_inh and lambda_exc:
#
# For prediction, the lambda parameters should be thought of (approximately) as
# changes in RT, in msec.  Thus, lambda_inh is a positive number, representing
# an increase in RT produced by inhibition from incongruent activation, and
# lambda_exc is a negative number, representing a decrease in RT produced by
# facilitation from congruent activation.
#
# During fitting, the system will find best estimates for lambda_inh and
# (optionally) lambda_exc. Should either of these values have the opposite
# of its expected positive/negative sign, it indicates that the data suggest
# the effect of the term was opposite to its usual direction. 

library(gamlss)
library(pracma)

source("asr_fit.R")
source("asr_distribution.R")
source("asr_predictions.R")

#################################################################################
# 1. Fitting an ASR model (i.e. estimating model parameters) from reaction time data
#################################################################################

#################################################################################
# asr_fit_demo: Date for three real subjects.
#################################################################################
asr_fit_demo <- function()
{
  # Allow some time for large data sets...
  print("Working...please wait")
  
  # Load some sample data to which the model can be fit
  test_data_df <- read.csv("asr_sample_data.csv", header = TRUE)
  
  # Fit each subject individually. Here, we fit subject 2
  subject_02 <- test_data_df[test_data_df$SubNo == 2,]
  
  # Set up separate vectors of the congruent and incongruent reaction times
  subject_02_con <- subject_02[subject_02$Congru == 1,]
  subject_02_inc <- subject_02[subject_02$Congru == 2,]
  
  # Initialise soa based on the experimental design from which
  # the data are derived
  soa <- 0
  
  # Sample subject 2. LambdaExc fixed
  print("Fixed lambda_exc - parameter estimates")
  
  # Function asr-fit will display the best model estimates (i.e. for starting
  # parameters that result in the smallest error) The output object of asr_fit
  # contains estimate and error values for all starting parameters as well as
  # the best estimates. User can do further work with those data, as desired.
  fixed_exc <- asr_fit(subject_02_con$RT, subject_02_inc$RT, soa)
  
  # Sample subject 2. LambdaExc allowed to vary
  print("Lambda_exc allowed to vary - parameter estimates")
  varying_exc <- asr_fit(subject_02_con$RT, subject_02_inc$RT, soa, estimate_excitation = TRUE)

}

#################################################################################
# 2. Generating predicted reaction time distributions and delta plots from 
#    specific values of the ASR model parameters
#################################################################################
# asr_predictions <- function(tau_A, tau_B, mu_C, sigma_C, lambda_inh, lambda_exc, soa,
#                             percentiles = seq(0.005,0.995,0.01),
#                             show_plots = TRUE)
##############################################################################
asr_predictions_demo <- function()
{
  tau_A <- 108
  tau_B <- 91
  mu_C <- 423
  sigma_C <- 47.5
  lambda_exc <- -10
  lambda_inh <- 67
  soa <- 25
  
  asr_pred_df <- asr_predictions(tau_A = tau_A, tau_B = tau_B, mu_C = mu_C, sigma_Cab = sigma_C, sigma_Cba = sigma_C, lambda_inh = lambda_inh, lambda_exc = lambda_exc, soa = soa)
  head(asr_pred_df)
}

