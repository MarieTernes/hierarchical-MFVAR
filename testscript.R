#####################################################################################################################
#-------------------- Test script for obtaining the hierarchical estimator & coincident indicator ------------------#
# "Hierarchical Regularizers for Mixed-Frequency Vector Autoregressions" by Alain Hecq, Marie Ternes and Ines Wilms #
#####################################################################################################################

cat("\014")       # Clear Console
graphics.off()    # Close graphs
rm(list=ls())     # Clean memory

# Load libraries
library(Matrix)
library(RSpectra)
library(spcov)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

# Source functions
source("sparseMFVAR_Functions.R")
source("inputs_Functions.R")
source("coincident_indicator.R")
sourceCpp("FunctionsRcpp.cpp")

#########################################################################
#-------------- Data example & Data pre-processing ---------------------#
#########################################################################
load("example_quarterly.RData")
load("example_monthly.RData")
data = make_data_matrices(data_quarterly = example_quarterly, data_monthly = example_monthly)
Y_VAR = data$Y_VAR
k_q = data$k_q
k_m = data$k_m
k_w = data$k_w
k_d = data$k_d
K = data$K

#############################################################################
#-------------- Obtaining the hierarchical regularizer ---------------------#
#############################################################################
MFHIER <- sparseMFVAR(Y_VAR = Y_VAR, k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, 
                      penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER",
                      inputsProx = NULL,l_lambda_beta = 10, standardize = TRUE, epsilon = 1e-3, 
                      max_iter = 400)
MFHIER$betas  

# Get optimal tuning parameter via time series cross validation 
# rolling window and one-step-ahead MSFE as a cross-validation score
CV <- MFVAR_cv(Y_VAR = Y_VAR, k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, 
               penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER",
               lambdas_beta = MFHIER$lambdas_beta, cvcut = 20, standardize = TRUE,
               epsilon = 1e-3, max_iter = 400)
CV$lambda_opt
CV$lambda_opt_oneSE
which(MFHIER$lambdas_beta == CV$lambda_opt)
which(MFHIER$lambdas_beta == CV$lambda_opt_oneSE)

#######################################################################################################
#------------------------------------- out-of-sample-exercise ----------------------------------------#
# To perform a "pseudo" out-of-sample exercise it is possible to use MFHIER as explained above. 
# If inputsProx = NULL (default) it is calculated within the function sparseMFVAR
# However, since the object inputsProx does not change across the samples, it is possible to increase
# the speed of the code and give inputsProx as an input to MFHIER.
# The object inputsProx can be calculated using the function inputs_sparseMFVAR
#######################################################################################################
inputsProx = inputs_sparseMFVAR(k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, 
                                penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER")


#################################################################################################
#--------------------------- Construction of the coincident indicator --------------------------#
#------------------------ for a variable in the Y_VAR time series matrix -----------------------#
#################################################################################################
ind = coincident_indicator(Y_VAR = Y_VAR, betas = MFHIER$betas, 
                          k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d,  
                          l_lambda_sigma = 10, col_coinci = 1)
ind$maxcor
ind$coincInd
ind$series_coincInd

# Plot standardized variables vs standardized coincident indicator
plot(ind$y, type ='l', col ='black', ylab = "", lwd = 1.7)
lines(ind$coincInd, col = 'blue', lwd = 1.7)
legend("bottomleft", bty = "n", legend = c("GDP growth", "Coincident indicator"), lty = c(1),col = c("black", "blue"))



  