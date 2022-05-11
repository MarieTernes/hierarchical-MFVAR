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
library(corrplot)
library(igraph)

# Source functions
source("sparseMFVAR_Functions_R1.R")
source("inputs_function_R1.R")
source("coincident_indicator_R1.R")
sourceCpp("FunctionsRcpp_R1.cpp")

#########################################################################
#-------------- Data example & Data pre-processing ---------------------#
#########################################################################
load("mfvar_small.RData")
data = mfvar_small
Y_VAR = data$Y_VAR[1:125,]
k_q = data$k_q
k_m = data$k_m
k_w = data$k_w
k_d = data$k_d
K = data$K
p = 1

#############################################################################
#-------------- Obtaining the hierarchical regularizer ---------------------#
#############################################################################
inputsProx = inputs_sparseMFVAR(k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p,
                                penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER")

# If inputsProx = NULL (default) it is calculated within the function sparseMFVAR
# However, for out-of-sample exercises since the object inputsProx does not change across the samples, it is possible to increase the speed of the code and give inputsProx as an input to MFHIER.
MFHIER <- sparseMFVAR(Y_VAR = Y_VAR, k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p, 
                      penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER", 
                      inputsProx = inputsProx, l_lambda_beta = 10, standardize = TRUE, epsilon = 1e-3, 
                      max_iter = 400)

# Get optimal tuning parameter via time series cross validation 
# rolling window and one-step-ahead MSFE as a cross-validation score
CV <- MFVAR_cv(Y_VAR = Y_VAR, k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p,
               penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER",
               inputsProx = inputsProx,
               lambdas_beta = MFHIER$lambdas_beta, cvcut = 20, standardize = TRUE,
               epsilon = 1e-3, max_iter = 400)
cv_opt = which(MFHIER$lambdas_beta == CV$lambda_opt)


#################################################################################################
#--------------------------- Construction of the coincident indicator --------------------------#
#------------------------ for a variable in the Y_VAR time series matrix -----------------------#
#################################################################################################
ind = coincident_indicator(Y_VAR = Y_VAR, betas = MFHIER$betas, 
                           k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p, 
                           l_lambda_sigma = 10, col_coinci = 1)

# Plot variables vs standardized coincident indicator
GDP = ts(ind$y, start = c(1987,4), end = c(2018,4), frequency = 4)
indicator = ts(ind$coincInd, start = c(1987,4), end = c(2018,4), frequency = 4)

# Small
plot(GDP, type ='l', main ="Coincident indicator from small MF-VAR", col ='black', ylab = "", lwd = 1.7)
lines(indicator, col = c('#2b8cbe'), lwd = 1.9)
legend("bottomleft", bty = "n", legend = c("GDP growth", "Coincident indicator"), lty = c(1), lwd = c(1.7,1.7), col = c("black", '#2b8cbe'))

  
