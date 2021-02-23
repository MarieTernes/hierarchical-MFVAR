library(spcov)

coincident_indicator <- function(Y_VAR, betas, k_q, k_m, k_w, k_d, l_lambda_sigma, col_coinci, col_excl = NULL){
  # Input:
  # Y_VAR: TxK matrix of time series; matrix can be constructed from make_data_matrices function-> change this maybe later
  # betas: K^2 x l_lambda_beta matrix of estimated autoregressive coefficients of the MFVAR across different lambdas
  # k_q, k_m, k_w, k_d: number of quarterly, monthly, weekly and daily series
  # l_lambda_sigma: scalar, specifies how many values the tuning parameter grid 
  #                 for regularization of variance covariance matrix should contain
  # col_coinci: scalar that indicates the column of the low-frequency series in Y_VAR of 
  #             which one wishes to construct the coincident indicator
  # col_excl: scalar or vector (optional), indicates the column(s) of the series in Y_VAR that one wishes to exclude
  #           when constructing the coincident indicator, e.g. other low-frequency variables
  
  # Function: Investigate which high-frequency variables nowcast a low-frequency variable and 
  #           construct a coincident indicator from those
  
  # Output:
  # y: Standardized variable of interest for which the coincident indicator is constructed
  # coincInd: standardized coincident indicator 
  # maxcor: correlation achieved between y and coincInd
  # series_coincInd: list of variables selected for the construction of coincInd
  # Sigmas: Array with (ncol(betas) x l_lambda_sigma) regularized variance-covariance matrices
  # lambdas_sigma: (ncol(betas) x l_lambda_sigma) tuning parameter grid
  # l_index_beta_opt, l_index_sigma_opt: Optimal index of tuning parameter 
  #                                      according to maximum correlation between y and coincInd
  # correlation: (ncol(betas) x l_lambda_sigma) matrix that gives correlation 
  #              between y and all coincInd across the tuning parameter matrix
  # X_FPC_raw: Array with (ncol(betas) x l_lambda_sigma) unstandardized coincident indicators
  # PVE: Array with (ncol(betas) x l_lambda_sigma) proportion of variance explained
  
  p = 1
  
  Y_VAR = scale(Y_VAR)
  
  MFVARdata <- MFVARmodel(Y_VAR = Y_VAR, p = p)
  Y = MFVARdata$Y
  X = MFVARdata$X
  N = MFVARdata$N # T-p
  K = MFVARdata$K # total number of regressors
  l_lambda_beta = ncol(betas)
  
  varcov <- sparseCov(Y, X, betas, k_q, k_m, k_w, k_d, l_lambda_sigma, col_coinci)
  
  coincInd <- construct_coincInd(Y_VAR, varcov$Sigmas, k_q, k_m, k_w, k_d, col_coinci, col_excl)
  
  out <- list("y" = coincInd$y, "coincInd" = coincInd$x, "maxcor" = coincInd$maxcor, "series_coincInd" = coincInd$selected_opt,
              "Sigmas" = varcov$Sigmas, "lambdas_sigma" = varcov$lambdas_sigma, 
              "l_index_beta_opt" =  coincInd$l_index_beta_opt, "l_index_sigma_opt" = coincInd$l_index_sigma_opt,
              "correlation" = coincInd$correlation,"X_FPC_raw" = coincInd$X_FPC, "PVE" = coincInd$PVE, 
              "series_selected" = coincInd$selected_list)

  return(out)
}

sparseCov <- function(Y, X, betas, k_q, k_m, k_w, k_d, l_lambda_sigma, var_col){
  K = ncol(X)
  N = nrow(X)
  yhat <- apply(betas, 2, yhats_function, Xdata = X, K = K)
  P = offdiag_penalty_mat(k_q, k_m, k_w, k_d)
  l_lambda_beta = ncol(betas)
  lambda_sigma_index = round(seq(from = 1, to = (K-k_q), length = l_lambda_sigma))
  if(!is.null(colnames(Y))){ # time series names
    varnames <- colnames(Y)
  } else {
    varnames <- NULL
  }
  
  Sigmas = array(rep(NA, K^2*l_lambda_beta*l_lambda_sigma), c(K, K,l_lambda_beta,l_lambda_sigma))
  dimnames(Sigmas) = list(varnames,
                          varnames, 
                          paste0("lambda_beta", 1:l_lambda_beta),  paste0("lambda_sigma", 1:l_lambda_sigma))
  lambdas_sigma =  matrix(NA, l_lambda_beta, l_lambda_sigma)
  rownames(lambdas_sigma) = c(paste0("lambda_beta", 1:l_lambda_beta))
  colnames(lambdas_sigma) = c(paste0("lambda_sigma", 1:l_lambda_sigma))
  
  

  for(i in 1:l_lambda_beta){
    U <- Y - matrix(yhat[,i], N, K)
    varcov = cov(U)
    sol_path = sort(abs(varcov[((k_q+1):K),var_col]), decreasing = T)
    lambda_sigma_seq = sol_path[lambda_sigma_index]
    lambda_sigma_seq[l_lambda_sigma] =  lambda_sigma_seq[l_lambda_sigma]*0.99  #1% less so it's below min
    lambdas_sigma[i,] = lambda_sigma_seq
    
    for(j in 1:l_lambda_sigma){
      ADMM = ProxADMM(A = varcov, del = 0.005, lam = lambda_sigma_seq[j], P = P, maxiters = 200)
      Sigmas[,,i,j] = ADMM$Z
    }
  }
  return(list("Sigmas"=Sigmas, "lambdas_sigma"=lambdas_sigma))
}  

# Coincident indicator
construct_coincInd <- function(Y_VAR, Sigmas, k_q, k_m, k_w, k_d, col_coinci, col_excl){
  l_lambda_beta = dim(Sigmas)[3]
  l_lambda_sigma = dim(Sigmas)[4]
  K = ncol(Y_VAR)
  n = nrow(Y_VAR)
  #p = 1
  
  correlation = matrix(NA,l_lambda_beta, l_lambda_sigma) #correlation matrix
  rownames(correlation) = c(paste0("lambda_beta", 1:l_lambda_beta))
  colnames(correlation) = c(paste0("lambda_sigma", 1:l_lambda_sigma))
  X_FPC = array(rep(NA), c(n,l_lambda_beta,l_lambda_sigma), dimnames = list(paste0("t", 1:n), paste0("lambda_beta", 1:l_lambda_beta),  paste0("lambda_sigma", 1:l_lambda_sigma)))
  PVE = array(rep(NA), c(K-k_q,l_lambda_beta,l_lambda_sigma), dimnames = list(paste0("k", 1:(K-k_q)), paste0("lambda_beta", 1:l_lambda_beta),  paste0("lambda_sigma", 1:l_lambda_sigma)))
  
  posdef <- apply(Sigmas, c(3,4), min_eigen)
  selected_list = replicate(n=l_lambda_beta, expr=list())
  
  Y_VAR_std = scale(Y_VAR)
  
  for(i in 1:l_lambda_beta){
    for(j in 1:l_lambda_sigma){
      index = which(Sigmas[-c(col_coinci, col_excl),col_coinci,i,j] != 0)
      selected_list[[i]][[j]] = index
      if(length(index)>1 & posdef[i,j]>0){
        PCA = coincInd_aux(X_std = Y_VAR_std[,-c(col_coinci, col_excl)],  index_nonzero = index) #data_standardized[,-c(1:3)]
        X_FPC[,i,j] = PCA$x # First principal component score vector
        PVE[1:length(index),i,j] = PCA$pve  # Proportion of Variance Explained
        correlation[i,j] = cor(scale(X_FPC[,i,j]), Y_VAR_std[,col_coinci]) # Correlation of between series at interest and coincident indicator
      }
      if(length(index)==1 & posdef[i,j]>0){
        X_FPC[,i,j] = Y_VAR_std[,index]
        correlation[i,j] = cor(Y_VAR_std[,index], Y_VAR_std[,1])
      }
    }
  }
  
  
  lambdas_maxcorr = which(abs((correlation)) == max(abs((correlation)), na.rm = TRUE), arr.ind = TRUE)
  l_index_beta_opt = lambdas_maxcorr[1,1]
  l_index_sigma_opt = lambdas_maxcorr[1,2]
  maxcor = correlation[l_index_beta_opt,  l_index_sigma_opt]
  
  if(maxcor > 0){
    X_FPC_opt = scale(X_FPC[, l_index_beta_opt, l_index_sigma_opt])
  }else{
    X_FPC_opt = -scale(X_FPC[, l_index_beta_opt, l_index_sigma_opt])
  }
  
  PVE_opt = PVE[,l_index_beta_opt, l_index_sigma_opt]
  selected_list_opt = selected_list[[l_index_beta_opt]][[l_index_sigma_opt]]
  
  out <- list("y" =  Y_VAR_std[,col_coinci], "x" = X_FPC_opt,"maxcor" = abs(maxcor),"selected_opt" = selected_list_opt,
              "l_index_beta_opt" =  l_index_beta_opt, "l_index_sigma_opt" = l_index_sigma_opt,
              "correlation" = correlation,"X_FPC_raw" = X_FPC, "PVE" = PVE, 
              "selected_list" = selected_list)
    
  return(out)
  
}

# Penalty matrix of the same dimension as Sigma,  
# Should be used to penalize only off-diagonal elements
offdiag_penalty_mat <- function(k_q = 1, k_m = 1, k_w = 0, k_d = 0){
  m1 = 3
  m2 = 12
  m3 = 60
  K = k_q+k_m*m1+k_w*m2+k_d*m3
  P = matrix(1,K,K)
  
  if(k_q != 0){
    for(i in 1:k_q){
      P[i,i] = 0 
    }
  }
  if(k_m != 0){  
    for(i in 1:k_m){
      P[(k_q+1+m1*(i-1)):(k_q+m1*i),(k_q+1+m1*(i-1)):(k_q+m1*i)] = 0
    }
  }
  if(k_w != 0){  
    for(i in 1:k_w){
      P[(k_q+m1*k_m+1+m2*(i-1)):(k_q+m1*k_m+m2*i),(k_q+m1*k_m+1+m2*(i-1)):(k_q+m1*k_m+m2*i)] = 0
    }
  }
  if(k_d != 0){  
    for(i in 1:k_d){
      P[(k_q+m1*k_m+m2*k_w+1+m3*(i-1)):(k_q+m1*k_m+m2*k_w+m3*i),(k_q+m1*k_m+m2*k_w+1+m3*(i-1)):(k_q+m1*k_m+m2*k_w+m3*i)] = 0
    }
  }
  return(P)
}

# Auxiliary function for construct_coincInd
coincInd_aux <- function(X_std, index_nonzero){
  # Input
  # X_std: standardized TxK time series matrix
  # index_nonzero: vector that indicated columns of X that should be included in the PCA 
  X_selected = X_std[,index_nonzero]
  cor_mat = cor(X_selected) # correlation matrix of selected columns of X
  PC = eigen(cor_mat) # eigenvalues & eigenvectors
  PVE = PC$values/ sum(PC$values) # Proportion of Variance Explained
  FPC = PC$vectors[,1] # First principal component (first eigenvector)
  X_new = X_selected%*%FPC # First principal component score vector
  return(list("x"= X_new, "pve"= PVE, "fev" = FPC))
}

# Check if min eigenvalue of varcov is strictly positive
# Then varcov is positive definite
min_eigen <- function(varcov){
  ev <- eigen(varcov)
  V <- min(ev$values)
  return(V)
}



