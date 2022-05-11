# R functions 
sparseMFVAR <- function(Y_VAR, k_q, k_m, k_w, k_d, p = 1,
                        penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER",
                        inputsProx = NULL,
                        #data_quarterly, data_monthly = NULL, data_weekly = NULL, data_daily = NULL, 
                        l_lambda_beta = 10, standardize = TRUE, epsilon = 1e-3, max_iter = 400){
  # Inputs:
  # Y_VAR: TxK matrix of time series; matrix can be constructed from make_data_matrices function-> change this maybe later
  # k_q: Number of quarterly series
  # k_m: # of monthly series
  # k_w: # of weekly series
  # k_d: # of daily series
  # p: autoregressive lag order of the MF-VAR (p = 1 default)
  # penalty_H_on_L: Penalty High Frequency on Low Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # penalty_L_on_H: Penalty Low Frequency on High Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # penalty_own_on_own: Penalty High/Low Frequency on High/Low Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # inputsProx: output object from inputs_sparseMFVAR function, if not supplied program computes it
  # l_lambda_beta: scalar, specifies how many values the tuning parameter grid should contain
  # standardize: whether to standardize the data (default = TRUE)
  # epsilon: a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
  # max_iter: scaler, specifies maximum number of iterations in proximal gradient algorithm
  
  # Function: Sparse Estimation of the mixed-frequency Vector AutoRegressive (MFVAR) Model
  
  # Output: 
  # Y_VAR: TxK matrix of time series (used for estimation, standardized if standarized = TRUE)
  # Y_VAR_raw: original TxK matrix of time series (unstandardized)
  # betas: Matrix of estimated autoregressive coefficients of the MFVAR across different lambdas
  #        coefficients in column 1 correspond to the sparsest solution, in the last column to the most dense solution
  # lambdas_betas: tuning parameter grid
  # K: Number of time series
  # k_q = Number of quarterly series
  # k_m = # of monthly series
  # k_w = # of weekly series
  # k_d = # of daily series
  # series_names: names of times series in Y_VAR
  # iter: vector containing iterations until convergence for each lambda 
  
  # Preliminaries
  # data = make_data_matrices(data_quarterly, data_monthly, data_weekly, data_daily)
  # Y_VAR = data$Y_VAR
  if(is.null(inputsProx)){
    inputsProx = inputs_sparseMFVAR(k_q, k_m, k_w, k_d, p,
                                    penalty_H_on_L = penalty_H_on_L, penalty_L_on_H = penalty_L_on_H, penalty_own_on_own = penalty_own_on_own)
    
  }
  # inputsProx = inputs_sparseMFVAR(data$k_q, data$k_m, data$k_w, data$k_d,
  #                                 penalty_H_on_L = penalty_H_on_L, penalty_L_on_H = penalty_L_on_H, penalty_own_on_own = penalty_own_on_own)
  
  
  # - 1 to adjust for the different indexing in Rcpp
  B_matrix_1 = inputsProx$B_matrix_1 
  B_matrix_2 = inputsProx$B_matrix_2
  index = inputsProx$index-1
  penalty_matrix = inputsProx$penalty_matrix-1
  columns_matrix = inputsProx$columns_matrix-1
  rows_L1 = inputsProx$rows_L1-1
  rows_HIER = inputsProx$rows_HIER-1
  indices_B_L1 = inputsProx$indices_B_L1-1 
  indices_B_HIER = inputsProx$indices_B_HIER-1
  column_L1 = inputsProx$column_L1-1 
  column_HIER = inputsProx$column_HIER-1
  K = inputsProx$K
  
  
  # standardize data 
  if(standardize == TRUE){
    Y_VAR_raw = Y_VAR
    Y_VAR = scale(Y_VAR)
  }else{
    Y_VAR_raw = Y_VAR
  }
  
  # Construct necessary data matrices
  MFVARdata <- MFVARmodel(Y_VAR = Y_VAR, p = p)
  Y = MFVARdata$Y
  X = MFVARdata$X
  y = MFVARdata$y 
  Xkron = MFVARdata$Xkron
  y_Xkron = MFVARdata$y_Xkron
  XtranspX_kron = MFVARdata$XtranspX_kron
  s = MFVARdata$s
  
  # Construct tuning parameter grid
  lambda_beta_seq = lambdagrid_beta(epsilon2 = .01, 
                                    Y, X, y_Xkron, XtranspX_kron,
                                    l_lambda_beta,
                                    B_matrix_1, B_matrix_2, index,
                                    penalty_matrix,
                                    columns_matrix,
                                    rows_L1, rows_HIER,
                                    indices_B_L1, indices_B_HIER,
                                    column_L1, column_HIER,
                                    K, p, s, epsilon = 10^-1, max_iter2 = 50)
  
  # MFVAR estimation across entire tuning parameter sequence (using warm-starts)
  PLS_results = PLS_Cpp(y_Xkron = y_Xkron, XtranspX_kron = XtranspX_kron,
                        B_matrix_1 = B_matrix_1, B_matrix_2 = B_matrix_2, index = index,
                        penalty_matrix = penalty_matrix,
                        columns_matrix = columns_matrix,
                        rows_L1 = rows_L1, rows_HIER = rows_HIER,
                        indices_B_L1 = indices_B_L1, indices_B_HIER = indices_B_HIER,
                        column_L1 = column_L1, column_HIER = column_HIER,
                        K = K, p = p, s = s, lambda = lambda_beta_seq, epsilon = epsilon, max_iter = max_iter)
  
  betas_PLS = PLS_results$betas_PLS
  iter_PLS = PLS_results$iter_PLS
  
  
  out <- list("Y_VAR" = Y_VAR,"Y_VAR_raw" = Y_VAR_raw, 
              "betas" = betas_PLS, "lambdas_beta"= lambda_beta_seq,
              "K" = K, "k_q" = k_q, "k_m" = k_m, "k_w" = k_w, "k_d" = k_d, 
              "iter"= iter_PLS)
  return(out)
  
}


# CV with rolling window and one-step-ahead mean-squared forecast error (MSFE) as a cross-validation score
MFVAR_cv <- function(Y_VAR, k_q, k_m, k_w, k_d, p = 1,
                     penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER",
                     inputsProx = NULL, 
                     lambdas_beta, cvcut, standardize = TRUE,
                     epsilon = 1e-3, max_iter = 400){
 
  
  # Inputs:
  # Y_VAR: TxK matrix of time series; matrix can be constructed from make_data_matrices function-> change this maybe later
  # k_q, k_m, k_w, k_d: number of quarterly, monthly, weekly and daily series
  # p: autoregressive lag order of the MF-VAR (p = 1 default)
  # penalty_H_on_L: Penalty High Frequency on Low Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # penalty_L_on_H: Penalty Low Frequency on High Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # penalty_own_on_own: Penalty High/Low Frequency on High/Low Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty)
  # lambdas_beta: vector, tuning parameter grid 
  # cvcut: number of observations used for forecast evaluation in the time series cross-validation procedure. The remainder is used for model estimation.
  # standardize: whether to standardize the data (default = TRUE)
  # epsilon: a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
  # max_iter: scaler, specifies maximum number of iterations in proximal gradient algorithm
  
  # Function: Time series cross-validation with rolling window and one-step-ahead mean-squared forecast error (MSFE) as a cross-validation score
  
  # Output:
  # lambdas_beta: tuning parameter grid
  # MSFE_avg: MSFE cross-validation scores for each value of the sparsity parameter in the considered grid
  # MSFE_all: MSFE cross-validation full output
  # lambda_opt: Optimal value of the sparsity parameter as selected by the time-series cross-validation procedure
  # lambda_optSE: Optimal value of the sparsity parameter as selected by the time-series cross-validation procedure and after applying the one-standard-error rule

  
  t <- nrow(Y_VAR)
  N1 <- t-cvcut
  lambda_beta_seq = lambdas_beta
  MSFEmatrix = matrix(NA, nrow = cvcut, ncol = length(lambdas_beta))
  rownames(MSFEmatrix) = paste0("n", 1:cvcut)
  colnames(MSFEmatrix) = paste0("lambdas_beta", 1:length(lambdas_beta))
  
  if(is.null(inputsProx)){
    inputsProx = inputs_sparseMFVAR(k_q, k_m, k_w, k_d, p,
                                    penalty_H_on_L = penalty_H_on_L, penalty_L_on_H = penalty_L_on_H, penalty_own_on_own = penalty_own_on_own)
  }
  
  # - 1 to adjust for the different indexing in Rcpp
  B_matrix_1 = inputsProx$B_matrix_1 
  B_matrix_2 = inputsProx$B_matrix_2
  index = inputsProx$index-1
  penalty_matrix = inputsProx$penalty_matrix-1
  columns_matrix = inputsProx$columns_matrix-1
  rows_L1 = inputsProx$rows_L1-1
  rows_HIER = inputsProx$rows_HIER-1
  indices_B_L1 = inputsProx$indices_B_L1-1 
  indices_B_HIER = inputsProx$indices_B_HIER-1
  column_L1 = inputsProx$column_L1-1 
  column_HIER = inputsProx$column_HIER-1
  K = inputsProx$K
  
  for(n in 1:cvcut){
    # standardize data 
    if(standardize == TRUE){
      Y_VARaux = scale(Y_VAR[(n:(N1+n-1)),])
      Y_VARaux_mean = colMeans(Y_VAR[(n:(N1+n-1)),])
      Y_VARaux_sd = apply(Y_VAR[(n:(N1+n-1)),], 2, sd)
      
      MFVARdata <- MFVARmodel(Y_VAR = Y_VARaux, p = p)
      Ytrain = MFVARdata$Y
      Xtrain = MFVARdata$X
      y = MFVARdata$y 
      Xkron = MFVARdata$Xkron
      y_Xkron = MFVARdata$y_Xkron
      XtranspX_kron = MFVARdata$XtranspX_kron
      s = MFVARdata$s
      
      Ytest = (Y_VAR[N1+n,]-Y_VARaux_mean)/Y_VARaux_sd
      Xtest = embed(Y_VARaux, dimension = p + 1)[nrow(Ytrain),c(1:(K*p))] 
    }else{
      Y_VARaux = Y_VAR[(n:(N1+n-1)),]
  
      MFVARdata <- MFVARmodel(Y_VAR = Y_VARaux, p = p)
      Ytrain = MFVARdata$Y
      Xtrain = MFVARdatal$X
      y = MFVARdata$y 
      Xkron = MFVARdata$Xkron
      y_Xkron = MFVARdata$y_Xkron
      XtranspX_kron = MFVARdata$XtranspX_kron
      s = MFVARdata$s
      
      Ytest = Y_VAR[N1+n,]
      Xtest = embed(Y_VARaux, dimension = p + 1)[nrow(Ytrain),c(1:(K*p))]  
    }
    
    # MFVAR estimation across entire tuning parameter sequence (using warm-starts)
    PLS_results = PLS_Cpp(y_Xkron = y_Xkron, XtranspX_kron = XtranspX_kron,
                          B_matrix_1 = B_matrix_1, B_matrix_2 = B_matrix_2, index = index,
                          penalty_matrix = penalty_matrix,
                          columns_matrix = columns_matrix,
                          rows_L1 = rows_L1, rows_HIER = rows_HIER,
                          indices_B_L1 = indices_B_L1, indices_B_HIER = indices_B_HIER,
                          column_L1 = column_L1, column_HIER = column_HIER,
                          K = K, p = p, s = s, lambda = lambda_beta_seq, epsilon = epsilon, max_iter = max_iter)
    
    # CV score
    yhat_out_PLS <- apply(PLS_results$betas_PLS, 2, yhats_function, Xdata = Xtest, K = K, p = p)
    SFE_out_PLS <- apply(yhat_out_PLS, 2, SFE, ydata = c(Ytest)) # matrix K x l_lambda_beta 
    MSFEmatrix[n,] = colMeans(SFE_out_PLS) 
  }
  
  CV_score = colMeans(MSFEmatrix)

  min_CV = which.min(CV_score)
  se_min_CV = sd(MSFEmatrix[,min_CV])/sqrt(cvcut)
  SE_rule = CV_score[min_CV]+se_min_CV
  sparse_CV = which(CV_score <= SE_rule, arr.ind = TRUE)[1]

  lambda_opt <- lambda_beta_seq[min_CV]
  lambda_opt_oneSE <- lambda_beta_seq[sparse_CV]
  
  # Output
  out <- list("lambdas_beta"=lambda_beta_seq,  
              "lambda_opt"=lambda_opt,"lambda_opt_oneSE"=lambda_opt_oneSE,
              "MSFE_avg"=CV_score, "MSFE_all"=MSFEmatrix)
  return(out)
}



MFVARmodel<-function(Y_VAR, p = 1){
  # Preliminaries
  K <- ncol(Y_VAR) # Number of Endogenous variables
  
  if(!is.null(colnames(Y_VAR))){ # time series names
    varnames <- colnames(Y_VAR)
  } else {
    varnames <- NULL
  }
  
  # Lagged predictor matrices
  data_lags = embed(Y_VAR, dimension = p+1)
  Y = data_lags[, 1:K] 
  X = data_lags[, -c(1:K)]
  colnames(Y) = varnames
  colnames(X) = rep(varnames, p)
  
  # Create Xkron (dim: TKxK^2)
  I = Diagonal(K)
  Xkron <- kronecker(I, X) # sparse dgCMatrix
  # Create y vector (dim: TK x 1)
  y = c(Y) # Y1 - Y2 - .... - YN of length (n-p)*K x 1
  
  # Save values for calculation of gradient (saves recalculating them each time)
  y_Xkron = as.matrix(t(y)%*%Xkron) # 1 x K^2
  XtranspX_kron = kronecker(I, crossprod(X)) # K^2 x K^2 
  
  # step size for proximal gradient algorithm
  Xsvd = svds(X, k = 1, nu = 0, nv = 0)$d
  s = (Xsvd)^(-2)
  
  N = nrow(Y) # N = T-p
  
  out<-list("Y"=Y, "X"=X,
            "y" = y, "Xkron" = Xkron, "y_Xkron" = y_Xkron, "XtranspX_kron" = XtranspX_kron,
            "K"=K, "N" = N,"p"=p, "s" = s)
  return(out)
}




# Construction of tuning parameter grid
lambdagrid_beta <- function(epsilon2 = .01, 
                            Y, X, y_Xkron, XtranspX_kron,
                            l_lambda_beta,
                            B_matrix_1, B_matrix_2, index,
                            penalty_matrix,
                            columns_matrix,
                            rows_L1, rows_HIER,
                            indices_B_L1, indices_B_HIER,
                            column_L1, column_HIER,
                            K, p, s, epsilon = 10^-1, max_iter2 = 50){

  mat = t(X)%*%Y
  frob_norm = rep(NA,K)
  #Frobenius norm between your design matrix X and each of the response variables and take the largest value of those.
  for(k in 1:K){
    frob_norm[k] = norm(cbind(mat[,k]), type = "F")
  }
  lambda_max_start = max(frob_norm)*2
  lambda_max = find_lambda_max(lambda_max_start = lambda_max_start, epsilon2 = epsilon2, 
                               y_Xkron = y_Xkron, XtranspX_kron = XtranspX_kron,
                               B_matrix_1 = B_matrix_1, B_matrix_2 = B_matrix_2, index = index,
                               penalty_matrix = penalty_matrix,
                               columns_matrix = columns_matrix,
                               rows_L1 = rows_L1, rows_HIER = rows_HIER,
                               indices_B_L1 = indices_B_L1, indices_B_HIER = indices_B_HIER,
                               column_L1 = column_L1, column_HIER = column_HIER,
                               K = K, p = p, s = s, epsilon = epsilon, max_iter = max_iter2)
  lambda_min = lambda_max/10^5
  lambda_beta_seq <- c(exp(seq(log(lambda_max),log(lambda_min), length = l_lambda_beta)))
  lambda_beta_seq[l_lambda_beta] = 0 #only possible if K^2 < N
  
  return(lambda_beta_seq)
}

# Find the lambda that sets all parameters equal to 0
find_lambda_max <- function(lambda_max_start, epsilon2 = 1, y_Xkron, XtranspX_kron,
                            B_matrix_1, B_matrix_2, index, 
                            penalty_matrix,
                            columns_matrix, 
                            rows_L1, rows_HIER,
                            indices_B_L1, indices_B_HIER, 
                            column_L1, column_HIER,
                            K, p, s,
                            epsilon = 10^-3, max_iter = 100){
  # Inputs
  # lambda_max_start: starting value for lambda to try out 
  # epsilon2: a small positive numeric value giving the tolerance for convergence between the difference of lambda_right (sparser) and lambda_left (denser)
  
  lambda_left = 0
  lambda_right = lambda_max_start
  
  while(lambda_right - lambda_left > epsilon2){
    result <- algorithm1_Cpp_sparse(y_Xkron = y_Xkron, XtranspX_kron = XtranspX_kron,
                                    B_matrix_1 = B_matrix_1, B_matrix_2 = B_matrix_2, index = index,
                                    penalty_matrix = penalty_matrix,
                                    columns_matrix = columns_matrix, 
                                    rows_L1 = rows_L1, rows_HIER = rows_HIER,
                                    indices_B_L1 = indices_B_L1, indices_B_HIER = indices_B_HIER, 
                                    column_L1 = column_L1, column_HIER = column_HIER,
                                    K = K, p = p, s = s, 
                                    lambda = lambda_right, epsilon = epsilon, max_iter = max_iter)
    if(sum(result$B) == 0){
      lambda_works = lambda_right
      lambda_right = (lambda_left+lambda_right)/2
    }else{
      lambda_left = lambda_right
      lambda_right = lambda_right*1.5
    }
  }
  return(lambda_works)
}

# Fitted values 
yhats_function <- function(beta_vec, Xdata, K, p){
  beta_matrix <- matrix(beta_vec, K*p, K)
  yhat <- Xdata%*%beta_matrix
  return(c(yhat)) # returns vectorized version
}

# Squared forecast error
SFE <- function(ydata, yhat){
  (ydata-yhat)^2
}

# Standard error of the mean
sem <- function(x){
  se <- sd(x)/sqrt(length(x))
  return(se)
}

MSEbetas <- function(beta, betahat){
  mean((betahat - beta)^2) 
}

# False Positive/Non-zero Rate
FPR  <- function(beta, betahat){
  # Positive = non-zero parameter
  # false Positive : zero estimated parameter (but actually it's actually non-zero/positive) !!false negative
  
  # Positive = non-zero parameter
  # Negative = zero parameter
  # false Negative : zero estimated parameter (but actually it's actually non-zero/positive) 
  
  beta_logical <- (beta!=0) + (betahat==0) #false positive , false negative
  FPRvalue <- sum(beta_logical==2)/sum(beta!=0)
  return(FPRvalue)
}

# False Negative / Zero Rate
FNR  <- function(beta, betahat){
  # Negative = a zero parameter
  # false Negative : non-zero estimated parameter (but actually it's actually zero)
  
  # Positive = non-zero parameter
  # Negative = zero parameter
  # false Positive : non-zero estimated parameter (but actually it's actually zero/negative)
  
  beta_logical <- (beta==0) + (betahat!=0) #false neg
  FNRvalue <- sum(beta_logical==2)/sum(beta==0)
  return(FNRvalue)
}

# Matthews correlation coefficient
mcc_fct <- function(beta, betahat){
  preds = (betahat!=0) #non-zero = TRUE, zero = FALSE
  actual = (beta!=0)
  
  coef = mcc(preds, actual)
  return(coef)
}

# Function to check for stability
stable_check <- function(B){
  #library(plotrix)
  eigen = eigen(B)$values
  #plot(1, type="n", xlab="", ylab="", xlim=c(-1, 1), ylim=c(-1, 1), asp=1)
  #draw.circle(0,0,1)
  #lines(x = eigen, y = rep(0,length(eigen)), type = "p")
  
  check = max(abs(eigen)) #sum(abs(eigen) < 1)
  if(check < 1){
    is_stable = TRUE
  } else{
    is_stable = FALSE}
  
  return(list(eigen = eigen, is_stable = is_stable))
}

simMFVAR <- function(n, BVAR, mu, sigma, p, burn_in, standardize = TRUE){
  innov = mvrnorm(n+burn_in+p, mu = mu, Sigma = sigma)
  start_value = rbind(innov[1:p,])
  innov = innov[-c(1:p),]
  sim = VAR.sim(BVAR, n = n+burn_in, lag = p, include = c("none"), starting = start_value, innov = innov, returnStarting = T) 
  if(standardize==TRUE){
    sim_data = scale(sim[-c(1:burn_in),])
  }else{
    sim_data = sim[-c(1:burn_in),]
  }
  return(sim_data)
}


  
