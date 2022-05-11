# Construction of matrix that gives column that element belongs to
columns_function <- function(k_q = 1, k_m = 1, k_w = 0, k_d = 0, p = 1){
  
  # k_q, k_m, k_w, k_d: number of quarterly, monthly, weekly and daily series
  # p: number of lags
  m1 = 3
  m2 = 12
  m3 = 60
  
  K = k_q+k_m*m1+k_w*m2+k_d*m3
  num_groups = (k_q+k_m+k_w+k_d)^2
  multi = K^2*c(0:(p-1))
  
  # Case 1: quarterly & monthly 
  if(k_q != 0 & k_m != 0 & k_w == 0 & k_d == 0){
    max_size_group = m1^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1:p] = (K*(l-1)+j) + multi
        i = i+1
      }
    }
    
    # pi
    for(l in 1:k_q){
      for(j in 1:k_m){
        z = c((K*(l-1)+k_q+2*j+(j-2)):(K*(l-1)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # theta
    for(l in 1:k_m){
      for(j in 1:k_q){
        z = c(K*(k_q+(l-1)*3)+ j, K*(k_q+(l-1)*3+1)+ j, K*(k_q+(l-1)*3+2)+ j)
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # rho
    for(l in 1:k_m){
      for(j in 1:k_m){
        z = c((K*(k_q+(l-1)*3)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3)+k_q+2*j+(j-2)+2),
              (K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)+2),
              (K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    rownames(columns_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2))
  }
  
  # Case 2: quarterly & weekly 
  if(k_q != 0 & k_m == 0 & k_w != 0 & k_d == 0){
    max_size_group = m2^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1:p] = K*(l-1)+j + multi
        i = i+1
      }
    }
    
    #alpha
    for(l in 1:k_q){
      for(j in 1:k_w){
        z = c((K*(l-1)+k_q+m1*k_m+m2*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # mu
    for(l in 1:k_w){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+(l-1)*12+ c(0:11)) +j 
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # gamma
    for(l in 1:k_w){
      for(j in 1:k_w){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+12) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1 
      }
    }
    rownames(columns_matrix) = c(rep("phi",k_q^2), rep("alpha", k_q*k_w), rep("mu", k_q*k_w), rep("gamma", k_w^2))
  }
  
  # Case 3: quarterly & daily
  if(k_q != 0 & k_m == 0 & k_w == 0 & k_d != 0){
    max_size_group = m3^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1:p] = K*(l-1)+j + multi
        i = i+1
      }
    }
    
    # epsilon
    for(l in 1:k_q){
      for(j in 1:k_d){
        z = c((K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # eta
    for(l in 1:k_d){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ j
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # kappa
    b = c()
    for(l in 1:k_d){
      for(j in 1:k_d){                                                                
        for(m in 1:m3){
          a = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*k_w+m3*(j-1)+m
          b = c(b,a)
        }                                                                   
        z = sort(b)
        columns_matrix[i,1:(m3^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
        b = c()
      }
    }
    
    rownames(columns_matrix) = c(rep("phi",k_q^2), rep("epsilon", k_q*k_d), rep("eta", k_d*k_q), rep("kappa",k_d^2))
    
  }
  
  # Case 4: quarterly & monthly & weekly
  if(k_q != 0 & k_m != 0 & k_w != 0 & k_d == 0){
    max_size_group = m2^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1:p] = K*(l-1)+j + multi
        i = i+1
      }
    }
    
    # pi
    for(l in 1:k_q){
      for(j in 1:k_m){
        z = c((K*(l-1)+k_q+2*j+(j-2)):(K*(l-1)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # theta
    for(l in 1:k_m){
      for(j in 1:k_q){
        z = c(K*(k_q+(l-1)*3)+ j, K*(k_q+(l-1)*3+1)+ j, K*(k_q+(l-1)*3+2)+ j)
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # rho
    for(l in 1:k_m){
      for(j in 1:k_m){
        z = c((K*(k_q+(l-1)*3)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3)+k_q+2*j+(j-2)+2),
              (K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)+2),
              (K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    #alpha
    for(l in 1:k_q){
      for(j in 1:k_w){
        z = c((K*(l-1)+k_q+m1*k_m+m2*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # beta
    for(l in 1:k_m){
      for(j in 1:k_w){
        z = c((K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*(j-1)+1):(K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*(j-1)+1):(K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+(l-1)*3+2)+k_q++m1*k_m+m2*(j-1)+1):(K*(k_q+(l-1)*3+2)+k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m1*m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # mu
    for(l in 1:k_w){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+(l-1)*12+ c(0:11)) +j 
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1
      }
    }
    
    # lambda
    for(l in 1:k_w){
      for(j in 1:k_m){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+2*j+(j-2)+2) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+2*j+(j-2)+2), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+2*j+(j-2)+2), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1*m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1  
      }
    }
    
    # gamma
    for(l in 1:k_w){
      for(j in 1:k_w){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+12) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1 
      }
    }
    
    rownames(columns_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),
                                 rep("alpha", k_q*k_w), rep("beta", k_m*k_w), rep("mu", k_q*k_w), rep("lambda", k_m*k_w), rep("gamma", k_w^2))
  }
  
  # Case 5: quarterly & monthly & daily
  if(k_q != 0 & k_m != 0 & k_w == 0 & k_d != 0){
    max_size_group = m3^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1:p] = K*(l-1)+j + multi
        i = i+1
      }
    }
    
    # pi
    for(l in 1:k_q){
      for(j in 1:k_m){
        z = c((K*(l-1)+k_q+2*j+(j-2)):(K*(l-1)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # theta
    for(l in 1:k_m){
      for(j in 1:k_q){
        z = c(K*(k_q+(l-1)*3)+ j, K*(k_q+(l-1)*3+1)+ j, K*(k_q+(l-1)*3+2)+ j)
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # rho
    for(l in 1:k_m){
      for(j in 1:k_m){
        z = c((K*(k_q+(l-1)*3)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3)+k_q+2*j+(j-2)+2),
               (K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)+2),
               (K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # epsilon
    for(l in 1:k_q){
      for(j in 1:k_d){
        z = c((K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # delta
    for(l in 1:k_m){
      for(j in 1:k_d){
        z =  c((K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
               (K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
               (K*(k_q+(l-1)*3+2)+k_q++m1*k_m+m2*k_w+m3*(j-1)+1):(K*(k_q+(l-1)*3+2)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m1*m3*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1
      }
    }
    
    # eta
    for(l in 1:k_d){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ j
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # tau
    for(l in 1:k_d){
      for(j in 1:k_m){
        help_1 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60 + c(0:59))+ k_q+ 2*j+(j-2)
        help_2 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+ c(0:59))+ k_q+2*j+(j-2)+1
        help_3 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+ c(0:59))+ k_q+2*j+(j-2)+2
        
        z = sort(c(help_1,help_2,help_3))
        columns_matrix[i,1:(m1*m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # kappa
    b = c()
    for(l in 1:k_d){
      for(j in 1:k_d){                                                                
        for(m in 1:m3){
          a = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*k_w+m3*(j-1)+m
          b = c(b,a)
        }                                                                   
        
        z = sort(b)
        columns_matrix[i,1:(m3^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
        b = c()
      }
    }
    
    rownames(columns_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),    
                                 rep("epsilon", k_q*k_d), rep("delta", k_m*k_d),rep("eta", k_d*k_q),
                                 rep("tau", k_m*k_d), rep("kappa",k_d^2))
    
  }
  
  # Case 6: quarterly & weekly & daily 
  if(k_q != 0 & k_m == 0 & k_w != 0 & k_d != 0){
    max_size_group = m3^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1*p] = K*(l-1)+j + multi
        i = i+1
      }
    }
    
    #alpha
    for(l in 1:k_q){
      for(j in 1:k_w){
        z = c((K*(l-1)+k_q+m1*k_m+m2*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # mu
    for(l in 1:k_w){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+(l-1)*12+ c(0:11)) +j 
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1
      }
    }
    
    # gamma
    for(l in 1:k_w){
      for(j in 1:k_w){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+12) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1 
      }
    }
    
    # epsilon
    for(l in 1:k_q){
      for(j in 1:k_d){
        z = c((K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # xi
    for(l in 1:k_w){
      for(j in 1:k_d){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m2*m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1 
      }
    }
    
    # eta
    for(l in 1:k_d){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ j
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # nu
    for(l in 1:k_d){
      for(j in 1:k_w){
        help_1 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +1
        help_2 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +2
        help_3 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +3
        help_4 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +4
        help_5 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +5
        help_6 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +6
        help_7 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +7
        help_8 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +8
        help_9 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +9
        help_10 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +10
        help_11 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +11
        help_12 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +12
        
        z = sort(c(help_1,help_2,help_3, help_4, help_5, help_6, help_7, help_8, help_9, help_10, help_11, help_12))
        columns_matrix[i,1:(m2*m3*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1 
      }
    }
    
    # kappa
    b = c()
    for(l in 1:k_d){
      for(j in 1:k_d){                                                                
        for(m in 1:m3){
          a = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*k_w+m3*(j-1)+m
          b = c(b,a)
        }                                                                   
        
        z = sort(b)
        columns_matrix[i,1:(m3^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
        b = c()
      }
    }
    
    rownames(columns_matrix) = c(rep("phi",k_q^2),
                                 rep("alpha", k_q*k_w), rep("mu", k_q*k_w), rep("gamma", k_w^2), 
                                 rep("epsilon", k_q*k_d), rep("xi", k_w*k_d), rep("eta", k_d*k_q),
                                 rep("nu", k_w*k_d), rep("kappa",k_d^2))
    
  }
  
  # Case 7: quarterly & monthly & weekly & daily 
  if(k_q != 0 & k_m != 0 & k_w != 0 & k_d != 0){
    max_size_group = m3^2*p
    
    columns_matrix = matrix(NA, num_groups, max_size_group)
    
    # phi
    i = 1
    
    for(l in 1:k_q){
      for(j in 1:k_q){
        columns_matrix[i,1:p] = K*(l-1)+j + multi
        i = i+1
      }
    }
    
    # pi
    for(l in 1:k_q){
      for(j in 1:k_m){
        z = c((K*(l-1)+k_q+2*j+(j-2)):(K*(l-1)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # theta
    for(l in 1:k_m){
      for(j in 1:k_q){
        z = c(K*(k_q+(l-1)*3)+ j, K*(k_q+(l-1)*3+1)+ j, K*(k_q+(l-1)*3+2)+ j)
        columns_matrix[i,1:(m1*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # rho
    for(l in 1:k_m){
      for(j in 1:k_m){
        z = c((K*(k_q+(l-1)*3)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3)+k_q+2*j+(j-2)+2),
              (K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+1)+k_q+2*j+(j-2)+2),
              (K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)):(K*(k_q+(l-1)*3+2)+k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    #alpha
    for(l in 1:k_q){
      for(j in 1:k_w){
        z = c((K*(l-1)+k_q+m1*k_m+m2*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # beta
    for(l in 1:k_m){
      for(j in 1:k_w){
        z = c((K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*(j-1)+1):(K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*(j-1)+1):(K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+(l-1)*3+2)+k_q++m1*k_m+m2*(j-1)+1):(K*(k_q+(l-1)*3+2)+k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m1*m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # mu
    for(l in 1:k_w){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+(l-1)*12+ c(0:11)) +j 
        columns_matrix[i,1:(m2*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1
      }
    }
    
    # lambda
    for(l in 1:k_w){
      for(j in 1:k_m){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+2*j+(j-2)+2) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+2*j+(j-2)+2), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+2*j+(j-2)+2), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+2*j+(j-2)+2),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+2*j+(j-2)) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+2*j+(j-2)+2))
        columns_matrix[i,1:(m1*m2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1  
      }
    }
    
    # gamma
    for(l in 1:k_w){
      for(j in 1:k_w){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*(j-1)+12) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*(j-1)+12), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*(j-1)+12),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*(j-1)+12))
        columns_matrix[i,1:(m2^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1 
      }
    }
    
    # epsilon
    for(l in 1:k_q){
      for(j in 1:k_d){
        z = c((K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(l-1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # delta
    for(l in 1:k_m){
      for(j in 1:k_d){
        z = c((K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(k_q+(l-1)*3)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+1):(K*(k_q+(l-1)*3+1)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+(l-1)*3+2)+k_q++m1*k_m+m2*k_w+m3*(j-1)+1):(K*(k_q+(l-1)*3+2)+k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m1*m3*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1
      }
    }
    
    # xi
    for(l in 1:k_w){
      for(j in 1:k_d){
        z = c((K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60) , 
              (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+1)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60), 
              (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+2)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+3)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60), 
              (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+4)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+5)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+6)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+7)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+8)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+9)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+10)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60),
              (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+1) : (K*(k_q+m1*k_m+(l-1)*12+11)+ k_q+m1*k_m+m2*k_w+m3*(j-1)+60))
        columns_matrix[i,1:(m2*m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1 
      }
    }
    
    # eta
    for(l in 1:k_d){
      for(j in 1:k_q){
        z = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ j
        columns_matrix[i,1:(m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # tau
    for(l in 1:k_d){
      for(j in 1:k_m){
        help_1 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60 + c(0:59))+ k_q+ 2*j+(j-2)
        help_2 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+ c(0:59))+ k_q+2*j+(j-2)+1
        help_3 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+ c(0:59))+ k_q+2*j+(j-2)+2
        
        z = sort(c(help_1,help_2,help_3))
        columns_matrix[i,1:(m1*m3*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
      }
    }
    
    # nu
    for(l in 1:k_d){
      for(j in 1:k_w){
        help_1 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +1
        help_2 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +2
        help_3 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +3
        help_4 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +4
        help_5 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +5
        help_6 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +6
        help_7 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +7
        help_8 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +8
        help_9 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +9
        help_10 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +10
        help_11 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +11
        help_12 = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*(j-1) +12
        
        z = sort(c(help_1,help_2,help_3, help_4, help_5, help_6, help_7, help_8, help_9, help_10, help_11, help_12))
        columns_matrix[i,1:(m2*m3*p)] = c(outer(z, multi, FUN = "+"))
        
        i = i+1 
      }
    }
    
    # kappa
    b = c()
    for(l in 1:k_d){
      for(j in 1:k_d){                                                                
        for(m in 1:m3){
          a = K*(k_q+m1*k_m+m2*k_w+(l-1)*60+c(0:59))+ k_q+m1*k_m+m2*k_w+m3*(j-1)+m
          b = c(b,a)
        }                                                                   
        
        z = sort(b)
        columns_matrix[i,1:(m3^2*p)] = c(outer(z, multi, FUN = "+"))
        i = i+1
        b = c()
      }
    }
    
    rownames(columns_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),
                                 rep("alpha", k_q*k_w), rep("beta", k_m*k_w), rep("mu", k_q*k_w), rep("lambda", k_m*k_w), rep("gamma", k_w^2), 
                                 rep("epsilon", k_q*k_d), rep("delta", k_m*k_d), rep("xi", k_w*k_d), rep("eta", k_d*k_q),
                                 rep("tau", k_m*k_d), rep("nu", k_w*k_d), rep("kappa",k_d^2))
    
  }
  
  return(columns_matrix)
}


# Construction of B matrix 
construct_B_matrix <- function(k_q = 1, k_m = 1, k_w = 0, k_d = 0, p = 1){
  
  # k_q, k_m, k_w, k_d: number of quarterly, monthly, weekly and daily series
  m1 = 3
  m2 = 12
  m3 = 60
  
  K = k_q+k_m*m1+k_w*m2+k_d*m3
  num_groups = (k_q+k_m+k_w+k_d)^2
  
  # Case 1: quarterly & monthly 
  if(k_q != 0 & k_m != 0 & k_w == 0 & k_d == 0){
    max_size_group = m1^2*p
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # pi
    for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # theta
    for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # rho
    for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1^2*p)] = 0 
    }
    rownames(B_matrix_1) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2))
  }
  
  # Case 2: quarterly & weekly 
  if(k_q != 0 & k_m == 0 & k_w != 0 & k_d == 0){
    max_size_group = m2^2*p
    
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # alpha 
    for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # mu
    for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # gamma
    for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2^2*p)] = 0 
    }
    rownames(B_matrix_1) = c(rep("phi",k_q^2), rep("alpha", k_q*k_w), rep("mu", k_q*k_w), rep("gamma", k_w^2))
  }
  
  # Case 3: quarterly & daily
  if(k_q != 0 & k_m == 0 & k_w == 0 & k_d != 0){
    max_size_group = m3^2*p
    
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2*p)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # epsilon
    for(i in (k_q^2+1):(k_q*k_d+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # eta
    for(i in (k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # kappa
    for(i in (2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)){
      B_matrix_1[i,1:(m3^2*p)] = 0 
    }
    
    rownames(B_matrix_1) = c(rep("phi",k_q^2), rep("epsilon", k_q*k_d), rep("eta", k_d*k_q), rep("kappa",k_d^2))
    
  }
  
  # Case 4: quarterly & monthly & weekly
  if(k_q != 0 & k_m != 0 & k_w != 0 & k_d == 0){
    max_size_group = m2^2*p
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # pi
    for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # theta
    for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # rho
    for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1^2*p)] = 0 
    }
    
    # alpha 
    for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # beta
    for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m2*p)] = 0 
    }
    
    # mu
    for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # lambda
    for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m2*p)] = 0 
    }
    
    # gamma
    for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2^2*p)] = 0 
    }
    rownames(B_matrix_1) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),
                             rep("alpha", k_q*k_w), rep("beta", k_m*k_w), rep("mu", k_q*k_w), rep("lambda", k_m*k_w), rep("gamma", k_w^2))
  }
  
  # Case 5: quarterly & monthly & daily
  if(k_q != 0 & k_m != 0 & k_w == 0 & k_d != 0){
    max_size_group = m3^2*p
    
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # pi
    for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # theta
    for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # rho
    for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1^2*p)] = 0 
    }
    
    # epsilon
    for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # delta
    for(i in (k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m3*p)] = 0 
    }
    
    # eta
    for(i in (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # tau 
    for(i in (2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m3*p)] = 0 
    }
    
    # kappa
    for(i in (2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m3^2*p)] = 0 
    }
    
    rownames(B_matrix_1) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),    
                             rep("epsilon", k_q*k_d), rep("delta", k_m*k_d),rep("eta", k_d*k_q),
                             rep("tau", k_m*k_d), rep("kappa",k_d^2))
    
  }
  
  # Case 6: quarterly & weekly & daily 
  if(k_q != 0 & k_m == 0 & k_w != 0 & k_d != 0){
    max_size_group = m3^2*p
    
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # alpha 
    for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # mu
    for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # gamma
    for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2^2*p)] = 0 
    }
    
    # epsilon
    for(i in (k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # xi
    for(i in (k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2*m3*p)] = 0 
    }
    
    # eta
    for(i in (k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # nu
    for(i in (2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m2*m3*p)] = 0 
    }
    
    # kappa
    for(i in (2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
      B_matrix_1[i,1:(m3^2*p)] = 0 
    }
    
    rownames(B_matrix_1) = c(rep("phi",k_q^2),
                             rep("alpha", k_q*k_w), rep("mu", k_q*k_w), rep("gamma", k_w^2), 
                             rep("epsilon", k_q*k_d), rep("xi", k_w*k_d), rep("eta", k_d*k_q),
                             rep("nu", k_w*k_d), rep("kappa",k_d^2))
    
  }
  
  # Case 7: quarterly & monthly & weekly & daily 
  if(k_q != 0 & k_m != 0 & k_w != 0 & k_d != 0){
    max_size_group = m3^2*p
    
    B_matrix_1 = matrix(NaN, nrow = num_groups, ncol = max_size_group)
    
    # phi
    for(i in 1:(k_q^2)){
      B_matrix_1[i,1:p] = 0 
    }
    
    # pi
    for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # theta
    for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*p)] = 0 
    }
    
    # rho
    for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1^2*p)] = 0 
    }
    
    # alpha 
    for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # beta
    for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m2*p)] = 0 
    }
    
    # mu
    for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2*p)] = 0 
    }
    
    # lambda
    for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m2*p)] = 0 
    }
    
    # gamma
    for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2^2*p)] = 0 
    }
    
    # epsilon
    for(i in (k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # delta
    for(i in (k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m3*p)] = 0 
    }
    
    # xi
    for(i in (k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2*m3*p)] = 0 
    }
    
    # eta
    for(i in (k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m3*p)] = 0 
    }
    
    # tau 
    for(i in (2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m1*m3*p)] = 0 
    }
    
    # nu
    for(i in (2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m2*m3*p)] = 0 
    }
    
    # kappa
    for(i in (2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
      B_matrix_1[i,1:(m3^2*p)] = 0 
    }
    
    rownames(B_matrix_1) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),
                             rep("alpha", k_q*k_w), rep("beta", k_m*k_w), rep("mu", k_q*k_w), rep("lambda", k_m*k_w), rep("gamma", k_w^2), 
                             rep("epsilon", k_q*k_d), rep("delta", k_m*k_d), rep("xi", k_w*k_d), rep("eta", k_d*k_q),
                             rep("tau", k_m*k_d), rep("nu", k_w*k_d), rep("kappa",k_d^2))
    
  }
  
  return(B_matrix_1)
}


# Hierarchy Penalty function
# 1 lowest penalty ... x highest penalty 
hierarchy_indices = function(penalty, row = 1, column = 1, p = 1){
  # penalty: form of penalty on nested group (horizontal = 2, vertical/matrix form = 3)
  # row: length of regressor group (m0 = 1, m1 = 3, m2 = 12 or m3 = 60)
  # column: length of dependent variable group (m0 = 1, m1 = 3, m2 = 12 or m3 = 60)
  
  # horizontal penalty (left to right) 
  # concerns groups:  pi, alpha, epsilon
  if(penalty == 2){
    indices = c(1:(column*p))
  }
  
  # vertical penalty if column = 1 (concerns groups: phi, theta, mu, eta)
  # symmetric matrix penalty if column = row (rho, gamma, kappa)
  # matrix penalty otherwise (concerns groups: lambda, tau, nu, beta, delta, xi)
  if(penalty == 3){
    mat <- matrix(NA, row, column)
    for(i in 1:row){
      mat[row-(i-1), ] <- seq(from = i, by = 1, length = column)
    }
    indices_help = c(t(mat))
    indices = c(outer(indices_help, max(indices_help)*c(0:(p-1)), FUN = "+" ))
  }
  
  # add another option regarding the matrix penalty with option 2 (care about different frequencies)
  # would affect beta, lambda, delta, xi, tau, nu
  # penalty = 4 (this would have to be changed in penalty_function when hierarchy_indices is being called)
  return(indices)
}


# Construction of general penalty indicator (L1 = 1, HIER = 2)
# and hierarchy penalty matrix function

penalty_function = function(k_q = 1, k_m = 1, k_w = 0, k_d = 0, p = 1,
                            penalty_H_on_L = "L1", penalty_L_on_H = "L1", penalty_own_on_own = "L1"){
  
  # k_q, k_m, k_w, k_d: number of quarterly, monthly, weekly and daily series
  # penalty_H_on_L: Penalty High Frequency on Low Frequency ("L1" or "HIER")
  # penalty_L_on_H: Penalty Low Frequency on High Frequency ("L1" or "HIER")
  # penalty_own_on_own: Penalty High/Low Frequency on High/Low Frequency ("L1" or "HIER")
  m1 = 3
  m2 = 12
  m3 = 60
  
  K = k_q+k_m*m1+k_w*m2+k_d*m3
  num_groups = (k_q+k_m+k_w+k_d)^2
  
  # Case 1: quarterly & monthly 
  if(k_q != 0 & k_m != 0 & k_w == 0 & k_d == 0){
    max_size_group = m1^2*p
    
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
    }
    
  }
  
  # Case 2: quarterly & weekly 
  if(k_q != 0 & k_m == 0 & k_w != 0 & k_d == 0){
    max_size_group = m2^2*p
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2), rep("alpha", k_q*k_w), rep("mu", k_q*k_w), rep("gamma", k_w^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[ (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      
      for(i in  (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[ (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in  (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[ (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in  (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
    }
  }
  
  # Case 3: quarterly & daily
  if(k_q != 0 & k_m == 0 & k_w == 0 & k_d != 0){
    max_size_group = m3^2*p
    
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2), rep("epsilon", k_q*k_d), rep("eta", k_d*k_q), rep("kappa",k_d^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)] = 3  #eta
      
      for(i in (k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_q*k_d+k_q^2)] = 2 #epsilon
      
      for(i in (k_q^2+1):(k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_q*k_d+k_q^2)] = 2 #epsilon
      penalty_indicator[(2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_q*k_d+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)] = 3  #eta
      
      for(i in (k_q^2+1):(k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_q*k_d+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, , p = p) #phi
      }
      for(i in (k_q^2+1):(k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_q^2+1):(2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_q^2+1):(k_d^2+2*k_q*k_d+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    
  }
  
  # Case 4: quarterly & monthly & weekly
  if(k_q != 0 & k_m != 0 & k_w != 0 & k_d == 0){
    max_size_group = m2^2*p
    
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),
                                 rep("alpha", k_q*k_w), rep("beta", k_m*k_w), rep("mu", k_q*k_w), rep("lambda", k_m*k_w), rep("gamma", k_w^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
    }
    
  }
  
  # Case 5: quarterly & monthly & daily
  if(k_q != 0 & k_m != 0 & k_w == 0 & k_d != 0){
    max_size_group = m3^2*p
    
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),    
                                 rep("epsilon", k_q*k_d), rep("delta", k_m*k_d),rep("eta", k_d*k_q),
                                 rep("tau", k_m*k_d), rep("kappa",k_d^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      penalty_indicator[(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
      for(i in (2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      penalty_indicator[(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1): (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (k_m*k_d+k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_m*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
      for(i in (2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_m*k_d+2*k_q*k_d+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
  }
  
  # Case 6: quarterly & weekly & daily 
  if(k_q != 0 & k_m == 0 & k_w != 0 & k_d != 0){
    max_size_group = m3^2*p
    
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2),
                                 rep("alpha", k_q*k_w), rep("mu", k_q*k_w), rep("gamma", k_w^2), 
                                 rep("epsilon", k_q*k_d), rep("xi", k_w*k_d), rep("eta", k_d*k_q),
                                 rep("nu", k_w*k_d), rep("kappa",k_d^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      penalty_indicator[(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #nu
      
      
      for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[(k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #xi
      
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups)
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      penalty_indicator[(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #nu
      penalty_indicator[(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
      for(i in (2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      penalty_indicator[(k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #xi
      penalty_indicator[(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      for(i in (2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      penalty_indicator[(k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #xi
      penalty_indicator[(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #nu
      
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      for(i in (k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_q*k_w+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)] = 3 #gamma
      penalty_indicator[(k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #xi
      penalty_indicator[(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3  #nu
      penalty_indicator[(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_q^2+1):(2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_q^2+1):(k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (k_w^2+2*k_q*k_w+k_q^2+1):(k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      for(i in (k_w*k_d+k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_q*k_d+k_w*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
      for(i in (2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_q*k_d+k_w^2+2*k_q*k_w+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
  }
  
  # Case 7: quarterly & monthly & weekly & daily 
  if(k_q != 0 & k_m != 0 & k_w != 0 & k_d != 0){
    max_size_group = m3^2*p
    
    penalty_matrix = matrix(NA, nrow = num_groups, ncol = max_size_group)
    rownames(penalty_matrix) = c(rep("phi",k_q^2), rep("pi", k_m*k_q), rep("theta", k_m*k_q), rep("rho", k_m^2),
                                 rep("alpha", k_q*k_w), rep("beta", k_m*k_w), rep("mu", k_q*k_w), rep("lambda", k_m*k_w), rep("gamma", k_w^2), 
                                 rep("epsilon", k_q*k_d), rep("delta", k_m*k_d), rep("xi", k_w*k_d), rep("eta", k_d*k_q),
                                 rep("tau", k_m*k_d), rep("nu", k_w*k_d), rep("kappa",k_d^2))
    
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      penalty_indicator[(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      penalty_indicator[(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      penalty_indicator[(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #nu
      
      
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
      for(i in (k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
      for(i in (2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #xi
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      
    }
    if(penalty_H_on_L == "L1" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      penalty_indicator[(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      penalty_indicator[(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #nu
      penalty_indicator[(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
      for(i in (2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
      for(i in (2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "L1" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      penalty_indicator[(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #xi
      penalty_indicator[(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      for(i in (2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "L1"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      penalty_indicator[(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #xi
      penalty_indicator[(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #tau
      penalty_indicator[(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #nu
      
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
      for(i in (k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      for(i in (k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) 
      }
      for(i in (2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
      
    }
    if(penalty_H_on_L == "HIER" & penalty_L_on_H == "HIER" & penalty_own_on_own == "HIER"){
      penalty_indicator = rep(1,num_groups) 
      penalty_indicator[1:(k_q^2)] = 3 #phi
      penalty_indicator[(k_q^2+1):(k_m*k_q+k_q^2)] = 2   #pi
      penalty_indicator[(k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)] = 3   #theta
      penalty_indicator[(2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)] = 3 #rho
      penalty_indicator[(k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #alpha
      penalty_indicator[(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #beta
      penalty_indicator[(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #mu
      penalty_indicator[(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #lambda
      penalty_indicator[(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #gamma
      penalty_indicator[(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 2 #epsilon
      penalty_indicator[(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #delta
      penalty_indicator[(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #xi
      penalty_indicator[(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #eta
      penalty_indicator[(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  
      penalty_indicator[(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3  #nu
      penalty_indicator[(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)] = 3 #kappa
      
      for(i in 1:(k_q^2)){
        penalty_matrix[i, 1:p] = hierarchy_indices(penalty = 3, row = 1, column = 1, p = p) #phi
      }
      for(i in (k_q^2+1):(k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 2, row = 1, column = 3, p = p) #pi
      }
      for(i in (k_m*k_q+k_q^2+1):(2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*p)] = hierarchy_indices(penalty = 3, row = 3, column = 1, p = p) #theta
      }
      for(i in (2*k_m*k_q+k_q^2+1):(k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1^2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 3, p = p) #rho
      }
      for(i in (k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 2, row = 1, column = 12, p = p) #alpha
      }
      for(i in (k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 3, column = 12, p = p) #beta
      }
      for(i in (k_m*k_w+k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 1, p = p) #mu
      }
      for(i in (2*k_q*k_w+k_m*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 3, p = p) #lambda
      }
      for(i in (2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2^2*p)] = hierarchy_indices(penalty = 3, row = 12, column = 12, p = p) #gamma
      }
      for(i in (k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 2, row = 1, column = 60, p = p) #epsilon
      }
      for(i in (k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 3, column = 60, p = p) #delta
      }
      for(i in (k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 12, column = 60, p = p) #xi
      }
      for(i in (k_w*k_d+k_m*k_d+k_q*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 1, p = p) #eta
      }
      for(i in (2*k_d*k_q+k_w*k_d+k_m*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m1*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 3, p = p) #tau
      }
      for(i in (2*k_m*k_d+2*k_d*k_q+k_w*k_d+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m2*m3*p)] = hierarchy_indices(penalty = 3, row = 60, column = 12, p = p) #nu
      }
      for(i in (2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2+1):(k_d^2+2*k_w*k_d+2*k_m*k_d+2*k_d*k_q+k_w^2+2*k_m*k_w+2*k_q*k_w+k_m^2+2*k_m*k_q+k_q^2)){
        penalty_matrix[i, 1:(m3^2*p)] = hierarchy_indices(penalty = 3, row = 60, column = 60, p = p) #kappa
      }
    }
  }
  
  return(list("penalty_indicator" = penalty_indicator, "penalty_matrix" = penalty_matrix))
}


# Construct all objects necessary as input for sparseMFVAR
inputs_sparseMFVAR <- function(k_q, k_m, k_w, k_d, p,
                               penalty_H_on_L = "HIER", penalty_L_on_H = "HIER", penalty_own_on_own = "HIER"){
  
  # Inputs:
  # k_q, k_m, k_w, k_d: number of quarterly, monthly, weekly and daily series
  # penalty_H_on_L: Penalty High Frequency on Low Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # penalty_L_on_H: Penalty Low Frequency on High Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  # penalty_own_on_own: Penalty High/Low Frequency on High/Low Frequency ("L1" (lasso penalty) or "HIER" (hierarchical penalty))
  
  # Function: Returns object which is necessary input for sparseMFVAR function
  
  # Output: 
  # B_matrix_1, B_matrix_2: matrices in each row we save parameter estiamtes of a certain parameter group
  # index: indices of all non NA values in B_matrix_1, B_matrix_2
  # penalty_matrix: Hierarchy penalty matrix (if parameter group has "HIER" as penalty)
  # columns_matrix: Matrix that gives the column that a parameter belongs to in the lagged coefficient matrix (Xkron = kronecker(I, X)) 
  # rows_L1, rows_HIER:  Rows (in B_matrix) which have L1 penalty or HIER penalty
  # indices_B_L1, indices_B_HIER: Indices of the specific values in rows that have either L1 or HIER penalty
  # column_L1, column_HIER: Columns that parameter belongs to (split into either L1 or HIER penalty)
  # K: number of time series
  
  
  m1 = 3
  m2 = 12
  m3 = 60
  K = k_q+k_m*m1+k_w*m2+k_d*m3 # number of time series
  
  # Construction of B matrix 
  B_matrix_1 = construct_B_matrix(k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p)
  B_matrix_2 = B_matrix_1
  index = which(!is.nan(B_matrix_2)) #indices of all non NA values
  
  # Construction of general penalty indicator (L1 = 1, HIER = 2 or 3)
  # and hierarchy penalty matrix (if penalty indicator = 2 or 3)
  penalty = penalty_function(k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p,
                             penalty_H_on_L = penalty_H_on_L, penalty_L_on_H = penalty_L_on_H, penalty_own_on_own = penalty_own_on_own)
  penalty_indicator = penalty$penalty_indicator
  penalty_matrix = penalty$penalty_matrix
  
  # Construction of matrix that gives the column that a parameter belongs to in the Xkron matrix 
  columns_matrix = columns_function(k_q = k_q, k_m = k_m, k_w = k_w, k_d = k_d, p = p)
  
  # Rows (in B_matrix) which have L1 penalty or HIER penalty
  rows_L1 = which(penalty_indicator == 1) 
  rows_HIER = which(penalty_indicator != 1)
  
  # Indices of the specific values in rows that have either L1 or HIER penalty
  indices_B_L1 <- which(!is.nan(B_matrix_2[rows_L1, ])) 
  indices_B_HIER <- which(!is.nan(B_matrix_2[rows_HIER, ])) 
  
  # Columns that parameter belongs to (split into either L1 or HIER penalty)
  column_L1 <- columns_matrix[rows_L1,][indices_B_L1]
  column_HIER <- columns_matrix[rows_HIER,][indices_B_HIER]
  
  
  
  out <- list(
    #"y" = y, "Xkron" = Xkron, "y_Xkron" = y_Xkron, "XtranspX_kron" = XtranspX_kron,
    "B_matrix_1" = B_matrix_1, "B_matrix_2" = B_matrix_2, "index" = index,
    "penalty_matrix" = penalty_matrix,
    "columns_matrix" = columns_matrix, 
    "rows_L1" = rows_L1, "rows_HIER" = rows_HIER,
    "indices_B_L1" = indices_B_L1, "indices_B_HIER" = indices_B_HIER, 
    "column_L1" = column_L1, "column_HIER" = column_HIER,
    "K"= K) 
  #"s" = s)
  
  return(out)
}

# Data preprocessing
make_data_matrices <- function(data_quarterly, data_monthly = NULL, data_weekly = NULL, data_daily = NULL){
  # Input
  # data_quarterly: Quarterly data in format cbind(LF1, LF2,..., LFk_q), where k_q = # of quarterly series
  # data_monthly: Monthly data in format cbind(HF1, HF2,..., HFk_m), where k_m = # of monthly series
  # data_weekly: Weekly data in format cbind(HHF1, HHF2,..., HHFk_w), where k_w = # of weekly series
  # data_daily: daily data in format cbind(HHHF1, HHHF2,..., HHHFk_d), where k_d = # of daily series
  
  # Function: Preprocesses data of different frequencies such that one gets a suitable matrix of time series
  
  # Output:
  # Y_VAR: TxK matrix of time series with columns
  #        Q_var1, Q_var2,..., Q_vark_q,
  #        M3_var1(over t), M2_var1, M1var1,..., M3_vark_m(over t), M2_vark_m,...,M1vark_m,
  #        W12_var1(over t), W11_var1,...,W1var1,..., W12_vark_w(over t), W11_vark_w,...,W1vark_w,
  #        D60_var1(over t), D59_var1,...,D58var1,..., D1_vark_1(over t), D60_vark_d(over t), D59_var_k_d,..., D1_var_k_d(over t)
  # series_names: names of time series
  # K: total number of time series
  # k_q = # of quarterly series
  # k_m = # of monthly series
  # k_w = # of weekly series
  # k_d = # of daily series
  
  
  k_q <- ncol(data_quarterly) #k_q = # of quarterly series
  m1 <- 3 #relation quarterly to monthly
  m2 <- 12 #relation quarterly to weekly
  m3 <- 60 #relation quarterly to daily 
  k_m <- 0 #k_m = # of monthly series
  k_w <- 0 #k_w = # of weekly series
  k_d <- 0 #k_d = # of daily series
  n <- nrow(data_quarterly)
  
  if(!is.null(colnames(data_quarterly))){ # time series names
    seriesQ_names <- colnames(data_quarterly)
  }else{
    seriesQ_names <- NULL
  }
  
  Y_Q = matrix(data_quarterly, n, k_q) ##Q_var1, Q_var2,..., Q_vark_q
  colnames(Y_Q) = seriesQ_names
  Y_VAR = Y_Q
  
  
  if(!is.null(data_monthly)){
    k_m <- ncol(data_monthly)
    if(!is.null(colnames(data_monthly))){ # time series names
      seriesM_names <- colnames(data_monthly)
    }else{
      seriesM_names <- NULL
    }
    
    Y_M = matrix(NA, n, m1*k_m)
    colnames(Y_M) = paste(rep(seriesM_names, each = m1), rep(c("M3","M2","M1"), k_m))
    
    j = 1
    for(i in 1:n){
      Y_M[i,] = data_monthly[(j*m1):((j-1)*m1+1),] #M3_var1(over t), M2_var1, M1var1,..., M3_vark_m(over t), M2_vark_m,...,M1vark_m
      j = j+1
    }
    Y_VAR = cbind(Y_VAR, Y_M)
  }
  
  if(!is.null(data_weekly)){
    k_w <- ncol(data_weekly)
    if(!is.null(colnames(data_weekly))){ # time series names
      seriesW_names <- colnames(data_weekly)
    }else{
      seriesW_names <- NULL
    }
    
    Y_W = matrix(NA, n, m2*k_w)
    colnames(Y_W) = paste(rep(seriesW_names, each = m2), rep(c(paste0("W", m2:1)), k_w))
    
    j = 1
    for(i in 1:n){
      Y_W[i,] = data_weekly[(j*m2):((j-1)*m2+1),] #W12_var1(over t), W11_var1,...,W1var1,..., W12_vark_w(over t), W11_vark_w,...,W1vark_w
      j = j+1
    }
    Y_VAR = cbind(Y_VAR, Y_W)
  }
  
  if(!is.null(data_daily)){
    k_d <- ncol(data_daily)
    if(!is.null(colnames(data_daily))){ # time series names
      seriesD_names <- colnames(data_daily)
    }else{
      seriesD_names <- NULL
    }
    
    Y_D = matrix(NA, n, m3*k_d)
    colnames(Y_D) = paste(rep(seriesD_names, each = m3), rep(c(paste0("D", m3:1)), k_d))
    
    j = 1
    for(i in 1:n){
      Y_D[i,] = data_daily[(j*m3):((j-1)*m3+1),] #D60_var1(over t), D59_var1,...,D58var1,..., D1_vark_1(over t), D60_vark_d(over t), D59_var_k_d,..., D1_var_k_d(over t)
      j = j+1
    }
    Y_VAR = cbind(Y_VAR, Y_D)
  }
  
  K = k_q+k_m*m1+k_w*m2+k_d*m3 # total number of regressors
  
  if(!is.null(colnames(Y_VAR))){ # time series names
    series_names <- colnames(Y_VAR)
  }else{
    series_names <- NULL
  }
  
  # output
  out <- list("Y_VAR" = Y_VAR, "series_names" = series_names,
              "K" = K, "k_q" = k_q, "k_m" = k_m, "k_w" = k_w, "k_d" = k_d)
  return(out)
}







