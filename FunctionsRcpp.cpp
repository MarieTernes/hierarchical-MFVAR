// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;


arma::vec arma_sort(arma::vec x, arma::vec y){
  // Input
  // x: vector 
  // y: vector consisting of indices according to which x is sorted
  
  // Function: x is sorted in term of y
  return x(arma::sort_index(y));
}


arma::vec softgroup(const arma::vec& u, const double& lambda){
  // Input
  // u: vector (of dimension K)
  // lambda: scalar, tuning parameter
  
  // Function: Groupwise soft-thresholding
  
  // Output
  // sg: vector after groupwise soft-thresholding
  
  arma::vec sg = std::max(0.0, 1 - lambda/arma::norm(u,2))*u;
  return(sg);
}

arma::vec softnestedgroup(const arma::vec& u, const double& lambda, const arma::vec& group_indices){
  // Input
  // u: vector (of dimension K)
  // lambda: scalar, tuning parameter
  // group_indices: vector that indicates to which nested group an element in u belongs
  
  // Function: Nested group soft-thresholding
  
  // Output
  // sg: vector after nested group soft-thresholding
  
  arma::vec r = u;
  int total_length = r.size();
  arma::vec r_gh = zeros(total_length);
  arma::uvec indices_included_groups;
  int group_length = 0; 
  int w = 0; 
  
  // Iterate through each nested group
  // In each step do groupwise soft thresholding
  for(int h=0; h <= max(group_indices); ++h){
    if(h==0){
      r_gh = r;
      r_gh = softgroup(r_gh, lambda);
      r = r_gh;
    }
    if(h!=0){
      indices_included_groups = find(group_indices >= h); 
      
      // Take only elements of variable r that are still allowed to be updated
      r_gh = r.elem(indices_included_groups);
      // weight w
      group_length = r_gh.size();
      w = ((total_length - group_length) + 1); 
      // Do weighted groupwise soft thresholding 
      r_gh = softgroup(r_gh, w * lambda); 
      // Update the values in r (and keep unchanged for future iterations) which belong to group h 
      r.elem(indices_included_groups) = r_gh;
    }
  }
  return(r);
}

arma::vec softelem(const arma::vec& u, const double& lambda){
  // Input
  // u : vector
  // lambda : scalar, tuning parameter
  
  // Function : Elementwise soft-thresholding of a vector
  
  // Output
  // r : vector after elementwise soft-thresholding
  
  arma::vec r = zeros(u.size());
  
  for(int h=0; h < u.size(); ++h){
    r(h) = ((u(h) > 0) - (u(h) < 0)) * std::max(0.0, std::abs(u(h)) - lambda);
  }
  return(r);
}


arma::mat gradient_fct(const arma::mat& y_Xkron, const arma::mat& XtranspX_kron, const arma::vec& B){
  // Input 
  // y_Xkron: matrix dimension 1 x K^2
  // XtranspX_kron: matrix dimension K^2 x K^2 
  // B: vector of dimension K^2
  
  // Function: Calculates the negative gradient of (simplified) convex and differentiable objective function
  
  arma::mat gradient = -(y_Xkron - B.t() * XtranspX_kron);
  return(gradient);
}


// [[Rcpp::export]]
Rcpp::List algorithm1_Cpp_sparse(const arma::mat& y_Xkron, const arma::sp_mat& XtranspX_kron,
                                 arma::mat B_matrix_1, arma::mat B_matrix_2, const arma::uvec& index,
                                 const arma::mat& penalty_matrix,
                                 const arma::mat& columns_matrix, 
                                 const arma::uvec& rows_L1, const arma::uvec& rows_HIER,
                                 const arma::uvec& indices_B_L1, const arma::uvec& indices_B_HIER, 
                                 const arma::vec& column_L1, const arma::vec& column_HIER,
                                 const int& K, const double& s,
                                 const double& lambda, const double& epsilon, const int& max_iter){
  
  // y_Xkron:  1 x K^2 vector (calculated from t(y)%*%Xkron)
  // XtranspX_kron: K^2 x K^2 matrix (calculated from XtranspX_kron = kronecker(I, crossprod(X)) 
  // B_matrix_1, B_matrix_2: matrices (copies) to save estimated parameter values 
  // penalty_matrix: hierarchy penalty matrix (if penalty indicator = 2 or 3)
  // columns_matrix: matrix that gives the column that a parameter belongs to in the Xkron matrix 
  // rows_L1, rows_HIER: Rows (in B_matrix) which have L1 penalty or HIER penalty
  // indices_B_L1, indices_B_HIER: Indices of the specific values in rows that have either L1 or HIER penalty
  // column_L1, column_HIER: Columns that parameter belongs to (split into either L1 or HIER penalty)
  // K: number of explanatory variables
  // s: step size for gradient
  // lambda: tuning parameter
  // epsilon: convergence criterion
  // max_iter: maximum number of iterations allowed 
  
  double r = 2.0; // iteration r 
  bool convergence = false;
  arma::mat B_old =  B_matrix_2; // make copy to compare convergence later 
  arma::vec b_new = zeros(sqrt(K));
  arma::vec b_old = zeros(sqrt(K));
  
  int row_update = 0;
  const int nrows = B_matrix_1.n_rows;
  const arma::vec& rows_total = linspace<vec>(0, nrows-1, nrows);
  arma::uvec rows_others;
  
  arma::uvec indices_B_update;  
  arma::uvec indices_B_others; 

  arma::uvec column_L1_index;
  arma::mat column_submat_update;
  arma::vec column_update;
  arma::uvec column_update_index;
  arma::mat column_submat_others;
  arma::vec column_others;
  
  arma::mat B_submat_1_update;
  arma::mat B_submat_2_update;
  arma::mat B_submat_2_others;
  
  arma::mat penalty_submat_update;
  
  arma::vec B_hat;
  arma::vec B_others;
  arma::vec B;
  arma::vec B_hat_update;
  arma::vec B_hat_update_sorted;
  
  arma::mat gradient_full;
  arma::vec gradient;
  arma::vec penalty_index;
  arma::vec update_step;
  
  
  while(convergence == false && r < max_iter){
    r = r + 1; // update iteration
    B_old = B_matrix_2; // make copy to compare convergence later 
    
    /////////////////////////////////////////////
    //-------- Update step L1 penalty ---------//
    /////////////////////////////////////////////
    
    if(rows_L1.size() != 0){
      //submatrices
      B_submat_1_update = B_matrix_1.rows(rows_L1);
      B_submat_2_update = B_matrix_2.rows(rows_L1);
      B_submat_2_others = B_matrix_2.rows(rows_HIER);
      
      //B update
      B_hat = B_submat_2_update.elem(indices_B_L1) + ((r - 2)/(r + 1)) * ( B_submat_2_update.elem(indices_B_L1)- B_submat_1_update.elem(indices_B_L1)) ;
      B_others = B_submat_2_others.elem(indices_B_HIER);
      
      B = arma_sort(join_cols(B_hat,B_others), join_cols(column_L1,column_HIER)); // order the elements in beta again
      
      // Gradient
      gradient_full = -(y_Xkron - B.t() * XtranspX_kron);
      column_L1_index = conv_to<uvec>::from( column_L1);
      gradient = gradient_full.cols(column_L1_index).t();
      
      //Update elementwise soft thresholding
      update_step = B_hat - s*gradient;
      B_hat_update = softelem(update_step, s*lambda);
      
      // Save new values
      B_hat_update_sorted = arma_sort(join_cols(B_hat_update,B_others), join_cols(column_L1,column_HIER));
      for(int i = 0; i < rows_L1.size(); ++i){
        row_update = rows_L1(i);
        indices_B_update = arma::find_finite(B_matrix_2.row(row_update));
        column_submat_update = columns_matrix.row(row_update);
        column_update = column_submat_update.elem(indices_B_update);
        column_update_index = conv_to<uvec>::from( column_update);
        
        B_matrix_1(row_update, span(0, indices_B_update.size() -1) ) = B_matrix_2(row_update, span(0, indices_B_update.size() -1) );
        
        //Save new updated values 
        B_matrix_2(row_update, span(0, indices_B_update.size() -1) )  = B_hat_update_sorted.elem(column_update_index).t();
      }
      
    }
    /////////////////////////////////////////////
    //------ Update step HIER penalty ---------//
    /////////////////////////////////////////////
    
    if(rows_HIER.size() != 0){
      for(int i=0; i < rows_HIER.size(); ++i){
        row_update = rows_HIER(i);
        rows_others = find(rows_total != row_update);
        
        indices_B_update = arma::find_finite(B_matrix_2.row(row_update));
        indices_B_others = arma::find_finite(B_matrix_2.rows(rows_others));
        
        //submatrices 
        column_submat_update = columns_matrix.row(row_update);
        column_submat_others = columns_matrix.rows(rows_others);
        
        column_update = column_submat_update.elem(indices_B_update);
        column_others =  column_submat_others.elem(indices_B_others);
        
        B_submat_1_update = B_matrix_1.row(row_update);
        B_submat_2_update = B_matrix_2.row(row_update);
        B_submat_2_others = B_matrix_2.rows(rows_others);
        
        penalty_submat_update = penalty_matrix.row(row_update);
        
        // B update 
        B_hat = B_submat_2_update.elem(indices_B_update) + ((r - 2)/(r + 1))*( B_submat_2_update.elem(indices_B_update)- B_submat_1_update.elem(indices_B_update)) ;
        
        B_matrix_1(row_update, span(0, indices_B_update.size() -1) ) = B_matrix_2(row_update, span(0, indices_B_update.size() -1) );
        B_others = B_submat_2_others.elem(indices_B_others);
        
        B = arma_sort(join_cols(B_hat,B_others), join_cols(column_update,column_others)); // order the elements in beta again
        
        // Gradient
        gradient_full = -(y_Xkron - B.t() * XtranspX_kron);
        column_update_index = conv_to<uvec>::from( column_update);
        gradient = gradient_full.cols(column_update_index).t();
        
        // Penalty
        penalty_index = penalty_submat_update.elem(indices_B_update);
        
        //Update nestedgroup soft thresholding
        update_step = B_hat - s*gradient;
        B_hat_update = softnestedgroup(update_step, s*lambda, penalty_index);
        
        //Save new updated values 
        B_matrix_2(row_update, span(0, indices_B_update.size() -1) )  = B_hat_update.t();
      }
    }
    //Get all non NA values of old and updated B_matrix 
    b_old = B_old.elem(index);
    b_new = B_matrix_2.elem(index);
    
    arma::vec d = abs(b_old-b_new);
    double c = d.max();
   
    //Check if estimator has reached convergence
    if(c <= epsilon){
      convergence = true;
    }
  }
  
  // Reorganize B_new = B_matrix_2[index] so that in the end function returns the right VAR matrix structure
  arma::vec columns = columns_matrix.elem(index);
  arma::vec B_final = arma_sort(b_new, columns);
  
  Rcpp::List results = Rcpp::List::create(
                       Rcpp::Named("B") = B_final,
                       Rcpp::Named("iter") = r);
  return(results);
}


// [[Rcpp::export]]
Rcpp::List PLS_Cpp(const arma::mat& y_Xkron, const arma::sp_mat& XtranspX_kron,
                    arma::mat B_matrix_1, arma::mat B_matrix_2, const arma::uvec& index,
                    const arma::mat& penalty_matrix,
                    const arma::mat& columns_matrix, 
                    const arma::uvec& rows_L1, const arma::uvec& rows_HIER,
                    const arma::uvec& indices_B_L1, const arma::uvec& indices_B_HIER, 
                    const arma::vec& column_L1, const arma::vec& column_HIER,
                    const int& K, const double& s,
                    const arma::vec& lambda_beta, const double& epsilon, const int& max_iter){
  
  arma::mat betas_PLS = zeros(pow(K,2), lambda_beta.size());
  arma::vec iter_PLS = zeros(lambda_beta.size());
  
  for(int i = 0; i < lambda_beta.size(); ++i){
    Rcpp::List list_algorithm1Cpp_PLS = algorithm1_Cpp_sparse(y_Xkron, XtranspX_kron,
                                                              B_matrix_1, B_matrix_2,
                                                              index, penalty_matrix,
                                                              columns_matrix,
                                                              rows_L1, rows_HIER,
                                                              indices_B_L1, indices_B_HIER,
                                                              column_L1, column_HIER,
                                                              K, s, lambda_beta(i), epsilon, max_iter);
    
    arma::vec B_vec = list_algorithm1Cpp_PLS["B"];
    int iter = list_algorithm1Cpp_PLS["iter"];
    betas_PLS.col(i) = B_vec;
    iter_PLS(i) = iter;
    
    // warm start for next, smaller lambda (previous solution is new starting value)
    arma::vec columns = columns_matrix.elem(index);  
    arma::uvec columns_index = conv_to<uvec>::from(columns);
    B_matrix_1.elem(index) = B_vec.elem(columns_index);
    B_matrix_2 = B_matrix_1;
  }
  Rcpp::List results = Rcpp::List::create(
    Rcpp::Named("betas_PLS") = betas_PLS,
    Rcpp::Named("iter_PLS") = iter_PLS);
  return(results);
}

