HdS = function(design_obj,fitted_beta){
  # cdm_obj is an object returned by the create_design_matrix function
  # fitted_beta is a vector of the coeff obtained after the fitting of a model
  
  # getting all the info we need 
  pmain = design_obj$pmain
  pinterac = pmain*(pmain-1)/2
  pquadra = pmain
  main = design_obj$xtilde[,1:pmain]
  means_main = design_obj$means[1:pmain]
  sds_main = design_obj$sds[1:pmain]
  interac = design_obj$xtilde[,(pmain+1):(pmain+pinterac)]
  means_interac = design_obj$means[(pmain+1):(pmain+pinterac)]
  sds_interac = design_obj$sds[(pmain+1):(pmain+pinterac)]
  quadra = design_obj$xtilde[,(pmain+pinterac+1):(2*pmain+pinterac)]
  means_quadra = design_obj$means[(pmain+pinterac+1):(2*pmain+pinterac)]
  sds_quadra = design_obj$sds[(pmain+pinterac+1):(2*pmain+pinterac)]

  fitted_beta_main = fitted_beta[1:pmain]
  fitted_beta_interac = fitted_beta[(pmain+1):(pmain+pinterac)]
  fitted_beta_quadra = fitted_beta[(pmain+pinterac+1):(pmain+pinterac+pquadra)]
  
  tmp_names = make_interac_names(pmain)
  
  # 1. Hirarcical Descaling for the quadratic effects
  HdS_quadra = function(j){return(fitted_beta_quadra[j]/(sds_main[j])**2)}
  HdS_beta_quadra = unlist(sapply(1:pquadra,HdS_quadra))
  
  # 2. Hierarchical descaling for the interac effects
  HdS_interac = function(j){
    
    # The idea of the lines below
    # For a fixed j:
    #  - look for and store the indexes of the variable names starting by j ( for instance if pmain = 5 for j = 1 we look fore the indexes of the name 12, 13, 14 and 15)
    # - we create a temporary vector of the coeff of interac effect and we fill it with the coeff of the interac effects selected by the previous step
    # - for a fixed j we nedd to find all the k such that gamma_j_k exist (with the same example as before k = 2,3,4,5) in order to find the standard deviation of the associate covar X_k. We store the values in a temporary vector of sd
    # - Finally for a fixed j since j<k and k in [1,p] wa have two vectors of size p-j: we compute the element wise quotient and get a vector of size p-j
    # - we returna matrix of two columns ans p-j lines. The first column is the original idexes of the coeff, the second one is the descaled coeff
    tmp_idx = which(tmp_names[,1]==j)
    tmp_beta_interac = fitted_beta_interac[tmp_idx]
    tmp_sds_main = sds_main[(j+1):(pmain)]
    
    HdS_beta_interac = tmp_beta_interac/(tmp_sds_main*sds_main[j])
    
    return(unname(cbind(tmp_idx,HdS_beta_interac)))
  }
  # we apply the above function to all j in 1:(p -1) precisely because in gamma_j_k j needs to be strictly lower than k. 
  # We stock the result in a list of p-1 elements.
  # Each element is a 2 column matrix: we concatenate those p-1 elements into one matrix with the 
  # do.call function then reorder it by its first column and finally return the reordered second column.
  HdS_beta_interac = do.call(rbind,(sapply(1:(pmain-1),HdS_interac)))
  HdS_beta_interac = HdS_beta_interac[order(HdS_beta_interac[,1]),2]
  
  # 3. Hierarchical Descaling for the main effects
  HdS_main = function(j){
    
    # computing the sum (see the formulas in chen et al 2020)
    
    tmp_idx = which(tmp_names[,1]==j | tmp_names[,2] == j)
    #print(tmp_idx)
    tmp_beta_interac = fitted_beta_interac[tmp_idx]
    #print(tmp_beta_interac)
    tmp_means_main = means_main[-j]
    tmp_sds_main = sds_main[-j]
    # To avoid a loop on k we write the sum of the formula as a dot product
    sum = (tmp_beta_interac %*% (tmp_means_main/(tmp_sds_main*sds_main[j])))
    
    HdS_beta_main = fitted_beta_main[j]/sds_main[j] - (2*fitted_beta_quadra[j]*means_main[j])/(sds_main[j]**2) - sum
    return(HdS_beta_main)
    
    
  }
  HdS_beta_main = unlist(sapply(1:pmain,HdS_main))
  
  # we combine the three vectors obtained and return it in its original order: main, interac and quadratic effect. The vector returned is now comparable (in term of indexes) with the one given in argument to the function HdS()
  HdS_beta = c(HdS_beta_main,HdS_beta_interac,HdS_beta_quadra)
  names(HdS_beta) = colnames(xtilde)
  return(HdS_beta)
  
}
