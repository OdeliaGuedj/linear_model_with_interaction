# 1. Functions to evaluate the performances of an algorithm
eval_performance = function(true_beta,fitted_beta){
  TP = length(which((fitted_beta != 0)&(true_beta != 0)))
  TN = length(which((fitted_beta == 0)&(true_beta == 0)))
  FP = length(which((fitted_beta != 0)&(true_beta == 0)))
  FN = length(which((fitted_beta == 0)&(true_beta != 0)))
  sensitivity_ = TP/(TP+FN)
  specificity_ = TN/(TN+FP)
  roc_ = pROC::roc(num_to_binary(true_beta),num_to_binary(fitted_beta), quiet = T)
  auc_ = pROC::auc(roc_)
  return(data.frame(sensitivity = sensitivity_, specificity = specificity_, auc = auc_))
}

compute_MSH = function(true_beta,fitted_beta,main_effects,pmain){
  fitted_beta_main = fitted_beta[1:pmain]
  fitted_beta_interac = fitted_beta[(pmain+1):(pmain*(pmain+1)/2)]
  fitted_beta_quadra = fitted_beta[(pmain*(pmain+1)/2 + 1 ):length(fitted_beta)]
  interac_names = make_interac_names(pmain)
  # we select the indexes of the selected interaction effect and their names
  selected_interac_idx = which(fitted_beta_interac != 0)
  selected_interac_names = c(interac_names[selected_interac_idx,1],interac_names[selected_interac_idx,2])
  # the parent of XjXk are Xj and Xk : we find them with the interaction names and remove the duplicates (we do the same with the quadratic effects)
  parent_of_selected_interac = selected_interac_names[!duplicated(selected_interac_names)]
  parent_of_selected_quadra = which(fitted_beta_quadra != 0)
  # theoretical total number of non-zero parent main-effects (parent of selected order2)
  parent_of_selected_order2 = union(parent_of_selected_interac,parent_of_selected_quadra)
  # correctly identified non‐zero parent main effects (selected parent of selected order2)
  true_parent_of_selected_order2 = intersect(parent_of_selected_order2, which(fitted_beta_main != 0))
  num = length(true_parent_of_selected_order2)
  den = length(parent_of_selected_order2)
  MSH = ifelse(length(parent_of_selected_order2) == 0, 1, num/den)
  return(MSH)
}

compute_RMSE = function(y_test, y_pred) {
  return(sqrt(mean((y_pred - y_test)^2)))}

compute_main_coverage = function(true_beta,fitted_beta,pmain){
  true_main_coeff = true_beta[1:pmain]
  fitted_main_coeff = fitted_beta[1:pmain]
  names(true_main_coeff) = 1:pmain
  names(fitted_main_coeff) = 1:pmain
  activ_true_main_coeff = unname(which(true_main_coeff != 0))
  activ_fitted_main_coeff = unname(which(fitted_main_coeff != 0))
  if((sum(activ_true_main_coeff %in% activ_fitted_main_coeff) == length(activ_true_main_coeff))&(sum(activ_true_main_coeff %in% activ_fitted_main_coeff) <= length(activ_fitted_main_coeff)))
    main.cov = 1
  else
    main.cov = 0
  return(main.cov)
}

compute_order2_coverage = function(true_beta,fitted_beta,pmain){
  true_order2_coeff = true_beta[(pmain+1):length(true_beta)]
  fitted_order2_coeff = fitted_beta[(pmain+1):length(fitted_beta)]
  names(true_order2_coeff) = c(make_interac_names(pmain)[,3],make_quadra_names(pmain))
  names(fitted_order2_coeff) = names(true_order2_coeff)
  activ_true_order2_coeff = unname(which(true_order2_coeff != 0))
  activ_fitted_order2_coeff = unname(which(fitted_order2_coeff != 0))
  if((sum(activ_true_order2_coeff %in% activ_fitted_order2_coeff) == length(activ_true_order2_coeff))&(sum(activ_true_order2_coeff %in% activ_fitted_order2_coeff) <= length(activ_fitted_order2_coeff)))
    order2.cov = 1
  else
    order2.cov = 0
  return(order2.cov)
}

compute_main_exact_select = function(true_beta,fitted_beta,pmain){
  true_main_coeff = true_beta[1:pmain]
  fitted_main_coeff = fitted_beta[1:pmain]
  names(true_main_coeff) = 1:pmain
  names(fitted_main_coeff) = 1:pmain
  activ_true_main_coeff = unname(which(true_main_coeff != 0))
  activ_fitted_main_coeff = unname(which(fitted_main_coeff != 0))
  if(setequal(activ_true_main_coeff, activ_fitted_main_coeff))  
    main.exact.select = 1
  else
    main.exact.select = 0
  return(main.exact.select)
}

compute_order2_exact_select = function(true_beta,fitted_beta,pmain){
  true_order2_coeff = true_beta[(pmain+1):length(true_beta)]
  fitted_order2_coeff = fitted_beta[(pmain+1):length(fitted_beta)]
  names(true_order2_coeff) = c(make_interac_names(pmain)[,3],make_quadra_names(pmain))
  names(fitted_order2_coeff) = names(true_order2_coeff)
  activ_true_order2_coeff = unname(which(true_order2_coeff != 0))
  activ_fitted_order2_coeff = unname(which(fitted_order2_coeff != 0))
  if(setequal(activ_true_order2_coeff, activ_fitted_order2_coeff))
    order2.exact.select = 1
  else
    order2.exact.select = 0
  return(order2.exact.select)
}


compute_model_size = function(fitted_beta){
  return(length(which(fitted_beta!=0)))
}
# 2. Intermediate functions
num_to_binary = function(vec){return(as.vector(ifelse(vec == 0,0,1)))}

make_interac_names = function(pmain){
  num_names = combinations(n = pmain, r = 2, v = 1:pmain, repeats.allowed = F)
  num_names = cbind(num_names,apply(num_names,1,function(x)paste(x[1],x[2])))
  return(num_names)
}

make_quadra_names = function(pmain){
  num_names = sapply(1:pmain, function(x) paste(x,x))
  return(num_names)
}