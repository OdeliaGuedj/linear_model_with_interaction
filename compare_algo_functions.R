compare_algo_quadraGaussian = function(design_obj, output_obj){
  
  nb = length(output_obj$true_beta[which(output_obj$true_beta[1:design_obj$pmain] !=0) ])
  true_beta = output_obj$true_beta
  pmain = design_obj$pmain
  n = nrow(design_obj$xtilde)
  
  RMSE = matrix(nrow = 1, ncol = 9); MSH = matrix(nrow = 1, ncol = 9); 
  sensitivity = matrix(nrow = 1, ncol = 9); specificity = matrix(nrow = 1, ncol = 9); 
  AUC = matrix(nrow = 1, ncol = 9); main_coverage = matrix(nrow = 1, ncol = 9); 
  order2_coverage = matrix(nrow = 1, ncol = 9); main_exact_select = matrix(nrow = 1, ncol = 9); 
  order2_exact_select = matrix(nrow = 1, ncol = 9); model_size = matrix(nrow = 1, ncol = 9);
  times_ = matrix(nrow = 1, ncol = 9)
  
  algo_names = c("ORACLE", "All Pairs LASSO", "HdS LASSO", "RAMP", "HierNet", "HierNetProj", "FAMILY", "PIE", "SPRINTR")
  colnames(RMSE) = algo_names; colnames(MSH) = algo_names; colnames(sensitivity) = algo_names;
  colnames(specificity) = algo_names; colnames(AUC) = algo_names; colnames(main_coverage) = algo_names;
  colnames(order2_coverage) = algo_names; colnames(main_exact_select) = algo_names; 
  colnames(order2_exact_select) = algo_names; colnames(model_size) = algo_names; colnames(times_) = algo_names;

  trainIndex = caret::createDataPartition(1:n, p=0.8, list=FALSE)
  x_train = as.data.frame(design_obj$xtilde[trainIndex,1:pmain]); x_test = as.data.frame(design_obj$xtilde[-trainIndex,1:pmain]);
  xtilde_train = as.data.frame(design_obj$xtilde[trainIndex,]); xtilde_test = as.data.frame(design_obj$xtilde[-trainIndex,]);
  xtilde_Hs_train = as.data.frame(design_obj$xtilde_Hs[trainIndex,]); xtilde_Hs_test = as.data.frame(design_obj$xtilde_Hs[-trainIndex,]);
  y_train = output_obj$Y[trainIndex]; y_test = output_obj$Y[-trainIndex];
  
  interac_names= make_interac_names(pmain)
  quadra_names= make_quadra_names(pmain)
  

  #Oracle
  print("ORACLE")
  start_time = Sys.time()
  oracle_fit = lm(as.matrix(xtilde_train) %*% true_beta ~., data = xtilde_train)
  oracle_beta = oracle_fit$coefficients[-1]
  oracle_predict  = predict(oracle_fit, newx = xtilde_test)
  end_time = Sys.time()
  times_[1,"ORACLE"] = as.numeric(end_time-start_time)
  RMSE[1,"ORACLE"] =  compute_RMSE(y_test,oracle_predict)
  MSH[1,"ORACLE"] =  compute_MSH(true_beta = true_beta ,fitted_beta = oracle_beta, main_effects = nb, pmain = design_obj$pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = oracle_beta)
  sensitivity[1,"ORACLE"] = performance$sensitivity
  specificity[1,"ORACLE"] = performance$specificity
  AUC[1,"ORACLE"] = performance$auc
  main_coverage[1,"ORACLE"] = compute_main_coverage(true_beta = true_beta, fitted_beta = oracle_beta,pmain = pmain)
  order2_coverage[1,"ORACLE"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = oracle_beta,pmain = pmain)
  main_exact_select[1,"ORACLE"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = oracle_beta,pmain = pmain)
  order2_exact_select[1,"ORACLE"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = oracle_beta,pmain = pmain)
  model_size[1,"ORACLE"] = compute_model_size(fitted_beta = oracle_beta)
  
  #All pairs LASSO
  print("All pairs LASSO")
  start_time = Sys.time()
  lambda_min_allPairsLASSO = glmnet::cv.glmnet(as.matrix(xtilde_train),y_train,alpha=1,standardize =T)$lambda.min
  all_pair_lasso_fit = glmnet::glmnet(as.matrix(xtilde_train),y_train,alpha=1,standardize =T, lambda = lambda_min_allPairsLASSO)
  all_pair_lasso_beta = as.vector(all_pair_lasso_fit$beta)
  all_pair_lasso_predict = predict(all_pair_lasso_fit, s = lambda_min_allPairsLASSO, newx = as.matrix(xtilde_test))
  end_time = Sys.time()
  times_[1,"All Pairs LASSO"] = as.numeric(end_time-start_time)
  RMSE[1,"All Pairs LASSO"] =  compute_RMSE(y_test,all_pair_lasso_predict)
  MSH[1,"All Pairs LASSO"] =  compute_MSH(true_beta = true_beta ,fitted_beta = all_pair_lasso_beta, main_effects = nb, pmain = pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = all_pair_lasso_beta)
  sensitivity[1,"All Pairs LASSO"] = performance$sensitivity
  specificity[1,"All Pairs LASSO"] = performance$specificity
  AUC[1,"All Pairs LASSO"] = performance$auc
  main_coverage[1,"All Pairs LASSO"] = compute_main_coverage(true_beta = true_beta, fitted_beta = all_pair_lasso_beta,pmain = pmain)
  order2_coverage[1,"All Pairs LASSO"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = all_pair_lasso_beta,pmain = pmain)
  main_exact_select[1,"All Pairs LASSO"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = all_pair_lasso_beta,pmain = pmain)
  order2_exact_select[1,"All Pairs LASSO"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = all_pair_lasso_beta,pmain = pmain)
  model_size[1,"All Pairs LASSO"] = compute_model_size(fitted_beta = all_pair_lasso_beta)
  
  #Hds LASSO
  print("Hds LASSO Chen et al 2020")
  start_time = Sys.time()
  lambda_min_HdsLASSO = find_best_lambda(x_train = xtilde_Hs_train,
                                         y_train = y_train,
                                         x_valid = xtilde_Hs_test[1:floor(0.5*n),],
                                         y_valid = y_test[1:floor(0.5*n)],
                                         method = 2)$lambda_min_2
  HdS_lasso_fit = glmnet::glmnet(as.matrix(xtilde_Hs_train),y_train,alpha=1,standardize =F, lambda = lambda_min_HdsLASSO)
  HdS_lasso_beta = HdS(simu_$data,HdS_lasso_fit$beta)
  HdS_lasso_fit$beta = as.matrix(HdS_lasso_beta)
  HdS_lasso_predict = predict(HdS_lasso_fit, s = lambda_min_HdsLASSO, newx = as.matrix(x_H_s_test))
  end_time = Sys.time()
  times_[1,"HdS LASSO"] = as.numeric(end_time-start_time)
  RMSE[1,"HdS LASSO"] =  compute_RMSE(y_test,HdS_lasso_predict)
  MSH[1,"HdS LASSO"] =  compute_MSH(true_beta = true_beta ,fitted_beta = HdS_lasso_beta, main_effects = nb, pmain = pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = HdS_lasso_beta)
  sensitivity[1,"HdS LASSO"] = performance$sensitivity
  specificity[1,"HdS LASSO"] = performance$specificity
  AUC[1,"HdS LASSO"] = performance$auc
  main_coverage[1,"HdS LASSO"] = compute_main_coverage(true_beta = true_beta, fitted_beta = HdS_lasso_beta,pmain = pmain)
  order2_coverage[1,"HdS LASSO"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = HdS_lasso_beta,pmain = pmain)
  main_exact_select[1,"HdS LASSO"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = HdS_lasso_beta,pmain = pmain)
  order2_exact_select[1,"HdS LASSO"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = HdS_lasso_beta,pmain = pmain)
  model_size[1,"HdS LASSO"] = compute_model_size(fitted_beta = HdS_lasso_beta)
  
  
  #RAMP
  print("RAMP Hao et al 2016")
  start_time = Sys.time()
  ramp_fit = RAMP::RAMP(main_train,y_train)
  ramp_beta = rep(0,get_col_number(pmain))
  names(ramp_beta) = names(true_beta)
  ramp_beta[ramp_fit$mainInd] = ramp_fit$beta.m
  ramp_beta[trimws(gsub("X"," ",ramp_fit$interInd),"l")] = ramp_fit$beta.i
  ramp_predict = predict(ramp_fit, newdata = main_test)
  end_time = Sys.time()
  times_[1,"RAMP"] = as.numeric(end_time-start_time)
  RMSE[1,"RAMP"] =  compute_RMSE(y_test,ramp_predict)
  MSH[1,"RAMP"] =  compute_MSH(true_beta = true_beta ,fitted_beta = ramp_beta, main_effects = nb, pmain = pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = ramp_beta)
  sensitivity[1,"RAMP"] = performance$sensitivity
  specificity[1,"RAMP"] = performance$specificity
  AUC[1,"RAMP"] = performance$auc
  main_coverage[1,"RAMP"] = compute_main_coverage(true_beta = true_beta, fitted_beta = ramp_beta,pmain = pmain)
  order2_coverage[1,"RAMP"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = ramp_beta,pmain = pmain)
  main_exact_select[1,"RAMP"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = ramp_beta,pmain = pmain)
  order2_exact_select[1,"RAMP"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = ramp_beta,pmain = pmain)
  model_size[1,"RAMP"] = compute_model_size(fitted_beta = ramp_beta)
  
  
  #HierNet
  print("HierNet Bien et al 2013")
  start_time = Sys.time()
  lambda_min_HierNet = hierNet::hierNet.cv(hierNet::hierNet.path(scale(main_train,T,T),y_train),scale(main_train,T,T),y_train, trace =0)$lamhat.1se
  hiernet_fit = hierNet::hierNet(scale(main_train,T,T),y_train,lam=lambda_min_HierNet,strong = strong_hiernet, trace = 0)
  hiernet_beta = c(hiernet_fit$bp - hiernet_fit$bn,
                   sapply(1:dim(interac_names)[1],function(i) hiernet_fit$th[interac_names[i,1],interac_names[i,2]]) ,
                   diag(hiernet_fit$th))
  names(hiernet_beta) = c(1:pmain, interac_names[,3], quadra_names)
  hiernet_predict = predict(hiernet_fit,scale(main_test,T,T))
  end_time = Sys.time()
  times_[1,"HierNet"] = as.numeric(end_time-start_time)
  RMSE[1,"HierNet"] =  compute_RMSE(y_test,hiernet_predict)
  MSH[1,"HierNet"] =  compute_MSH(true_beta = true_beta ,fitted_beta = hiernet_beta, main_effects = nb, pmain = pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = hiernet_beta)
  sensitivity[1,"HierNet"] = performance$sensitivity
  specificity[1,"HierNet"] = performance$specificity
  AUC[1,"HierNet"] = performance$auc
  main_coverage[1,"HierNet"] = compute_main_coverage(true_beta = true_beta, fitted_beta = hiernet_beta,pmain = pmain)
  order2_coverage[1,"HierNet"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = hiernet_beta,pmain = pmain)
  main_exact_select[1,"HierNet"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = hiernet_beta,pmain = pmain)
  order2_exact_select[1,"HierNet"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = hiernet_beta,pmain = pmain)
  model_size[1,"HierNet"] = compute_model_size(fitted_beta = hiernet_beta)
  
  
  #Family
  print("FAMILY Harris et al 2016")
  W_test = c()
  W_train = c()
  # !!!!!!!! BOUCLES Comment faire mieux ??
  main_s_train = scale(main_train, scale = F)
  main_s_test = scale(main_test, scale = F)
  for(i in 1:(pmain+1)){
    for(j in 1:(pmain+1)){
      W_test =  cbind(W_test,cbind(1,main_s_test)[,j]*cbind(1,main_s_test)[,i])
      W_train = cbind(W_train,cbind(1,main_s_train)[,j]*cbind(1,main_s_train)[,i])
    }
  }
  B = matrix(0,ncol = pmain+1,nrow = pmain+1)
  rownames(B)<- c("inter" , 1:(nrow(B)-1))
  colnames(B)<- c("inter" , 1:(nrow(B)-1))
  B[apply(sapply(sample_, function(i) i== colnames(B)), 2, function(i) which(i == T)),] = rep(1,(pmain+1))
  B[,apply(sapply(sample_, function(i) i== colnames(B)), 2, function(i) which(i == T))] = rep(1,(pmain+1))
  Y_train_FAMILY = as.vector(W_train%*%as.vector(B)+rnorm(nrow(W_train),sd = noise_sd))
  Y_test_FAMILY = as.vector(W_test%*%as.vector(B)+rnorm(nrow(W_test),sd = noise_sd))
  
  alphas = c(0.01,0.5,0.99)
  lambdas = seq(0.1,1,length = 50)
  start_time = Sys.time()
  family_fit =  FAMILY(main_s_train, main_s_train, Y_train_FAMILY, lambdas ,alphas, quad = T,iter=500, verbose = TRUE )
  family_all_pred=  predict(family_fit, as.matrix(main_s_test), as.matrix(main_s_test))
  mse_family = apply(family_all_pred,c(2,3), "-" ,Y_test_FAMILY)
  mse_family =  apply(mse_family^2,c(2,3),sum)
  im = which(mse_family==min(mse_family),TRUE)
  family_predict = family_all_pred[,im[,1],im[,2]]
  end_time = Sys.time()
  times_[1,"FAMILY"] = as.numeric(end_time-start_time)
  beta_family = coef(family_fit, XequalZ =  T)[[im[2]]][[im[1]]]
  mains = beta_family$mains
  
  quadra = beta_family$interacts[which(beta_family$interacts[,1]==beta_family$interacts[,2]),]
  quadra = cbind(paste(quadra[,1],quadra[,2]),quadra[,3])
  quadra2 = data.frame(names = quadra_names, values = rep(0,length(quadra_names)))
  quadra2$values[quadra2$names %in% quadra[,1]] = quadra[,2]
  quadra = quadra2
  interac = beta_family$interacts[which((beta_family$interacts[,1]!=beta_family$interacts[,2])&(beta_family$interacts[,1]<beta_family$interacts[,2])),]
  interac = cbind(paste(interac[,1],interac[,2]),interac[,3])
  interac2 = data.frame(names = interac_names[,3], values = rep(0,dim(interac_names)[1]))
  interac2$values[interac2$names %in% interac[,1]] = interac[,2]
  interac = interac2
  family_beta = as.numeric(c(mains[,2], interac[,2],quadra[,2]))
  names(family_beta) = c(mains[,1],interac[,1],quadra[,1])
  RMSE[1,"FAMILY"] =  compute_RMSE(y_test,family_predict)
  MSH[1,"FAMILY"] =  compute_MSH(true_beta = true_beta ,fitted_beta = family_beta, main_effects = nb, pmain = pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = family_beta)
  sensitivity[1,"FAMILY"] = performance$sensitivity
  specificity[1,"FAMILY"] = performance$specificity
  AUC[1,"FAMILY"] = performance$auc
  main_coverage[1,"FAMILY"] = compute_main_coverage(true_beta = true_beta, fitted_beta = family_beta,pmain = pmain)
  order2_coverage[1,"FAMILY"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = family_beta,pmain = pmain)
  main_exact_select[1,"FAMILY"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = family_beta,pmain = pmain)
  order2_exact_select[1,"FAMILY"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = family_beta,pmain = pmain)
  model_size[1,"FAMILY"] = compute_model_size(fitted_beta = family_beta)
  
  #PIE
  print("PIE Wang et al 2019")
  start_time = Sys.time()
  main_beta = as.vector(coef(cv.glmnet(as.matrix(main_train),y_train),s="lambda.min"))[-1]
  names(main_beta) = 1:pmain
  coeff_PIE = PIE(main_train,y_train) # la cv est faites à l'intereieur par defaut avec une grille de 50 lambdas
  colnames(coeff_PIE) = 1:pmain
  rownames(coeff_PIE) = 1:pmain
  quadra_pie = diag(coeff_PIE)
  names(quadra_pie) = make_quadra_names(pmain)
  coeff_PIE[lower.tri(coeff_PIE, diag = T)] = NA
  coeff_PIE = as.data.frame(as.matrix(coeff_PIE))
  coeff_PIE$r = rownames(coeff_PIE)
  interac_pie = melt(coeff_PIE, na.rm = T, id.vars = "r", value.name = "coeff", variable.name = "c")
  interac_pie$r = as.numeric(as.character(interac_pie$r))
  interac_pie$c = as.numeric(as.character(interac_pie$c))
  interac_pie_names = paste(interac_pie$r,interac_pie$c)
  interac_pie = interac_pie$coeff
  names(interac_pie) = interac_pie_names
  pie_beta = c(main_beta,interac_pie,quadra_pie)
  pie_predict = as.matrix(x_test)%*%pie_beta
  end_time = Sys.time()
  times_[1,"PIE"] = as.numeric(end_time-start_time)
  RMSE[1,"PIE"] =  compute_RMSE(y_test,pie_predict)
  MSH[1,"PIE"] =  compute_MSH(true_beta = true_beta ,fitted_beta = pie_beta, main_effects = nb, pmain = pmain)
  performance = eval_performance(true_beta = true_beta, fitted_beta = pie_beta)
  sensitivity[1,"PIE"] = performance$sensitivity
  specificity[1,"PIE"] = performance$specificity
  AUC[1,"PIE"] = performance$auc
  main_coverage[1,"PIE"] = compute_main_coverage(true_beta = true_beta, fitted_beta = pie_beta,pmain = pmain)
  order2_coverage[1,"PIE"] = compute_order2_coverage(true_beta = true_beta, fitted_beta = pie_beta,pmain = pmain)
  main_exact_select[1,"PIE"] = compute_main_exact_select(true_beta = true_beta, fitted_beta = pie_beta,pmain = pmain)
  order2_exact_select[1,"PIE"] = compute_order2_exact_select(true_beta = true_beta, fitted_beta = pie_beta,pmain = pmain)
  model_size[1,"PIE"] = compute_model_size(fitted_beta = pie_beta)
  print("fini")
  
  betas = cbind(true_beta, oracle_beta, all_pair_lasso_beta, HdS_lasso_beta, ramp_beta, hiernet_beta, family_beta, pie_beta)
  
  return(list(RMSE = RMSE,
              MSH = MSH,
              sens = sensitivity,
              spec = specificity,
              auc = AUC,
              main_coverage = main_coverage,
              order2_coverage = order2_coverage,
              main_exact_select = main_exact_select,
              order2_select_exact = order2_exact_select,
              model_size = model_size,
              times_ = times_,
              betas = betas))
}
