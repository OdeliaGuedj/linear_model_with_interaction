# 1. Simu the design matrix 

simu_main_var = function(n = 1000, pmain = 10, rho = 0, sigma = 1, compute_mean_sd = T){
  create_cor = function(i){return(rho**(i-1)*sigma^2)}
  vec = sapply(c(1:pmain), FUN = create_cor)
  Sigma = lqmm::make.positive.definite(toeplitz(vec))
  main = mvtnorm::rmvnorm(n, mean = rep(0,pmain) , sigma = Sigma)
  colnames(main) = 1:pmain
  if(compute_mean_sd){
    means_main = apply(main,2,mean)
    sds_main = apply(main,2,sd)}
  else
  {
    means_main = NULL
    sds_main = NULL}
  return(list(X = as.data.frame(main), means_X = means_main, sds_X = sds_main, Sigma = Sigma))
}

generate_order2_from_main = function(X){
  pmain = ncol(X)
  quadra = sapply(1:pmain, function(i){X[,i] * X[,i]})
  colnames(quadra) = make_quadra_names(pmain)
  if(pmain == 1){
    interac = NULL
    means_interac = NULL
    sds_interac = NULL
    means_quadra = mean(quadra)
    sds_quadra = sd(quadra)
  }
  else if(pmain > 1){
      names_interac = make_interac_names(pmain)
      interac = sapply(1:nrow(names_interac),FUN = function(i){X[,names_interac[i,1]] * X[,names_interac[i,2]]})
      colnames(interac) = names_interac[,3]
      means_interac = apply(interac,2,mean)
      sds_interac = apply(interac,2,sd)
      means_quadra = apply(quadra,2,mean)
      sds_quadra = apply(quadra,2,sd)
  }
  
  return(list(Z = as.matrix(cbind(interac,quadra)), means_Z = c(means_interac,means_quadra), sds_Z= c(sds_interac, sds_quadra)))
}

simu_design_matrix = function(n=1000, pmain = 10, rho = 0,sigma = 1, H_scale = F, ordinary_scale = F){
  main_obj = simu_main_var(n, pmain, rho, sigma, compute_mean_sd = T)
  X = main_obj$X
  order2_obj = generate_order2_from_main(X)
  xtilde = as.matrix(cbind(X,order2_obj$Z))
  if(ordinary_scale){xtilde_s = as.matrix(scale(xtilde))} else {xtilde_s = NULL}
  if(H_scale == T){
    X_s = scale(X)
    order2_obj_s = generate_order2_from_main(X_s)
    xtilde_Hs = as.matrix(cbind(X_s,Z_s = order2_obj_s$Z))
  }
  else{xtilde_Hs = NULL}

  return(list(xtilde=xtilde, xtilde_Hs=xtilde_Hs, xtilde_s = xtilde_s, pmain=pmain, means = c(main_obj$means_X, order2_obj$means_Z),
              sds = c(main_obj$sds_X, order2_obj$sds_Z), rho = rho))
}

compute_design_proj_matrix = function(xtilde, pmain){
  xtilde_proj = xtilde
  for (j in ((pmain+1):(pmain + (pmain) * (pmain+1)/2))){
    model_j = lm(xtilde[,j]~xtilde[,1:pmain])
    xtilde_proj[,j] = xtilde_proj[,j] - predict(model_j)
  }
  return(xtilde_proj)
}

# 2. Simu a quadratic model with strong or weak hierarchy
simu_quadraGauss_output = function(design, pmain, nb = 3, heredity = "strong",
                                   SNR=10,noise_sd = NULL, main_idx_sample = NULL, 
                                   true_beta = NULL){
  if(is.null(true_beta)){
    true_beta = rep(0, ncol(design))
    names(true_beta) = colnames(design)
    if(heredity == "strong"){
      if(is.null(main_idx_sample)) {non_0_main_idx = sort(sample(1:pmain,nb))}
      non_0_quadra_idx = sapply(non_0_main_idx, function(x) paste(x,x))
      non_0_interac_idx = apply(combinations(n=nb,r=2,repeats.allowed = F, v = non_0_main_idx),1,
                                function(i) paste(i[1],i[2]))
      true_beta[c(non_0_main_idx,non_0_interac_idx,non_0_quadra_idx)] = 1
    }
    else if(heredity == "weak"){
      if(is.null(main_idx_sample)) non_0_main_idx = sort(sample(1:pmain,nb))
      tmp = combinations(n=pmain,r=2,repeats.allowed = T, v = 1:pmain)
      tmp = tmp[tmp[,1] %in% non_0_main_idx |  tmp[,2] %in% non_0_main_idx,]
      tmp = tmp[sample(nb*(nb+1)/2),]
      non_0_order2_idx = apply(tmp,1,function(i) paste(i[1],i[2]))
      true_beta[c(non_0_main_idx,non_0_order2_idx)] = 1
    }
  }
  if(is.null(noise_sd)){noise_sd = set_sigma_from_SNR(SNR = SNR, design = design, beta = true_beta)}
  else{SNR = compute_snr(noise_sd = noise_sd, design = design, beta = true_beta)}
  epsilon = rnorm(n = nrow(design), mean = 0, sd = noise_sd)
  Y = as.matrix(design)%*%true_beta + epsilon
  return(list(Y = Y,true_beta = true_beta, SNR = as.numeric(SNR), noise_sd = as.numeric(noise_sd)))
}
 
simu_quadraGaussian = function(n = 1000, pmain = 10, rho = 0, sigma = 1, heredity = "strong", 
                         SNR = 10, noise_sd = NULL, H_scale = F, ordinary_scale = F, nb = 3, 
                         true_beta = NULL, main_idx_sample = NULL){
  design_obj = simu_design_matrix(n = n, pmain = pmain, rho = rho, sigma = sigma, H_scale = H_scale, ordinary_scale = ordinary_scale)
  output_obj = simu_quadraGauss_output(design = design_obj$xtilde, pmain = design_obj$pmain, nb = nb, heredity = heredity, SNR=SNR, noise_sd = noise_sd, true_beta = true_beta, main_idx_sample = main_idx_sample)
  return(list(design_obj = design_obj, output_obj = output_obj))
} 
  
# 3. Intermediate functions 
make_interac_names = function(pmain){
  num_names = combinations(n = pmain, r = 2, v = 1:pmain, repeats.allowed = F)
  num_names = cbind(num_names,apply(num_names,1,function(x)paste(x[1],x[2])))
  return(num_names)
}

make_quadra_names = function(pmain){
  num_names = sapply(1:pmain, function(x) paste(x,x))
  return(num_names)
}
get_col_number = function(p){return((p*(p-1))/2 + 2*p)}
set_sigma_from_SNR = function(SNR,design,beta){
  return(sqrt(var(design%*%beta))/SNR)}

compute_snr = function(noise_sd,design,beta){
  return(var(design%*%beta)/(noise_sd)**2)
}

create_partition = function(p_train = 0.7, p_test = 0.15, p_valid = 0.15, n){
  splitSample = sample(1:3, size=n, prob=c(p_train,p_test,p_valid), replace = TRUE)
  return(splitSample)
}