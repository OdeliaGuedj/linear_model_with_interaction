withdraw_proj = function(xtilde, pmain){
  p2 = pmain + (pmain) * (pmain+1)/2
  xtilde2 = xtilde
  for (j in ((pmain+1):(pmain + (pmain) * (pmain+1)/2))){
    model_j = lm(xtilde[,j]~as.matrix(xtilde[,1:pmain]))
    xtilde2[,j] = xtilde[,j] - predict(model_j)
  }
  return(xtilde2)
}


HierNetProj = function(x, y){
  x = as.matrix(x)
  pmain = ncol(x)
  y = as.numeric(y)
  interac_names = make_interac_names(pmain)
  z = as.matrix(generate_order2_from_main(x)$Z)
  xtilde = cbind(x, z)
  z_proj = withdraw_proj(xtilde,pmain)[,-c(1:pmain)]
  xtilde_proj = cbind(x, z_proj)
  lambda_min = hierNet::hierNet.cv(
    hierNet::hierNet.path(scale(x,T,T),y, zz = z_proj),scale(x,T,T),y, trace =0)$lamhat.1se
  fit = hierNet::hierNet(scale(x,T,T),y,lam=lambda_min,strong = T, zz = z_proj, trace = 0)
  theta = c(sapply(1:dim(interac_names)[1],function(i) fit$th[as.numeric(interac_names[i,1]),as.numeric(interac_names[i,2])]) ,
            diag(fit$th))
  gamma = fit$bp - fit$bn + solve(t(x)%*%x)%*%t(x) %*% z %*% theta
  phi = c(gamma,theta)
  names(phi) = c(1:pmain,interac_names[,3], make_quadra_names(pmain))
  return(list(lambda_min = lambda_min, fit = fit, phi = phi))
}
