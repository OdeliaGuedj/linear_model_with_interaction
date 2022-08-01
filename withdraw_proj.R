withdraw_proj = function(xtilde, pmain){
  p2 = pmain + (pmain) * (pmain+1)/2
  xtilde2 = xtilde
  for (j in ((pmain+1):(pmain + (pmain) * (pmain+1)/2))){
    model_j = lm(xtilde[,j]~as.matrix(xtilde[,1:pmain]))
    xtilde2[,j] = xtilde[,j] - predict(model_j)
  }
  return(xtilde2)
}