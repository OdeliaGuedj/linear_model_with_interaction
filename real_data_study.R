# 1. Load/shape the data and remove NA
source('~/Dropbox/Thèse/Thèse Odelia/Code/simu_functions.R', echo=TRUE)
source('~/Dropbox/Thèse/Thèse Odelia/Code/eval_algo_performance_functions.R', echo=T)
source('~/Dropbox/Thèse/Thèse Odelia/Code/withdraw_proj.R', echo=TRUE)

###############################################
## DO NOT RUN THE MATRIX ARE ALREADY STRORED ##
###############################################

data = read.table("communities.data",sep =",", na.strings = "?")
x = data[,6:127]
na_count = cbind(1:ncol(x), sapply(x, function(j) sum(length(which(is.na(j))))))
na_count
x = x[,-c(97:113,117:120,122)]
which(is.na(x$V31))
x[131,"V31"] = mean(x$V31, na.rm = T)
pmain = ncol(x)
n = nrow(x)
var_names = cbind(colnames(x),1:pmain)
colnames(x) = 1:pmain
xtilde = cbind(x,generate_order2_from_main(x)$Z)
y = data[,128]
xtilde_proj = withdraw_proj(xtilde = xtilde, pmain = pmain)

write.csv(x,"x.csv",row.names = F)
write.csv(y, "y.csv", row.names = F)
write.csv(xtilde,"xtilde.csv", row.names = F)
write.csv(xtilde_proj, "xtilde_proj.csv", row.names = F)

# 2. Data partition

x = read.csv("x.csv")
y = read.csv("y.csv")
xtilde = read.csv("xtilde.csv")
xtilde_proj = read.csv("xtilde_proj.csv")

pmain = ncol(x)
n = nrow(x)
interac_names = make_interac_names(pmain)
quadra_names = make_quadra_names(pmain)
colnames(x) = 1:pmain
colnames(y) = NULL
colnames(xtilde) = c(1:pmain,interac_names[,3],quadra_names)
colnames(xtilde_proj) = colnames(xtilde)

trainIndex = caret::createDataPartition(1:n, p=0.8, list=FALSE)
x_train = as.data.frame(xtilde[trainIndex,1:pmain]); x_test = as.data.frame(xtilde[-trainIndex,1:pmain])
xtilde_train = as.data.frame(xtilde[trainIndex,]); xtilde_test = as.data.frame(xtilde[-trainIndex,])
xtilde_proj_train = as.data.frame(xtilde_proj[trainIndex,]); xtilde_proj_test = as.data.frame(xtilde_proj[-trainIndex,])
y_train = y[trainIndex]; y_test = y[-trainIndex]

# 3. Using algorithms

## 3.1 All pairs LASSO

start_time_apl = Sys.time()
lambda_min_allPairsLASSO = glmnet::cv.glmnet(as.matrix(xtilde_train),y_train,alpha=1,standardize =T)$lambda.min
all_pair_lasso_fit = glmnet::glmnet(as.matrix(xtilde_train),y_train,alpha=1,standardize =T, lambda = lambda_min_allPairsLASSO)
all_pair_lasso_beta = as.vector(all_pair_lasso_fit$beta)
names(all_pair_lasso_beta) = colnames(xtilde)
all_pair_lasso_predict = predict(all_pair_lasso_fit, s = lambda_min_allPairsLASSO, newx = as.matrix(xtilde_test))
end_time_apl = Sys.time()

real_data_APL =  list("lambdaMin" = lambda_min_allPairsLASSO, 
                 "fit" = all_pair_lasso_fit, 
                 "beta" = all_pair_lasso_beta, 
                 "y_hat" = all_pair_lasso_predict, 
                 "time" = as.numeric(end_time_apl-start_time_apl), 
                 "rmse" = compute_RMSE(y_test, all_pair_lasso_predict),
                 "main" = length(which(all_pair_lasso_beta[1:pmain]!=0)),
                 "order2" = length(which(all_pair_lasso_beta[-c(1:pmain)]!=0)),
                 "model_size" = compute_model_size(fitted_beta = all_pair_lasso_beta)
                 )
save(real_data_APL,file = "real_data_APL.RData")

## 3.2 HierNet

start_time_hiernet = Sys.time()
lambda_min_HierNet = hierNet::hierNet.cv(hierNet::hierNet.path(scale(x_train,T,T),y_train),scale(x_train,T,T),y_train, trace =0)$lamhat.1se
hiernet_fit = hierNet::hierNet(scale(x_train,T,T),y_train,lam=lambda_min_HierNet,strong = T, trace = 0)
hiernet_beta = c(hiernet_fit$bp - hiernet_fit$bn,
                 sapply(1:dim(interac_names)[1],function(i) hiernet_fit$th[as.numeric(interac_names[i,1]),as.numeric(interac_names[i,2])]) ,
                 diag(hiernet_fit$th))
names(hiernet_beta) = c(1:pmain, interac_names[,3], quadra_names)
hiernet_predict = predict(hiernet_fit,scale(x_test,T,T))
end_time_hiernet = Sys.time()

real_data_HierNet =  list("lambdaMin" = lambda_min_HierNet, 
                      "fit" = hiernet_fit, 
                      "beta" = hiernet_beta, 
                      "y_hat" = hiernet_predict, 
                      "time" = as.numeric(end_time_hiernet-start_time_hiernet), 
                      "rmse" = compute_RMSE(y_test, hiernet_predict),
                      "main" = length(which(hiernet_beta[1:pmain]!=0)),
                      "order2" = length(which(hiernet_beta[-c(1:pmain)]!=0)),
                      "model_size" = compute_model_size(fitted_beta = hiernet_beta)
)
save(real_data_HierNet,file = "real_data_HierNet.RData")

## 3.3 HierNetProj

start_time_hiernetproj = Sys.time()
lambda_min_HierNetProj = hierNet::hierNet.cv(
  hierNet::hierNet.path(scale(x_train,T,T),y_train, 
                        zz = as.matrix(xtilde_proj_train[,-c(1:pmain)]))
  ,scale(x_train,T,T),y_train, trace =0)$lamhat.1se
HierNetProj_fit = hierNet::hierNet(scale(x_train,T,T),y_train,lam=lambda_min_HierNetProj,strong = T,
                                   zz = as.matrix(xtilde_proj_train[,-c(1:pmain)]), trace = 0)
theta = c(sapply(1:dim(interac_names)[1],function(i) HierNetProj_fit$th[as.numeric(interac_names[i,1]),as.numeric(interac_names[i,2])]) ,
  diag(HierNetProj_fit$th))
gamma = HierNetProj_fit$bp - HierNetProj_fit$bn + 
  solve(t(as.matrix(x_train))%*%as.matrix(x_train))%*%t(as.matrix(x_train)) %*% as.matrix(xtilde_train[,-c(1:pmain)]) %*% theta
HierNetProj_beta = c(gamma, theta)
names(HierNetProj_beta) = c(1:pmain, interac_names[,3], quadra_names)
HierNetProj_predict = predict(HierNetProj_fit,scale(x_test,T,T))
end_time_hiernetproj = Sys.time()

mean((hiernet_predict - HierNetProj_predict)^2)

real_data_HierNetProj =  list("lambdaMin" = lambda_min_HierNetProj, 
                          "fit" = HierNetProj_fit, 
                          "beta" = HierNetProj_beta, 
                          "y_hat" = HierNetProj_predict, 
                          "time" = as.numeric(end_time_hiernetproj-start_time_hiernetproj), 
                          "rmse" = compute_RMSE(y_test, HierNetProj_predict),
                          "main" = length(which(HierNetProj_beta[1:pmain]!=0)),
                          "order2" = length(which(HierNetProj_beta[-c(1:pmain)]!=0)),
                          "model_size" = compute_model_size(fitted_beta = HierNetProj_beta)
)
save(real_data_HierNetProj, file = "real_data_HierNetProj.RData")

# 4. Results

tab_results_real_data = rbind(round(cbind(real_data_APL$rmse,real_data_HierNet$rmse,real_data_HierNetProj$rmse),2),
                              round(cbind(real_data_APL$time,real_data_HierNet$time,real_data_HierNetProj$time),2),
                              as.character(cbind(real_data_APL$main,real_data_HierNet$main,real_data_HierNetProj$main)),
                              as.character(cbind(real_data_APL$order2,real_data_HierNet$order2,real_data_HierNetProj$order2))
)
row.names(tab_results_real_data) = c("RMSE", "Time (sec)", "# main", "# order 2")
colnames(tab_results_real_data) = c("APL", "HierNet", "HierNetProj")

save(tab_results_real_data, file = "tab_results_real_data.RData")
kable(tab_results_real_data, "latex")
