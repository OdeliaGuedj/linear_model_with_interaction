# Linear Model with order 2 interactions
This repository contains the reproducible code of our paper *Oracle inequality for the Hierarchical LASSO (O.Guedj, A.Guilloux and M-L.Taupin)* whose preprint will be soon availaible on arxiv.

The files:
- simu_functions.R contains 12 functions needed to simulate a Gaussian model allowing oreder interactions and preserving the strong or weak hierarchy principle.
- eval_algo_performance_functions.R contains 11 functions including 8 computing indicators of the performance of an algorithm (both in prediction and selection) such as the AUC, the RMSE, ...

## Data simulation: an example

```r
source('./simu_functions.R', echo=TRUE)
n = 1000
pmain = 20
heredity = "strong"
data = simu_quadraGaussian(n=n,pmain=pmain,heredity=heredity)
xtilde = data$design_obj$xtilde
y = data$output_obj$Y
```
