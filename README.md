# Linear Model with order 2 interactions
This repository contains the reproducible code of our paper *Fast Oracle inequality for penalized quadratic regression (O.Guedj, A.Guilloux and M-L.Taupin)* whose preprint will be soon availaible on arxiv.

The files:
- simu_functions.R contains functions needed to simulate a Gaussian model allowing order 2 interactions and preserving the strong or weak hierarchy principle.
- eval_algo_performance_functions.R contains functions computing indicators of the performance of an algorithm (both in prediction and selection) such as the AUC, the RMSE, ...
- data_simulations.R contains a parallelized version of the main simulation function: simu_quadraGaussian(), and the code used to simulate all of our settings. Be aware that the function _mcapply()_ was used and it does'nt not work on Windows !
- HierarchicalDescalingChenetal2020.R contains the implementation of the method of Chen et al (2020): *An Easy-to-Implement Hierarchical Standardization for Variable Selection Under Strong Heredity Constraint. Journal of statistical theory and practice, 14(3), 1-32*
- compare_algo_functions.R is the file containing a function comparing several algorithms performing variable selection and prediction in a quadratic linear regression.
- HierNerProj.source.R contains the source code of our algorithm _HierNetProj_.

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
