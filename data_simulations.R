library(parallel)
library(gtools)
num_cor = detectCores() -1
MC = 50

PAR.simu_quadraGaussian = function(i,n = 1000, pmain = 10, rho = 0, sigma = 1, heredity = "strong", 
                                   SNR = 10, noise_sd = NULL, H_scale = T, ordinary_scale = T, nb = pmain/2, 
                                   true_beta = NULL, main_idx_sample = NULL) {
  simu_quadraGaussian(n = n, pmain = pmain, rho = rho, sigma = sigma, heredity = heredity, 
                      SNR = SNR, noise_sd = noise_sd, H_scale = H_scale, ordinary_scale = ordinary_scale, nb = nb, 
                      true_beta = true_beta, main_idx_sample = main_idx_sample)}

for(n in c(200,600,1000,10000)){
  for(pmain in c(50,100,150)){
    for(rho in c(0,0.2,0.5,0.8)){
      for(snr in c(1,10,100,1000)){
        name = paste0('data_n',n,"_pmain",pmain,"_rho",rho,"_snr",snr)
        data = list()
        print(name)
        data = mclapply(1:MC, mc.cores = num_cor, PAR.simu_quadraGaussian,
                        n = n, pmain = pmain, rho = rho, SNR = snr, H_scale = T,ordinary_scale = T, nb = pmain/2)
        save(data, file = file.path(name))
        assign(name, get(load(name)))
      }
    }
  }
}