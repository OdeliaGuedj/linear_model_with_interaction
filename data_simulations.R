for(n in c(200,600,1000,10000)){
  for(pmain in c(50, 100, 150)){
    for(rho in c(0,0.2,0.5,0.8)){
      for(snr in c(1,10,100,1000)){
        name = paste0('data_n',n,"_pmain",pmain,"_rho",rho,"_snr",snr)
        data = simu_quadraGaussian(n = n,pmain = pmain,rho = rho,SNR = snr,H_scale = T,ordinary_scale = T,nb=pmain/2) 
        save(data, file = file.path(name))
        assign(name, get(load(name)))
        rm(data)
      }
    }
  }
}