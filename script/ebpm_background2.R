
library(ebpm)

## Y: J ,K
## s: K
## a, b : L
## mu : J
## Pi: L, K
ebpm_background_workhorse <- function(Y, s, a, b, mu, Pi,
                                      grids = list(al = 10^seq(-10,10), bl = 10^seq(-10,10)),
                                      maxiter = 10, inner_iter = 0){
  K = ncol(Y)
  J = nrow(Y)
  Y_rs = rowSums(Y)
  progress = c()
  for(iter in 1: maxiter){
    ## update mu
    #mu <- update_mu(Y = Y, a = a, b = b, s = s, mu = mu, Pi = Pi, maxiter = inner_iter)
    if(iter == 1){mu <- mu}
    else{mu <- Y_rs/s_v}

    ## update Pi
    ll = 0
    s_v = replicate(J, 0)
    for(k in 1:K){
      g_init = gammamix(pi = Pi[,k], shape = grids$al, scale = 1/grids$bl)
      fit = ebpm_gamma_mixture(x = Y[,k], s = s[k] * mu, shape = grids$al, scale = 1/grids$bl, g_init = g_init)
      Pi[,k] = fit$fitted_g$pi
      ll = ll + fit$log_likelihood
      s_v = s_v + s[k] * fit$posterior$mean
    }
    progress = c(progress, ll)
  }
  return(list(Pi = Pi, mu = mu, progress = progress))
}




update_mu <- function(Y, a, b, s, mu, Pi, maxiter){
  if(maxiter == 0){return(mu)}
}
