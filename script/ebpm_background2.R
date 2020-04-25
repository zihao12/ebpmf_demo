
library(ebpm)

ebpm_background <- function(Y, s,
                            mu = NULL, Pi = NULL, grids = NULL,
                            maxiter = 10){
  if(is.null(grids)){
    al = 10^seq(-10,10)
    grids = list(al = al, bl = al)
  }

  L = length(grids$al)
  K = ncol(Y)

  if(is.null(mu)){
    mu = rowMeans(Y)/mean(s)
  }

  if(is.null(Pi)){
    Pi = matrix(runif(n = L * K), ncol = K)
    Pi = t(t(Pi)/colSums(Pi))
  }

  fit = ebpm_background_workhorse(Y = Y, s = s, mu = mu, Pi = Pi, grids = grids, maxiter = maxiter)
  return(fit)
}


## Y: J ,K
## s: K
## a, b : L
## mu : J
## Pi: L, K
ebpm_background_workhorse <- function(Y, s, mu, Pi, grids, maxiter){
  K = ncol(Y)
  J = nrow(Y)
  Y_rs = rowSums(Y)
  progress = c()
  pos.mean = matrix(, nrow = J, ncol = K)
  pos.mean_log = matrix(, nrow = J, ncol = K)
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
      pos.mean[,k] = fit$posterior$mean
      pos.mean_log[,k] = fit$posterior$mean_log
      s_v = s_v + s[k] * fit$posterior$mean
    }
    progress = c(progress, ll)
  }
  return(list(Pi = Pi, mu = mu,
              posterior = list(mean = pos.mean, mean_log = pos.mean_log),
              progress = progress))
}




update_mu <- function(Y, a, b, s, mu, Pi, maxiter){
  if(maxiter == 0){return(mu)}
}
