rm(list = ls())
source("nb_means.R")
set.seed(123)

exper_nb_means <- function(K = 3, p = 999, maxiter = 10, seed = 123, verbose =  FALSE){
  set.seed(seed)
  eps = 1e-3
  signal = 10
  s = replicate(K, 1)
  mu = 10 + 100*runif(p)

  Phi = runif(n = p*K) ## now Phi is just index
  mask_small = Phi < 0.8
  Phi[mask_small] = 1e-4
  id_tmp = rbinom(sum(!mask_small),1, 0.5)
  Phi[!mask_small] = id_tmp * 1e+4 + (1- id_tmp) * 1e+8
  Phi = matrix(Phi, nrow = p)
  A = 1/Phi

  ## simulate data from the model
  V = matrix(rgamma(n = p*K, shape = A, rate = A), nrow = p)
  Theta = (mu %o% s) * V
  Y = matrix(rpois(n = p * K, lambda = Theta), nrow = p)

  ll_oracle = loglikelihood.nb_means(Y, s, mu, A)
  #print(sprintf("ll_oracle : %f", ll_oracle))

  mu0 = runif(p)
  A0 = matrix(runif(p * K), nrow = p, ncol = K)
  runtime <- system.time(
    fit_init_random1 <- nb_means(Y, s, mu = NULL, A = NULL, maxiter = maxiter, verbose = verbose,
                              control = list(gradient = TRUE, hessian = FALSE), seed = 123)
  )
  fit_init_random1[["runtime"]] = runtime[[3]]


  runtime <- system.time(
    fit_init_random2 <- nb_means(Y, s, mu = NULL, A = NULL, maxiter = maxiter, verbose = verbose,
                              control = list(gradient = TRUE, hessian = FALSE), seed = 1234)
  )
  fit_init_random2[["runtime"]] = runtime[[3]]

  runtime <- system.time(
    fit_init_truth <- nb_means(Y, s, mu, A, maxiter = maxiter, verbose = verbose,
                              control = list(gradient = TRUE, hessian = FALSE))
  )
  fit_init_truth[["runtime"]] = runtime[[3]]

  return(list(Y = Y, Theta = Theta, mu = mu, Phi = Phi, ll_oracle = ll_oracle,
              fit_init_random1 = fit_init_random1, fit_init_random2 = fit_init_random2, fit_init_truth = fit_init_truth))
}

exper <- exper_nb_means(K = 10, p = 99, maxiter = 1000, verbose = TRUE)


