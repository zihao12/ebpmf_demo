##
rm(list = ls())
source("nb_means.R")
set.seed(123)

K = 3
p = 999
eps = 1e-3
signal = 10
s = replicate(K, 1)
mu = runif(p) 
Phi = matrix(eps * runif(p * K), nrow = p, ncol = K)
Phi[1:(p/3),1] = 1 + signal * runif(p/3)
Phi[((p/3)+1):(2*p/3), 2] = 1 + signal * runif(p/3)
Phi[((2*p/3)+1):p,3] = 1 + signal * runif(p/3)
A = 1/Phi

## simulate data from the model
V = matrix(rgamma(n = p*K, shape = A, rate = A), nrow = p)
Theta = (mu %o% s) * V
Y = matrix(rpois(n = p * K, lambda = Theta), nrow = p)

fit = mle.nb_means.workhorse(Y, s, mu, A, maxiter = 100, verbose = TRUE,
                            control = list(gradient = TRUE, hessian = FALSE))
plot(fit$progress)