##
rm(list = ls())
setwd("~/Desktop/git/ebpmf_demo/script/")
source("ebpm_background2.R")
source("nb_means.R")

set.seed(123)

K = 50
p = 999
eps = 1e-3
signal = 10
s = replicate(K, 1)
# mu = runif(p)
# Phi = matrix(eps * runif(p * K), nrow = p, ncol = K)
# Phi[1:(p/3),1] = 1 + signal * runif(p/3)
# Phi[((p/3)+1):(2*p/3), 2] = 1 + signal * runif(p/3)
# Phi[((2*p/3)+1):p,3] = 1 + signal * runif(p/3)
# A = 1/Phi

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

# fit1 = mle.nb_means.workhorse(Y, s, mu, A, maxiter = 100, verbose = TRUE,
#                             control = list(method = "descend",gradient = TRUE, hessian = FALSE))

## Y, s, a, b, mu, Pi

a = 10^seq(-10,10)
b = a
Pi = matrix(runif(n = K * length(a)), ncol = K)
Pi = t(t(Pi) / colSums(Pi))

fit = mle.nb_means.workhorse(Y, s, mu, A, maxiter = 100, verbose = TRUE,
                              control = list(method = "grid",gradient = TRUE, hessian = FALSE))
fit2 = ebpm_background_workhorse(Y = Y, s = s, a = a, b = b, mu = mu, Pi = Pi)
plot(fit2$progress)
#points(fit1$progress, col = "red")
