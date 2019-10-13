## ------------------------------------------------------------------------
library(ebpmf)
library(gtools)
library(NNLM)


## ------------------------------------------------------------------------
sim_mgamma <- function(dist){
  pi = dist$pi
  a = dist$a
  b = dist$b
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}

## simulate a poisson mean problem
## to do:
simulate_pm  <-  function(n, p, dl, df, K,scale_b = 10, seed = 123){
  set.seed(seed)
  ## simulate L
  a = replicate(dl,1)
  b = 10*runif(dl)
  pi <- rdirichlet(1,rep(1/dl, dl))
  gl = list(pi = pi, a = a, b= b)
  L = matrix(replicate(n*K, sim_mgamma(gl)), ncol = K)
  ## simulate F
  a = replicate(df,1)
  b = 10*runif(df)
  pi <- rdirichlet(1,rep(1/df, df))
  gf = list(pi = pi, a = a, b= b)
  F = matrix(replicate(p*K, sim_mgamma(gf)), ncol = K)
  ## simulate X
  lam = L %*% t(F)
  X = matrix(rpois(n*p, lam), nrow = n)
  Y = matrix(rpois(n*p, lam), nrow = n)
  ## prepare output
  g = list(gl = gl, gf = gf)
  out = list(X = X, Y = Y, L = L, F = F, g = g)
  return(out)
}

n = 50
p = 100
K = 2
dl = 10
df = 10
scale_b = 5
sim = simulate_pm(n, p, dl, df, K, scale_b = scale_b, seed =12)



## ----warning  = F, message=F, results='hide'-----------------------------
methods = c(); runtimes = c(); ll_trains = c(); ll_vals = c(); RMSEs = c()
## ebpmf_exponential_mixture
start = proc.time()
out = ebpmf::ebpmf_exponential_mixture(sim$X, K, m = 2, maxiter.out = 100)
runtime = proc.time() - start
lam_fit = out$qg$qls_mean %*% t(out$qg$qfs_mean)
ll_train = sum(dpois(sim$X, lambda = lam_fit, log = T))
ll_val   = sum(dpois(sim$Y, lambda = lam_fit, log = T))
RMSE     = mean((lam_fit - (sim$L %*% t(sim$F)))^2) 

methods = c(methods, "ebpmf_exponential_mixture")
runtimes = c(runtimes, runtime[[3]])
ll_trains = c(ll_trains, ll_train)
ll_vals   = c(ll_vals, ll_val)
RMSEs = c(RMSEs, RMSE)

## nnmf
W0 = out$qg$qls_mean
H0 = t(out$qg$qfs_mean)
start = proc.time()
out = NNLM::nnmf(sim$X, K,init = list(W0 = W0, H0 = H0), loss = "mkl", method = "lee", max.iter = 100, rel.tol = -1)
runtime = proc.time() - start
lam_fit = out$W %*% out$H
ll_train = sum(dpois(sim$X, lambda = lam_fit, log = T))
ll_val   = sum(dpois(sim$Y, lambda = lam_fit, log = T))
RMSE     = mean((lam_fit - (sim$L %*% t(sim$F)))^2) 

methods = c(methods, "NNMF")
runtimes = c(runtimes, runtime[[3]])
ll_trains = c(ll_trains, ll_train)
ll_vals   = c(ll_vals, ll_val)
RMSEs = c(RMSEs, RMSE)

## ebpmf_point_gamma
start = proc.time()
out = ebpmf::ebpmf_point_gamma(sim$X, K,maxiter.out = 100)
runtime = proc.time() - start
lam_fit = out$qg$qls_mean %*% t(out$qg$qfs_mean)
ll_train = sum(dpois(sim$X, lambda = lam_fit, log = T))
ll_val   = sum(dpois(sim$Y, lambda = lam_fit, log = T))
RMSE     = mean((lam_fit - (sim$L %*% t(sim$F)))^2) 
methods = c(methods, "ebpmf_point_gamma")
runtimes = c(runtimes, runtime[[3]])
ll_trains = c(ll_trains, ll_train)
ll_vals   = c(ll_vals, ll_val)
RMSEs = c(RMSEs, RMSE)

df <- data.frame(method = methods, runtime = runtimes, ll_train = ll_trains, ll_val = ll_vals, RMSE = RMSEs)


## ------------------------------------------------------------------------
df


## ------------------------------------------------------------------------
X = read.csv("../10xgenomics/cd14_monocytes/filtered_matrices_mex/hg19/Y.csv")
Y = read.csv("../10xgenomics/cd14_monocytes/filtered_matrices_mex/hg19/Yhat.csv")
real = list(X = as.matrix(X), Y = as.matrix(Y))
K = 2
maxiter.out = 20


## ----warning  = F, message=F, results='hide'-----------------------------
methods = c(); runtimes = c(); ll_trains = c(); ll_vals = c(); 
## ebpmf_exponential_mixture
start = proc.time()
out = ebpmf::ebpmf_exponential_mixture(real$X, K, m = 2, maxiter.out = maxiter.out)
runtime = proc.time() - start
lam_fit = out$qg$qls_mean %*% t(out$qg$qfs_mean)
ll_train = sum(dpois(real$X, lambda = lam_fit, log = T))
ll_val   = sum(dpois(real$Y, lambda = lam_fit, log = T))

methods = c(methods, "ebpmf_exponential_mixture")
runtimes = c(runtimes, runtime[[3]])
ll_trains = c(ll_trains, ll_train)
ll_vals   = c(ll_vals, ll_val)

## nnmf
W0 = out$qg$qls_mean
H0 = t(out$qg$qfs_mean)
start = proc.time()
out = NNLM::nnmf(real$X, K,init = list(W0 = W0, H0 = H0), loss = "mkl", method = "lee", max.iter = maxiter.out, rel.tol = -1)
runtime = proc.time() - start
lam_fit = out$W %*% out$H
ll_train = sum(dpois(real$X, lambda = lam_fit, log = T))
ll_val   = sum(dpois(real$Y, lambda = lam_fit, log = T))

methods = c(methods, "NNMF")
runtimes = c(runtimes, runtime[[3]])
ll_trains = c(ll_trains, ll_train)
ll_vals   = c(ll_vals, ll_val)

## ebpmf_point_gamma
start = proc.time()
out = ebpmf::ebpmf_point_gamma(real$X, K,maxiter.out = maxiter.out)
runtime = proc.time() - start
lam_fit = out$qg$qls_mean %*% t(out$qg$qfs_mean)
ll_train = sum(dpois(real$X, lambda = lam_fit, log = T))
ll_val   = sum(dpois(real$Y, lambda = lam_fit, log = T))
methods = c(methods, "ebpmf_point_gamma")
runtimes = c(runtimes, runtime[[3]])
ll_trains = c(ll_trains, ll_train)
ll_vals   = c(ll_vals, ll_val)

df <- data.frame(method = methods, runtime = runtimes, ll_train = ll_trains, ll_val = ll_vals, RMSE = RMSEs)


## ------------------------------------------------------------------------
df

