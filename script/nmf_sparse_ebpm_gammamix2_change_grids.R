## script for running `ebpm_gamma_mixture2` on dataset from https://zihao12.github.io/ebpmf_demo/nmf_sparse_data_prep
## here I change grids every D iterations

rm(list = ls())
set.seed(123)
devtools::load_all("../../ebpm")
library("ebpmf.alpha")
library(parallel)

## specify function used
maxiter = 5000
D = maxiter ### should be divided by maxiter
func_name = sprintf("ebpm_gamma_mixture2_change_grids_per_%d_iter", D)
ebpm_func = ebpm::ebpm_gamma_mixture2
pm_control = NULL

print(func_name)
ncpu = parallel::detectCores()
print(sprintf("%d cpus available", ncpu))

## setup for fitting
verbose = TRUE
tol = -1e+10

## load data
data = readRDS("../data/nmf_sparse_data.Rds")
X = data$X
k = ncol(data$L)
init = data$init

#browser()

## fit


start = proc.time()
M = as.integer(maxiter/D)
ELBO = c()
for(i in 1:M){
  fit <- ebpmf.alpha::ebpmf(X, K = k,pm_func = ebpm_func,
                            pm_control = pm_control, init = init,
														maxiter = D, verbose = verbose, tol = tol)
  ELBO = c(ELBO, fit$ELBO)
  ## need to remove previous g()
  qg = fit$qg
  qg$gls = replicate(k, list(NULL))
  qg$gfs = replicate(k, list(NULL))
  init = list(qg = qg)
}
t1 = proc.time() - start
fit[["runtime"]] = t1
fit[["ELBO"]] = ELBO

## save results
saveRDS(fit, sprintf("../data/nmf_sparse_%s.Rds", func_name))

print(sessionInfo())
