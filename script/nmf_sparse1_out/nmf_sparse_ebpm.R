## script for running `ebpm` on dataset from https://zihao12.github.io/ebpmf_demo/nmf_sparse_data_prep

rm(list = ls())
set.seed(123)
devtools::load_all("../../ebpm")
library("ebpmf.alpha")
library(parallel)

## specify function used
args = commandArgs(trailingOnly=TRUE)
func_name = args[1]
#func_name = "ebpm_two_gamma_fast5"

pm_control = NULL
if(func_name == "ebpm_two_gamma"){
	ebpm_func = ebpm::ebpm_two_gamma
	pm_control = list(n_iter = 10)
}
if(func_name == "ebpm_two_gamma_fast5"){
	ebpm_func = ebpm::ebpm_two_gamma_fast5
	pm_control = list(n_iter = 10)
}
if(func_name == "ebpm_gamma_mixture2"){ebpm_func = ebpm::ebpm_gamma_mixture2}
if(func_name == "ebpm_point_gamma"){ebpm_func = ebpm::ebpm_point_gamma}


print(func_name)
ncpu = parallel::detectCores()
print(sprintf("%d cpus available", ncpu))

## setup for fitting
verbose = TRUE
maxiter = 1000
tol = -1e+10

## load data
data = readRDS("../data/nmf_sparse_data.Rds")
X = data$X
k = ncol(data$L)
init = data$init

#browser()

## fit
t1 <- system.time(
  fit <- ebpmf.alpha::ebpmf(X, K = k,pm_func = ebpm_func,
                            pm_control = pm_control, init = init,
														maxiter = maxiter, verbose = verbose, tol = tol)
)
fit[["runtime"]] = t1

## save results
saveRDS(fit, sprintf("../data/nmf_sparse_%s.Rds", func_name))

print(sessionInfo())
