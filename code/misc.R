sample_point_gamma_one <- function(point_gamma_){
  if(rbinom(1,1, point_gamma_$pi0) == 1){
    return(0)
  }else{
    return(rgamma(1,shape = point_gamma_$shape, scale = point_gamma_$scale))
  }
}

sample_point_gamma  <- function(n, point_gamma_, seed = 123){
  set.seed(seed)
  out = replicate(n, sample_point_gamma_one(point_gamma_))
  return(out)
}

## simulate a poisson mean problem from mixture of exponential
sample_expmix  <-  function(n,gammamix_, seed = 123){
  set.seed(seed)
  a = gammamix_$shape
  b = 1/gammamix_$scale
  pi = gammamix_$pi
  lam = replicate(n, sim_mgamma(a, b, pi))
  return(lam)
}
sim_mgamma <- function(a,b,pi){
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}

# Scale each column of A so that the entries in each column sum to 1;
# i.e., colSums(scale.cols(A)) should return a vector of ones.
scale.cols <- function (A)
  apply(A,2,function (x) x/sum(x))

# Convert the parameters (factors & loadings) for the Poisson model to
# the factors and loadings for the multinomial model. The return value
# "s" gives the Poisson rates for generating the "document" sizes.
poisson2multinom <- function (F, L) {
  L <- t(t(L) * colSums(F))
  s <- rowSums(L)
  L <- L / s
  F <- scale.cols(F)
  return(list(F = F,L = L,s = s))
}


KL <- function(true,est){
  sum(ifelse(true==0,0,true * log(true/est)) + est - true)
}

JS  <- function(true,est){
  0.5*(KL(true, est) + KL(est, true))
}

RMSE <- function(true, est){
  sqrt(mean((true - est)^2))
}

log_lik <- function(X, lam){
  return(sum(dpois(x = X, lambda = lam , log = T)))
}
