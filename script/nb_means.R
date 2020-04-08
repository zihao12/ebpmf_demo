## solve a nb_means problem
library(stats)

## model
## Y_jk ~ NB(1/phi_jk, 1/(1 + s_k mu_j phi_jk))
## In the code I use A_jk = 1/Phi_jk instead

## algorithm
## Y_jk ~ Pois(s_k mu_j v_jk)
## v_jk ~ Gamma(a_jk, a_jk), a_jk = 1/phi_jk



nb_means <- function(Y, s, mu = NULL, A = NULL, maxiter = 10, verbose = FALSE,
                     control = list(method = "descend", gradient = TRUE, hessian = FALSE), seed = 123){
  ## initialization
  init = init_nb_means(Y = Y,s = s,mu = mu, A = A, seed = seed)
  mu = init$mu
  A = init$A
  ## compute MLE
  fit = mle.nb_means.workhorse(Y = Y, s = s, mu = mu, A  = A, maxiter = maxiter, verbose = verbose, control = control)
  ## compute posterior
  posterior = posterior_nb_means(Y = Y, s = s, g = fit$g)
  return( list(fitted_g = fit$g, posterior = posterior, log_likelihood = fit$log_likelihood, progress = fit$progress) )
}

init_nb_means <- function(Y, s, mu, A, seed = 123){
  set.seed(seed)
  if(is.null(mu)){
    mu = rowMeans(Y)/mean(s)
  }
  if(is.null(A)){
    J = nrow(Y)
    K = ncol(Y)
    A = matrix(runif(n = J * K), ncol = K)
  }
  return( list(mu = mu, A = A) )
}


posterior_nb_means <- function(Y, s, g){
  A = g$A
  mu = g$mu
  V.pos = (A + Y)/(A + mu %o% s)
  V_log.pos = psigamma(A + Y) - log(A + mu %o% s)
  lam_pm = mu * V.pos
  lam_log_pm =  log(mu) + V_log.pos
  posterior = list(mean = lam_pm, mean_log = lam_log_pm)
  return(posterior)
}

## Input
## Y: count matrix of dimension [J, K]
## s: non-negative vector of length K
## mu: initiliazation of mean vector, length J
## A: initialization of 1/variance matrix, dimension [J, K]

## Output
## list(mu = mu*, A = A*)


mle.nb_means.workhorse <- function(Y, s, mu, A, maxiter = 10, verbose = FALSE,
                                  control = list(method = "descend", gradient = TRUE, hessian = FALSE)){
  J = nrow(Y)
  K = ncol(Y)
  Y_rsum = rowSums(Y)
  mask_zero = Y_rsum == 0
  progress = c(loglikelihood.nb_means(Y = Y, s = s, mu = mu, A = A))
  if(verbose){print("iter 			loglik\n")}

  for(iter in 1:maxiter){
    ## E-step
    V.pos = (A + Y)/(A + mu %o% s)
    V_log.pos = psigamma(A + Y) - log(A + mu %o% s)
    ## M-step
    ### update mu
    mu = Y_rsum/rowSums(t(s * t(V.pos))) ## TODO: fix the NAN issue
    #mu = Y_rsum / rowSums(V.pos %*% diag(s))
    mu[mask_zero] = 0
    ### update A
    A = update_A(A, V_log.pos, V.pos, control = control)
    ll = loglikelihood.nb_means(Y = Y, s = s, mu = mu, A = A)
    progress <- c(progress, ll)
    if(verbose){print(sprintf("%d 		%f\n", iter, ll))}
  }
  g = list(mu = mu, A = A)
  log_likelihood = progress[iter]
  return(list(g = g, progress = progress, log_likelihood = log_likelihood))
}


## maximize J(a) = sum( c*a + a*log(a) - lgamma(a) )

update_A <- function(A, V_log.pos, V.pos, control){
  K = ncol(A)
  J = nrow(A)
  if(control$method == "descend"){
    for(k in 1:K){
      #if(iter == 153 && k == 6){browser()}
      tmp = try(update_a(a = A[,k], c = V_log.pos[,k] - V.pos[,k],
                         gradient = control$gradient, hessian = control$hessian))
      if(class(tmp) == "try-error"){
        A[, k] = update_a(a = A[,k], c = V_log.pos[,k] - V.pos[,k],
                          gradient = FALSE, hessian = FALSE)
      }
      else{A[, k] = tmp}
    }
  }

  if(control$method == "grid"){
    #browser()
    ## for each a_jk, I use a grid search. Parallelize over J
    ## TODO: parallelize over J*K elements (vectorize them...)
    grid_range = 10^seq(-10, 10)
    grids = t(matrix(replicate(J, grid_range), nrow = length(grid_range)))
    for(k in 1:K){
      # vals = Ja.val(a = grids, c = V_log.pos[,k] - V.pos[,k]) ## to check!!
      # idx = max.col(vals) ## the max value per row
      lls <- - Ja.val.element(a = grids, c = V_log.pos[,k] - V.pos[,k])
      #idx = max.col(lls) ## this function is buggy???
      idx = apply(lls, 1, which.max)
      A[,k] = grid_range[idx]
    }
  }
  return(A)
}


update_a <- function(a, c, gradient, hessian){
  fn_params = list(c = c, gradient = gradient, hessian = hessian)
  opt = do.call(nlm, c(list(obj.nb_means, a),fn_params, iterlim = 5, check.analyticals = FALSE))
  return(opt$estimate)
}

loglikelihood.nb_means <- function(Y, s, mu, A){
  Phi = 1 / A
  ll = lgamma(Y + A) - lgamma(A) - (Y  + A) * log(1 + (mu %o% s)*Phi ) + Y * (mu + log(Phi))
  return(sum(ll))
}

obj.nb_means <- function(a, c, gradient = TRUE, hessian =   FALSE){
  out = Ja.val(a, c)
  if(gradient){attr(out, 'gradient') <- Ja.grad(a, c)}
  if(hessian){attr(out, 'hessian') <- Ja.hess(a)}
  return(out)
}

Ja.val.element <- function(a, c){
  - ( c*a + a*log(a) - lgamma(a)  )
}

Ja.val <- function(a, c){
  - sum( c*a + a*log(a) - lgamma(a) )
}

Ja.grad <- function(a, c){
  - (  log(a) - psigamma(a, deriv = 0) + c + 1  )
}

Ja.hess <- function(a){
  diag(- ( 1/a - psigamma(a, deriv = 1) ))
}
