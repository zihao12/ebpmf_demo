## solve a nb_means problem


## model
## Y_jk ~ NB(1/phi_jk, 1/(1 + s_k mu_j phi_jk))

## Input
## Y: count matrix of dimension [J, K]
## s: non-negative vector of length K
## mu: initiliazation of mean vector, length J
## phi: initialization of variance matrix, dimension [J, K]

## Output
## list(mu = mu*, phi = phi*)
mle.nb_means.workhorse <- function(Y, s, mu, Phi, maxiter = 10){
  J = nrow(Y)
  K = ncol(Y)

  for(iter in 1:maxiter){
    mu = update_mu(Y = Y, s = s, mu = mu, Phi = Phi)
    for(k in 1:K){
      Phi[, k] = update_phi(y = Y[,k], s = s[k], mu = mu, phi = Phi[,k])
    }
  }

  return(list(mu = mu, Phi = Phi))
}

update_mu <- function(Y,s, mu, Phi){
  ## TODO
}

update_phi <- function(y, s, mu, phi){
  ## TODO
}