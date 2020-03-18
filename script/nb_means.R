## solve a nb_means problem
library(CVXR)

## model
## Y_jk ~ NB(1/phi_jk, 1/(1 + s_k mu_j phi_jk))

## algorithm
## Y_jk ~ Pois(s_k mu_j v_jk)
## v_jk ~ Gamma(a_jk, a_jk), a_jk = 1/phi_jk

## Input
## Y: count matrix of dimension [J, K]
## s: non-negative vector of length K
## mu: initiliazation of mean vector, length J
## A: initialization of 1/variance matrix, dimension [J, K]

## Output
## list(mu = mu*, A = A*)
mle.nb_means.workhorse <- function(Y, s, mu, A, maxiter = 10){
  J = nrow(Y)
  K = ncol(Y)

  for(iter in 1:maxiter){
    print(iter)
    #if(iter == 2){browser()}
    ## E-step
    V.pos = (A + Y)/(A + mu %o% s)
    V_log.pos = psigamma(A + Y) - log(A + mu %o% s)
    ## M-step
    ### update mu
    mu = rowSums(Y)/rowSums(t(s * t(V.pos))) ## TODO: fix the NAN issue
    ### update A
    for(k in 1:K){
      A[, k] = update_a(a = A[,k], c = V_log.pos[,k] - V.pos[,k])
    }
    ll = loglikelihood.nb_means(Y = Y, s = s, mu = mu, A = A)
    print(ll)
  }
  
  return(list(mu = mu, A = A))
}


## maximize J(a) = sum( c*a + a*log(a) - lgamma(a) )

update_a <- function(a, c){
  #browser()
  n = length(a)
  x <- Variable(n)
  x@value = a ## initialize
  constraints <- list(x >= 0)
  objective <- Ja.val(x, c)
  prob <- Problem(Minimize(objective), constraints)
  result <- solve(prob)
  a = result$getValue(x)
  return(a)
}

loglikelihood.nb_means <- function(Y, s, mu, A){
  
  Phi = 1 / A
  ll = lgamma(Y + A) - lgamma(A) - (Y  + A) * log(1 + (mu %o% s)*Phi ) + Y * (mu + log(Phi))
  return(sum(ll))
}


Ja.val <- function(a, c){
  browser()
  - sum( c*a + a*log(a) - log(gamma(a@value)) ) ## CVXR doesn't know lgamma?
}

# Ja.grad <- function(a, c){
#   - (  log(a) - psigamma(a, deriv = 0) + c + 1  )
# }

# Ja.hess <- function(a){
#   - ( 1/a - psigamma(a, deriv = 1) )
# }