####################################################
## Main functions for Poisson random effects example
## Plus data generation
####################################################

library(mcmcse)
library(foreach)
library(doParallel)

set.seed(8024248)
ni_s <- 5
I <- 50
c <- 10
sigma_eta <- 3
one_mat <- rep(1, I)
data <- matrix(0, nrow = I, ncol = ni_s)
identity_mat <- diag(1, I, I)
mu <- rnorm(1, 0, c)
eta_vec <- rnorm(I, mu, sigma_eta)
tol_nr <- 1e-8

# data generation

for (j in 1:ni_s)       
{
  data[,j] <- rpois(I, exp(eta_vec))
}
data

# log of true target i.e. log(p(eta, mu|y))

log_p <- function(eta, mu, data)   
{
  densval <- sum(eta^2)/(2*sigma_eta^2) - mu*sum(eta)/(sigma_eta^2) + 
                      ni_s*sum(exp(eta)) - sum(eta*apply(data, 1, sum)) +
                              (mu^2/2)*(I/(sigma_eta^2) + 1/(c^2))
  return(-densval)
}

log_plam <- function(eta, mu, data, lambda, y)   
{
  densval <- sum(eta^2)/(2*sigma_eta^2) - mu*sum(eta)/(sigma_eta^2) + 
    ni_s*sum(exp(eta)) - sum(eta*apply(data, 1, sum)) +
    (mu^2/2)*(I/(sigma_eta^2) + 1/(c^2)) + sum((y-c(eta, mu))^2)/(2*lambda)
  return(-densval)
}

grad_logp <- function(mu, sigma, eta)  # function evaluates gradient of log target 
{
  term_exp_eta <- (ni_s)*exp(eta)
  term_y <- apply(data, 1, sum)
  term_mu_sigma <- (mu/(sigma^2))
  grad_vec_eta <- - eta/(sigma^2) - term_exp_eta + term_y + term_mu_sigma
  grad_mu <- - mu*(I/(sigma^2) + 1/(c^2)) + sum(eta)/(sigma^2)
  grad_value <- c(grad_vec_eta, grad_mu)
  return(grad_value)
}

true_hessian <- function(sigma, eta)  # function evaluates hessian of eta's
{
  term_eta_diag <- 1/(sigma^2) + (ni_s)*exp(eta)
  return(-term_eta_diag)
}

proxfunc <- function(eta, mu, lambda, eta_initial, mu_initial, sigma)
{
  # For starting values
  iter <- 0
  eta_next <- eta_initial
  mu_next <- mu_initial

  grad_vec <- grad_logp(mu_next, sigma, eta_next) - (c(eta_next, mu_next) - c(eta, mu))/lambda
  while (sqrt(sum(grad_vec^2)) > tol_nr) 
  {
   # For eta's
    iter <- iter + 1

    # calculating hessian
    eta_hessian <- true_hessian(sigma, eta_next) - 1/lambda
    eta_grad <- grad_vec[1:I]    

    # defining the update 
    eta_next <- eta_next - eta_grad/eta_hessian

    mu.num <- (sum(eta_next)/sigma^2 + mu/lambda)
    mu.den <- I/(sigma^2) + 1/(c^2) + 1/lambda
    mu_next <- mu.num/mu.den
    grad_vec <- grad_logp(mu_next, sigma, eta_next) - (c(eta_next, mu_next) - c(eta, mu))/lambda
  }
  optima <- c(eta_next, mu_next)
  return(optima)
}

#gradient of log target^lambda

 grad_logplam <- function(eta, mu, lambda, eta_initial, mu_initial, sigma) 
  {
  term_prox <- proxfunc(eta, mu, lambda, eta_initial, mu_initial, sigma)
  term <- c(eta, mu)
  ans <-  (term-term_prox)/lambda
  return(-ans)
 }
 
 # pxbarker and mybarker proposal
 bark.prop <- function(eta, mu, lambda, eta_initial, mu_initial, sigma, delta)
 {
   aux_var <- rnorm(length(eta)+1, 0, 1)
   z <- sqrt(delta)*aux_var
   denom_prod <- z*grad_logplam(eta, mu, lambda, eta_initial, mu_initial, sigma)
   prob <- 1 / (1 + exp(- denom_prod))
   unifs <- runif(length(eta)+1)
   prop <- (c(eta,mu) + z)*(unifs <= prob) + (c(eta,mu) - z)*(unifs > prob)
   return(prop)
 }
 
 # true barker proposal
 bark.prop_true <- function(eta, mu, sigma, delta)
 {
   aux_var <- rnorm(length(eta)+1, 0, 1)
   z <- sqrt(delta)*aux_var
   denom_prod <- z*grad_logp(mu, sigma, eta)
   prob <- 1 / (1 + exp(- denom_prod))
   unifs <- runif(length(eta)+1)
   prop <- (c(eta,mu) + z)*(unifs <= prob) + (c(eta,mu) - z)*(unifs > prob)
   return(prop)
 }
 
 # barker log density
 log_bark.dens <- function(curr_point, prop_point, grad_curr_point, delta)
 {
   rw_dens <- sum(dnorm(prop_point-curr_point, 0, sqrt(delta), log = TRUE))
   exp_term <- - (grad_curr_point*(prop_point-curr_point))
   denom <- sum(log1p(exp(exp_term)))
   dens_val <- rw_dens - denom
   return(dens_val)
 }
 
 # to evaluate asymptotic covariance matrix
 asymp_covmat_fn <- function(chain, weights)
 {
   wts_mean <- mean(weights)
   num <- chain*weights
   sum_mat <- apply(num, 2, sum)
   is_est <- sum_mat / sum(weights)
   input_mat <- cbind(num, weights)  # input samples for mcse
   Sigma_mat <- mcse.multi(input_mat)$cov  # estimated covariance matrix of the tuple
   kappa_eta_mat <- cbind(diag(1/wts_mean, ncol(chain)), -is_est/wts_mean) # derivative of kappa matrix
   asymp_covmat <- (kappa_eta_mat %*% Sigma_mat) %*% t(kappa_eta_mat) # IS asymptotic variance
   return(asymp_covmat)
 }
 
##### MYMALA samples function

mymala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.mym <- matrix(0, nrow = iter, ncol = I+1)
  wts_is_est <- numeric(length = iter)
  
  # starting values
  samp_current <- c(eta_start, mu_start)
  samp.mym[1,] <- samp_current
  prox_val.curr <- proxfunc(eta_start, mu_start, lambda, 
                                eta_start, mu_start, sigma)
  targ_val.curr <- log_plam(prox_val.curr[1:I], prox_val.curr[I+1], data, lambda,
                             samp_current)
  
  # weights calculation
  psi_val <- - log_p(eta_start, mu_start, data)
  psi_lambda_val <- - targ_val.curr
  wts_is_est[1] <-  psi_lambda_val -  psi_val
  
   accept <- 0
  for (i in 2:iter) 
  {
    # proposal step
     prop.mean <- samp_current + 
           (delta / 2)*grad_logplam(samp_current[1:I], samp_current[I+1],
                                lambda, eta_start, mu_start, sigma)
     samp_next <- rnorm(I + 1, prop.mean, sqrt(delta))   
    
    # calculating prox values
    prox_val.next <- proxfunc(samp_next[1:I], samp_next[I+1],
                               lambda, eta_start, mu_start, sigma)
    
    targ_val.next <- log_plam(prox_val.next[1:I], prox_val.next[I+1],data, lambda,
                              samp_next)
    
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                       (delta / 2)*grad_logplam(samp_next[1:I], samp_next[I+1],
                         lambda, eta_start, mu_start, sigma),sqrt(delta), log = TRUE))
   
     q.curr_to_next <- sum(dnorm(samp_next, prop.mean,  sqrt(delta), log = TRUE))
    
    mh.ratio <- targ_val.next + q.next_to_curr - (targ_val.curr + q.curr_to_next)  # mh  ratio
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.mym[i,] <- samp_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      
      # weights
      psi_val <- - log_p(samp_next[1:I], samp_next[I+1], data)
      psi_lambda_val <- - targ_val.curr
      wts_is_est[i] <- psi_lambda_val - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.mym[i,] <- samp_current
      wts_is_est[i] <- wts_is_est[i-1]
    }
    samp_current <- samp.mym[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  object <- list(samp.mym, wts_is_est)
  return(object)
}

##### PxMALA samples function

px.mala <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.pxm <- matrix(0, nrow = iter, ncol = I+1)
  
  #  starting values
  samp_current <- c(eta_start, mu_start)
  samp.pxm[1,] <- samp_current
  U_sampcurr <- log_p(samp_current[1:I], samp_current[I+1], data)
  accept <- 0
  
  for (i in 2:iter)
  {
    # proposal step
    prop.mean <- samp_current + 
           (delta / 2)*grad_logplam(samp_current[1:I], samp_current[I+1],
                               lambda, eta_start, mu_start, sigma)
    samp_next <- rnorm(length(samp_current), prop.mean,  sqrt(delta))   
    
    U_sampnext <- log_p(samp_next[1:I], samp_next[I+1], data)
    
    q.next_to_curr <- sum(dnorm(samp_current, samp_next + 
                       (delta / 2)*grad_logplam(samp_next[1:I], samp_next[I+1],
                         lambda, eta_start, mu_start, sigma),sqrt(delta), log = TRUE))
    
    q.curr_to_next <- sum(dnorm(samp_next, prop.mean,  sqrt(delta), log = TRUE))
   
    mh.ratio <- U_sampnext + q.next_to_curr - (U_sampcurr + q.curr_to_next)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.pxm[i,] <- samp_next
      U_sampcurr <- U_sampnext
      accept <- accept + 1
    }
    else
    {
      samp.pxm[i,] <- samp_current
    }
    samp_current <- samp.pxm[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(accept/iter)
  return(samp.pxm)
}

#### MyBarker samples function  

mybarker <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.bark <- matrix(0, nrow = iter, ncol = I+1)
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  samp_current <- c(eta_start, mu_start)
  samp.bark[1,] <- samp_current
  prox_val.curr <- proxfunc(eta_start, mu_start, lambda, eta_start, mu_start, sigma)
  targ_val.curr <- log_plam(prox_val.curr[1:I], prox_val.curr[I+1], data, lambda,
                            c(eta_start, mu_start))
  
  # weights calculation
  psi_val <- - log_p(eta_start, mu_start, data)
  psi_lambda_val <- - targ_val.curr
  wts_is_est[1] <-  psi_lambda_val -  psi_val
  
  accept <- 0
  # For barker  
  for (i in 2:iter) 
  {
    # proposal step
    samp_next <- bark.prop(samp_current[1:I], samp_current[I+1], lambda, 
                                  eta_start, mu_start, sigma, delta)
    
    prox_val.next <- prox_val.next <- proxfunc(samp_next[1:I], samp_next[I+1],
                                               lambda, eta_start, mu_start, sigma)
    targ_val.next <- log_plam(prox_val.next[1:I], prox_val.next[I+1],data, lambda,
                              c(samp_next[1:I], samp_next[I+1]))
    
    grad_samp_curr <- grad_logplam(samp_current[1:I], samp_current[I+1], lambda,
                                   eta_start, mu_start, sigma)
    grad_samp_next <- grad_logplam(samp_next[1:I], samp_next[I+1], lambda,
                                   eta_start, mu_start, sigma)
    
    mh.ratio <- targ_val.next + log_bark.dens(samp_next, samp_current, grad_samp_next, delta) - 
                  targ_val.curr - log_bark.dens(samp_current, samp_next, grad_samp_curr, delta)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- samp_next
      prox_val.curr <- prox_val.next
      targ_val.curr <- targ_val.next
      
      # weights
      psi_val <- - log_p(samp_next[1:I], samp_next[I+1], data)
      psi_lambda_val <- - targ_val.curr
      wts_is_est[i] <- psi_lambda_val - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- samp_current
      wts_is_est[i] <- wts_is_est[i-1]
    }
    samp_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      # print(cat(i, j))
      # print(j)
      print(i)
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, wts_is_est, acc_rate)
  return(object)
}

## PxBarker samples

px.barker <- function(eta_start, mu_start, lambda, sigma, iter, delta, data)
{
  samp.bark <- matrix(0, nrow = iter, ncol = I+1)
  
  #  starting values
  samp_current <- c(eta_start, mu_start)
  samp.bark[1,] <- samp_current
  U_sampcurr <- log_p(samp_current[1:I], samp_current[I+1], data)
  accept <- 0
  
  # For barker
  accept <- 0
  for (i in 2:iter)
  {
    # proposal step
    samp_next <- bark.prop(samp_current[1:I], samp_current[I+1], lambda, 
                                        eta_start, mu_start, sigma, delta)
    
    U_sampnext <- log_p(samp_next[1:I], samp_next[I+1], data)
    
    grad_samp_curr <- grad_logplam(samp_current[1:I], samp_current[I+1], lambda,
                                   eta_start, mu_start, sigma)
    grad_samp_next <- grad_logplam(samp_next[1:I], samp_next[I+1], lambda,
                                   eta_start, mu_start, sigma)
    
    mh.ratio <- U_sampnext + log_bark.dens(samp_next, samp_current, grad_samp_next, delta) - 
                  U_sampcurr - log_bark.dens(samp_current, samp_next, grad_samp_curr, delta)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- samp_next
      U_sampcurr <- U_sampnext
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- samp_current
    }
    samp_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, acc_rate)
  return(object)
}

## Barker samples

barker <- function(eta_start, mu_start, sigma, iter, delta, data)
{
  samp.bark <- matrix(0, nrow = iter, ncol = I+1)
  
  #  starting values
  samp_current <- c(eta_start, mu_start)
  samp.bark[1,] <- samp_current
  U_sampcurr <- log_p(samp_current[1:I], samp_current[I+1], data)
  accept <- 0
  
  # For barker
  accept <- 0
  for (i in 2:iter)
  {
    # proposal step
    samp_next <- bark.prop_true(samp_current[1:I], samp_current[I+1], sigma, delta)
    
    U_sampnext <- log_p(samp_next[1:I], samp_next[I+1], data)
    
    grad_samp_curr <- grad_logp(samp_current[I+1], sigma, samp_current[1:I])
    grad_samp_next <- grad_logp(samp_next[I+1], sigma, samp_next[1:I])
    
    mh.ratio <- U_sampnext + log_bark.dens(samp_next, samp_current, grad_samp_next, delta) - 
      U_sampcurr - log_bark.dens(samp_current, samp_next, grad_samp_curr, delta)
    
    if(log(runif(1)) <= mh.ratio)
    {
      samp.bark[i,] <- samp_next
      U_sampcurr <- U_sampnext
      accept <- accept + 1
    }
    else
    {
      samp.bark[i,] <- samp_current
    }
    samp_current <- samp.bark[i,]
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.bark, acc_rate)
  return(object)
}

##  myhmc samples

myhmc <- function(eta_start, mu_start,lambda, sigma, iter, data, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = I+1)
  wts_is_est <- numeric(length = iter)
  
  # starting value computations
  samp <- c(eta_start, mu_start)
  proxval_curr <- proxfunc(eta_start, mu_start, lambda, eta_start, mu_start, sigma)
  samp.hmc[1,] <- samp
  
  # weights
  psi_val <- - log_p(eta_start, mu_start, data)
  psi_lambda_val <- - log_plam(proxval_curr[1:I], proxval_curr[I+1], data, lambda, 
                               c(eta_start, mu_start))
  wts_is_est[1] <- psi_lambda_val - psi_val
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*(I+1)), nrow = iter, ncol = I+1)
  accept <- 0
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logplam(samp[1:I], samp[I+1],
                          lambda, eta_start, mu_start, sigma)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logplam(samp[1:I], samp[I+1],
                            lambda, eta_start, mu_start, sigma)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    #  proximal values
    proxval_prop <-  proxfunc(samp[1:I], samp[I+1], lambda, eta_start, mu_start, sigma)
    
    U_curr <- - log_plam(proxval_curr[1:I], proxval_curr[I+1], data, lambda,
                         q_current)
    U_prop <- - log_plam(proxval_prop[1:I], proxval_prop[I+1], data, lambda, 
                         samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      proxval_curr <- proxval_prop
      
      # weights
      psi_val <- - log_p(samp[1:I], samp[I+1], data)
      wts_is_est[i] <- U_prop - psi_val
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      psi_val <- - log_p(q_current[1:I], q_current[I+1], data)
      wts_is_est[i] <- U_curr - psi_val
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, wts_is_est, acc_rate)
  return(object)
}

## pxhmc samples

pxhmc <- function(eta_start, mu_start,lambda, sigma, iter, data, eps_hmc, L)
{
  samp.hmc <- matrix(0, nrow = iter, ncol = I+1)
  
  # starting value computations
  samp <- c(eta_start, mu_start)
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*(I+1)), nrow = iter, ncol = I+1)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logplam(samp[1:I], samp[I+1],
                          lambda, eta_start, mu_start, sigma)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logplam(samp[1:I], samp[I+1],
                            lambda, eta_start, mu_start, sigma)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- - log_p(q_current[1:I], q_current[I+1], data)
    U_prop <- - log_p(samp[1:I], samp[I+1], data)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))
    }
  }
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}
