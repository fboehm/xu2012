#' Update the mean vector mu (after Xu et al., 2012)
#'
#' @param mu a mean vector
#' @param trinary a trinary matrix
#' @param tau prior sd of mu
#' @param sd_jump sd of jump distribution
#'
#' @export
update_mu <- function(mu, trinary, y, alpha, tau = 1, sd_jump = 0.1){
  imax <- nrow(trinary)
  tmax <- ncol(trinary)
  # make a proposal beta
  eps <- rnorm(n = imax, mean = 0, sd = sd_jump)
  mu_pro <- mu + eps
  # eval posterior at each of 1. mu and 2. mu_pro
  ## prior ratio calcs
  logprior_ratio <- sum(dnorm(mu_pro,
                       mean = 0, sd = tau,
                       log = TRUE)
  - dnorm(mu_pro, mean = 0,
          sd = tau, log = TRUE))
  ####################################
  ## lik ratio calcs
  loglik_mat_mu <- matrix(data = 0, nrow = imax, ncol = tmax)
  loglik_mat_mu_pro <- matrix(data = 0, nrow = imax, ncol = tmax)
  # all entries of the matrices are 0, initially.
  # we then update those entries that have trinary indicator
  # equal to zero
  for (i in 1:imax){
    for (t in 1:tmax){
      if (trinary[i, t] == 0)
        # subtract mu[i] and alpha[t]
        loglik_mat_mu[i, t] <- dnorm(y[i, t] - mu[i] - alpha[t],
                                        mean = 0, sd = sigma[i],
                                        log = TRUE)
      loglik_mat_mu_pro[i, t] <- dnorm(y[i, t] - mu_pro[i] - alpha[t],
                                           mean = 0, sd = sigma[i],
                                           log = TRUE)
    }
  }
  # take sum of entries in each matrix
  loglik_ratio <- sum(loglik_mat_mu_pro) - sum(loglik_mat_mu)
  # add logprior_ratio and loglik_ratio; exponentiate
  acc_ratio <- exp(loglik_ratio + logprior_ratio)
  u <- runif(n = 1, min = 0, max = 1)
  if (u < acc_ratio){
    out <- mu_pro
  }
  else {
    out <- mu
  }
  return(out)
}
