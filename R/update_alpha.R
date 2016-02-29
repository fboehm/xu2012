#' Update the alpha vector, the vector of sample-specific random effects
#'
#' @param alpha vector of sample-specific random effects
#' @param y data matrix
#' @param trinary trinary indicator matrix
#' @export
update_alpha <- function(alpha, y, trinary){
  tmax <- length(alpha)
  imax <- nrow(y)
  eps <- rnorm(n = tmax, mean = 0, sd = 0.1) # is 0.1 a reasonable sd??
  eps[tmax] <- - sum(eps[1:(tmax - 1)])
  alpha_prop <- alpha + eps
  logprior_ratio <- sum(dnorm(alpha_prop, mean = 0, sd = 0.1, log = TRUE) -
    dnorm(alpha, mean = 0, sd = 0.1, log = TRUE))
  # to calculate the log lik ratio, we need to
  # consider only those y's for which the corresponding
  # trinary indicator has value zero...
  # because, then, the lik is not equal to one (it's normal)
  loglik_mat_alpha <- matrix(data = 0, nrow = imax, ncol = tmax)
  loglik_mat_alpha_prop <- matrix(data = 0, nrow = imax, ncol = tmax)
  # all entries of the matrices are 0, initially.
  # we then update those entries that have trinary indicator
  # equal to zero
  for (i in 1:imax){
    for (t in 1:tmax){
      if (trinary[i, t] == 0)
      # subtract mu[i] and alpha[t]
      loglik_mat_alpha[i, t] <- dnorm(y[i, t] - mu[i] - alpha[t],
                                      mean = 0, sd = sigma[i],
                                      log = TRUE)
      loglik_mat_alpha_prop[i, t] <- dnorm(y[i, t] - mu[i] - alpha_prop[t],
                                      mean = 0, sd = sigma[i],
                                      log = TRUE)
    }
  }
  # take sum of entries in each matrix
  loglik_ratio <- sum(loglik_mat_alpha_prop) - sum(loglik_mat_alpha)
  acc_ratio <- exp(loglik_ratio + logprior_ratio)
  u <- runif(n = 1, min = 0, max = 1)
  if (u < acc_ratio){
    out <- alpha_prop
  }
  else {
    out <- alpha
  }
  return(out)
}
