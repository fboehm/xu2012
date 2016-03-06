#' Update the alpha vector, the vector of sample-specific random effects
#'
#' @param alpha vector of sample-specific random effects
#' @param y data matrix
#' @param trinary trinary indicator matrix
#' @param mu mean vector
#' @param sigma standard deviations vector
#' @param kminus kminus vector for uniform distributions
#' @param kplus kplus vector for uniform distributions
#' @param tau hyperparameter for prior variance
#' @export
update_alpha <- function(alpha, y, trinary, mu, sigma, kminus, kplus, tau = 1){
  tmax <- length(alpha)
  out <- numeric(length = tmax)
  imax <- nrow(y)
  lik <- numeric(length = imax)
  # calculate likelihood for each row
  for (t in 1:(tmax-1)){
    lik <- numeric(length = imax)
    # calculate likelihood for each component of alpha,ie, each value of little t
    for (i in 1:imax){
      if (trinary[i, t] == 0){
        lik[i] <- dnorm(y[i, t] - mu[i], sd = sigma[i])
      }
      if (trinary[i, t] == 1){
        lik[i] <- (alpha[t] > y[i, t] - mu[i] - kplus[i]) *
          (alpha[t] < y[i, t] - mu[i])
      }
      if (trinary[i, t] == -1){
        lik[i] <- (alpha[t] > y[i, t] - mu[i]) *
          (alpha[t] < y[i, t] - mu[i] + kminus[i])
      }
    }
    lik <- prod(lik)


  }
  return(out)
}
