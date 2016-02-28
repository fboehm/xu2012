#' Update trinary matrix
#'
#' @param binary a binary matrix of zeros and ones
#' @param y a data matrix
#' @param pp a vector of probabilities
#' @param mu a vector of means
#' @param sigma a vector of standard deviations
#' @param kplus a vector of kplus values
#' @export


update_trinary <- function(binary, y, pp, mu, sigma, kplus){
  imax <- nrow(binary)
  tmax <- ncol(binary)
  trinary <- matrix(nrow = imax, ncol = tmax)
  for (i in 1:imax){
    for (t in 1:tmax){
      if (binary[i, t] == 0) trinary[i, t] <- -1
      if (binary[i, t] == 1){
        p2 <- dnorm(y[i, t], mean = alpha[t] + mu[i], sd = sigma[i])
        p1 <- (y[i, t] < alpha[t] + mu[i] + kplus[i] &
                 y[i, t] > alpha[t] + mu[i] ) / kplus[i]
        prob <- pp[i] * p1 / (pp[i] * p1 + (1 - pp[i]) * p2)
        trinary[i, t] <- rbinom(n = 1, size = 1, prob = prob)
      }
    }
  }
  return(trinary)
}
