#' Update sigma in Gibbs sampler
#'
#' @param y data matrix
#' @param trinary trinary matrix
#' @param alpha alpha vector of sample random effects
#' @param mu mu vector of means
#' @param shape prior shape for prior on inverse sigma^2
#' @param rate prior rate for inverse sigma^2
#' @export
update_sigma <- function(y, trinary, alpha, mu, shape = 2, rate = 0.1){
  # determine imax, tmax
  imax <- nrow(y)
  out <- numeric(length = imax)
  for (i in 1:imax){
    ind <- trinary[i, ] == 0
    # determine posterior shape parameter value
    shape_post <- shape + sum(ind) / 2
    # determine posterior rate parameter value
    rate_post <- rate + sum(((y[i, ] - alpha - mu[i])[trinary == 0])^2) / 2
    foo <- rgamma(n = 1, shape = shape_post, rate = rate_post)
    out[i] <- sqrt(1 / foo)
  }
  return(out)
}
