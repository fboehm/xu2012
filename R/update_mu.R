#' Update the mean vector mu (after Xu et al., 2012)
#'
#' @param sigma a standard deviations vector
#' @param trinary a trinary matrix
#' @param tau prior sd of mu
#'
#' @export
update_mu <- function(trinary, y, alpha, tau = 1){
  imax <- nrow(trinary)
  tmax <- ncol(trinary)
  out <- numeric(length = imax)
  for (i in 1:imax){
    indic <- trinary[i,] == 0
    post_var <- 1 / (1 / tau + sum(indic) / sigma[i] ^ 2)
    post_sd <- sqrt(post_var)
    post_mean <- (sum(y[i, indic] - alpha[indic]) / (sigma[i]) ^ 2) * (post_var)
    # determine upper and lower bounds of mu
    mu_lower = max((y[i,] - alpha)[trinary[i,] == -1],(y[i,] - alpha - kplus[i])[trinary[i,] == 1], -Inf)
    mu_upper = min((y[i,] - alpha + kminus[i])[trinary[i,] == -1],(y[i,] - alpha)[trinary[i,] == 1], Inf)
    ## note how we define the sets over which we're taking min and max.
    if(min(mu_lower, mu_upper) > post_mean + 5 * post_sd){
      out[i] <- min(mu_lower, mu_upper)
    }
    else if(max(mu_lower, mu_upper) < post_mean - 5 * post_sd){
      out[i] <- max(mu_lower, mu_upper)
    }
    else{
      out[i] = mc2d::rtrunc("rnorm", n = 1, linf = min(mu_lower, mu_upper), lsup=max(mu_lower, mu_upper), mean = post_mean, sd = post_sd)
    }
  }
  return(out)
}
