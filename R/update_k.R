#' Update kminus or kplus vector
#'
#' @param k a vector of kminus values
#' @param trinary a trinary matrix
#' @param y data matrix
#' @param indic indicator of whether to update kminus (-1) or kplus (1)
#' @param shape prior shape parameter for inv gamma
#' @param scale prior scale parameter for inv gamma
#' @export
update_k <- function(k, indic = 1, trinary, y, shape = 10, scale = 0.02, sd_jump = 0.1){
  imax <- nrow(trinary)
  tmax <- ncol(trinary)
  # make a proposal kminus
  eps <- rnorm(n = imax, mean = 0, sd = sd_jump)
  k_pro <- kminus + eps
  # prior_ratio
  prior_ratio <- MCMCpack::dinvgamma(k_pro, shape = shape, scale = scale) /
    MCMCpack::dinvgamma(k, shape = shape, scale = scale)
  # lik ratio calcs
  loglik_mat_pro <- matrix(data = 0, nrow = imax, ncol = tmax)
  loglik_mat <- matrix(data = 0, nrow = imax, ncol = tmax)
  # all entries of the matrices are 0, initially.
  # we then update those entries that have trinary indicator
  # equal to zero
  for (i in 1:imax){
    for (t in 1:tmax){
      if (trinary[i, t] == indic & (y[i, t] - mu[i] - alpha[t]) * indic > indic * k[i]){
        # subtract mu[i] and alpha[t]
        loglik_mat[i, t] <- dnorm(y[i, t] - mu[i] - alpha[t],
                                     mean = 0, sd = sigma[i],
                                     log = TRUE)
        loglik_mat_pro[i, t] <- dnorm(y[i, t] - mu[i] - alpha[t],
                                       mean = 0, sd = sigma[i],
                                       log = TRUE)
      }
      if (trinary[i, t] == 0 & (y[i, t] - mu[i] - alpha[t]) * indic < indic * k[i]){
        loglik_mat[i, t] <- - log(k[i])
        loglik_mat_pro[i, t] <- - log(k_pro[i]) # equals log(1 / k_pro[i])
      }
    } # end loop over t
  } # end loop over i
  # take sum of entries in each matrix

  loglik_ratio <- sum(loglik_mat_mu_pro) - sum(loglik_mat_mu)

}

