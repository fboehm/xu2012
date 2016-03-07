#' Update kminus or kplus vector
#'
#' @param trinary a trinary matrix
#' @param y data matrix
#' @param indic indicator of whether to update kminus (-1) or kplus (1)
#' @param shape prior shape parameter for gamma
#' @param rate prior rate parameter for gamma
#' @param alpha alpha vector of sample random effects
#' @param mu mean vector
#' @export
update_k <- function(indic = -1, trinary, y,
                     shape = 10, rate = 50){
  # define imax
  imax <- nrow(y)
  out <- numeric(length = imax)
  for (i in 1:imax){
    ind <- trinary[i, ] == indic
    bb <- max(indic * (y[i, ] - mu[i] - alpha)[ind], 0)
    lower <- -Inf
    upper <- 1 / bb
    #calculate posterior parameter values for gamma
    shape_post <- shape + sum(ind)
    foo <- rtrunc("rgamma", n = 1, linf = lower,
           lsup = upper, shape = shape_post,
           rate = rate)
    out[i] <- 1 / foo
  }
  return(out)
}
