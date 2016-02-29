#' Update the bernoulli parameter pp vector
#'
#' @param trinary trinary matrix
#' @export
update_pp <- function(trinary){
  imax <- nrow(trinary)
  out <- numeric(length = imax)
  for (i in 1:imax){
    n0 <- sum(trinary[i, ] == 0)
    n1 <- sum(trinary[i, ] == 1)
    out[i] <- rbeta(n = 1, shape1 = 1 + n1, shape2 = 1 + n0)
  }
  return(out)
}
