#' Sample a trinary matrix
#'
#' @param binary a binary matrix
#' @param pp a probability vector
#' @export


sample_trinary <- function(binary, pp){
  imax <- nrow(binary)
  tmax <- ncol(binary)
  trinary <- matrix(nrow = imax, ncol = tmax)
  for (i in 1:imax){
    for (t in 1:tmax){
      trinary[i, t] <- -(binary[i, t] == 0) + (binary[i, t] == 1) * rbinom(n = 1, size = 1, prob = 1 - pp[i])
    }
  }
  return(trinary)
}
