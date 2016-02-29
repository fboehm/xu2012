#' Draw a binary matrix
#'
#' @param beta a beta matrix
#' @param G a graph matrix
#' @param tmax number of vectors to draw
#' @export
sample_binary <- function(beta, G, tmax = 300){
  imax <- 3
  n_vecs <- 2 ^ imax
  ints <- 0:(n_vecs - 1)
  bins <- int_to_bin(ints)
  probs_unnorm <- apply(FUN = calc_prob_binary_unnorm,
                        X = bins, MARGIN = 2,
                        betamatrix = beta)
  probs <- probs_unnorm / sum(probs_unnorm) #normalized probabilities
  sampled <- sample(1:n_vecs, size = tmax, replace = TRUE, prob = probs)
  return(bins[, sampled])
}

#' Calculate unnormalized probability for binary vector
#'
#' @param binary_vec a vector of zeros and ones
#' @param betamatrix a matrix of beta coefficients
#' @export
calc_prob_binary_unnorm <- function(binary_vec, betamatrix){
  term1 <- diag(betamatrix) %*% binary_vec
  matrix_prod1 <- (binary_vec - boehm::expit(diag(betamatrix))) %*%
    t(binary_vec - boehm::expit(diag(betamatrix)))
  matrix_prod <- betamatrix * matrix_prod1
  term2 <- sum(matrix_prod[lower.tri(matrix_prod)])
  return(exp(as.numeric(term1 + term2)))
}


#' Converts integer to vector of zeros and ones
#'
#' @param x a numeric vector of positive integers
#' @return a matrix in which each column corresponds to one integer of x vector
#' @export
int_to_bin <- function(x){
  char <- R.utils::intToBin(x)
  foo <- stringr::str_split(char, pattern = "")
  out <- sapply(FUN = as.numeric, X = foo)
  return(out)
}
