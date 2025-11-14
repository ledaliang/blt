#' Compute Hamming Distance Matrix
#'
#' Computes the pairwise Hamming distances between rows of a binary matrix.
#' The function returns a symmetric distance matrix where entry \code{(i, j)}
#' is the number of positions where rows \code{i} and \code{j} differ.
#'
#' @param data A numeric or logical matrix where each row is a sequence.
#'
#' @return A symmetric numeric distance matrix.
#'
#' @examples
#' mat <- matrix(c(0,1,0, 1,1,0, 0,0,1), nrow = 3, byrow = TRUE)
#' hamming_dist_matrix(mat)
#'
#' @export
hamming_dist_matrix <- function(data) {
  ## Data in the form of lumpedD
  n <- nrow(data)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist_matrix[i,j] <- sum(data[i,] != data[j,])
      dist_matrix[j,i] <- dist_matrix[i,j]
    }
  }
  return(dist_matrix)
}


#' Flatten Upper Triangular Matrix by Diagonals
#'
#' Converts the upper triangular part of a square matrix into a vector,
#' traversing along successive diagonals starting from the main diagonal.
#'
#' For example, the matrix:
#' \preformatted{
#' [theta1, theta12, theta13]
#' [0     , theta2 , theta23]
#' [0     , 0      , theta3 ]
#' }
#' becomes:
#' \code{c(theta1, theta2, theta3, theta12, theta23, theta13)}.
#'
#' @param M A square matrix.
#' @return A numeric vector containing the upper triangular elements of \code{M},
#'         ordered by diagonal.
#'
#' @examples
#' M <- matrix(c(1, 2, 3,
#'               0, 4, 5,
#'               0, 0, 6), nrow = 3, byrow = TRUE)
#' flatten_upper(M)
#' # Output: 1 4 6 2 5 3
#'
#' @export
flatten_upper <- function(M) {
  # Flattens theta upper triangular theta matrix into theta vector.
  # Input:
  # theta = [theta1, theta12, theta13]
  #         [0     , theta2 , theta23]
  #         [0     , 0      , theta3 ]
  # Output:
  # theta = [theta1, theta2, theta3, theta12, theta23, theta13]
  S <- nrow(M)
  out <- c()
  for (k in 0:(S-1)) {
    for (i in 1:(S-k)) {
      j <- i + k
      out <- c(out, M[i, j])
    }
  }
  return(out)
}



#' Reconstruct Upper Triangular Matrix from Flattened Vector
#'
#' Converts a vector created by \code{flatten_upper()} back into a square
#' upper triangular matrix (including main diagonal). Zeros are filled in
#' the lower triangular part.
#'
#' @param v A numeric vector containing the upper triangular elements, ordered
#'          by diagonal (as returned by \code{flatten_upper()}).
#'
#' @return A square matrix \code{M} such that \code{flatten_upper(M) == v}.
#'
#' @examples
#' theta <- c(1,2,3,4,5,6)
#' unflatten_upper(theta)
#' #      [,1] [,2] [,3]
#' # [1,]    1    4    6
#' # [2,]    0    2    5
#' # [3,]    0    0    3
#'
#' @export
unflatten_upper <- function(v) {
  # Solve quadratic to find matrix size S from length of v
  len <- length(v)
  S <- (-1 + sqrt(1 + 8 * len)) / 2
  if (S != floor(S)) stop("Vector length is not valid for triangular matrix.")

  M <- matrix(0, nrow = S, ncol = S)

  idx <- 1
  for (k in 0:(S-1)) {
    for (i in 1:(S-k)) {
      j <- i + k
      M[i, j] <- v[idx]
      idx <- idx + 1
    }
  }
  return(M)
}


