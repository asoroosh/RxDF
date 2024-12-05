#' Generate a Time Series with First-Order Autoregressive Correlation
#'
#' This function generates a time series with first-order autoregressive (AR(1)) correlation.
#'
#' @param ar1 Numeric. The first-order autoregressive coefficient. Must be between -1 and 1 for stationarity.
#' @param ndpr Integer. The number of time points in the generated time series. Must be a positive integer.
#' @return A numeric vector of length `ndpr`, representing the generated time series.
#' @examples
#' # Generate an AR(1) time series with coefficient 0.8 and 100 time points
#' ts <- GenTsAC1(0.8, 100)
#' plot(ts, type = "l", main = "Generated AR(1) Time Series", ylab = "Value", xlab = "Time")
#'
#' # Generate a white noise series with AR(1) coefficient 0
#' ts_white_noise <- GenTsAC1(0, 50)
#'
#' @export
GenTsAC1 <- function(ar1, ndpr) {
  if (!is.numeric(ar1) || ar1 < -1 || ar1 > 1) {
    stop("The 'ar1' parameter must be a numeric value between -1 and 1.")
  }
  
  if (!is.numeric(ndpr) || ndpr <= 0 || ndpr != round(ndpr)) {
    stop("The 'ndpr' parameter must be a positive integer.")
  }
  
  # Initialize the time series
  ts = numeric(ndpr)
  ts[1] = rnorm(1)
  for (t in 2:ndpr) {
    ts[t] = ts[t - 1] * ar1 + rnorm(1)
  }
  return(ts)
}


#' Generate Autocorrelated Time Series from Matrix Normal Distribution
#'
#' This function generates time series with specified autocorrelation and cross-correlation structures
#' by sampling from a Matrix Normal Distribution using Cholesky decomposition. It handles non-PSD matrices
#' by converting them to the nearest positive semi-definite (PSD) matrix.
#'
#' @param mu Numeric vector. Expected values, of length equal to the number of time series required.
#' @param sigR Numeric matrix or scalar. Covariance matrix between time series. If two univariate time series
#'   are required, a scalar indicating the off-diagonal elements suffices.
#' @param sigC Numeric vector, list of vectors, or matrix. Covariance of columns (serial correlation). 
#'   - If a vector, it is converted to a Toeplitz matrix for universal AC structure across all time series.
#'   - If a list of vectors, each vector specifies the autocorrelation structure for a specific time series.
#'   - If a matrix, it is treated as a universal AC structure for all time series.
#' @param ndp Integer. Number of data points (length of time series).
#' @param verboseflag Logical. If TRUE, warnings about non-PSD matrices will be printed. Default is TRUE.
#' @return A numeric matrix where each row represents a time series.
#' @examples
#' # Example 1: Single autocorrelation structure for all time series
#' ts1 <- corrautocorr(c(0, 0), 0.6, c(0.6, 0.3, 0.2), 1000)
#' cor(ts1)
#'
#' # Example 2: Different autocorrelation structures for each time series
#' sigC_list <- list(c(0.6, 0.3, 0.2), c(0.5, 0.25, 0.1))
#' ts2 <- corrautocorr(c(0, 0), 0.6, sigC_list, 1000)
#' cor(ts2)
#'
#' @references
#' Higham, N. J. (1988). Computing the nearest symmetric positive semi-definite matrix. Linear Algebra and its Applications.
#' Soroosh Afyouni, Stephen M. Smith, & Thomas E. Nichols (2018). Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation.
#' @export
corrautocorr <- function(mu, sigR, sigC, ndp, verboseflag = TRUE) {
  # Handle `sigC` as a vector, list of vectors, or matrix
  if (is.atomic(sigC) && !is.list(sigC)) {
    # sigC is a vector
    sigC <- make_toeplitz(sigX = sigC, ndp = ndp)
  } else if (is.list(sigC)) {
    # sigC is a list of vectors
    if (length(sigC) != length(mu)) {
      stop("When `sigC` is a list, it must have the same length as `mu`.")
    }
    sigC <- lapply(sigC, function(vec) make_toeplitz(sigX = vec, ndp = ndp))
  }
  
  # Handle `sigR` as a scalar
  if (is_scalar(sigR)) {
    if (length(mu) != 2) {
      stop("sigR must be a square matrix when more than 2 time series are required.")
    }
    sigR <- matrix(sigR, nrow = 2, ncol = 2)
    diag(sigR) <- 1
  }
  
  # Ensure `sigR` is PSD
  CsigR <- tryCatch(
    chol(sigR),
    error = function(e) {
      if (verboseflag) warning("sigR is not positive semi-definite. Using nearest PSD matrix.")
      chol(nearestSPD(sigR))
    }
  )
  
  # Generate random data and apply transformations
  z <- matrix(rnorm(length(mu) * ndp), nrow = length(mu), ncol = ndp)
  MVDisk <- t(CsigR) %*% z
  
  # Handle `sigC` as a list of Toeplitz matrices
  if (is.list(sigC)) {
    t <- matrix(NA, nrow = length(mu), ncol = ndp)
    for (i in seq_along(mu)) {
      CsigC <- tryCatch(
        chol(sigC[[i]]),
        error = function(e) chol(nearestSPD(sigC[[i]]))
      )
      t[i, ] <- mu[i] + MVDisk[i, ] %*% t(CsigC)
    }
  } else {
    # Ensure `sigC` is PSD
    CsigC <- tryCatch(
      chol(sigC),
      error = function(e) {
        if (verboseflag) warning("sigC is not positive semi-definite. Using nearest PSD matrix.")
        chol(nearestSPD(sigC))
      }
    )
    t <- matrix(mu, nrow = length(mu), ncol = ndp, byrow = TRUE) + MVDisk %*% t(CsigC)
  }
  
  return(t)
}





#' Convert a Matrix to the Nearest Positive Semi-Definite Matrix
#'
#' This function converts a square matrix to the nearest positive semi-definite (PSD) matrix
#' using the algorithm described by Higham (1988). The input matrix is symmetrized and adjusted
#' iteratively to ensure all eigenvalues are non-negative.
#'
#' @param A Numeric matrix. A square matrix to be converted to its nearest PSD matrix.
#' @return A numeric matrix that is positive semi-definite and as close as possible to the input matrix `A`.
#' @examples
#' # Example of converting a non-PSD matrix to PSD
#' mat <- matrix(c(4, -2, 2, -2, 3, -1, 2, -1, 2), nrow = 3, byrow = TRUE)
#' psd_mat <- nearestSPD(mat)
#' eigen(psd_mat)$values  # All eigenvalues should now be non-negative
#'
#' @references
#' Higham, N. J. (1988). Computing the nearest symmetric positive semi-definite matrix. 
#' Linear Algebra and its Applications, 103, 103-118.
#' @export
nearestSPD <- function(A) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("Input must be a square matrix.")
  }
  
  # Step 1: Symmetrize the matrix
  B <- (A + t(A)) / 2
  
  # Step 2: Compute the symmetric polar factor
  eig <- eigen(B)
  H <- eig$vectors %*% diag(pmax(eig$values, 0)) %*% t(eig$vectors)
  
  # Step 3: Adjust the matrix
  Ahat <- (B + H) / 2
  Ahat <- (Ahat + t(Ahat)) / 2  # Ensure symmetry
  
  # Step 4: Ensure all eigenvalues are non-negative
  p <- 1
  while (p != 0) {
    eigenvalues <- eigen(Ahat)$values
    if (min(eigenvalues) >= 0) {
      p <- 0
    } else {
      Ahat <- Ahat + diag(abs(min(eigenvalues)) * 1e-8, nrow(Ahat))
    }
  }
  
  return(Ahat)
}


#' Check if a Value is a Scalar
#'
#' This function checks whether the input is a numeric scalar. 
#' A scalar is defined as a single numeric value (i.e., a vector of length 1).
#'
#' @param x Any R object. The object to be checked.
#' @return A logical value: `TRUE` if the input is a numeric scalar, otherwise `FALSE`.
#' @examples
#' # Check scalar values
#' is_scalar(5)          # TRUE
#' is_scalar(c(1, 2, 3)) # FALSE
#' is_scalar("text")     # FALSE
#' is_scalar(NA)         # FALSE
#' @export
is_scalar <- function(x) {
  is.numeric(x) && length(x) == 1
}

#' Create a Toeplitz Matrix with Specified Off-Diagonal Entries
#'
#' This function generates a Toeplitz matrix of a specified size, where the main diagonal
#' is filled with 1s, and the off-diagonal entries are specified by a given vector. If the
#' matrix size exceeds the length of the vector, the remaining entries are set to 0.
#'
#' @param sigX Numeric vector. The off-diagonal entries to fill the Toeplitz matrix.
#' @param ndp Integer. The desired size of the Toeplitz matrix (number of rows and columns).
#' @return A numeric Toeplitz matrix of dimensions `ndp x ndp` with the specified structure:
#' - The main diagonal is filled with 1s.
#' - The first sub- and super-diagonals are filled with the first element of `sigX`, 
#'   the second sub- and super-diagonals with the second element of `sigX`, and so on.
#' - Any remaining diagonals are filled with 0 if `sigX` is shorter than `ndp - 1`.
#' @examples
#' # Example 1: Create a 5x5 Toeplitz matrix with specified off-diagonal values
#' sigX <- c(0.6, 0.3, 0.2)
#' ndp <- 5
#' make_toeplitz(sigX, ndp)
#'
#' # Example 2: Create a larger Toeplitz matrix (10x10) with the same off-diagonal values
#' make_toeplitz(sigX, 10)
#'
#' @export
make_toeplitz <- function(sigX, ndp) {
  diag_values = c(1, sigX)  
  full_values = c(diag_values, rep(0, ndp - length(diag_values)))  # Extend with zeros
  toeplitz_matrix <- toeplitz(full_values[1:ndp])  # Ensure size matches ndp
  return(toeplitz_matrix)
}

