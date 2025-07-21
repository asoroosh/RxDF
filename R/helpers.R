
#' Simulate a Single AR(1) Time Series
#'
#' Generates a single autoregressive time series of order 1 (AR(1)) with
#' specified length and autoregressive coefficient. The series is initialized
#' with a standard normal value and evolves as:
#' \deqn{Y_t = \phi Y_{t-1} + \epsilon_t, \quad \epsilon_t \sim \mathcal{N}(0, 1)}
#'
#' @param Tt Integer. Number of time points in the time series.
#' @param phi Numeric. AR(1) coefficient (typically between -1 and 1 for stationarity).
#'
#' @return A numeric vector of length `Tt` representing the simulated AR(1) time series.
#'
#' @examples
#' sim_ar1(Tt = 100, phi = 0.8)
#'
#' @seealso \code{\link{sim_many_ts_same_ar1}} for simulating multiple AR(1) series.
#'
#' @export
sim_ar1 <- function(Tt, phi) {
  y <- numeric(Tt)
  y[1] <- rnorm(1)
  for (t in 2:Tt) {
    y[t] <- phi * y[t - 1] + rnorm(1)
  }
  return(y)
}

#' Simulate Multiple AR(1) Time Series with Shared Coefficient
#'
#' Generates multiple independent realisations of an AR(1) (autoregressive of order 1)
#' time series, each of length `Tt` and with the same autoregressive coefficient `phi`.
#' Each series is initialized with a standard normal value and evolves according to:
#' \deqn{Y_t = \phi Y_{t-1} + \epsilon_t, \quad \epsilon_t \sim \mathcal{N}(0, 1)}
#'
#' @param Tt Integer. Number of time points in each simulated time series.
#' @param phi Numeric. AR(1) coefficient (typically between -1 and 1 for stationarity).
#' @param n_series Integer. Number of independent AR(1) series to simulate.
#'
#' @return A numeric matrix of dimensions `Tt` x `n_series`, where each column
#'         corresponds to one simulated AR(1) time series.
#'
#' @examples
#' sim_many_ts_same_ar1(Tt = 100, phi = 0.7, n_series = 50)
#'
#' @seealso \code{\link{sim_ar1}} for simulating a single AR(1) process.
#'
#' @export
sim_many_ts_same_ar1 <- function(Tt, phi, n_series) {
  return(replicate(n_series, sim_ar1(Tt, phi)))
}

#' Ensure Correct Dimensionality of Time Series Data. It ensures matrics are
#' number of data points x number of time series.
#'
#' This utility function checks and standardizes the input `Y` to ensure it has
#' the expected time series length `Tt` as rows. It accepts either a vector
#' (single time series) or a matrix (multiple time series).
#'
#' @param Y A numeric vector (length Tt) or matrix (with Tt rows or Tt columns).
#' @param Tt Integer. The expected length of each time series (i.e., the number of time points).
#'
#' @return A numeric matrix with `Tt` rows. If `Y` was a vector, it is converted
#'         to a matrix with one column. If `Y` is a matrix with `Tt` columns but
#'         not `Tt` rows, it is transposed.
#'
#' @examples
#' check_dim(1:100, 100)              # returns 100x1 matrix
#' check_dim(matrix(rnorm(300), 100), 100)  # returns 100x3 matrix
#'
#' @export
check_dim <- function(Y,Tt){
  if (is.vector(Y)) {
    if (length(Y) != Tt) stop("Length of vector Y must equal Tt.")
    Y <- matrix(Y, ncol = 1)
  }

  if (!is.null(dim(Y)) && length(dim(Y)) == 2) {
    if (!(Tt %in% dim(Y))) stop("Dimensions are wrong!")
    if (ncol(Y) != Tt) Y <- t(Y)
  }
  return(Y)
}

#' Demean a matrix along rows or columns
#'
#' Removes the mean from each row or column of a matrix `Y`.
#' This is useful for centering time series prior to further processing (e.g., autocorrelation, FFT).
#' The dimension along which demeaning is performed is controlled by the `dim` argument.
#'
#' @param Y Numeric matrix. The data matrix to demean.
#' Rows or columns represent separate time series depending on the `dim` argument.
#' @param dim Integer. Direction along which to demean:
#'   - `dim = 1`: demean columns (subtract mean of each column),
#'   - `dim = 2`: demean rows (subtract mean of each row, default).
#'
#' @return Numeric matrix of the same dimensions as `Y` with means removed along the specified dimension.
#'
#' @details
#' If `dim = 1`, each column is centered so that its mean is zero.
#' If `dim = 2`, each row is centered so that its mean is zero.
#' An error is raised if `dim` is not 1 or 2.
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(12), nrow = 3, ncol = 4)
#' demean_ts(Y, dim = 1)  # Demean columns
#' demean_ts(Y, dim = 2)  # Demean rows
#'
#' @export
demean_ts <- function(Y, dim = 2)
if (dim == 1) {
  # Demean columns (each column is a timeseries)
  Y <- scale(Y, center = TRUE, scale = FALSE)
} else if (dim == 2) {
  # Demean rows (each row is a timeseries)
  row_means <- rowMeans(Y)
  Y <- Y - row_means
} else {
  stop("Invalid dim: must be 1 or 2")
}

#' Estimate AR(1) Coefficients Using a Rough Moment-Based Estimator
#'
#' Computes rough estimates of AR(1) coefficients for one or more time series using
#' a simplified moment-based formula:
#' \deqn{\hat{\phi} = \frac{\sum_{t=2}^{T} Y_{t-1} Y_t}{\sum_{t=1}^{T} Y_t^2}}
#' This method is fast and suitable for simulations but does not correspond to
#' the ordinary least squares (OLS) estimator. Use with caution for inference.
#'
#' @param Y A numeric vector (single time series) or matrix (multiple series), where
#'          each column is a time series of length `Tt`.
#' @param Tt Integer. The expected length (number of time points) of each time series.
#'
#' @return A numeric vector of estimated AR(1) coefficients, one per column of `Y`.
#'
#' @details If `Y` is a vector, it is first converted into a matrix with one column.
#'          If `Y` is not in column-wise format, an attempt is made to transpose it
#'          so that each column corresponds to a time series of length `Tt`.
#'
#' @examples
#' Y <- sim_many_ts_same_ar1(100, 0.6, 10)
#' est_rough_ar1(Y, 100)
#'
#' @seealso \code{\link{check_dim}}, \code{\link{sim_ar1}}, \code{\link{sim_many_ts_same_ar1}}
#'
#' @export
est_rough_ar1 <- function(Y, Tt){
  Y <- check_dim(Y, Tt)
  Y <- demean_ts(Y, dim = 2)
  arone <- apply(Y, 1, function(g) {
    sum(g[1:(Tt - 1)] * g[2:Tt]) / sum(g^2)
  })
  return(arone)
}

#' Compute Two-Tailed P-Values from a Correlation Matrix via Fisher Z-Transform
#'
#' Applies the Fisher Z-transformation to a correlation matrix and computes
#' two-tailed p-values for testing the null hypothesis that each correlation
#' coefficient is zero.
#'
#' @param netmat A numeric matrix representing pairwise correlation coefficients
#'               (typically symmetric, with values between -1 and 1).
#' @param R2Zcrt Numeric. A correction factor to scale the Fisher-transformed values
#'               (e.g., \eqn{\sqrt{df - 3}} where `df` is the degrees of freedom).
#'
#' @return A matrix of the same dimensions as `netmat`, containing two-tailed
#'         p-values for each correlation.
#'
#' @details The Fisher Z-transformation is applied as \eqn{Z = \text{atanh}(r) \cdot R2Zcrt},
#' where \eqn{r} is the raw correlation. Two-tailed p-values are then computed as:
#' \deqn{p = 2 \cdot \Phi(-|Z|)}
#' where \eqn{\Phi} is the standard normal cumulative distribution function.
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(100*10), 100, 10)
#' corr_mat <- cor(mat)
#' pvals <- fisher_pval_matrix(corr_mat, sqrt(97))  # df = 100 - 3
#'
#' @export
fisher_pval_matrix <- function(netmat, R2Zcrt) {
  Z <- atanh(netmat) * R2Zcrt
  P <- 2 * pnorm(-abs(Z))
  return(list(Z = Z, P = P))
}


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
GenTsAC1 <- function(ar1, Tt) {
  if (!is.numeric(ar1) || ar1 < -1 || ar1 > 1) {
    stop("The 'ar1' parameter must be a numeric value between -1 and 1.")
  }

  if (!is.numeric(Tt) || Tt <= 0 || Tt != round(Tt)) {
    stop("The 'Tt' parameter must be a positive integer.")
  }

  # Initialize the time series
  ts = numeric(Tt)
  ts[1] = rnorm(1)
  for (t in 2:Tt) {
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
corrautocorr <- function(mu,
                         sigR,
                         sigC,
                         Tt,
                         verboseflag = TRUE) {

  # Handle `sigC` as a vector, list of vectors, or matrix
  if (is.atomic(sigC) && !is.list(sigC)) {

    # sigC is a vector
    sigC <- make_toeplitz(sigX = sigC, Tt = Tt)
  } else if (is.list(sigC)) {

    # sigC is a list of vectors
    if (length(sigC) != length(mu)) {
      stop("When `sigC` is a list, it must have the same length as `mu`.")
    }
    sigC <- lapply(sigC, function(vec) make_toeplitz(sigX = vec, Tt = Tt))

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
  z <- matrix(rnorm(length(mu) * Tt), nrow = length(mu), ncol = Tt)
  MVDisk <- t(CsigR) %*% z

  # Handle `sigC` as a list of Toeplitz matrices
  if (is.list(sigC)) {
    t <- matrix(NA, nrow = length(mu), ncol = Tt)
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

    t <- matrix(mu,
                nrow = length(mu),
                ncol = Tt,
                byrow = TRUE) + MVDisk %*% t(CsigC)
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
#' make_toeplitz(sigX, Tt)
#'
#' # Example 2: Create a larger Toeplitz matrix (10x10) with the same off-diagonal values
#' make_toeplitz(sigX, 10)
#'
#' @export
make_toeplitz <- function(sigX, Tt) {
  diag_values = c(1, sigX)
  full_values = c(diag_values, rep(0, Tt - length(diag_values)))  # Extend with zeros
  toeplitz_matrix <- toeplitz(full_values[1:Tt])  # Ensure size matches Tt
  return(toeplitz_matrix)
}

#' Calculate the Autocorrelation Length
#'
#' This function calculates the autocorrelation length, which measures the
#' degree of autocorrelation in a time series. It uses the squared values
#' of the autocorrelation function (ACF) up to a specified number of lags.
#' The implementation is adapted from Straatsma et al. (2016).
#'
#' @param ts Numeric vector. The time series for which the autocorrelation length is calculated.
#' @param T Integer. The maximum number of lags to consider in the autocorrelation function.
#' @return A single numeric value representing the autocorrelation length.
#' @references
#' Straatsma, T. P., Berendsen, H. J. C., & Stam, A. J. (2016).
#' Estimation of statistical errors in molecular simulation calculations.
#' http://doi.org/10.1080/00268978600100071
#'
#' Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols.
#' "Effective Degrees of Freedom of the Pearson's Correlation Coefficient
#' under Serial Correlation." bioRxiv (2018): 453795.
#' @examples
#' Y <- rnorm(100)  # Example time series
#' T <- 50           # Number of lags
#' autocorr_length <- AutoCorrLength(Y, T)
#' print(autocorr_length)
#'
#' @export
AutoCorrLength <- function(Y, nlag) {
  # Ensure the input is valid
  if (!is.numeric(Y)) stop("The input time series must be numeric.")
  if (!is.numeric(nlag) || nlag <= 0 || nlag != round(nlag)) stop("T must be a positive integer.")

  acf_values <- acf(Y,
                    lag.max = nlag,
                    plot = FALSE)$acf[-1]
  CorrLeng <- sum(acf_values^2)
  return(CorrLeng)
}


#' Pairwise Timeseries Sum Matrix
#'
#' Computes a 3D array containing the element-wise sums of all unique pairs of timeseries
#' in a matrix `Y0`. For each timepoint, this function calculates the sum of each pair
#' of timeseries values at that timepoint and stores them in a symmetric 2D matrix slice.
#' The diagonal elements are set to zero.
#'
#' @param Y0 Numeric matrix. A matrix of size `Tt × n` or `n × Tt`, where `Tt` is the number
#' of timepoints and `n` is the number of timeseries. The function transposes `Y0` internally
#' if needed to ensure rows represent timepoints.
#' @param Tt Integer. Number of timepoints (`Tt` should match one dimension of `Y0`).
#'
#' @return A 3D array of dimensions `n × n × Tt`. For each timepoint `t`, the slice `SM0[ , , t]`
#' contains the element-wise sums of all pairs of timeseries at timepoint `t`.
#' The matrix is symmetric for each timepoint, and diagonal elements are set to zero.
#'
#' @details
#' - `SM0[i, j, t]` gives `Y0[t, i] + Y0[t, j]` for `i ≠ j`.
#' - `SM0[i, i, t]` is explicitly set to 0 for all `i` and `t`.
#' - Intended for fast computation of pairwise sums across all timeseries and timepoints.
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(5 * 100), nrow = 100, ncol = 5)
#' SM0 <- SumMat(Y, 100)
#' dim(SM0)  # 5 x 5 x 100
#'
#' @export
SumMat <- function(Y0, Tt) {
  if (!(Tt %in% dim(Y0))) stop("Dimensions are wrong!")
  if (nrow(Y0) != Tt) Y0 <- t(Y0)

  nn <- ncol(Y0)
  SM0 <- array(0, dim = c(nn, nn, Tt))
  for (t in 1:Tt) {
    yt <- Y0[t, ]
    SM0[, , t] <- outer(yt, rep(1, nn)) + outer(rep(1, nn), yt)
  }
  for (i in 1:nn) SM0[i, i, ] <- 0

  return(SM0)
}

#' Pairwise Timeseries Product Matrix
#'
#' Computes a 3D array containing the element-wise products of all pairs of timeseries
#' in a matrix `Y0`. For each timepoint, this function calculates the product of each pair
#' of timeseries values at that timepoint and stores them in a symmetric 2D matrix slice.
#' The diagonal elements are set to zero.
#'
#' @param Y0 Numeric matrix. A matrix of size `Tt × n` or `n × Tt`, where `Tt` is the number
#' of timepoints and `n` is the number of timeseries. The function transposes `Y0` internally
#' if needed to ensure rows represent timepoints.
#' @param Tt Integer. Number of timepoints (`Tt` should match one dimension of `Y0`).
#'
#' @return A 3D array of dimensions `n × n × Tt`. For each timepoint `t`, the slice `PM0[ , , t]`
#' contains the element-wise products of all pairs of timeseries at timepoint `t`.
#' The matrix is symmetric for each timepoint, and diagonal elements are set to zero.
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(5 * 100), nrow = 100, ncol = 5)
#' PM0 <- ProdMat(Y, 100)
#' dim(PM0)  # 5 x 5 x 100
#'
#' @export
ProdMat <- function(Y0, Tt) {
  if (!(Tt %in% dim(Y0))) stop("Dimensions are wrong!")
  if (nrow(Y0) != Tt) Y0 <- t(Y0)

  nn <- ncol(Y0)
  PM0 <- array(0, dim = c(nn, nn, Tt))

  for (t in 1:Tt) {
    yt <- Y0[t, ]
    PM0[, , t] <- yt %o% yt
  }

  for (i in 1:nn) PM0[i, i, ] <- 0

  return(PM0)
}

#' Apply Tukey Taper to Rows of a Matrix
#'
#' Applies a Tukey taper window to each row of a matrix `acs`, independently, over the first `M` columns (lags).
#' After column `M`, all values are set to zero. This is typically used to smooth autocorrelation or cross-correlation
#' sequences to reduce variance from high lags.
#'
#' @param acs Numeric matrix. A matrix of size `I × Tt`, where `I` is the number of time series (rows)
#' and `Tt` is the number of lags (columns). Each row is treated as an independent sequence.
#' @param Tt Integer. Number of lags (columns). Used to check consistency and ensure correct orientation.
#' @param M Integer. Length of the taper window. The Tukey taper will be applied over columns 1 to `M`.
#' Columns `M+1` to `Tt` are set to zero.
#'
#' @return A matrix of the same dimension as `acs` (`I × Tt`), where each row has been multiplied by
#' the Tukey taper in columns 1 to `M`, and zeros elsewhere.
#'
#' @details
#' The Tukey taper weights applied are defined as:
#' \deqn{
#'   w(m) = \frac{1 + \cos(m \pi / M)}{2}, \quad m = 1, 2, \dots, M
#' }
#'
#' The taper smoothly down-weights larger lags, which is useful when estimating quantities like
#' effective degrees of freedom in correlated data.
#'
#'
#' @export
tukey_taper_me <- function(acs, Tt, M = sqrt(Tt)) {
  acs <- check_dim(acs,Tt)

  tt_ts <- matrix(0, nrow = nrow(acs), ncol = Tt)

  M <- round(M)
  w <- (1 + cos((1:M) * pi / M)) / 2
  tt_ts[,1:M] <- sweep(acs[, 1:M, drop = FALSE], 2, w, "*")

  return(tt_ts)
}

#' Curb taper autocorrelation sequences
#'
#' Performs a hard truncation ("curb taper") of autocorrelation or cross-correlation sequences
#' by zeroing out all lag values beyond lag `M` for each timeseries.
#'
#' @param acs Numeric matrix. A matrix of size `I × Tt`, where `I` is the number of timeseries
#' and `Tt` is the number of lags (columns).
#' @param Tt Integer. Number of lags (columns), used for consistency check.
#' @param M Integer. Maximum lag to retain. All columns from `M+1` to `Tt` will be set to zero.
#'
#' @return A matrix of same dimensions as `acs` (`I × Tt`), with values zeroed beyond lag `M`.
#'
#' @details
#' This function implements the "curb taper" described by Anderson (1984) by simply truncating
#' autocorrelation sequences at lag `M`. No gradual tapering is applied.
curb_taper_me <- function(acs, Tt, M) {
  acs <- check_dim(acs,Tt)

  M <- round(M)

  if (Tt == M){
    stop('The curbing window is full length of timeseries. Tt = M.')
  }

  ct_ts <- acs
  ct_ts[,(M + 1):Tt] <- 0

  return(ct_ts)
}


#' Find breakpoint where autocorrelation becomes insignificant
#'
#' Determines the lag index at which the absolute autocorrelation sequence first falls
#' below an approximate 95% confidence bound for white noise. This identifies where
#' autocorrelation ceases to be statistically significant and can serve as a heuristic
#' cutoff for truncating or tapering.
#'
#' @param acs Numeric vector of length `Tt`. The autocorrelation sequence (lags 1 to `Tt`).
#' @param Tt Integer. Total number of lags. Used for consistency check and orientation.
#'
#' @return Integer scalar `where2stop` indicating the largest lag before autocorrelation
#' becomes insignificant:
#' - If all autocorrelations are significant ⇒ returns `Tt`.
#' - If autocorrelation is insignificant from the start ⇒ returns 0.
#' - Otherwise ⇒ returns index of last significant lag before first insignificance.
#'
#' @details
#' The threshold for significance is defined as:
#' \deqn{
#'   \text{bound} = \frac{1.96}{\sqrt{Tt}}
#' }
#' corresponding to the 95% confidence bound for white noise autocorrelations.
#'
#' The function checks when `|acs| <= bound` and reports the first lag where this occurs,
#' then returns `where2stop` as the lag immediately preceding that point.
#'
#' @examples
#' set.seed(123)
#' acs <- acf(rnorm(100), plot = FALSE, lag.max = 30)$acf[-1]  # Exclude lag 0
#' find_break_point(acs, 30)
#'
#' @export
find_break_point <- function(acs, Tt) {
  acs <- check_dim(acs,Tt)

  # 95% confidence bound for white noise
  bnd <- qnorm(0.975) / sqrt(Tt)

  # Find first lag where abs(acs) <= bnd
  first_below <- which(abs(acs) <= bnd)[1]

  if (is.na(first_below)) {
    where2stop <- Tt  # Never drops below threshold
  } else if (first_below == 1) {
    where2stop <- 0   # Below threshold from start
  } else {
    where2stop <- first_below - 1  # Stop right before first below-threshold lag
  }

  return(where2stop)
}


#' Shrink autocorrelation sequence beyond significance threshold
#'
#' Applies an adaptive truncation ("shrinkage") procedure to an autocorrelation sequence.
#' It identifies the first lag where autocorrelations fall below a 95% confidence bound
#' and sets all subsequent lags to zero.
#'
#' @param acs Numeric vector or matrix representing a single autocorrelation sequence.
#' Length must equal `Tt`. If input is a matrix, it must be 1 row or column.
#' @param Tt Integer. Number of lags (length of `acs`).
#'
#' @return Numeric vector of length `Tt` containing the shrunk autocorrelation sequence.
#'
#' @details
#' - If all autocorrelations are below the confidence bound from the beginning, returns all zeros.
#' - Uses `find_break_point()` internally to determine truncation point.
#' - Calls `curb_taper_me()` internally for hard truncation at the break point.
#'
#' @seealso \code{\link{find_break_point}}, \code{\link{curb_taper_me}}
#'
shrink_me <- function(acs, Tt) {
  acs <- check_dim(acs,Tt)
  where2stop <- find_break_point(acs, Tt)

  if (where2stop == 0) {
    srnkd_ts <- rep(0, Tt)
  } else {
    srnkd_ts <- curb_taper_me(acs, Tt, where2stop)
  }

  return(srnkd_ts)
}



