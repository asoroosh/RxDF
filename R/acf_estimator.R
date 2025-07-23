#' Fast Autocovariance and Autocorrelation Calculation Using FFT
#'
#' Computes autocovariance and autocorrelation functions efficiently for multiple time series using the Fast Fourier Transform (FFT).
#' Optionally applies zero-padding for performance optimization and unbiased lag-specific normalization to match R's `acf()` behavior.
#'
#' @param Y Numeric matrix. A matrix of size IxT where I is the number of time series (rows) and T is their length (columns).
#'          Each row is treated as an individual time series.
#' @param Tt Integer. The number of datapoints in the timeseries. Is always used for sanity checking.
#' @param lg Optional integer vector of lags to return (not used in current implementation).
#' @param pad Logical. If TRUE (default), zero-pads time series to next power of two for improved FFT performance.
#' @param bias_correction Logical. If TRUE (default), applies unbiased lag-dependent normalization (scales by N / (N - k) at lag k) to match R's `acf(..., type = "correlation")`.
#'
#' @return A list with:
#'   - `ACOV`: Autocovariance matrix of size IxTt (unbiased if `bias_correction = TRUE`).
#'   - `ACOR`: Autocorrelation matrix of size IxTt (normalized by variance at lag 0).
#'
#' @details
#' This function uses the Wiener–Khinchin theorem to compute autocovariance/autocorrelation via the FFT.
#' The FFT implementation assumes periodic extension; proper bias correction ensures equivalence to the sample autocovariance
#' definition used in R's `acf()` function.
#'
#' If `pad = TRUE`, time series are explicitly zero-padded to length `2^nextpow2(2*T - 1)` to avoid wrap-around effects and optimize FFT speed.
#' If `bias_correction = TRUE`, autocovariance estimates at lag k are rescaled by N / (N - k) to correct for the reduction in number of terms at larger lags.
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(500), nrow = 5, ncol = 100)  # 5 time series of length 100
#' result <- acf_fft(Y, Tt = 50, pad = TRUE, bias_correction = TRUE)
#' print(result$acor)
#' print(result$acov)
#'
#' @export

acf_fft <- function(Y,
                   Tt,
                   lag = NULL,
                   pad = TRUE,
                   bias_correction = TRUE) {


  Y <- check_dim(Y, Tt)
  Y <- demean_ts(Y, dim = 2)

  # zero pad
  if (pad){
    nfft <- 2^.nextpow2(2 * Tt - 1)
    Y_padded <- .pad_ts(Y, Tt, nfft = nfft)
  } else {
    Y_padded <- Y
    nfft <- Tt
  }

  # Rows = timeseries
  yfft <- (t(mvfft(t(Y_padded))))
  yspec <- yfft * Conj(yfft)

  acov_biased <- t(Re((mvfft(t(yspec), inverse = TRUE))))/ (nfft*nfft)
  acov_biased <- acov_biased[,1:Tt, drop = FALSE]

  if (bias_correction){
    lag_seq <- 0:(Tt-1)
    scale_vec <- matrix((1-lag_seq/Tt),nrow = 1, ncol = Tt)
    acov <- sweep(acov_biased, 2, scale_vec, "*")
  } else {
    acov <- acov_biased
  }

  # Or can alternatively be rowSums(abs(Y)^2)
  var0 <- acov[, 1, drop = FALSE]
  acor <- sweep(acov, 1, var0, "/")

  if (!is.null(lag)){
    acov <- acov[, 1:lag]
    acor <- acor[, 1:lag]
  }

  # Compute 95% confidence intervals for autocorrelations
  bnd <- (sqrt(2) * qnorm(0.975)) / sqrt(Tt)  # Assuming normality
  CI <- c(-bnd, bnd)

  # Return results as a list
  return(list(acor = acor,
              CI = CI,
              acov = acov))
}


#' Cross-Autocorrelation Function via Fast Fourier Transform (FFT)
#'
#' Computes full-lag cross-correlation sequences between all pairs of timeseries in `Y`
#' using the Fast Fourier Transform (FFT), efficiently handling multi-dimensional data.
#' The function supports optional zero-padding for performance and returns results
#' normalized to unit variance at lag 0 (Pearson correlation at lag 0).
#'
#' @param Y Numeric matrix. A matrix of size `I × Tt` where `I` is the number of timeseries (rows)
#' and `Tt` is the number of datapoints (columns). Each row is treated as a separate timeseries.
#' @param Tt Integer. Number of datapoints (columns) in `Y`, used for consistency check.
#' @param lag Optional integer. Maximum lag to compute. If `NULL` (default), computes up to `Tt - 1`.
#' If `lag = 0`, returns only the lag 0 correlation matrix (Pearson's correlation).
#' @param pad Logical (default `TRUE`). If `TRUE`, zero-pads each timeseries to the next power-of-two
#' length sufficient for efficient FFT and to avoid circular wrap-around effects.
#'
#' @return A list with:
#' \describe{
#'   \item{xC}{3D array of dimension `I × I × L`, where `L = 2 * (lag - 1) + 1`.
#'   Each slice `xC[, , k]` gives the cross-correlation matrix at lag `lidx[k]`.
#'   Symmetric property ensured: `xC[i, j, k] = xC[j, i, k]`. Diagonal entries contain autocorrelations.}
#'   \item{lidx}{Integer vector of lag indices corresponding to the third dimension of `xC`,
#'   ranging from `-(lag - 1)` to `lag - 1`.}
#' }
#'
#' @details
#' The cross-correlation at lag 0 (`lidx == 0`) exactly matches Pearson's correlation matrix.
#' Autocorrelations (diagonal entries) are normalized such that `xC[i, i, lag = 0] = 1`.
#'
#' If `pad = TRUE`, time series are zero-padded to length `2^nextpow2(2 * Tt - 1)` before FFT.
#'
#' The output `xC` is fully symmetric across rows and columns for all lags.
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(4 * 100), nrow = 4)
#' result <- xacf_fft(Y, Tt = 100)
#' dim(result$xC)  # should be 4 x 4 x 199 for default lags
#' result$lidx[which.max(result$lidx == 0)]  # lag 0
#'
#' # Verify symmetry at lag 0
#' lag0_idx <- which(result$lidx == 0)
#' all.equal(result$xC[, , lag0_idx], t(result$xC[, , lag0_idx]))
#'
#' # Check that diagonal autocorrelations at lag 0 equal 1
#' diag(result$xC[, , lag0_idx])
#'
#' @export
xacf_fft <- function(Y, Tt,
                     lag = NULL,
                     pad = TRUE) {

  Y <- check_dim(Y, Tt)
  Y <- demean_ts(Y, dim = 2)
  I <- nrow(Y)

  # Handle max lag argument
  if (!is.null(lag)) {
    mxL <- lag
    if (mxL == 0) mxL <- 1  # Handle lag = 0 as Pearson's correlation case
  } else {
    mxL <- Tt
  }

  # Determine FFT length: next power of two >= 2*T - 1
  # zero pad
  if (pad){
    nfft <- 2^.nextpow2(2 * Tt - 1)
    Y_padded <- .pad_ts(Y, Tt, nfft = nfft)
  } else {
    Y_padded <- Y
    nfft <- Tt
  }

  # Compute FFT along rows
  Yf <- t(mvfft(t(Y_padded)))

  # Initialize output
  mxLcc <- (mxL - 1) * 2 + 1
  xC <- array(0, dim = c(I, I, mxLcc))

  # Lag index vector
  lidx <- seq(-(mxL - 1), mxL - 1)

  # Upper triangle pairs only for efficiency
  idx <- which(upper.tri(matrix(1, I, I), diag = T), arr.ind = TRUE)

  for (k in seq_len(nrow(idx))) {
    i <- as.numeric(idx[k, 1])
    j <- as.numeric(idx[k, 2])

    # Cross-correlation via inverse FFT
    cxy <- fft(Yf[i, ] * Conj(Yf[j, ]), inverse = TRUE)/nfft
    cxy <- Re(cxy)

    # Normalize length and wrap-around
    if (mxL > 1){
      cxy <- c(cxy[(nfft - mxL + 2):nfft], cxy[1:mxL])
    } else {
      cxy <- cxy[1:mxL]
      }

    # Normalize by sqrt(Var_i * Var_j)
    norm_factor <- sqrt(sum(abs(Y[i, ])^2) * sum(abs(Y[j, ])^2))
    cxy <- cxy / norm_factor

    # Store result
    xC[i, j, ] <- cxy
  }

  if (I > 1) xC <- .symmetrise_slices(xC)

  return(list(xC = xC, lidx = lidx))
}


.symmetrise_slices <- function(x) {
  for (k in 1:dim(x)[3]) {
    x[, , k] <- x[, , k] + t(x[, , k]) - diag(diag(x[, , k]))
  }
  return(x)
}

.pad_ts <- function(Y,Tt, nfft = NULL){
  if (is.null(nfft)){
    nfft <- 2^.nextpow2(2 * Tt - 1)
    }

  pad_width <- nfft - Tt
  Y_padded <- cbind(Y, matrix(0, nrow(Y), pad_width))
  return(Y_padded)
}


.nextpow2 <- function(x) {
  ceiling(log2(x))
}



#' AR(1) Monte Carlo Estimation of Unbiased Variance
#'
#' This function estimates the AR(1) coefficient and computes the variance correction
#' factor for a time series dataset using a Monte Carlo simulation. It is adapted
#' from the FSLnets toolbox and designed for estimating global correction factors.
#'
#' @param Y Numeric matrix. A time series dataset where rows represent time points
#'   and columns represent variables (e.g., nodes or subjects).
#' @param T Integer. The number of time points in the time series.
#' @return A list containing:
#'   - `V`: Estimated variance of the null data correlations.
#'   - `Z`: Fisher-transformed correlation matrix scaled by the correction factor.
#'   - `P`: Two-tailed p-values for the transformed correlations.
#'   - `R2Zcrt`: Correction factor derived from the standard deviation of null data.
#'   - `arone`: Estimated AR(1) coefficients for each column in the input data.
#' @references
#' Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols.
#' "Effective Degrees of Freedom of the Pearson's Correlation Coefficient
#' under Serial Correlation." bioRxiv (2018): 453795.
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(500), nrow = 100, ncol = 5)  # Example time series dataset
#' T <- nrow(Y)
#' result <- AR1MC(Y, T)
#' print(result$V)
#'
#' @export
AR1MC <- function(Y, Tt, n_sim = 1000) {

  Y <- check_dim(Y, Tt)
  Y <- demean_ts(Y, dim = 2)

  arone <- est_rough_ar1(Y, Tt)
  arone0 <- median(arone)
  netmat <- cor(Y)

  # Create null data using the estimated AR(1) coefficient
  Yr <- sim_many_ts_same_ar1(Tt, arone0, n_sim)

  # Compute null correlations
  Yr_corr <- cor(Yr)

  # Extract upper triangular elements (excluding diagonal)
  IDX <- which(upper.tri(matrix(1, n_sim, n_sim)))
  Yrc <- Yr_corr[IDX]

  # Compute variance and correction factor
  V <- var(Yrc)
  R2Zcrt <- 1 / sd(atanh(Yrc))

  f2p <- fisher_pval_matrix(netmat, R2Zcrt)

  # Return results as a list
  return(list(V = V,
              Z = f2p$Z,
              P = f2p$P,
              R2Zcrt = R2Zcrt,
              arone = arone))
}
