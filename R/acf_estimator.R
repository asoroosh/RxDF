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
                   lg = NULL,
                   pad = TRUE,
                   bias_correction = TRUE) {


  Y <- check_dim(Y, Tt)
  Y <- demean_ts(Y, dim = 2)

  # zero pad
  if (pad){
    Y_padded <- .pad_ts(Y, Tt)
  } else {
    Y_padded <- Y
  }

  # Rows = timeseries
  yfft <- (t(mvfft(t(Y_padded))))
  yspec <- yfft * Conj(yfft)

  acov_biased <- t(Re((mvfft(t(yspec), inverse = TRUE))))/ (Tt*Tt)
  acov_biased <- acov_biased[,1:Tt]

  if (bias_correction){
    lag_seq <- 0:(Tt-1)
    acov <- acov_biased*(1-lag_seq/Tt)
  } else {
    acov <- acov_biased
  }

  # Or can alternatively be rowSums(abs(Y)^2)
  var0 <- acov[, 1]
  acor <- acov / var0

  if (!is.null(lg)){
    acov <- acov[, 1:lg]
    acor <- acor[, 1:lg]
  }

  # Compute 95% confidence intervals for autocorrelations
  bnd <- (sqrt(2) * qnorm(0.975)) / sqrt(Tt)  # Assuming normality
  CI <- c(-bnd, bnd)

  # Return results as a list
  return(list(acor = acor,
              CI = CI,
              acov = acov))
}


xacf_fft <- function(Y, Tt, lag = NULL, pad = TRUE) {

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
  }

  # Compute FFT along rows
  Yf <- t(mvfft(t(Y_padded)))
  mYf <- Mod(Yf)

  # Initialize output
  mxLcc <- (mxL - 1) * 2 + 1
  xC <- array(0, dim = c(I, I, mxLcc))

  # Lag index vector
  lidx <- seq(-(mxL - 1), mxL - 1)

  # Upper triangle pairs only for efficiency
  idx <- which(upper.tri(matrix(1, I, I)), arr.ind = TRUE)

  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]
    j <- idx[k, 2]

    # Cross-correlation via inverse FFT
    cxy <- Mod(fft(mYf[i, ] * Conj(mYf[j, ]), inverse = TRUE))
    cxy <- Re(cxy)

    # Normalize length and wrap-around
    cxy <- c(cxy[(nfft - mxL + 2):nfft], cxy[1:mxL])

    # Normalize by sqrt(Var_i * Var_j)
    norm_factor <- sqrt(sum(abs(Y[i, ])^2) * sum(abs(Y[j, ])^2))
    cxy <- cxy / norm_factor

    # Store result
    xC[i, j, ] <- cxy
  }

  # Symmetrize (lower triangle = transpose of upper triangle)
  xC <- xC + aperm(xC, c(2, 1, 3))

  # Autocorrelations (diagonal entries)
  for (i in 1:I) {
    cxx <- fft(mYf[i, ] * Conj(mYf[i, ]), inverse = TRUE)
    cxx <- Re(cxx)
    cxx <- c(cxx[(nfft - mxL + 2):nfft], cxx[1:mxL])
    cxx <- cxx / sum(abs(Y[i, ])^2)
    xC[i, i, ] <- cxx
  }

  return(list(xC = xC, lidx = lidx))
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

