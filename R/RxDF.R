#' Robust variance estimation for correlated time series network
#'
#' Computes a robust estimator of the variance of correlation coefficients
#' between multiple standardized time series, adjusting for serial correlation
#' and cross-correlation effects using an analytic formula.
#' Supports adaptive or fixed truncation strategies to curb influence from
#' high-lag autocorrelations.
#'
#' @param Y Numeric matrix of dimension `nn × Tt`, where rows represent time series
#' and columns represent timepoints.
#' @param Tt Integer. Number of timepoints (columns of `Y`).
#' @param truncate Either `"adaptive"` for adaptive truncation based on significance
#' of autocorrelations, or a numeric scalar specifying fixed truncation lag `M`.
#' Default is `"adaptive"`.
#' @param verbose Logical. If `TRUE` (default), prints progress and diagnostics.
#' @param TVflag Logical. If `TRUE` (default), enforces minimum variance equal to
#' textbook white noise variance when autocorrelation-adjusted variance falls below it.
#'
#' @return A list with elements:
#' \describe{
#'   \item{VarHatRho}{Matrix of robust variance estimates for all pairs of correlations
#'   (diagonal set to 0).}
#'   \item{Stat}{List containing:
#'     \describe{
#'       \item{z}{Matrix of Fisher-transformed z-statistics for each correlation.}
#'       \item{p}{Matrix of two-sided p-values corresponding to `z`.}
#'       \item{TV}{Matrix of textbook variance estimates under white noise assumption.}
#'       \item{EVE}{Indices of pairs where autocorrelation-adjusted variance was smaller
#'       than textbook variance (only if `TVflag = TRUE`).}
#'       \item{W2S}{Matrix of truncation lags used for each pair (only if `truncate = "adaptive"`).}
#'     }
#'   }
#' }
#'
#' @details
#' The method uses fast FFT-based computation of autocovariances and cross-covariances
#' and implements the analytic formula proposed in Afyouni et al. (2019) for estimating
#' variance inflation due to autocorrelation structure.
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(5 * 100), nrow = 5, ncol = 100)
#' result <- RxDF(Y, Tt = 100)
#' str(result)
#'
#' @export
RxDF <- function(Y,
                 Tt,
                 truncate = "adaptive",
                 verbose = TRUE,
                 TVflag = TRUE) {

  Y <- check_dim(Y, Tt)
  nn <- nrow(Y)

  # Standardize each row
  Y <- Y/apply(Y, 1, stats::sd)

  # Correlation matrix
  rho <- stats::cor(t(Y))
  diag(rho) <- 0

  # Autocorrelation matrix
  ac <- acf_fft(Y, Tt)$acor
  ac <- ac[, 2:(Tt-1)]  # remove lag-0 and final lag
  nLg <- Tt - 2

  # Cross-correlation
  xcf <- xacf_fft(Y, Tt)$xC
  xc_n <- xcf[, , 2:(Tt-1)]  # positive lags
  xc_p <- xcf[, , (Tt+1):(dim(xcf)[3] - 1)]  # negative lags

  W2S <- matrix(0, nn, nn)

  # Truncation
  if (!is.null(truncate)) {
    if (is.character(truncate) && truncate == "adaptive") {
      if (verbose) message("-- Adaptive truncation.")
      for (i in 1:nn) {
        for (j in 1:nn) {
          W2S[i, j] <- max(
            find_break_point(ac[i, ], nLg),
            find_break_point(ac[j, ], nLg)
          )
        }
      }
      for (i in 1:nn) {
        ac[i, ] <- shrink_me(ac[i, ], nLg)
        for (j in 1:nn) {
          xc_n[i, j, ] <- curb_taper_me(xc_n[i, j, ], nLg, W2S[i, j])
          xc_p[i, j, ] <- curb_taper_me(xc_p[i, j, ], nLg, W2S[i, j])
        }
      }
    } else if (is.numeric(truncate)) {
      M <- truncate
      if (verbose) message(paste0("-- Truncation with M = ", M))
      for (i in 1:nn) {
        ac[i, ] <- curb_taper_me(ac[i, ], nLg, M)
        for (j in 1:nn) {
          xc_n[i, j, ] <- curb_taper_me(xc_n[i, j, ], nLg, M)
          xc_p[i, j, ] <- curb_taper_me(xc_p[i, j, ], nLg, M)
        }
      }
    } else {
      stop("Available options are 'adaptive' truncation or numeric truncation level.")
    }
  }

  # Weight matrix
  wgt <- seq(nLg, 1)
  wgtm3 <- array(rep(wgt, each = nn * nn), dim = c(nn, nn, nLg))

  Tp <- Tt - 1

  SumMat_ac2 <- SumMat(ac^2, nLg)
  SumMat_ac <- SumMat(ac, nLg)
  ProdMat_ac <- ProdMat(ac, nLg)

  var_hat_rho <- (
    Tp * (1 - rho^2)^2 +
      rho^2 * apply(wgtm3 * (SumMat_ac2 + xc_p^2 + xc_n^2), c(1, 2), sum) -
      2 * rho * apply(wgtm3 * (SumMat_ac * (xc_p + xc_n)), c(1, 2), sum) +
      2 * apply(wgtm3 * (ProdMat_ac + (xc_p * xc_n)), c(1, 2), sum)
  ) / Tt^2

  TV <- (1 - rho^2)^2 / Tt
  Stat <- list(EVE = matrix(0, 0, 2))  # initialize empty matrix, to save memory allocations

  if (any(var_hat_rho < TV) && TVflag) {
    idx_ex <- which(var_hat_rho < TV, arr.ind = TRUE)
    var_hat_rho[idx_ex] <- TV[idx_ex]
    if (verbose) {
      message(paste0(nrow(idx_ex) - nn, " edges had variance smaller than textbook variance!"))
    }
    Stat$EVE <- idx_ex
  }
  Stat$TV <- TV

  diag(var_hat_rho) <- 0  # diagonal is "rubbish"

  rf <- .xdf_atanh(rho)
  sf <- var_hat_rho / (1 - rho^2)^2  # delta method
  rzf <- rf / sqrt(sf)
  diag(rzf) <- 0

  f_pval <- 2 * stats::pnorm(-abs(rzf))
  diag(f_pval) <- 0

  Stat$z <- rzf
  Stat$p <- f_pval
  Stat$W2S <- W2S
  Stat$rho <- rho

  return(list(stat = Stat,
              var_hat_rho = var_hat_rho))
}

.xdf_atanh <- function(rho) {
  if (is.matrix(rho)) {
    result <- apply(rho, c(1, 2), .safe_atanh)
    dimnames(result) <- dimnames(rho)
    diag(rho) <- 0 # Just to make sure we are going to follow the original behaviour
    return(result)
  } else if (is.vector(rho)) {
    return(sapply(rho, .safe_atanh))
  } else {
    return(.safe_atanh(rho))
  }
}

.safe_atanh <- function(x) {
  max_val <- .Machine$double.xmax
  if (x >= 1) {
    return(max_val)
  } else if (x <= -1) {
    return(-max_val)
  } else {
    return(atanh(x))
  }
}
