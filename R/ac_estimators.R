
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
