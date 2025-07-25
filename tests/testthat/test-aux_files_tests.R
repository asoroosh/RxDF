
######################################
# test-demean_ts
######################################

test_that("demean_ts correctly demeans rows when dim=2", {
  Y <- matrix(1:12, nrow = 3, byrow = TRUE)  # 3 rows, 4 columns
  Y_demeaned <- demean_ts(Y, dim = 2)

  # Manual row means
  row_means <- rowMeans(Y)

  # Check result matches expected demeaned matrix
  expect_equal(Y_demeaned, Y - row_means)

  # Check row means are approximately zero
  expect_true(all(abs(rowMeans(Y_demeaned)) < 1e-8))
})

test_that("demean_ts correctly demeans columns when dim=1", {
  Y <- matrix(1:12, nrow = 3, byrow = FALSE)  # 3 rows, 4 columns
  Y_demeaned <- demean_ts(Y, dim = 1)

  # Check result matches scale()
  expect_equal(Y_demeaned, scale(Y, center = TRUE, scale = FALSE))

  # Check column means are approximately zero
  expect_true(all(abs(colMeans(Y_demeaned)) < 1e-8))
})

test_that("demean_ts throws error for invalid dim", {
  Y <- matrix(1:12, nrow = 3)
  expect_error(demean_ts(Y, dim = 3), "Invalid dim")
})

######################################
# test-check_dim
######################################

test_that("check_dim works for valid vector input", {
  Y <- 1:10
  result <- check_dim(Y, Tt = 10)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 10)
  expect_equal(nrow(result), 1)
})

test_that("check_dim works for properly oriented matrix", {
  Y <- matrix(1:30, nrow = 3, ncol = 10)  # 10 col = timepoints
  result <- check_dim(Y, Tt = 10)

  expect_identical(result, Y)  # Should remain unchanged
})

test_that("check_dim transposes if needed", {
  Y <- matrix(1:30, nrow = 10, ncol = 3)  # wrong orientation
  result <- check_dim(Y, Tt = 10)

  expect_equal(dim(result), c(3, 10))
})

test_that("check_dim errors for incorrect vector length", {
  Y <- 1:9
  expect_error(check_dim(Y, Tt = 10), "Length of vector Y must equal Tt")
})

test_that("check_dim errors for matrix with no dimension matching Tt", {
  Y <- matrix(1:12, nrow = 4, ncol = 3)  # neither dim == 10
  expect_error(check_dim(Y, Tt = 10), "Dimensions are wrong")
})


######################################
# fisher_pval_matrix
######################################

test_that("fisher_pval_matrix works on identity correlation matrix", {
  netmat <- diag(1, 5)  # perfect self-correlations on diagonal
  R2Zcrt <- 1

  result <- fisher_pval_matrix(netmat, R2Zcrt)
  Z <- result$Z
  P <- result$P

  # Check shapes
  expect_equal(dim(Z), c(5, 5))
  expect_equal(dim(P), c(5, 5))

  # Diagonal Z should be Inf because atanh(1) = Inf ⇒ p=0
  expect_true(all(is.infinite(diag(Z))))
  expect_true(all(diag(P) == 0))
})

test_that("fisher_pval_matrix correctly handles zero correlation", {
  netmat <- matrix(0, 4, 4)
  R2Zcrt <- 1

  result <- fisher_pval_matrix(netmat, R2Zcrt)

  expect_true(all(result$Z == 0))
  expect_true(all(result$P == 1))
})

test_that("fisher_pval_matrix output is bounded and symmetric", {
  set.seed(42)
  mat <- matrix(rnorm(100*5), 100, 5)
  corr_mat <- cor(mat)
  R2Zcrt <- sqrt(97)

  result <- fisher_pval_matrix(corr_mat, R2Zcrt)

  P <- result$P
  Z <- result$Z

  # Check shape
  expect_equal(dim(P), dim(corr_mat))
  expect_equal(dim(Z), dim(corr_mat))

  # Check symmetry
  expect_equal(P, t(P))
  expect_equal(Z, t(Z))

  # Check all P-values between 0 and 1
  expect_true(all(P >= 0 & P <= 1))
})

######################################
# test-SumMat
######################################

test_that("SumMat returns correct shape", {
  set.seed(1)
  Tt <- 10
  n <- 3
  Y <- matrix(rnorm(Tt * n), nrow = Tt, ncol = n)

  SM0 <- SumMat(Y, Tt)

  expect_equal(dim(SM0), c(n, n, Tt))
})

test_that("SumMat produces symmetric slices with zero diagonals", {
  set.seed(2)
  Tt <- 5
  n <- 4
  Y <- matrix(rnorm(Tt * n), nrow = Tt, ncol = n)

  SM0 <- SumMat(Y, Tt)

  for (t in 1:Tt) {
    expect_equal(SM0[, , t], t(SM0[, , t]), tolerance = 1e-8)
    expect_true(all(diag(SM0[, , t]) == 0))
  }
})

test_that("SumMat computes correct pairwise sums for simple case", {
  # Deterministic input
  Y <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
  # Y = [1 2 3; 4 5 6] (Tt=2, n=3)

  SM0 <- SumMat(Y, Tt = 2)

  # Check first timepoint manually
  expect_equal(SM0[, , 1], matrix(c(
    0, 1+2, 1+3,
    2+1, 0, 2+3,
    3+1, 3+2, 0
  ), nrow = 3, byrow = TRUE))

  # Check second timepoint manually
  expect_equal(SM0[, , 2], matrix(c(
    0, 4+5, 4+6,
    5+4, 0, 5+6,
    6+4, 6+5, 0
  ), nrow = 3, byrow = TRUE))
})

test_that("SumMat transposes correctly if needed", {
  set.seed(3)
  Tt <- 8
  n <- 3
  Y <- matrix(rnorm(Tt * n), nrow = n, ncol = Tt)  # transposed orientation

  SM0 <- SumMat(Y, Tt)

  expect_equal(dim(SM0), c(n, n, Tt))
})

######################################
# test-ProdMat
######################################

test_that("ProdMat returns correct shape", {
  set.seed(1)
  Tt <- 10
  n <- 3
  Y <- matrix(rnorm(Tt * n), nrow = Tt, ncol = n)

  PM0 <- ProdMat(Y, Tt)

  expect_equal(dim(PM0), c(n, n, Tt))
})

test_that("ProdMat produces symmetric slices with zero diagonals", {
  set.seed(2)
  Tt <- 5
  n <- 4
  Y <- matrix(rnorm(Tt * n), nrow = Tt, ncol = n)

  PM0 <- ProdMat(Y, Tt)

  for (t in 1:Tt) {
    expect_equal(PM0[, , t], t(PM0[, , t]), tolerance = 1e-8)
    expect_true(all(diag(PM0[, , t]) == 0))
  }
})

test_that("ProdMat computes correct pairwise products for simple case", {
  # Deterministic input
  Y <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
  # Y = [1 2 3; 4 5 6] (Tt=2, n=3)

  PM0 <- ProdMat(Y, Tt = 2)

  # Check first timepoint manually
  expect_equal(PM0[, , 1], matrix(c(
    0, 1*2, 1*3,
    2*1, 0, 2*3,
    3*1, 3*2, 0
  ), nrow = 3, byrow = TRUE))

  # Check second timepoint manually
  expect_equal(PM0[, , 2], matrix(c(
    0, 4*5, 4*6,
    5*4, 0, 5*6,
    6*4, 6*5, 0
  ), nrow = 3, byrow = TRUE))
})

test_that("ProdMat transposes input if needed", {
  set.seed(3)
  Tt <- 8
  n <- 3
  Y <- matrix(rnorm(Tt * n), nrow = n, ncol = Tt)  # transposed orientation

  PM0 <- ProdMat(Y, Tt)

  expect_equal(dim(PM0), c(n, n, Tt))
})

######################################
# test-make_toeplitz.R
######################################

test_that("make_toeplitz returns correct matrix dimensions", {
  sigX <- c(0.5, 0.3)
  Tt <- 5
  mat <- make_toeplitz(sigX, Tt)
  expect_equal(dim(mat), c(5, 5))
})

test_that("make_toeplitz diagonal entries are 1", {
  sigX <- c(0.8, 0.4, 0.2)
  Tt <- 6
  mat <- make_toeplitz(sigX, Tt)
  expect_true(all(diag(mat) == 1))
})

test_that("make_toeplitz fills correct off-diagonals", {
  sigX <- c(0.9, 0.7)
  Tt <- 4
  mat <- make_toeplitz(sigX, Tt)
  expect_equal(mat[1, 2], 0.9)
  expect_equal(mat[1, 3], 0.7)
  expect_equal(mat[1, 4], 0)  # Beyond sigX length ⇒ 0
})

test_that("make_toeplitz works when sigX is shorter than needed", {
  sigX <- c(0.6)
  Tt <- 4
  mat <- make_toeplitz(sigX, Tt)
  expect_equal(mat[1, 2], 0.6)
  expect_equal(mat[1, 3], 0)  # Default 0
  expect_equal(mat[1, 4], 0)
})


