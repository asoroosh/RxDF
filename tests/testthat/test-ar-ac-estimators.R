# tests/testthat/test-ar-ac-estimators.R
library(testthat)

test_that("est_rough_ar1 estimates AR(1) correctly for simulated data", {
  set.seed(123)
  phi <- 0.6
  n <- 1000
  m <- 5
  # Simulate m AR(1) series of length n
  Y <- replicate(m, arima.sim(model = list(ar = phi), n = n))

  est <- est_rough_ar1(Y, n)

  expect_length(est, m)
  expect_true(all(abs(est - phi) < 0.2))  # allow tolerance due to randomness
})

test_that("est_rough_ar1 works for single vector input", {
  set.seed(456)
  phi <- 0.8
  n <- 1000
  y <- as.vector(arima.sim(model = list(ar = phi), n = n))

  est <- est_rough_ar1(y, n)

  expect_length(est, 1)
  expect_true(abs(est - phi) < 0.2)
})

test_that("est_rough_ar1 returns numeric vector of correct length", {
  Y <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
  result <- est_rough_ar1(Y, 100)

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("est_rough_ar1 handles row-oriented input gracefully", {
  Y <- matrix(rnorm(100 * 4), ncol = 100, nrow = 4)  # rows = series (wrong orientation)

  # Should auto-correct orientation using check_dim() internally
  result <- est_rough_ar1(Y, 100)

  expect_length(result, 4)
})

######################################
# test-GenTsAC1
######################################

test_that("GenTsAC1 returns correct output length and type", {
  set.seed(123)
  ar1 <- 0.8
  Tt <- 100
  ts <- GenTsAC1(ar1, Tt)

  expect_length(ts, Tt)
  expect_type(ts, "double")
})

test_that("GenTsAC1 generates white noise for ar1 = 0", {
  set.seed(456)
  Tt <- 1000
  ts <- GenTsAC1(0, Tt)

  acf1 <- acf(ts, plot = FALSE, lag.max = 1)$acf[2]
  expect_true(abs(acf1) < 0.1)
})

test_that("GenTsAC1 generates highly autocorrelated series for ar1 = 0.9", {
  set.seed(789)
  Tt <- 1000
  ts <- GenTsAC1(0.9, Tt)

  acf1 <- acf(ts, plot = FALSE, lag.max = 1)$acf[2]
  expect_true(acf1 > 0.75)
})

test_that("GenTsAC1 errors on invalid ar1 values", {
  expect_error(GenTsAC1(-1.5, 100), "ar1")
  expect_error(GenTsAC1(1.5, 100), "ar1")
  expect_error(GenTsAC1("bad", 100), "ar1")
})

test_that("GenTsAC1 errors on invalid Tt values", {
  expect_error(GenTsAC1(0.5, -10), "Tt")
  expect_error(GenTsAC1(0.5, 0), "Tt")
  expect_error(GenTsAC1(0.5, 10.5), "Tt")
  expect_error(GenTsAC1(0.5, "bad"), "Tt")
})

######################################
# test-tukey_taper_me
######################################

test_that("tukey_taper_me returns correct shape", {
  set.seed(1)
  acs <- matrix(rnorm(5 * 20), nrow = 5, ncol = 20)
  result <- tukey_taper_me(acs, Tt = 20, M = 10)

  expect_equal(dim(result), dim(acs))
})

test_that("tukey_taper_me applies correct weights for simple case", {
  acs <- matrix(1, nrow = 1, ncol = 5)  # simple row of all 1s
  M <- 3
  Tt <- 5
  result <- tukey_taper_me(acs, Tt, M)

  # Compute expected weights manually
  w <- (1 + cos((1:M) * pi / M)) / 2
  expected <- c(w, 0, 0)

  expect_equal(as.numeric(result), expected, tolerance = 1e-8)
})

test_that("tukey_taper_me zeroes columns > M", {
  acs <- matrix(rnorm(3 * 10), nrow = 3, ncol = 10)
  M <- 4
  result <- tukey_taper_me(acs, Tt = 10, M = M)

  expect_true(all(result[, (M+1):10] == 0))
})

test_that("tukey_taper_me works when M = Tt (no zero columns)", {
  acs <- matrix(1, nrow = 2, ncol = 6)
  result <- tukey_taper_me(acs, Tt = 6, M = 6)

  expect_true(all(result[, 1:6-1] != 0))
})

test_that("tukey_taper_me works when M = 1 (only first column nonzero)", {
  acs <- matrix(1, nrow = 2, ncol = 6)
  result <- tukey_taper_me(acs, Tt = 6, M = 1)

  expect_true(all(result[, 2:6] == 0))
})

######################################
# test-curb_taper_me
######################################

test_that("curb_taper_me preserves values for M = Tt", {
  acs <- matrix(rnorm(2 * 8), nrow = 2, ncol = 8)
  expect_error(curb_taper_me(acs, Tt = 8, M = 8), "The curbing window is full length of timeseries. Tt = M.")
})

test_that("curb_taper_me zeroes all except first column for M = 1", {
  acs <- matrix(rnorm(4 * 5), nrow = 4, ncol = 5)
  result <- curb_taper_me(acs, Tt = 5, M = 1)

  expect_true(all(result[, 2:5] == 0))
  expect_equal(result[, 1], acs[, 1])
})

test_that("curb_taper_me works for single-row input", {
  acs <- matrix(1:10, nrow = 1)
  result <- curb_taper_me(acs, Tt = 10, M = 5)

  expect_equal(dim(result), c(1, 10))
  expect_equal(result[1, 1:5], 1:5)
  expect_true(all(result[1, 6:10] == 0))
})


######################################
# test-shrink_me.R
######################################

test_that("shrink_me returns all zeros when all values below confidence bound", {
  Tt <- 5
  acs <- rep(0.001, Tt)  # very small values ⇒ below confidence bound
  result <- shrink_me(acs, Tt)

  expect_equal(result, rep(0, Tt))
})

test_that("shrink_me truncates correctly at known break point", {
  Tt <- 6
  acs <- c(0.9, 0.85, 0.01, 0.01, 0.01, 0.01)  # falls below after lag 2
  result <- shrink_me(acs, Tt)

  expect_true(all(result[1:2] != 0))
  expect_true(all(result[3:6] == 0))
})

test_that("shrink_me preserves length and type", {
  Tt <- 8
  acs <- c(0.6, 0.3, 0.2, 0.1, 0.05, 0.01, 0.01, 0.01)
  result <- shrink_me(acs, Tt)

  expect_length(result, Tt)
  expect_type(result, "double")
})


######################################
# find_break_point.R
######################################

test_that("find_break_point returns 0 when all values below threshold", {
  Tt <- 10
  bnd <- qnorm(0.975) / sqrt(Tt)
  acs <- rep(bnd / 2, Tt)  # All well below threshold

  result <- find_break_point(acs, Tt)
  expect_equal(result, 0)
})

test_that("find_break_point returns Tt when all values above threshold", {
  Tt <- 10
  bnd <- qnorm(0.975) / sqrt(Tt)
  acs <- rep(bnd * 2, Tt)  # All well above threshold

  result <- find_break_point(acs, Tt)
  expect_equal(result, Tt)
})

test_that("find_break_point returns correct index when crossing occurs", {
  Tt <- 10
  bnd <- qnorm(0.975) / sqrt(Tt)
  acs <- c(rep(bnd * 2, 5), rep(bnd / 2, 5))  # First 5 above, then below

  result <- find_break_point(acs, Tt)
  expect_equal(result, 5)
})

test_that("find_break_point returns correct index when crossing at lag 3", {
  Tt <- 6
  bnd <- qnorm(0.975) / sqrt(Tt)
  acs <- c(bnd * 2, bnd * 2, bnd / 2, bnd / 2, bnd / 2, bnd / 2)

  result <- find_break_point(acs, Tt)
  expect_equal(result, 2)
})

test_that("find_break_point works for random acf with known structure", {
  set.seed(123)
  x <- rnorm(100)
  ac <- acf(x, plot = FALSE, lag.max = 20)$acf[-1]  # Exclude lag 0

  result <- find_break_point(ac, Tt = 20)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 20)
})


