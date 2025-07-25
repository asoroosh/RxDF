![R-CMD-check](https://github.com/sorooshafyouni/RxDF/actions/workflows/R-CMD-check.yaml/badge.svg)

# RxDF

# Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation

## Highlights
* Autocorrelation biases the standard error of Pearson's correlation and breaks the variance-stabilising property of Fisher's transformation.
* Commonly used methods (see `mis` directory) to adjust correlation standard errors are themselves biased when true correlation is non-zero due to a confounding effect.
* We propose a “xDF” method to provide accurate estimates of the variance of Pearson’s correlation — before or after Fisher’s transformation — that considers auto-correlation of each time series as well as instantaneous and lagged cross-correlation.
* Accounting for the autocorrelation in resting-state functional connectivity considerably alters the graph theoretical description of the human connectome.
* Severity of resting-state fMRI autocorrelation varies systematically with region of interest size, and is heterogeneous over subjects.

## Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Simulation Examples](#simulation-examples)
* [Real-world Examples](#real-world-examples)
* [Autocorrelation Estimation](#acf)

## Introduction <a name="introduction"></a>
This repository implements the xDF method introduced in:

*Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. "Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation." bioRxiv (2018): 453795.*

*Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. "Effective degrees of freedom of the Pearson's correlation coefficient under autocorrelation." NeuroImage 199 (2019): 609-625.*

The `RxDF()` function can be used to:
* Estimate the variance of Pearson's correlation coefficients
* Calculate z-statistics maps for large-scale correlation matrices (e.g., functional connectivity)
* Estimate accurate p-values for correlation coefficients, accounting for autocorrelation

### Why does it matter? <a name="why"></a>

The dependence between pairs of time series is commonly quantified by Pearson's correlation. However, if the time series are themselves dependent (i.e., exhibit temporal autocorrelation), the effective degrees of freedom (EDF) are reduced, the standard error of the sample correlation coefficient is biased, and Fisher's transformation fails to stabilise the variance. 

## Installation <a name="installation"></a>
To install RxDF you can clone and install the package using `devtools`:

```r
devtools::install_github('https://github.com/asoroosh/RxDF.git')
```

## Simulation Examples <a name="simulation-examples"></a>
The following examples demonstrate how to use `RxDF()` with simulated data.

### Uncorrelated, White Time Series

Generate five time series of length 1000, with no correlation or autocorrelation:


```r
set.seed(1)
Y <- matrix(rnorm(1000 * 5), nrow = 5)
xDF_out <- RxDF::RxDF(Y, Tt = 1000)
```

Z-scores and p-values of pairwise correlations:

```
xDF_out$stat$z
           [,1]       [,2]       [,3]         [,4]         [,5]
[1,]  0.0000000 -0.4768183  1.1914714  0.101135382  0.873463857
[2,] -0.4768183  0.0000000 -1.4767050  0.113717428  0.211015471
[3,]  1.1914714 -1.4767050  0.0000000  1.575651760 -0.710217772
[4,]  0.1011354  0.1137174  1.5756518  0.000000000 -0.002800874
[5,]  0.8734639  0.2110155 -0.7102178 -0.002800874  0.000000000

```

```
> xDF_out$stat$p
          [,1]      [,2]      [,3]      [,4]      [,5]
[1,] 0.0000000 0.6334915 0.2334686 0.9194430 0.3824103
[2,] 0.6334915 0.0000000 0.1397547 0.9094618 0.8328752
[3,] 0.2334686 0.1397547 0.0000000 0.1151061 0.4775691
[4,] 0.9194430 0.9094618 0.1151061 0.0000000 0.9977652
[5,] 0.3824103 0.8328752 0.4775691 0.9977652 0.0000000
```

The estimated and theoretical variances should match exactly due to independence. Note that the diagonal for both p-values and z-scores are set to zero. 


xDF also generates the theoritical variance and estimated variance for each pairwise correlation. For theortical variance, you can print out `xDF_out$stat$TV`, 

```
> xDF_out$stat$TV
             [,1]         [,2]         [,3]         [,4]         [,5]
[1,] 0.0010000000 0.0009995454 0.0009971655 0.0009999795 0.0009984755
[2,] 0.0009995454 0.0010000000 0.0009956498 0.0009999741 0.0009999109
[3,] 0.0009971655 0.0009956498 0.0010000000 0.0009950490 0.0009989819
[4,] 0.0009999795 0.0009999741 0.0009950490 0.0010000000 0.0010000000
[5,] 0.0009984755 0.0009999109 0.0009989819 0.0010000000 0.0010000000
```
and for the estimated variance, you can print out `xDF_out$var_hat_rho`,

```
> xDF_out$var_hat_rho
             [,1]         [,2]         [,3]         [,4]         [,5]
[1,] 0.0000000000 0.0009995454 0.0009971655 0.0009999795 0.0009984755
[2,] 0.0009995454 0.0000000000 0.0009956498 0.0009999741 0.0009999109
[3,] 0.0009971655 0.0009956498 0.0000000000 0.0009950490 0.0010087944
[4,] 0.0009999795 0.0009999741 0.0009950490 0.0000000000 0.0010000000
[5,] 0.0009984755 0.0009999109 0.0010087944 0.0010000000 0.0000000000
```

Note that in this case the theoritical variance and the variance estimated by xDF are identical. That is because each time series is white and there is no correlation between time series.

### Correlated, White Time Series

Generate time series that are correlated with each other, but each time series is a white noice (i.e., each time series is not correlated with lagged version of itself)

```r
set.seed(1)
Y <- RxDF::corrautocorr(mu = c(0, 0), 
                  sigR = 0.5,
                  sigC = list(c(0), c(0)),
                  Tt = 1000)
```

Check actual correlation and autocorrelation structure:


```
cor(t(Y))
          [,1]      [,2]
[1,] 1.0000000 0.5214198
[2,] 0.5214198 1.0000000
```

```r
acf = RxDF::acf_fft(Y, 1000)
```

Check the autocorrelation, for both timeseries, for the first four lags:

```
acf$acor[1,1:4]
 [1]  1.00000000 -0.01555711  0.01374823 -0.02649370

acf$acor[2,1:4]
 [1]  1.00000000 -0.01475334 -0.00053736  0.00884265
```

Run xDF:

```r
xDF_out = RxDF::RxDF(Y,1000)
```

Review z-scores and p-values:

```
xDF_out$stat$z
         [,1]     [,2]
[1,]  0.00000 18.28706
[2,] 18.28706  0.00000

xDF_out$stat$p
             [,1]         [,2]
[1,] 0.000000e+00 1.049162e-74
[2,] 1.049162e-74 0.000000e+00
```


### Uncorrelated, Autocorrelated Time Series

Generate time series that are uncorrelated, but each time series is highly autocorrelated (i.e., the autocorrelation function of each time series is non-zero)

```r
set.seed(1)
Y <- RxDF::corrautocorr(mu = c(0, 0), 
                  sigR = 0,
                  sigC = list(c(0.6, 0.3, 0.2), c(0.5, 0.25, 0.1)),
                  Tt = 1000)
```

Review the correlation between the two simulated time series:

```
cor(t(Y))
            [,1]        [,2]
[1,]  1.00000000 -0.01662165
[2,] -0.01662165  1.00000000
```

Run acutocorrelation function: 

```r
acf = RxDF::acf_fft(Y, 1000)
```

Check the autocorrelation function, for both timeseries, for four first lags:


```
acf$acor[1,1:4]
 [1]  1.000000000  0.591724895  0.296869664  0.193454466

acf$acor[2,1:4]
 [1]  1.000000000  0.507665217  0.237327446  0.085485073 
``` 

Run xDF:

```r
xDF_out <- RxDF::RxDF(Y, Tt = 1000)
```

Review z-scores and p-values:


```
xDF_out$stat$z
           [,1]       [,2]
[1,]  0.0000000 -0.3947852
[2,] -0.3947852  0.0000000

xDF_out$stat$p
          [,1]      [,2]
[1,] 0.0000000 0.6930014
[2,] 0.6930014 0.0000000

```

### Correlated and Autocorrelated Time Series

Generate time series that are correlated and also each time series is highly autocorrelated. 

```r
set.seed(1)
Y <- RxDF::corrautocorr(mu = c(0, 0), 
                  sigR = 0.5,
                  sigC = list(c(0.6, 0.3, 0.2), c(0.5, 0.25, 0.1)),
                  Tt = 1000)
```



```
cor(t(Y))
          [,1]      [,2]
[1,] 1.0000000 0.4843995
[2,] 0.4843995 1.0000000
```

```r
acf = RxDF::acf_fft(Y, 1000)
```

```
acf$acor[1,1:4]
 [1]  1.000000000  0.591724895  0.296869664  0.193454466

acf$acor[2,1:4]
 [1]  1.00000000  0.49270660  0.25838924  0.12638006
``` 

Run xDF:

```r
xDF_out <- RxDF::RxDF(Y, Tt = 1000)
```

Review z-scores and p-values:


```
xDF_out$stat$z
         [,1]     [,2]
[1,]  0.00000 11.52464
[2,] 11.52464  0.00000

xDF_out$stat$p
             [,1]         [,2]
[1,] 0.000000e+00 9.912127e-31
[2,] 9.912127e-31 0.000000e+00

```

## Real-data Example for Pearson's Correlation <a name="Examples"></a>

To be populated.

## Rapid estimate of autocorrelation function using xDF <a name="acf"></a>

xDF offers two key functions for estimation auto- and cross-correlations of time series. The function leverage Wiener–Khinchin theorem to estimate the acf functions. 

Simulate 256 time series, each of length 1000:

```r
library(microbenchmark)
set.seed(1)
M <- 1000
N <- 256
Ym <- matrix(rnorm(M * N), nrow = M, ncol = N)
```

For a single time series, benchmark `RxDF()` against R stat's `acf()`:

```r
bench_single_timeseries <- microbenchmark(
  acf_fft = RxDF::acf_fft(Ym[,1], M),
  acf_native_single = acf(Ym[,1], 
                            plot = FALSE, 
                            demean = TRUE, 
                            lag.max = M - 1),
  times = 100
)
```

Check the run time for single timeseries using R stat's `acf()` and `RxDF()`:

```
print(bench_single_timeseries)


Unit: microseconds
              expr     min       lq     mean  median       uq     max neval
           acf_fft 119.187 134.9515 148.8509 141.737 152.1100 406.433   100
 acf_native_single 471.295 501.0815 520.9714 512.049 521.8275 994.742   100
```

Evaluate the run time for multiple time series:

```r
acf_native_multiple <- function(Ym, M) {
  apply(Ym, 1, function(x) {
    acf(x, 
        plot = FALSE, 
        demean = TRUE, 
        lag.max = M - 1)
  })
}

bench_multiple_timeseries <- microbenchmark(
  acf_fft = RxDF::acf_fft(Ym, M),
  acf_native_multiple = acf_native_multiple(Ym, M),
  times = 100
)
print(bench_multiple_timeseries)
```

Check the run time for single timeseries using R stat's `acf()` and `RxDF()`:

```
Unit: milliseconds
                expr      min       lq     mean   median       uq      max neval
             acf_fft 18.35447 23.89995 34.25953 27.03169 32.56761 103.1304   100
 acf_native_multiple 78.04174 82.92463 91.41345 85.03917 88.40299 181.8708   100
```





