# RxDF -- WIP WIP WIP

# Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation

## Highlights
* Autocorrelation biases the standard error of Pearson's correlation and breaks the variance-stabilising property of Fisher's transformation.
* Commonly used methods (see `mis` directory) to adjust correlation standard errors are themselves biased when true correlation is non-zero due to a confounding effect.
* We propose a “xDF” method to provide accurate estimates of the variance of Pearson’s correlation -- before or after Fisher’s transformation -- that considers auto-correlation of each time series as well as instantaneous and lagged cross-correlation.
* Accounting for the autocorrelation in resting-state functional connectivity considerably alters the graph theoretical description of human connectome.
* Severity of resting state fMRI autocorrelation varies systematically with region of interest size, and is heterogeneous over subjects.



## Table of contents
* [Introduction](#introduction)
* [Intallation](#Install)
* [Simulation Examples](#Example for Pearson's Correlation)
* [Real-world Examples](#Example for Pearson's Correlation)
* [Real-world Examples](#acf)


## Introduction <a name="introduction"></a>
Collection of scripts to implement the xDF method introduced in

*Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. "Effective Degrees of Freedom of the Pearson's Correlation Coefficient under Serial Correlation." bioRxiv (2018): 453795.*

*Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. "Effective degrees of freedom of the Pearson's correlation coefficient under autocorrelation." NeuroImage 199 (2019): 609-625.*

The `RxDF()` can be used to:
* Estimate the variance of Pearson's correlation coefficients
* Calculate the z-statistics maps of large scale correlation maps (e.g., functional connectivity)
* Estimate accurate p-values for such correlation coefficients


## Installation <a name="Install"></a>
To install RxDF you can clone and install the package using `devtools`;

```{r}
devtools::install_github('https://github.com/asoroosh/RxDF.git')
```

## Simulation Example for Pearson's Correlation <a name="Examples"></a>

Using simulated data, we demonstrate various scenarios of how to leverage xDF. To simulate correlated and/or autocorrelated data, you can leverage `RxDF::corrautocorr()` function. 

### Uncorrelated, White Time Series

Generate five time series that are uncorrelated and there is no autocorrelation in each timeseries. 

```
Y <- matrix(rnorm(1000 * 5), nrow = 5)
xDF_out <- RxDF::RxDF(Y, Tt = 1000)
```

You can explore the estimated z-scores by printing `xDF_out$stat$z`,

```
xDF_out$stat$z
           [,1]         [,2]       [,3]       [,4]         [,5]
[1,]  0.0000000 -0.596202179 -1.5305090 -0.1491905  0.124873392
[2,] -0.5962022  0.000000000 -1.2460740 -0.2263384 -0.002191868
[3,] -1.5305090 -1.246073982  0.0000000 -0.6070897 -0.174628166
[4,] -0.1491905 -0.226338431 -0.6070897  0.0000000  0.273672755
[5,]  0.1248734 -0.002191868 -0.1746282  0.2736728  0.000000000

```

and p-values by printing `xDF_out$stat$p`,

```
> xDF_out$stat$p
          [,1]      [,2]      [,3]      [,4]      [,5]
[1,] 0.0000000 0.5510402 0.1258908 0.8814033 0.9006238
[2,] 0.5510402 0.0000000 0.2127372 0.8209382 0.9982511
[3,] 0.1258908 0.2127372 0.0000000 0.5437914 0.8613718
[4,] 0.8814033 0.8209382 0.5437914 0.0000000 0.7843361
[5,] 0.9006238 0.9982511 0.8613718 0.7843361 0.0000000
```

Note that the diagonal for both p-values and z-scores are set to zero. 


xDF also generates the theoritical variance and estimated variance for each pairwise correlation. For theortical variance, you can print out `xDF_out$stat$TV`, 

```
> xDF_out$stat$TV
             [,1]         [,2]         [,3]         [,4]         [,5]
[1,] 0.0010000000 0.0009992894 0.0009953279 0.0009999555 0.0009999688
[2,] 0.0009992894 0.0010000000 0.0009969002 0.0009998975 0.0010000000
[3,] 0.0009953279 0.0009969002 0.0010000000 0.0009992632 0.0009999390
[4,] 0.0009999555 0.0009998975 0.0009992632 0.0010000000 0.0009998502
[5,] 0.0009999688 0.0010000000 0.0009999390 0.0009998502 0.0010000000
```
and for the estimated variance, you can print out `xDF_out$var_hat_rho`,

```
> xDF_out$var_hat_rho
             [,1]         [,2]         [,3]         [,4]         [,5]
[1,] 0.0000000000 0.0009992894 0.0009953279 0.0009999555 0.0009999688
[2,] 0.0009992894 0.0000000000 0.0009969002 0.0009998975 0.0010000000
[3,] 0.0009953279 0.0009969002 0.0000000000 0.0009992632 0.0009999390
[4,] 0.0009999555 0.0009998975 0.0009992632 0.0000000000 0.0009998502
[5,] 0.0009999688 0.0010000000 0.0009999390 0.0009998502 0.0000000000
```

Note that in this case the theoritical variance and the variance estimated by xDF are identical. That is because each time series is white and there is no correlation between time series.

### Correlated, White Time Series

Generate time series that are correlated with each other, but each time series is a white noice (i.e., each time series is not correlated with lagged version of itself)

```
set.seed(1)
Y <- RxDF::corrautocorr(mu = c(0, 0), 
                  sigR = 0.5,
                  sigC = list(c(0), c(0)),
                  Tt = 1000)
```

```
cor(t(Y))
          [,1]      [,2]
[1,] 1.0000000 0.5214198
[2,] 0.5214198 1.0000000
```

```
acf = RxDF::acf_fft(Y, 1000)

acf$acor[1,1:10]
 [1]  1.00000000 -0.01555711  0.01374823 -0.02649370

acf$acor[2,1:10]
 [1]  1.00000000 -0.01475334 -0.00053736  0.00884265
```

```
xDF_out = RxDF::RxDF(Y,1000)
-- Adaptive truncation.

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

```
set.seed(1)
Y <- RxDF::corrautocorr(mu = c(0, 0), 
                  sigR = 0,
                  sigC = list(c(0.6, 0.3, 0.2), c(0.5, 0.25, 0.1)),
                  Tt = 1000)
```

```
cor(t(Y))
            [,1]        [,2]
[1,]  1.00000000 -0.01662165
[2,] -0.01662165  1.00000000
```

```
acf = RxDF::acf_fft(Y, 1000)
acf$acor[1,1:10]
 [1]  1.000000000  0.591724895  0.296869664  0.193454466

acf$acor[2,1:10]
 [1]  1.000000000  0.507665217  0.237327446  0.085485073 
``` 


```
xDF_out <- RxDF::RxDF(Y, Tt = 1000)

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

```
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

```
acf = RxDF::acf_fft(Y, 1000)

acf$acor[1,1:10]
 [1]  1.000000000  0.591724895  0.296869664  0.193454466

acf$acor[2,1:10]
 [1]  1.00000000  0.49270660  0.25838924  0.12638006
``` 

```
xDF_out <- RxDF::RxDF(Y, Tt = 1000)

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

## Rapid estimate of autocorrelation function using xDF <a name="acf"></a>

xDF offers two key functions for estimation auto- and cross-correlations of time series. The function leverage Wiener–Khinchin theorem to estimate the acf functions. 

```
Y <- RxDF::corrautocorr(mu = c(0, 0), 
                  sigR = 0.5,
                  sigC = list(c(0.6, 0.3, 0.2), c(0.5, 0.25, 0.1)),
                  Tt = 1000)
```





