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
* [Intall](#Install)


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
