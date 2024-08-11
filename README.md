
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flashfmZero

FlashfmZero is a computationally efficient approach to simultaneously
fine-map signals in any number of uncorrelated quantitative traits
(i.e. zero correlation), as may result from latent traits estimated from
factor analysis using a varimax rotation.

For correlated traits, the original
[flashfm](https://www.nature.com/articles/s41467-021-26364-y)
multi-trait fine-mapping method should be used, which is available in
this package for convenience.

Flashfm and flashfmZero output trait-specific results,leveraging
information between traits; for each trait, credible sets, SNP marginal
posterior probabilities of causality (MPP), and multi-SNP model
posterior probabilities (PP) are output.

For more details, please see:

F Zhou, WJ Astle, AS Butterworth, JL Asimit. (2024). Improved genetic
discovery and fine-mapping resolution of high-dimensional traits through
multivariate analyses of latent factors. *Pre-print*.

Website available at: <https://jennasimit.github.io/flashfmZero/>

We have applied these methods to GWAS results from 99 raw blood cell
traits and their 25 latent factors in the [INTERVAL
cohort](https://doi.org/10.1186/1745-6215-15-363). Our analysis scripts
are available here:

<https://github.com/fz-cambridge/flashfmZero-INTERVAL-analysis>

## System Requirements

flashfmZero could be installed with ease on versions of R \> 4.2.1. If
installing on a Windows machine, Rtools must be installed. Installation
time is estimated as 2 minutes.

## Installation Guide

### Short version

``` r
# install.packages("devtools")
devtools::install_github("jennasimit/flashfmZero")
```

### Longer version (if above fails)

The following packages from CRAN and Bioconductor are required:

``` r
install.packages("parallel")
install.packages("Matrix")
install.packages("gtools")
install.packages("rlist")
```

NB: Must have a Java JDK installed in order to run R2BGLiMS (it is
pre-installed in the flashfmZero package, no need install it
separately). This is only needed if you need to run single-trait
fine-mapping using JAM. If single-trait fine-mapping results are
available, then it is not necessary to have Java JDK installed.

``` r
remotes::install_github("jennasimit/flashfmZero")
library(flashfmZero)
```
