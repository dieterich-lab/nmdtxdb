
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nmdtx

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

nmdtx is the front-end for the NMD transcriptome project.

## Installation

You can install the development version of nmdtx like so:

``` r
devtools::install_github("https://github.com/dieterich-lab/nmd-app")
```

## Dev

### Resolving Renv::restore()

renv::restore() will failed because there are dependencies of ggtranscript
it cannot find. This can be solved by installing these with 
```
renv::install('bioc::GenomicRanges')
renv::restore()
```
As described in `deploy/Dockerfile`.

In addition, shared libraries from a conda enviroment my interfere with R, 
remeber to deactivate the enviroment once executing the Renv::update().

## Deployment

-   TODO
