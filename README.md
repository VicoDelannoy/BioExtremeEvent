
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BioExtremeEvent <img src="man/figures/package-sticker.png" align="right" style="float:right; height:120px;"/>

<!-- badges: start -->

[![Website](https://github.com/VicoDelannoy/BioExtremeEvent/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/VicoDelannoy/BioExtremeEvent/actions/workflows/pkgdown.yaml)
[![R CMD
Check](https://github.com/VicoDelannoy/BioExtremeEvent/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/VicoDelannoy/BioExtremeEvent/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<p align="left">

• <a href="#overview">Overview</a><br> •
<a href="#features">Features</a><br> •
<a href="#installation">Installation</a><br> •
<a href="#get-started">Get started</a><br> •
<a href="#long-form-documentations">Long-form documentations</a><br> •
<a href="#citation">Citation</a><br> •
<a href="#contributing">Contributing</a><br> •
<a href="#acknowledgments">Acknowledgments</a><br> •
<a href="#references">References</a>

</p>

## Overview

The R package `BioExtremeEvent` identifies and characterises an extreme event in  
time and space for a given GPS point (e.g. a sampling site) or for every pixel 
in an area.


## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
## Install < remotes > package (if not already installed) ----
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

## Install < BioExtremeEvent > from GitHub ----
remotes::install_github("VicoMarbec/BioExtremeEvent")
```

Then you can attach the package `BioExtremeEvent`:

``` r
library("BioExtremeEvent")
```

## Get started

For an overview of the main features of `BioExtremeEvent`, please read
the [Get
started](https://VicoMarbec.github.io/BioExtremeEvent/articles/BioExtremeEvent.html)
vignette.

## Long-form documentations

`BioExtremeEvent` provides **{{ NUMBER OF VIGNETTES }}** vignettes to
learn more about the package:

- the [Get
  started](https://VicoMarbec.github.io/BioExtremeEvent/articles/BioExtremeEvent.html)
  vignette describes the core features of the package
- **{{ LIST ADDITIONAL VIGNETTES }}**

## Citation

Please cite `BioExtremeEvent` as:

> Delannoy Victoria, Loiseau Nicolas, Villéger Sébastien, Cabrol Nicolas, Fièvre
Céleste (`r format(Sys.Date(), "%Y")`) BioExtremeEvent: An R package
to **characterise extreme event**. R package version 0.0.900. 
<https://github.com/VicoMarbec/BioExtremeEvent/>

## Contributing

All types of contributions are encouraged and valued. For more
information, check out our [Contributor
Guidelines](https://github.com/VicoMarbec/BioExtremeEvent/blob/main/CONTRIBUTING.md).

Please note that the `BioExtremeEvent` project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Acknowledgments

**{{ OPTIONAL SECTION }}**

## References

The dataset used in the exemple comes from :
Embury, O., Merchant, C.J., Good, S.A., Rayner, N.A., Høyer, J.L., Atkinson, C.,
Block, T., Alerskans, E., Pearson, K.J., Worsfold, M., McCarroll, N.,
Donlon, C., (2024). Satellite-based time-series of sea-surface temperature since
1980 for climate applications. Sci Data 11, 326. 
doi: https://doi.org/10.1038/s41597-024-03147-w
