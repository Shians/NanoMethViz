# NanoMethViz

<!-- badges: start -->
<!-- badges: end -->

NanoMethViz is a toolkit for visualising methylation data from Oxford Nanopore sequencing.

## Installation

You can install NanoMethViz from GitHub with:

``` r
remotes::install_github("shians/NanoMethViz", build_vignettes = TRUE)
```

## Example

This package currently works with data from nanopolish and f5c, to import your data please see the following vignette

``` r
vignette("ImportingData", package = "NanoMethViz")
```

An introductory example for plotting can be found in the package vignette:

``` r
vignette("Introduction", package = "NanoMethViz")
```

An example of the visualisation for Peg3

![](img/peg3_spaghetti.png)

## License

This project is licensed under Apache License, Version 2.0.
