
# dyngen

[![CRAN
Status](https://www.r-pkg.org/badges/version/dyngen)](https://cran.r-project.org/package=dyngen)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/dyngen)](https://cran.r-project.org/package=dyngen)
[![DOI](https://img.shields.io/badge/doi-10.1101/2020.02.06.936971-green)](https://doi.org/10.1101/2020.02.06.936971)
![R-CMD-check](https://github.com/dynverse/dyngen/workflows/R-CMD-check/badge.svg)<br><img src="man/figures/logo.png" align="right" />

dyngen is a novel, multi-modal simulation engine for studying dynamic
cellular processes at single-cell resolution. dyngen is more flexible
than current single-cell simulation engines, and allows better method
development and benchmarking, thereby stimulating development and
testing of novel computational methods.

A preprint is available on
[bioRxiv](https://doi.org/10.1101/2020.02.06.936971). Run
`citation("dyngen")` to obtain the corresponding citation information.
All source code for reproducing the results in this manuscript are
available on [GitHub](https://github.com/dynverse/dyngen_manuscript).

## System requirements

This package is supported for Linux, but should also work on Mac OS X
and Windows. It has been tested with [Github
Actions](https://github.com/dynverse/dyngen/actions?query=workflow%3AR-CMD-check)
for R 3.5, 3.6 and 4.0 on the following systems Ubuntu, Windows Server
and Mac OS X.

## Installation

## Step 1: Install from CRAN or GitHub

dyngen is available on CRAN, so you can install it with the following
command.

``` r
install.packages("dyngen")
```

If you would like to install the development version of dyngen from
GitHub instead, run the following command. Use at your own risk!

``` r
install.packages("remotes")
remotes::install_github("dynverse/dyngen@devel", dependencies = TRUE)
```

## Step 2: Configure host system

It’s recommended to let dyngen know where it can cache downloaded files
and how many cores the host system has. If you don’t perform these
steps, running dyngen simulations will take a lot longer than it needs
to.

To do so, start editing your Rprofile by running the following commands:

``` r
install.packages("usethis")
usethis::edit_r_profile()
```

Inside your Rprofile, add the following lines:

``` r
options(Ncpus = 8L) # change this to the number of cores in your system
options(dyngen_download_cache_dir = "~/.cache/dyngen")
```

## Vignettes

To learn about how dyngen, check out the example vignette below.
Expected execution time for rerunning the code is about 5 minutes.

-   [Advanced: Constructing a custom
    backbone](vignettes/advanced_constructing_backbone.md):  
    `vignette("advanced_constructing_backbone", package="dyngen")`
-   [Advanced: Running dyngen from a docker
    container](vignettes/advanced_run_dyngen_from_docker.md):  
    `vignette("advanced_run_dyngen_from_docker", package="dyngen")`
-   [Advanced: Simulating batch
    effects](vignettes/advanced_simulating_batch_effects.md):  
    `vignette("advanced_simulating_batch_effects", package="dyngen")`
-   [Advanced: Simulating a knockout
    experiment](vignettes/advanced_simulating_knockouts.md):  
    `vignette("advanced_simulating_knockouts", package="dyngen")`
-   [Getting started](vignettes/getting_started.md):  
    `vignette("getting_started", package="dyngen")`
-   [Showcase different backbones](vignettes/showcase_backbones.md):  
    `vignette("showcase_backbones", package="dyngen")`

## Latest changes

Check out `news(package = "dyngen")` or [NEWS.md](NEWS.md) for a full
list of changes.

<!-- This section gets automatically generated from NEWS.md -->

### Recent changes in dyngen 0.4.1

#### NEW FEATURES

-   `as_anndata()`: Added a function for converting the dyngen output to
    an anndata object.

-   The default number of cores used can be set by adding
    `options(Ncpus = ...)` to your Rprofile.

-   The default cache folder for dyngen can be set by adding
    `options(dyngen_download_cache_dir = ...)` to your Rprofile.

-   Combine simular models with different outputs using the
    `combine_models()` function.

#### MINOR CHANGES

-   `initialise_model()`: Change defaults of `num_cores` and
    `download_cache_dir` to `getOption("Ncpus")` and
    `getOption("dyngen_download_cache_dir")` respectively.

-   `as_dyno()`: Rename `wrap_dataset()` to `as_dyno()`.

#### BUG FIX

-   `as_dyno()`: Fix `drop = FALSE` bug when only one cell is being
    sampled.

#### DOCUMENTATION

-   Extended the vignettes with more examples.
