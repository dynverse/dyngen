
# dyngen

[![CRAN
Status](https://www.r-pkg.org/badges/version/dyngen)](https://cran.r-project.org/package=dyngen)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/dyngen)](https://cran.r-project.org/package=dyngen)
[![DOI](https://img.shields.io/badge/doi-10.1101/2020.02.06.936971-green)](https://doi.org/10.1101/2020.02.06.936971)
![R-CMD-check](https://github.com/dynverse/dyngen/workflows/R-CMD-check/badge.svg)
[![Coverage
Status](https://codecov.io/gh/dynverse/dyngen/branch/master/graph/badge.svg)](https://codecov.io/gh/dynverse/dyngen?branch=master)
<br><img src="man/figures/logo.png" align="right" />

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

## Installation

dyngen should work straight out of the CRAN box by running
`install.packages("dyngen")`. Having said that, there are a few
recommended steps that will make dyngen work even better. Check the
installation vignette below for more information!

## Vignettes

To learn about how dyngen, check out the example vignette below.
Expected execution time for rerunning the code is about 5 minutes.

-   [Advanced: Constructing a custom
    backbone](vignettes/advanced_constructing_backbone.md)  
    `vignette("advanced_constructing_backbone", package="dyngen")`
-   [Advanced: Running dyngen from a docker
    container](vignettes/advanced_run_dyngen_from_docker.md)  
    `vignette("advanced_run_dyngen_from_docker", package="dyngen")`
-   [Advanced: Simulating batch
    effects](vignettes/advanced_simulating_batch_effects.md)  
    `vignette("advanced_simulating_batch_effects", package="dyngen")`
-   [Advanced: Simulating a knockout
    experiment](vignettes/advanced_simulating_knockouts.md)  
    `vignette("advanced_simulating_knockouts", package="dyngen")`
-   [Getting started](vignettes/getting_started.md)  
    `vignette("getting_started", package="dyngen")`
-   [Installation instructions](vignettes/installation.md)  
    `vignette("installation", package="dyngen")`
-   [On runtime and scalability](vignettes/scalability_and_runtime.md)  
    `vignette("scalability_and_runtime", package="dyngen")`
-   [Showcase different backbones](vignettes/showcase_backbones.md)  
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

-   Store the timings throughout the dyngen execution. Extract the
    timings from a model using `get_timings()`.

#### MINOR CHANGES

-   `initialise_model()`: Change defaults of `num_cores` and
    `download_cache_dir` to `getOption("Ncpus")` and
    `getOption("dyngen_download_cache_dir")` respectively.

-   `as_dyno()`: Rename `wrap_dataset()` to `as_dyno()`.

-   `generate_experiment()`: Drastically speed up sampling of moleciles.

#### BUG FIX

-   `as_dyno()`: Fix `drop = FALSE` bug when only one cell is being
    sampled.

#### DOCUMENTATION

-   Extended vignettes:
    -   Advanced: Simulating batch effects
    -   Advanced: Simulating a knockout experiment
    -   Advanced: Running dyngen from a docker container
    -   Advanced: Constructing a custom backbone
