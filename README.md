
# dyngen

![R-CMD-check](https://github.com/dynverse/dyngen/workflows/R-CMD-check/badge.svg)<br><img src="man/figures/logo.png" align="right" />

dyngen is a novel, multi-modal simulation engine for studying dynamic
cellular processes at single-cell resolution. dyngen is more flexible
than current single-cell simulation engines, and allows better method
development and benchmarking, thereby stimulating development and
testing of novel computational methods.

A preprint is available on
[bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.06.936971v2).
Run `citation("dyngen")` to obtain the corresponding citation
information. All source code for reproducing the results in this
manuscript are available on
[GitHub](https://github.com/dynverse/dyngen_manuscript).

## System requirements

This package is supported for Linux, but should also work on Mac OS X
and Windows. It has been tested with [Github
Actions](https://github.com/dynverse/dyngen/actions?query=workflow%3AR-CMD-check)
for R 3.5, 3.6 and 4.0 on the following systems:

  - Mac OS X: Catalina (10.15.5)
  - Linux: Ubuntu Xenial (16.04.6)
  - Windows: Windows Server 2019 (10.0.17763)

## Installation

You can install dyngen by running the following command. This should
take no more than 10 minutes to install.

``` r
install.packages("remotes")
remotes::install_github("dynverse/dyngen", dependencies = TRUE)
```

To build the vignettes upon installation, run the following command
instead. This might take a while depending on on the computer this is
run on (\>30min).

``` r
remotes::install_github("dynverse/dyngen", dependencies = TRUE, build_vignettes = TRUE)
```

## Vignettes

To learn about how dyngen, check out the example vignette below.
Expected execution time for rerunning the code is about 5 minutes.

  - [Step-by-step example run](vignettes/example.md):  
    `vignette("example", package="dyngen")`
  - [Showcase different backbones](vignettes/showcase_backbones.md):  
    `vignette("showcase_backbones", package="dyngen")`

## Docker

To run dyngen in a docker container, run the following command.

``` sh
docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true -v `pwd`:/home/rstudio/workdir dynverse/dyngen
```

More information on running and building the docker container is
available [here](https://github.com/dynverse/dyngen/tree/master/docker).

## Latest changes

Check out `news(package = "dyngen")` or [NEWS.md](NEWS.md) for a full
list of changes.

<!-- This section gets automatically generated from NEWS.md -->

### Recent changes in dyngen 0.3.1

#### MAJOR CHANGES:

  - `wrap_dataset()`: Outputted `$counts` now contains counts of both
    spliced and unspliced reads, whereas `$counts_unspliced` and
    `$counts_spliced` contains separated counts.

  - Added a docker container containing the necessary code to run a
    dyngen simulation.

#### MINOR CHANGES:

  - Added logo to package.

  - Clean up internal code, mostly to satisfy R CMD check.

#### DOCUMENTATION:

  - Added two vignettes.

  - Expanded the README.

### Recent changes in dyngen 0.3.0 (2020-04-06)

#### NEW FEATURES:

  - Implement knockdown / knockouts / overexpression experiments.

  - Implement better single-cell regulatory activity by determining the
    effect on propensity values after knocking out a transcription
    factor.

  - Implement adding noise to the kinetic params of individual
    simulations.

  - Kinetics (transcription rate, translation rate, decay rate, …) are
    based on Schwannhausser et al. 2011.

  - Changed many parameter names to better explain its purpose.

#### MINOR CHANGES:

  - Fix module naming of backbones derived from `backbone_branching()`.

  - Allow to plot labels in `plot_simulation_expression()`.

  - Improve `backbone_disconnected()` and `backbone_converging()`.

  - Rename required columns in `backbone()` input data.

  - Use `backbone_linear()` to make `backbone_cyclic()` randomised.

  - Added a decay rate for pre-mRNAs as well.

  - Kinetics: redefine the decay rates in terms of the half-life of
    these molecules.

  - Only compute dimred if desired.

  - Allow computing the propensity ratios as ground-truth for rna
    velocity.

#### BUG FIXES:

  - Implement fix for double positives in `bblego` backbones.

  - Fix graph plotting mixup of interaction effects (up/down).

  - Made a fix to the computation of `feature_info$max_protein`.
