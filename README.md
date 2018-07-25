dyngen
================

[![Build Status](https://travis-ci.org/dynverse/dyngen.svg)](https://travis-ci.org/dynverse/dyngen) [![codecov](https://codecov.io/gh/dynverse/dyngen/branch/master/graph/badge.svg)](https://codecov.io/gh/dynverse/dyngen)

A package to generating synthetic single-cell data starting from regulatory networks. The data is generated in several steps: ![generation\_overview](man/figures/README_generation_overview_v1.svg)

``` r
library(tidyverse)
library(dyngen)

params <- simple_params
options(ncores = 1)
```

Generating a regulatory network and a system of differential equations

``` r
model <- invoke(generate_model_from_modulenet, params$model)
```

    ## [1] "Generating main network"
    ## [1] "Sampling targets"
    ## [1] "Randomizing network"
    ## [1] "Generating system"

``` r
plot_net(model)
```

![](man/figures/README_model-1.png)

Simulating this system

``` r
simulation <- invoke(simulate_multiple, params$simulation, model$system)
plot_simulation_space_time(simulation)
```

![](man/figures/README_simulate-1.png)

Extracting the gold standard trajectory from the simulations

``` r
gs <- invoke(extract_goldstandard, params$gs, simulation, model)
```

    ## [1] "Preprocessing"
    ## [1] "Extracting milestone paths"
    ## [1] "Processing operations"
    ## [1] "Extracting references"
    ## [1] "Mapping simulations onto reference"
    ## [1] "Postprocessing"

``` r
plot_goldstandard_heatmap(simulation, gs)
```

![](man/figures/README_extract_goldstandard-1.png)![](man/figures/README_extract_goldstandard-2.png)

Simulating the single-cell RNA-seq experiment

``` r
experiment <- invoke(run_experiment, params$experiment, simulation, gs)
```

Normalising the data

``` r
normalisation <- invoke(dynnormaliser::normalise_filter_counts, params$normalisation, experiment$counts)
```

Wrapping in a dynwrap dataset

``` r
task <- wrap_dyngen_dataset("readme_dataset", params, model, simulation, gs, experiment, normalisation)
```
