# dyngen 1.0.0

This version mostly upgrades dyngen's ease-of-use, such as better vignettes, conversion functions for working with dyngen datasets in other packages, and more useful ways of specifying platform-specific parameters (i.e. number of cores and cache location). Perhaps more excitingly, the dyngen documentation is more readable online at [https://dyngen.dynverse.org](https://dyngen.dynverse.org)!

## NEW FEATURES

* `as_anndata()`: Added a function for converting the dyngen output to an anndata object.

* `as_sce()`: Added a function for converting the dyngen output to an SingleCellExperiment object.

* `as_seurat()`: Added a function for converting the dyngen output to a Seurat object.

* The default number of cores used can be set by adding `options(Ncpus = ...)` to your Rprofile.

* The default cache folder for dyngen can be set by adding `options(dyngen_download_cache_dir = ...)` to your Rprofile.

* Combine similar models with different outputs using the `combine_models()` function.

* Store the timings throughout the dyngen execution. Extract the timings from a model using `get_timings()`.

## MAJOR CHANGES

* `generate_experiment()`: Map count density of reference dataset to simulation expression before sampling molecules.

## MINOR CHANGES

* `initialise_model()`: Change defaults of `num_cores` and `download_cache_dir`
  to `getOption("Ncpus")` and `getOption("dyngen_download_cache_dir")` respectively.
  
* `as_dyno()`: Rename `wrap_dataset()` to `as_dyno()`.

* `generate_experiment()`: Drastically speed up sampling of molecules.

## BUG FIX

* `as_dyno()`: Fix `drop = FALSE` bug when only one cell is being sampled.

* Removed names from feature ids in feature info (`unname(model$feature_info$feature_id)`). Thanks @milanmlft!

## DOCUMENTATION

* Extended vignettes:
  - Advanced: Simulating batch effects
  - Advanced: Simulating a knockout experiment
  - Advanced: Running dyngen from a docker container
  - Advanced: Constructing a custom backbone
  - Advanced: Tweaking parameters
  - Advanced: Comparison of characteristic features between dyngen and reference datasets


## Test environments
* local Fedora 32 install (R 4.0)
* ubuntu 16.04 (with Github Actions; R 3.3, 3.4, 3.5, 3.6, release)
* Mac OS X (with Github Actions; R release)
* Windows (with Github Actions; R release)
* win-builder (oldrelease, release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes
