# dyngen 1.0.1

## NEW FEATURES

* `generate_experiment()`: Allow toggling on or off the mapping of the reference library size & CPM distribution.

* `wrap_dataset()`: Repurpose function to be able to choose the output format using a parameter (#28).

* `as_list()`: Added a function for converting the dyngen output to a simple list object.

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


# dyngen 0.4.0

## MAJOR CHANGES

* `as_dyno()`: Outputted `$counts` now contains counts of both spliced and unspliced reads, whereas
  `$counts_unspliced` and `$counts_spliced` contains separated counts.
  
* Added a docker container containing the necessary code to run a dyngen simulation.
  
## MINOR CHANGES

* Added logo to package.

* Clean up internal code, mostly to satisfy R CMD check.

## DOCUMENTATION

* Added two vignettes.

* Expanded the README.

# dyngen 0.3.0 (2020-04-06)

## NEW FEATURES

* Implement knockdown / knockouts / overexpression experiments.

* Implement better single-cell regulatory activity by determining
  the effect on propensity values after knocking out a transcription factor.
  
* Implement adding noise to the kinetic params of individual simulations.

* Kinetics (transcription rate, translation rate, decay rate, ...) are 
  based on Schwannhausser et al. 2011.

* Changed many parameter names to better explain its purpose.

## MINOR CHANGES

* Fix module naming of backbones derived from `backbone_branching()`.

* Allow to plot labels in `plot_simulation_expression()`.

* Improve `backbone_disconnected()` and `backbone_converging()`.

* Rename required columns in `backbone()` input data.

* Use `backbone_linear()` to make `backbone_cyclic()` randomised.

* Added a decay rate for pre-mRNAs as well.

* Kinetics: redefine the decay rates in terms of the half-life of these molecules.

* Only compute dimred if desired.

* Allow computing the propensity ratios as ground-truth for rna velocity.

## BUG FIXES

* Implement fix for double positives in `bblego` backbones.

* Fix graph plotting mixup of interaction effects (up/down).

* Made a fix to the computation of `feature_info$max_protein`.


# dyngen 0.2.1 (2019-07-17)

* MAJOR CHANGES: Custom backbones can be defined using backbone lego pieces. See `?bblego` for more information.

* MAJOR CHANGES: Splicing reactions have been reworked to better reflect biology.

# dyngen 0.2.0 (2019-07-12)

Complete rewrite from `dyngen` from the bottom up.
 
* OPTIMISATION: All aspects of the pipeline have been optimised towards execution time and end-user usability.

* OPTIMISATION: `dyngen` 0.2.0 uses `gillespie` 0.2.0, which has also been rewritten entirely in `Rcpp`,
  thereby improving the speed significantly.
  
* OPTIMISATION: The transcription factor propensity functions have been refactored to make it much more 
  computationally efficient.
  
* OPTIMISATION: Mapping a simulation to the gold standard is more automised and less error-prone.

* FEATURE: A splicing step has been added to the chain of reaction events.

# dyngen 0.1.0 (2017-04-27)

 * INITIAL RELEASE: a package for generating synthetic single-cell data from regulatory networks.
   Key features are:
   
   - The cells undergo a dynamic process throughout the simulation.
   - Many different trajectory types are supported.
   - `dyngen` 0.1.0 uses `gillespie` 0.1.0, a clone of `GillespieSSA` that is much less
     error-prone and more efficient than `GillespieSSA`.

# dyngen 0.0.1 (2016-04-04)

 * Just a bunch of scripts on a repository, which creates random networks using `igraph` and 
   generates simple single-cell expression data using `GillespieSSA`.
