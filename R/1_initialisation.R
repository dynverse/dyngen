#' @export
initialise_model <- function(
  num_cells,
  num_features,
  pct_tfs,
  modulenet = modulenet_linear(),
  tfgen_params = tfgen_random(),
  networkgen_params = networkgen_realnet_sampler(),
  simulation_setup_params = simulation_setup_default(),
  simulation_params = simulation_default(),
  goldstandard_params = goldstandard_default(),
  experiment_sampler = experiment_sampler_snapshot(),
  normalisation_params = normalisation_default(),
  verbose = FALSE,
  num_cores = 1
) {
  lst(
    numbers = lst(
      num_cells,
      num_features,
      pct_tfs,
      num_tfs = num_features * pct_tfs,
      num_targets = num_features * (1 - pct_tfs),
      num_modules = nrow(modulenet$module_info)
    ),
    modulenet,
    tfgen_params,
    networkgen_params,
    simulation_setup_params,
    simulation_params,
    goldstandard_params,
    experiment_sampler,
    normalisation_params,
    verbose,
    num_cores
  )
}