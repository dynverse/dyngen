#' @export
initialise_model <- function(
  modulenet = modulenet_linear(),
  platform = platform_simple(),
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
    modulenet,
    platform,
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