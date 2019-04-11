#' @export
initialise_model <- function(
  num_cells,
  num_features,
  pct_tfs,
  pct_hks,
  dist_metric,
  modulenet = modulenet_linear(),
  tfgen_params = tfgen_random(),
  networkgen_params = networkgen_realnet_sampler(),
  simulation_setup_params = simulation_setup_default(),
  simulation_params = simulation_default(),
  goldstandard_params = goldstandard_default(),
  experiment_sampler = experiment_sampler_snapshot(),
  verbose = FALSE,
  num_cores = 1
) {
  dist_metric <- match.arg(dist_metric)
  assert_that(pct_tfs + pct_hks <= 1)
  
  lst(
    numbers = lst(
      num_cells,
      num_features,
      pct_tfs,
      num_tfs = num_features * pct_tfs,
      num_hks = num_features * pct_hks,
      num_targets = num_features * (1 - pct_tfs - pct_hks),
      num_modules = nrow(modulenet$module_info)
    ),
    dist_metric,
    modulenet,
    tfgen_params,
    networkgen_params,
    simulation_setup_params,
    simulation_params,
    goldstandard_params,
    experiment_sampler,
    verbose,
    num_cores
  )
}
formals(initialise_model)$dist_metric <- dynutils::list_distance_metrics()