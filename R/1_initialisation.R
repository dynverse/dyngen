#' Initial settings for simulating a dyngen dataset
#' 
#' @param num_cells The number of cells to sample.
#' @param num_tfs The number of transcription factors to generate.
#' @export
initialise_model <- function(
  num_cells,
  num_tfs,
  num_targets,
  num_hks,
  distance_metric,
  backbone,
  tf_network_params = tf_network_random(),
  feature_network_params = feature_network_default(),
  kinetics_params = kinetics_default(),
  simulation_params = simulation_default(),
  gold_standard_params = gold_standard_default(),
  experiment_params = experiment_snapshot(),
  verbose = FALSE,
  num_cores = 1,
  download_cache_dir = NULL
) {
  distance_metric <- match.arg(distance_metric)
  
  lst(
    numbers = lst(
      num_cells,
      num_tfs,
      num_targets,
      num_hks,
      num_features = num_tfs + num_targets + num_hks,
      num_modules = nrow(backbone$module_info)
    ),
    distance_metric,
    backbone,
    tf_network_params,
    feature_network_params,
    kinetics_params,
    simulation_params,
    gold_standard_params,
    experiment_params,
    verbose,
    num_cores,
    download_cache_dir
  )
}
formals(initialise_model)$distance_metric <- dynutils::list_distance_methods()