#' Initial settings for simulating a dyngen dataset
#' 
#' @param num_cells The number of cells to sample.
#' @param num_tfs The number of transcription factors to generate.
#' @export
initialise_model <- function(
  num_cells = round(10 ^ runif(1, 2, 4)),
  num_tfs = round(10 ^ runif(1, 1, log10(200))),
  num_targets = round(10 ^ runif(1, 1, 4)),
  num_hks = round(10 ^ runif(1, 1, 3)),
  distance_metric = "pearson",
  backbone,
  tf_network_params = tf_network_random(),
  feature_network_params = feature_network_default(),
  kinetics_params = kinetics_default(),
  simulation_params = simulation_default(),
  gold_standard_params = gold_standard_default(),
  experiment_params = experiment_snapshot(),
  verbose = TRUE,
  download_cache_dir = NULL,
  num_cores = 8
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
    download_cache_dir,
    num_cores
  )
}
formals(initialise_model)$distance_metric <- dynutils::list_distance_methods()