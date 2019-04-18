#!/usr/bin/Rscript

library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 60,
    num_targets = 200,
    num_hks = 500,
    distance_metric = "pearson",
    backbone = backbone_bifurcating_cycle(),
    tf_network_params = tf_network_random(min_tfs_per_module = 1),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(total_time = 20, num_simulations = 8),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    num_cores = 8,
    download_cache_dir = "~/.cache/dyngen"
  )
dyngen:::complete_function(model, ".")
