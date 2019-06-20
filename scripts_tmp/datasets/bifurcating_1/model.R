#!/usr/bin/Rscript

library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 500, num_tfs = 100, num_targets = 200, num_hks = 200,
    backbone = backbone_bifurcating(),
    tf_network_params = tf_network_random(min_tfs_per_module = 3),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(total_time = 10, num_simulations = 40),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen"
  )
complete_function(model, ".")
