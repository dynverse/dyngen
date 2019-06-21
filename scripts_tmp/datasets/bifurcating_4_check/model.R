#!/usr/bin/Rscript

library(dyngen)

set.seed(5)
model <- 
  initialise_model(
    num_cells = 10000, num_tfs = 200, num_targets = 4900, num_hks = 4900,
    backbone = backbone_bifurcating(),
    tf_network_params = tf_network_random(min_tfs_per_module = 5),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_default(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(total_time = 10, num_simulations = 100),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen"
  )
complete_function(model, ".")
