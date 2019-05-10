#!/usr/bin/Rscript

library(dyngen)

set.seed(4)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 60,
    num_targets = 400,
    num_hks = 500,
    distance_metric = "pearson",
    backbone = backbone_trifurcating(),
    tf_network_params = tf_network_random(min_tfs_per_module = 3),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(total_time = 10, num_simulations = 32),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    num_cores = 1,
    download_cache_dir = "~/.cache/dyngen"
  )
dyngen:::complete_function(model, ".")
