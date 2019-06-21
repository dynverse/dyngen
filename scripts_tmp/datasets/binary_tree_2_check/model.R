#!/usr/bin/Rscript

library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 60,
    num_targets = 300,
    num_hks = 500,
    distance_metric = "pearson",
    backbone = backbone_binary_tree(),
    tf_network_params = tf_network_random(),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_default(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(total_time = 10, num_simulations = 32),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen"
  )
complete_function(model, ".")
