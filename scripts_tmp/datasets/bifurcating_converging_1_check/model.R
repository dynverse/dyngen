#!/usr/bin/Rscript

library(dyngen)

set.seed(10)
model <- 
  initialise_model(
    backbone = backbone_bifurcating_converging(),
    tf_network_params = tf_network_random(),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(total_time = 10, num_simulations = 32),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen"
  )
complete_function(model, ".")
