generate_model_from_modulenet <- dynutils::inherit_default_params(c(add_targets_realnet, modulenet_to_genenet, generate_random_tree), generate_model_from_modulenet)
run_experiment <- dynutils::inherit_default_params(c(take_experiment_cells, filter_expression, add_housekeeping_poisson), run_experiment)
