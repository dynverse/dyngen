generate_model_from_modulenet <- dynutils::inherit_default_params(list(add_targets_realnet, modulenet_to_genenet, generate_random_tree, generate_system), generate_model_from_modulenet)
generate_system <- dynutils::inherit_default_params(list(randomize_gene_kinetics), generate_system)

# run_experiment <- dynutils::inherit_default_params(c(take_experiment_cells, filter_expression, add_housekeeping_poisson), run_experiment)

extract_goldstandard <- dynutils::inherit_default_params(c(extract_references, get_milestone_paths, preprocess_simulation_for_gs), extract_goldstandard)
