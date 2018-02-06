context("Test data generation")

test_that("a full dataset can be generated", {
  params <- simple_params
  options(ncores = 1)
  
  model <- invoke(generate_model_from_modulenet, params$model)
  pdf(tempfile())
  plot_model(model)
  dev.off()
  
  simulation <- invoke(simulate_multiple, params$simulation, model$system)
  simulation <- preprocess_simulation_for_gs(simulation, model, params$gs$smooth_window)
  plot_simulation(simulation)
  
  gs <- invoke(extract_goldstandard, params$gs, simulation, model, preprocess=FALSE)
  pdf(tempfile())
  plot_goldstandard(simulation, gs)
  dev.off()
  
  experiment <- invoke(run_experiment, params$experiment, simulation, gs)
  pdf(tempfile())
  plot_experiment(experiment)
  dev.off()
  
  pdf(tempfile())
  normalisation <- invoke(dynnormaliser::normalise_filter_counts, params$normalisation, experiment$counts, verbose = TRUE)
  plot_normalisation(experiment, normalisation)
  dev.off()
  
  task <- wrap_task(params, model, simulation, gs, experiment, normalisation)
})
