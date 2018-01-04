context("Test data generation")

test_that("a full dataset can be generated", {
  params <- simple_params
  options(ncores = 1)
  
  model <- invoke(generate_model_from_modulenet, params$model)
  
  simulation <- invoke(simulate_multiple, params$simulation, model$system)
  
  simulation <- preprocess_simulation_for_gs(simulation, model, params$gs$smooth_window)
  gs <- invoke(extract_goldstandard, params$gs, simulation, model, preprocess=FALSE)
  
  experiment <- invoke(run_experiment, params$experiment, simulation, gs)
  
  pdf(tempfile())
  normalisation <- invoke(dynutils::normalise_filter_counts, params$normalisation, experiment$counts, verbose = TRUE)
  dev.off()
  
  task <- wrap_task(params, model, simulation, gs, experiment, normalisation)
})