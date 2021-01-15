context("generate_dataset")

test_that("generating a dataset with linear backbone", {
  backbone <- backbone_linear()
  
  init <- initialise_model(
    backbone = backbone,
    num_cells = 10,
    num_tfs = nrow(backbone$module_info),
    num_targets = 5,
    num_hks = 5,
    gold_standard_params = gold_standard_default(
      census_interval = 1,
      tau = .1
    ),
    simulation_params = simulation_default(
      ssa_algorithm = ssa_etl(tau = .1),
      experiment_params = simulation_type_wild_type(num_simulations = 1),
      census_interval = 1,
      compute_cellwise_grn = TRUE,
      compute_dimred = TRUE,
      compute_rna_velocity = TRUE
    ),
    verbose = FALSE
  )
  out <- generate_dataset(init, make_plots = TRUE)
  
  expect_is(out$plot, "ggplot")
  
  expect_true(dynwrap::is_wrapper_with_expression(out$dataset))
})


