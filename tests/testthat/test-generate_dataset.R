context("generate_dataset")

test_that("generating a dataset with linear backbone", {
  backbone <- backbone_linear()
  
  init <- initialise_model(
    backbone = backbone,
    num_cells = 500,
    num_tfs = 100,
    num_targets = 50,
    num_hks = 25,
    simulation_params = simulation_default(
      ssa_algorithm = ssa_etl(tau = .1),
      census_interval = 10
    ),
    verbose = FALSE
  )
  out <- generate_dataset(init, make_plots = TRUE)
  
  expect_is(out$plot, "ggplot")
  
  expect_true(dynwrap::is_wrapper_with_expression(out$dataset))
})


