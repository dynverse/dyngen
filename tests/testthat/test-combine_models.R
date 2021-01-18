context("combine_models")

test_that("combining multiple models works", {
  backbone <- backbone_linear()
  
  model_common <- initialise_model(
    backbone = backbone,
    num_cells = 10,
    num_tfs = nrow(backbone$module_info),
    num_targets = 5,
    num_hks = 5,
    gold_standard_params = gold_standard_default(
      census_interval = 10,
      tau = .5
    ),
    simulation_params = simulation_default(
      ssa_algorithm = ssa_etl(tau = .5),
      experiment_params = simulation_type_wild_type(num_simulations = 1),
      census_interval = 2,
      compute_cellwise_grn = FALSE, # todo: set this to true
      compute_dimred = TRUE,
      compute_rna_velocity = FALSE # todo: set this to true
    ),
    verbose = FALSE
  ) %>% 
    generate_tf_network() %>% 
    generate_feature_network()
  
  model_a <- 
    model_common %>% 
    generate_kinetics() %>% 
    generate_gold_standard() %>% 
    generate_cells()
  
  model_b <- 
    model_common %>% 
    generate_kinetics() %>% 
    generate_gold_standard() %>% 
    generate_cells()
  
  model_ab <-
    combine_models(list(left = model_a, right = model_b)) %>% 
    generate_experiment()
  
  expect_is(model_ab$simulations$meta, "data.frame")
  expect_true(any(grepl("left_", model_ab$simulations$meta$from)))
  expect_true(any(grepl("right_", model_ab$simulations$meta$from)))
})


