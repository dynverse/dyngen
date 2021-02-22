context("combine_models")

# backbone <- backbone_linear()
# 
# model_common <- initialise_model(
#   backbone = backbone,
#   num_cells = 10,
#   num_tfs = nrow(backbone$module_info),
#   num_targets = 5,
#   num_hks = 5,
#   gold_standard_params = gold_standard_default(
#     census_interval = 10,
#     tau = .5
#   ),
#   simulation_params = simulation_default(
#     ssa_algorithm = ssa_etl(tau = .5),
#     experiment_params = simulation_type_wild_type(num_simulations = 1),
#     census_interval = 2,
#     compute_cellwise_grn = TRUE,
#     compute_dimred = TRUE,
#     compute_rna_velocity = TRUE
#   ),
#   verbose = FALSE
# ) %>% 
#   generate_tf_network() %>% 
#   generate_feature_network()
# 
# model_a <- 
#   model_common %>% 
#   generate_kinetics() %>% 
#   generate_gold_standard() %>% 
#   generate_cells()

test_that("combine before experiment", {
  # reuse example_model for the sake of time
  model_a <- example_model
  model_a$experiment <- NULL
  
  model_aa <-
    combine_models(list(left = model_a, right = model_a)) %>% 
    generate_experiment()
  
  expect_is(model_aa$simulations$meta, "data.frame")
  expect_true(any(grepl("left_", model_aa$simulations$meta$from)))
  expect_true(any(grepl("right_", model_aa$simulations$meta$from)))
  
  expect_is(model_aa$experiment$cell_info, "data.frame")
  expect_true(any(grepl("left_", model_aa$experiment$cell_info$from)))
  expect_true(any(grepl("right_", model_aa$experiment$cell_info$from)))
  
  expect_is(model_aa$simulations$dimred, "matrix")
  expect_is(model_aa$experiment$cellwise_grn, "Matrix")
  expect_is(model_aa$experiment$rna_velocity, "Matrix")
})

test_that("combine after experiment", {
  # reuse example_model for the sake of time
  model_b <- example_model
  model_b$simulation_params$compute_cellwise_grn <- FALSE
  model_b$simulation_params$compute_rna_velocity <- FALSE
  model_b$experiment$cellwise_grn <- NULL
  model_b$experiment$rna_velocity <- NULL
  
  model_bb <- combine_models(list(left = model_b, right = model_b))
  
  expect_is(model_bb$simulations$meta, "data.frame")
  expect_true(any(grepl("left_", model_bb$simulations$meta$from)))
  expect_true(any(grepl("right_", model_bb$simulations$meta$from)))
  
  expect_is(model_bb$experiment$cell_info, "data.frame")
  expect_true(any(grepl("left_", model_bb$experiment$cell_info$from)))
  expect_true(any(grepl("right_", model_bb$experiment$cell_info$from)))
  
  expect_null(model_bb$experiment$cellwise_grn)
  expect_null(model_bb$experiment$rna_velocity)
})


