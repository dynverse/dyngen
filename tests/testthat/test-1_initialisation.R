context("1_initialisation")

backbone <- backbone(
  module_info = tibble(
    module_id = c("A", "B"),
    basal = c(0.2, 0.8),
    burn = c(TRUE, FALSE),
    independence = c(0.4, 0.5)
  ),
  module_network = tibble(
    from = c("A", "A"),
    to = c("A", "B"),
    effect = c(-1L, 1L),
    strength = c(10, 1),
    hill = c(2, 3)
  ),
  expression_patterns = tibble(
    from = c("burn", "start"),
    to = c("start", "end"),
    module_progression = c("+A", "+B"),
    start = c(TRUE, FALSE),
    burn = c(TRUE, FALSE),
    time = c(10, 10)
  )
)

tf_network_params <- tf_network_default()
feature_network_params <- feature_network_default()
kinetics_params <- kinetics_default()
gold_standard_params <- gold_standard_default()
simulation_params <- simulation_default()
experiment_params <- experiment_snapshot()

test_that("Testing normal use case of initialisation", {
  model <- initialise_model(
    backbone = backbone,
    num_cells = 100,
    num_tfs = 10,
    num_targets = 50,
    num_hks = 40,
    distance_metric = "euclidean",
    tf_network_params = tf_network_params,
    feature_network_params = feature_network_params,
    kinetics_params = kinetics_params,
    gold_standard_params = gold_standard_params,
    simulation_params = simulation_params,
    experiment_params = experiment_params,
    verbose = FALSE,
    download_cache_dir = "test_dir",
    num_cores = 2
  )
  
  simulation_params$burn_time <- 12
  simulation_params$total_time <- 24
  
  expect_equal(model$backbone, backbone)
  expect_equal(model$numbers, list(num_cells = 100, num_tfs = 10, num_targets = 50, num_hks = 40, num_features = 100, num_modules = 2))
  expect_equal(model$distance_metric, "euclidean")
  expect_equal(model$tf_network_params, tf_network_params)
  expect_equal(model$feature_network_params, feature_network_params)
  expect_equal(model$kinetics_params, kinetics_params)
  expect_equal(model$gold_standard_params, gold_standard_params)
  expect_equal(model$simulation_params, simulation_params)
  expect_equal(model$experiment_params, experiment_params)
  expect_false(model$verbose)
  expect_equal(model$download_cache_dir, "test_dir")
  expect_equal(model$num_cores, 2)
})