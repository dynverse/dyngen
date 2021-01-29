context("generate_dataset")

backbone <- backbone_linear()

init <- initialise_model(
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
    compute_cellwise_grn = TRUE,
    compute_dimred = TRUE,
    compute_rna_velocity = TRUE
  ),
  verbose = FALSE
)

out <- generate_dataset(init, make_plots = TRUE)

test_that("generating a dataset with linear backbone", {
  expect_is(out$plot, "ggplot")
  
  expect_true(dynwrap::is_wrapper_with_expression(out$dataset))

    # test converting to SCE
  sce <- as_SCE(out$model)
  expect_equal(dim(sce), dim(t(out$dataset$counts)))
})

# skip if not rcannood because anndata is probably not installed
skip_on_cran()
skip_if_not(Sys.info()[["user"]] == "rcannood")

test_that("generating converting to anndata", {
  ad <- as_anndata(out$model)
  expect_equal(dim(ad$X), dim(out$dataset$counts))
})


