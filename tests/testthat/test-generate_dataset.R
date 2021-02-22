context("generate_dataset")

backbone <- backbone_linear()
total_time <- simtime_from_backbone(backbone)

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
    experiment_params = simulation_type_wild_type(num_simulations = 2),
    census_interval = 20,
    compute_cellwise_grn = TRUE,
    compute_dimred = TRUE,
    compute_rna_velocity = TRUE
  ),
  verbose = FALSE
)

out <- generate_dataset(init, make_plots = TRUE)

test_that("output contains the desired objects", {
  dataset <- out$dataset
  expected_names <- c(
    "cell_ids", "feature_ids", "counts", "counts_spliced", "counts_unspliced",
    "counts_protein", "expression", "cell_info", "feature_info", "milestone_ids",
    "milestone_network", "progressions", "milestone_percentages", "dimred", 
    "dimred_segment_points", "dimred_segment_progressions", "regulatory_network",
    "regulatory_network_sc", "regulators", "targets", "rna_velocity"
  )
  expect_true(all(expected_names %in% names(dataset)))
  expect_true(all(sapply(expected_names, function(nam) !is.null(dataset[[nam]]))))
  
  # simple checks
  expect_is(dataset$counts, "Matrix")
  expect_equal(dataset$cell_ids, rownames(dataset$counts))
  expect_equal(dataset$feature_ids, colnames(dataset$counts))
})

# skip on cran because some dependencies might be missing
skip_on_cran()

test_that("generating a dataset with linear backbone", {
  expect_is(out$plot, "ggplot")
  
  # test converting to SCE
  sce <- as_sce(out$model)
  expect_equal(dim(sce), dim(t(out$dataset$counts)))
  
  # test converting to Seurat
  obj <- as_seurat(out$model)
  expect_equal(dim(obj), dim(t(out$dataset$counts)))
  
  # test converting to anndata
  ad <- as_anndata(out$model)
  expect_equal(dim(ad$X), dim(out$dataset$counts))
  
  # test converting to dyno
  dataset <- as_dyno(out$model)
  expect_equal(dim(dataset$counts), dim(out$dataset$counts))
})


