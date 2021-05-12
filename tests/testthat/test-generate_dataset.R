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

test_that("check plot", {
  expect_is(out$plot, "ggplot")
})
  

test_that("test converting to SCE", {
  skip_if(!requireNamespace("SingleCellExperiment", quietly = TRUE))
  sce <- as_sce(out$model)
  expect_equal(dim(sce), dim(t(out$dataset$counts)))
})

test_that("test converting to Seurat", {
  skip_if(!requireNamespace("Seurat", quietly = TRUE))
  obj <- as_seurat(out$model)
  expect_equal(dim(obj), dim(t(out$dataset$counts)))
})

test_that("test converting to anndata", {
  skip_if(!requireNamespace("anndata", quietly = TRUE))
  
  # check if python anndata is installed
  py_anndata_available <-
    tryCatch({
      anndata::AnnData()
      TRUE
    }, error = function(e) {
      FALSE
    })
  skip_if(!py_anndata_available)
  
  ad <- as_anndata(out$model)
  expect_equal(dim(ad$X), dim(out$dataset$counts))
})
  
test_that("test converting to dyno", {
  skip_if(!requireNamespace("dynwrap", quietly = TRUE))
  dataset <- as_dyno(out$model)
  expect_equal(dim(dataset$counts), dim(out$dataset$counts))
})


