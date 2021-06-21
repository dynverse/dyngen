#' Generate a dataset
#' 
#' This function contains the complete pipeline for generating a dataset
#' with \pkg{dyngen}. In order to have more control over how the dataset 
#' is generated, run each of the steps in this function separately.
#' 
#' @param model A dyngen initial model created with [initialise_model()].
#' @param output_dir If not `NULL`, then the generated model and dynwrap 
#'   dataset will be written to files in this directory.
#' @param make_plots Whether or not to generate an overview of the dataset.
#' 
#' @inheritParams wrap_dataset
#' 
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @importFrom ggplot2 ggsave
#' @importFrom methods is
#' 
#' @export
#' 
#' @return A list containing a dyngen model (`li$model`) and a dynwrap dataset (`li$dataset`).
#' 
#' @examples
#' model <- 
#'   initialise_model(
#'     backbone = backbone_bifurcating()
#'   )
#' \dontshow{
#' # actually use a smaller example 
#' # to reduce execution time during
#' # testing of the examples
#' model <- initialise_model(
#'   backbone = model$backbone,
#'   num_cells = 5,
#'   num_targets = 0,
#'   num_hks = 0,
#'   gold_standard_params = gold_standard_default(census_interval = 1, tau = 0.1),
#'   simulation_params = simulation_default(
#'     burn_time = 10,
#'     total_time = 10,
#'     census_interval = 1,
#'     ssa_algorithm = ssa_etl(tau = 0.1),
#'     experiment_params = simulation_type_wild_type(num_simulations = 1)
#'   )
#' )
#' }
#' \donttest{
#' out <- generate_dataset(model, format = "list")
#'   
#' model <- out$model
#' dataset <- out$dataset
#' 
#' # can also generate other dataset formats:
#' # out <- generate_dataset(model, format = "dyno")
#' # out <- generate_dataset(model, format = "sce")
#' # out <- generate_dataset(model, format = "seurat")
#' # out <- generate_dataset(model, format = "anndata")
#' }
generate_dataset <- function(
  model,
  format = c("list", "dyno", "sce", "seurat", "anndata", "none"),
  output_dir = NULL,
  make_plots = FALSE, 
  store_dimred = model$simulation_params$compute_dimred,
  store_cellwise_grn = model$simulation_params$compute_cellwise_grn,
  store_rna_velocity = model$simulation_params$compute_rna_velocity
) {
  assert_that(is(model, "dyngen::init"))
  format <- match.arg(format)
  
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network() %>% 
    generate_kinetics() %>% 
    generate_gold_standard() %>%  
    generate_cells() %>% 
    generate_experiment()
  
  if (model$verbose && format != "none") cat("Wrapping dataset as ", format, "\n", sep = "")
  dataset <- wrap_dataset(
    model, 
    store_dimred = store_dimred, 
    store_cellwise_grn = store_cellwise_grn,
    store_rna_velocity = store_rna_velocity,
    format = format
  )

    # write to file
  if (!is.null(output_dir)) {
    if (model$verbose) cat("Writing model to file\n")
    dir.create(dirname(output_dir), showWarnings = FALSE, recursive = FALSE)
    
    if (format == "anndata") {
      dataset$write_h5ad(paste0(output_dir, "dataset.h5ad"))
    } else if (format != "none") {
      saveRDS(dataset, paste0(output_dir, "dataset.rds"))
    }
    
    saveRDS(model, paste0(output_dir, "model.rds"))
  }
  
  if (make_plots) {
    if (model$verbose) cat("Making plots\n")
    g <- plot_summary(model)
    
    if (!is.null(output_dir)) {
      ggsave(paste0(output_dir, "plot.pdf"), g, width = 30, height = 25)
    }
  }
  
  if (is.null(output_dir)) {
    out <- list(model = model, dataset = dataset)
    if (make_plots) {
      out$plot <- g
    }
    out
  } else {
    invisible()
  }
}
