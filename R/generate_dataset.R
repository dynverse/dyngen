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
#' @inheritParams as_dyno
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
#' out <- generate_dataset(model)
#'   
#' model <- out$model
#' dataset <- out$dataset
#' 
#' # dynplot::plot_dimred(dataset)
#' }
generate_dataset <- function(
  model, 
  output_dir = NULL,
  make_plots = FALSE, 
  store_dimred = model$simulation_params$compute_dimred,
  store_cellwise_grn = model$simulation_params$compute_cellwise_grn,
  store_rna_velocity = model$simulation_params$compute_rna_velocity
) {
  assert_that(is(model, "dyngen::init"))
  
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network() %>% 
    generate_kinetics() %>% 
    generate_gold_standard() %>%  
    generate_cells() %>% 
    generate_experiment()
  
  if (model$verbose) cat("Wrapping dataset\n")
  dataset <-
    as_dyno(
      model, 
      store_dimred = store_dimred, 
      store_cellwise_grn = store_cellwise_grn,
      store_rna_velocity = store_rna_velocity
    )
  
  # write to file
  if (!is.null(output_dir)) {
    if (model$verbose) cat("Writing model to file\n")
    dir.create(dirname(output_dir), showWarnings = FALSE, recursive = FALSE)
    saveRDS(dataset, paste0(output_dir, "dataset.rds"))
    saveRDS(model, paste0(output_dir, "model.rds"))
  }
  
  if (make_plots) {
    if (model$verbose) cat("Making plots\n")
    # make plots :scream:
    g1 <- plot_backbone_statenet(model) + labs(title = "Backbone state network")
    g2 <- plot_backbone_modulenet(model) + labs(title = "Backbone module reg. net.")
    g3 <- plot_feature_network(model) + labs(title = "TF + target reg. net.")
    g4 <- plot_gold_simulations(model) + labs(title = "Gold + simulations")
    g5 <- plot_gold_mappings(model, do_facet = FALSE) + labs(title = "Simulations to gold mapping")
    g6 <- plot_simulations(model) + labs(title = "Simulation time")
    g7 <- plot_gold_expression(model, what = "mol_mrna") + labs(title = "Gold mRNA expression over time")
    g8 <- plot_simulation_expression(model, what = "mol_mrna") + labs(title = "Simulation 1 mRNA expression over time")
    g9 <- plot_experiment_dimred(model) + labs(title = "Dim. Red. of final dataset")
    
    
    g <- patchwork::wrap_plots(
      g1, g2, g3, 
      g4, g5, g6,
      g7, g8, g9,
      byrow = TRUE,
      ncol = 3,
      widths = rep(1, 3),
      heights = rep(1, 3)
    ) +
      patchwork::plot_annotation(tag_levels = "A") +
      patchwork::plot_layout(guides = "collect")
    
    if (!is.null(output_dir)) {
      ggsave(paste0(output_dir, "plot.pdf"), g, width = 30, height = 25)
    }
  }
  
  if (is.null(output_dir)) {
    out <- list(dataset = dataset, model = model)
    if (make_plots) {
      out$plot <- g
    }
    out
  } else {
    invisible()
  }
}
