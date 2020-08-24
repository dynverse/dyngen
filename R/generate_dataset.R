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
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom ggplot2 ggsave
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' out <- 
#'   initialise_model(
#'     backbone = backbone_bifurcating()
#'   ) %>%
#'   generate_dataset()
#'   
#' model <- out$model
#' dataset <- out$dataset
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
    wrap_dataset(
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
    g3 <- plot_feature_network(model, show_targets = FALSE) + labs(title = "TF reg. net.")
    g4 <- plot_feature_network(model) + labs(title = "TF + target reg. net.")
    g5 <- plot_gold_simulations(model) + labs(title = "Gold + simulations")
    g6 <- plot_gold_mappings(model, do_facet = FALSE) + labs(title = "Simulations to gold mapping")
    g7 <- plot_simulations(model) + labs(title = "Simulation time")
    g8 <- plot_gold_expression(model, what = "mol_mrna") + labs(title = "Gold mRNA expression over time")
    g9 <- plot_simulation_expression(model, what = "mol_mrna") + labs(title = "Simulation 1 mRNA expression over time")
    g <- patchwork::wrap_plots(
      g1, g2, g3, g4, g5, g6, g7, g8, g9,
      byrow = TRUE,
      ncol = 3,
      widths = rep(1, 3),
      heights = rep(1, 3)
    ) +
      patchwork::plot_annotation(tag_levels = "A")
    
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
