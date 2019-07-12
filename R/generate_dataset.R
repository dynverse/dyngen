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
#' @param store_grn Whether or not to store the gene regulatory network
#'   in the dynwrap object as well.
#' 
#' @export
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom ggplot2 ggsave
generate_dataset <- function(model, output_dir = NULL, make_plots = FALSE, store_grn = FALSE) {
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network() %>% 
    generate_kinetics() %>% 
    generate_gold_standard() %>%  
    generate_cells() %>% 
    generate_experiment()
  
  dataset <-
    wrap_dataset(model, store_grn = store_grn)
  
  # write to file
  if (!is.null(output_dir)) {
    dir.create(basename(output_dir), showWarnings = FALSE, recursive = FALSE)
    write_rds(dataset, paste0(output_dir, "dataset.rds"), compress = "gz")
    write_rds(model, paste0(output_dir, "model.rds"), compress = "gz")
  }
  
  if (make_plots) {
    # make plots :scream:
    g1 <- plot_backbone_statenet(model) + labs(title = "Backbone state network")
    g2 <- plot_backbone_modulenet(model) + labs(title = "Backbone module reg. net.")
    g3 <- plot_feature_network(model, show_targets = FALSE) + labs(title = "TF reg. net.")
    g4 <- plot_feature_network(model) + labs(title = "TF + target reg. net.")
    g5 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3") + labs(title = "Gold + simulations")
    g6 <- plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Set3") + labs(title = "Simulations to gold mapping")
    g7 <- plot_simulations(model) + labs(title = "Simulation time")
    g8 <- plot_gold_expression(model, what = "w") + labs(title = "Gold pre-mRNA expression over time")
    g9 <- plot_simulation_expression(model, what = "w") + labs(title = "Simulation 1 pre-mRNA expression over time")
    g <- patchwork::wrap_plots(
      g1, g2, g3, g4, g5, g6, g7, g8, g9,
      byrow = TRUE,
      ncol = 3
    ) +
      patchwork::plot_annotation(tag_levels = "A")
    
    if (!is.null(output_dir)) {
      ggsave(paste0(output_dir, "plot.pdf"), g, width = 30, height = 25)
    } else {
      print(g)
    }
  }
  
  if (is.null(output_dir)) {
    list(dataset = dataset, model = model)
  } else {
    invisible()
  }
}
