# @importFrom dynplot plot_dimred plot_default plot_heatmap
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
complete_function <- function(model, directory) {
  print(model$numbers)
  
  sink(paste0(directory, "/log.txt"))
  on.exit(sink())
  
  ## GENERATE FEATURE NETWORK
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network() %>% 
    generate_kinetics()
  
  g1 <- plot_backbone(model)
  g2 <- plot_feature_network(model, show_targets = FALSE)
  g3 <- plot_feature_network(model)
  g <- patchwork::wrap_plots(
    patchwork::wrap_plots(g1, g2, ncol = 1),
    g3, nrow = 1
  )  
  ggsave(paste0(directory, "/plot_1_network.pdf"), g, width = 20, height = 12)
  
  ## GENERATE GOLD STANDARD
  model <- model %>% 
    generate_gold_standard()
  
  g4 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
  g5 <- plot_gold_expression(model)
  ggsave(paste0(directory, "/plot_2_gold.pdf"), g4, width = 8, height = 8)
  ggsave(paste0(directory, "/plot_3_gold_expression.pdf"), g5, width = 12, height = 12)
  
  ## GENERATE CELLS
  model <- model %>%  
    generate_cells()
  
  g6 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
  g7 <- plot_gold_mappings(model) + scale_colour_brewer(palette = "Set3")
  g8 <- plot_simulations(model)
  g9 <- plot_simulation_expression(model)
  g <- patchwork::wrap_plots(g6, g8, nrow = 1)
  ggsave(paste0(directory, "/plot_4_simulations.pdf"), g, width = 12, height = 6)
  ggsave(paste0(directory, "/plot_5_gold_mappings.pdf"), g7, width = 12, height = 10)
  
  pdf(paste0(directory, "/plot_6_simulation_expression.pdf"), width = 12, height = 12)
  tryCatch({
    for (sim_i in seq_len(model$simulation_params$num_simulations)) {
      g <- 
        plot_simulation_expression(model, simulation_i = sim_i) +
        labs(title = paste0("Simulation ", sim_i))
      print(g)
    }
  }, finally = {
    dev.off()
  })
  
  ## GENERATE EXPERIMENT
  model <- model %>% 
    generate_experiment()
  
  ## CONVERT TO DYNWRAP
  traj <- 
    model %>% 
    wrap_dyngen_dataset()
  
  # g10 <- dynplot::plot_dimred(traj)
  # # g11 <- dynplot::plot_default(traj)
  # # g <- patchwork::wrap_plots(g10, g11, nrow = 1)
  # g <- g10
  # g12 <- dynplot::plot_heatmap(traj, features_oi = 100)
  # ggsave(paste0(directory, "/plot_7_dataset_dimred.pdf"), g, width = 12, height = 6)
  # ggsave(paste0(directory, "/plot_8_heatmap.pdf"), g12, width = 12, height = 8)
  
  ## SAVE R OBJECTS
  write_rds(traj, paste0(directory, "/out_dataset.rds"), compress = "xz")
  write_rds(model, paste0(directory, "/out_model.rds"), compress = "gz")
  
  # :scream:
  g <- patchwork::wrap_plots(
    g1, g2, g3, g4, g5, g6, g7, g8, g9, #g10, #g11,
    #g12,
    byrow = FALSE,
    ncol = 4
  ) +
      patchwork::plot_annotation(tag_levels = "A")
  ggsave(paste0(directory, "/plot_all.pdf"), g, width = 40, height = 30)
  
  return()
}
