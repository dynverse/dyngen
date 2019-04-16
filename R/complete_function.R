
complete_function <- function(model, directory) {
  sink(paste0(directory, "/log.txt"))
  on.exit(sink())
  
  ## GENERATE FEATURE NETWORK
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network()
  
  g1 <- plot_backbone(model)
  g2 <- plot_feature_network(model, tfs_only = TRUE)
  g3 <- plot_feature_network(model)
  g <- patchwork::wrap_plots(
    patchwork::wrap_plots(g1, g2, ncol = 1),
    g3, nrow = 1
  )  
  ggsave(paste0(directory, "/plot_1_network.pdf"), g, width = 20, height = 12)
  
  ## GENERATE GOLD STANDARD
  model <- model %>% 
    generate_kinetics() %>% 
    generate_gold_standard()
  
  g4 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
  ggsave(paste0(directory, "/plot_2_gold.pdf"), g4, width = 8, height = 8)
  
  ## GENERATE CELLS
  model <- model %>%  
    generate_cells() 
  
  g5 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
  g6 <- plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")
  g7 <- plot_simulations(model)
  g <- patchwork::wrap_plots(g5, g7, nrow = 1)
  ggsave(paste0(directory, "/plot_3_simulations.pdf"), g, width = 12, height = 6)
  ggsave(paste0(directory, "/plot_4_gold_mappings.pdf"), g6, width = 12, height = 10)
  
  ## GENERATE EXPERIMENT
  model <- model %>% 
    generate_experiment()
  
  ## CONVERT TO DYNWRAP
  traj <- 
    model %>% 
    wrap_dyngen_dataset()
  
  g8 <- dynplot::plot_dimred(traj)
  g9 <- dynplot::plot_default(traj)
  g <- patchwork::wrap_plots(g8, g9, nrow = 1)
  g10 <- dynplot::plot_heatmap(traj, features_oi = 100)
  ggsave(paste0(directory, "/plot_5_dataset_dimred.pdf"), g, width = 12, height = 6)
  ggsave(paste0(directory, "/plot_6_heatmap.pdf"), g10, width = 12, height = 8)
  
  ## SAVE R OBJECTS
  write_rds(traj, paste0(directory, "/out_dataset.rds"), compress = "xz")
  write_rds(model, paste0(directory, "/out_model.rds"), compress = "gz")
  
  # :scream:
  g <- patchwork::wrap_plots(
    g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
    byrow = FALSE,
    ncol = 4
  )
  ggsave(paste0(directory, "/plot_all.pdf"), g, width = 40, height = 30)
  
  return()
}
