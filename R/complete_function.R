
complete_function <- function(model, directory) {
  sink(paste0(directory, "/log.txt"))
  on.exit(sink())
  
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network() %>% 
    generate_kinetics()
  
  pdf(paste0(directory, "/plot_1_network.pdf"), width = 8, height = 8)
  plot_backbone(model)
  plot_feature_network(model)
  dev.off()
  
  model <- model %>% 
    generate_gold_standard()
  
  g <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
  ggsave(paste0(directory, "/plot_2_gold.pdf"), g, width = 8, height = 8)
  
  model <- model %>%  
    generate_cells() 
  
  g1 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
  g2 <- plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")
  g3 <- plot_simulations(model)
  g <- patchwork::wrap_plots(g1, g3, nrow = 1)
  ggsave(paste0(directory, "/plot_3_simulations.pdf"), g, width = 12, height = 6)
  ggsave(paste0(directory, "/plot_4_gold_mappings.pdf"), g2, width = 12, height = 10)
  
  model <- model %>% 
    generate_experiment()
  
  traj <- 
    model %>% 
    wrap_dyngen_dataset()
  
  g1 <- dynplot::plot_dimred(traj)
  g2 <- dynplot::plot_default(traj)
  g <- patchwork::wrap_plots(g1, g2, nrow = 1)
  g3 <- dynplot::plot_heatmap(traj, features_oi = 100)
  ggsave(paste0(directory, "/plot_5_dataset_dimred.pdf"), g, width = 12, height = 6)
  ggsave(paste0(directory, "/plot_6_heatmap.pdf"), g3, width = 12, height = 8)
  
  write_rds(traj, paste0(directory, "/out_dataset.rds"), compress = "xz")
  write_rds(model, paste0(directory, "/out_model.rds"), compress = "xz")
  
  return()
}
