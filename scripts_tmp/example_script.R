library(tidyverse)
library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 80,
    num_targets = 0,
    num_hks = 0,
    distance_metric = "pearson",
    backbone = backbone_disconnected(),
    tf_network_params = tf_network_random(min_tfs_per_module = 2),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(burn_time = 6, total_time = 10, num_simulations = 50, ssa_algorithm = ssa_em(noise_strength = 6)),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen"
  )

model <- model %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_kinetics()

plot_backbone(model)
plot_feature_network(model, show_targets = FALSE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)

model <- model %>% 
  generate_gold_standard()

plot_gold_simulations(model, mapping = aes(comp_1, comp_2)) + scale_colour_brewer(palette = "Dark2")
plot_gold_expression(model, "w")
plot_gold_expression(model)

# model$num_cores <- 1
model <- model %>% 
  generate_cells() %>%
  generate_experiment()

plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")
plot_simulations(model)
plot_simulation_expression(model, 1)
plot_simulation_expression(model, 2)
plot_simulation_expression(model, 3)
plot_simulation_expression(model, 6)

plot_simulation_expression(model, 5, "w") + facet_wrap(~module_group, ncol = 3) + geom_vline(aes(xintercept = 0))

traj <- 
  model %>% 
  wrap_dyngen_dataset()

library(dynplot)
g1 <- dynplot(traj) +
  geom_cell_point(color = "grey80") +
  new_scale_fillcolour() +
  geom_trajectory_segments(aes(colour = milestone_percentages), size = 2) +
  geom_milestone_label(aes(fill = milestone_id)) +
  scale_milestones_fillcolour() +
  geom_velocity_arrow(stat = stat_velocity_grid(grid_n = 20))

g2 <- dynplot(dataset) +
  geom_cell_point(color = "grey80") +
  new_scale_fillcolour() +
  geom_trajectory_segments(aes(colour = milestone_percentages), size = 2) +
  geom_milestone_label(aes(fill = milestone_id)) +
  scale_milestones_fillcolour() +
  geom_velocity_arrow(stat = stat_velocity_cells()) 

g1 <- dynplot::plot_dimred(traj)
g2 <- dynplot::plot_default(traj)
patchwork::wrap_plots(g1, g2, nrow = 1)
dynplot::plot_heatmap(traj, features_oi = 100)

dynplot::plot_dimred(traj)

# TODO: also make plot for simulations