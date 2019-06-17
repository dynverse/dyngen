library(tidyverse)
library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 300,
    num_tfs = 7,
    num_targets = 10,
    num_hks = 5,
    distance_metric = "pearson",
    backbone = backbone_bifurcating(),
    tf_network_params = tf_network_random(min_tfs_per_module = 1),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(ssa_algorithm = ssa_em(.01, 2)),
    # simulation_params = simulation_default(burn_time = 2, total_time = 10, num_simulations = 10, ssa_algorithm = ssa_em(h = .01, noise_strength = 6)),
    simulation_params = simulation_default(burn_time = 2, total_time = 10, num_simulations = 10, ssa_algorithm = ssa_direct()),
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
plot_gold_simulations_proj(model, mapping = aes(comp_1, comp_2)) + scale_colour_brewer(palette = "Dark2")
plot_gold_expression(model, "w")
plot_gold_expression(model, c("w", "x"))
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
plot_simulation_expression(model, 5)

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
  geom_velocity_arrow(stat = stat_velocity_grid(grid_n = 10))
g1

g2 <- dynplot(traj) +
  geom_cell_point(color = "grey80") +
  new_scale_fillcolour() +
  geom_trajectory_segments(aes(colour = milestone_percentages), size = 2) +
  geom_milestone_label(aes(fill = milestone_id)) +
  scale_milestones_fillcolour() +
  geom_velocity_arrow(stat = stat_velocity_cells()) 

patchwork::wrap_plots(g1, g2, nrow = 1)



cell_layout <- layout_onedim(traj)
feature_modules <- get_features(traj)
# feature_modules <- model$feature_info %>% filter(is_hk) %>% transmute(feature_id, module_ix = 1)
feature_layout <- layout_modules(traj, feature_modules = feature_modules, cell_layout = cell_layout)
layout <- layout_heatmap(traj, feature_layout = feature_layout)

dynplot(traj, layout = layout) +
  geom_trajectory_segments(aes(color = milestone_percentages)) +
  geom_trajectory_connection() +
  geom_milestone_label(aes(fill = milestone_id, hjust = as.integer(type == "end"))) +
  scale_milestones_fillcolour() +
  new_scale_fillcolour() +
  geom_expression_raster() +
  scale_expression_fillcolour() +
  new_scale_fillcolour() +
  geom_tile(aes(x = x, y = 1))


feature_modules <- model$feature_info %>% filter(is_hk) %>% transmute(feature_id, module_ix = 1)
feature_layout <- layout_modules(traj, feature_modules = feature_modules, cell_layout = cell_layout)
layout <- layout_heatmap(traj, feature_layout = feature_layout)

dynplot(traj, layout = layout) +
  geom_trajectory_segments(aes(color = milestone_percentages)) +
  geom_trajectory_connection() +
  geom_milestone_label(aes(fill = milestone_id, hjust = as.integer(type == "end"))) +
  scale_milestones_fillcolour() +
  new_scale_fillcolour() +
  geom_expression_raster() +
  scale_expression_fillcolour() +
  new_scale_fillcolour() +
  geom_tile(aes(x = x, y = 1))

library(dyno)
out <- infer_trajectory(traj, ti_scorpius())


dynplot(out %>% dynwrap::add_dimred(dimred = dyndimred::dimred_landmark_mds(dynwrap::get_expression(traj), ndim = 2, distance_method = "spearman"))) +
  geom_cell_point(color = "grey80") +
  new_scale_fillcolour() +
  geom_trajectory_segments(aes(colour = milestone_percentages), size = 2) +
  geom_milestone_label(aes(fill = milestone_id)) +
  scale_milestones_fillcolour() #+
  # geom_velocity_arrow(stat = stat_velocity_grid(grid_n = 20))

