library(tidyverse)
library(dyngen)
library(fastgssa)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 10000,
    num_tfs = 100,
    num_targets = 4900,
    num_hks = 0,
    distance_metric = "pearson",
    backbone = backbone_bifurcating_loop(),
    tf_network_params = tf_network_random(min_tfs_per_module = 3),
    feature_network_params = feature_network_default(target_resampling = 5000),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(ssa_algorithm = ssa_em(tau = .005, noise_strength = 0)),
    simulation_params = simulation_default(burn_time = 2, total_time = 20, num_simulations = 24, ssa_algorithm = ssa_etl(.05), census_interval = .05, use_vector_optimisation = TRUE),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 1
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

# model$simulation_params$ssa_algorithm <- ssa_em(.01, 6)
# model$simulation_params$ssa_algorithm <- ssa_etl(.01)
# model$simulation_params$use_vector_optimisation <- FALSE
# model <- model %>% generate_cells(qsub = FALSE)
handle <- generate_cells(model, qsub = TRUE)
class(handle) <- c("qsub::qsub_config", "list")
model <- generate_cells(model, qsub = handle)

plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")
patchwork::wrap_plots(
  plot_simulations(model) + coord_equal(),
  plot_gold_simulations(model, mapping = aes(comp_1, comp_2)) + scale_colour_brewer(palette = "Dark2") + coord_equal(),
  nrow = 1
)
plot_simulation_expression(model, 1)
plot_simulation_expression(model, 2)
plot_simulation_expression(model, 3)
plot_simulation_expression(model, 5)

plot_simulation_expression(model, 5, "w") + facet_wrap(~module_group, ncol = 3) + geom_vline(aes(xintercept = 0))

model <- model %>%
  generate_experiment()

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
out <- infer_trajectory(traj, ti_slingshot())


dynplot(out %>% dynwrap::add_dimred(dimred = dyndimred::dimred_landmark_mds(dynwrap::get_expression(traj), ndim = 2, distance_method = "spearman"))) +
  geom_cell_point(color = "grey80") +
  new_scale_fillcolour() +
  geom_trajectory_segments(aes(colour = milestone_percentages), size = 2) +
  geom_milestone_label(aes(fill = milestone_id)) +
  scale_milestones_fillcolour() #+
  # geom_velocity_arrow(stat = stat_velocity_grid(grid_n = 20))

