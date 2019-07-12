library(tidyverse)
library(dyngen)
library(fastgssa)

set.seed(1)
# backbone <- backbone_bifurcating()

backbone <- bblego(
  bblego_start("A", num_modules = 2, type = "simple"),
  bblego_linear("A", "B", num_modules = 5),
  bblego_branching("B", c("C", "D"), num_modules = 4, type = "simple"),
  bblego_end("C", num_modules = 4),
  bblego_end("D", num_modules = 4)
)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 100,
    num_targets = 0,
    num_hks = 0,
    distance_metric = "pearson",
    backbone = backbone,
    tf_network_params = tf_network_default(min_tfs_per_module = 3, sample_num_regulators = function() 2),
    feature_network_params = feature_network_default(target_resampling = 5000),
    kinetics_params = kinetics_default(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(
      burn_time = simtime_from_backbone(backbone, TRUE),
      total_time = simtime_from_backbone(backbone)
    ),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8
  )

model <- model %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_kinetics()

plot_backbone_modulenet(model)
plot_feature_network(model, show_targets = FALSE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)

model <- model %>% 
  generate_gold_standard()

plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
plot_gold_expression(model, "w")
plot_gold_expression(model, c("w", "x"))
plot_gold_expression(model)
  
model <- model %>% generate_cells()
# handle <- generate_cells(model, qsub = TRUE)
# model <- generate_cells(model, qsub = handle)

plot_gold_mappings(model, do_facet = T) + scale_colour_brewer(palette = "Dark2")
plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
patchwork::wrap_plots(
  plot_simulations(model) + coord_equal(),
  plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2") + coord_equal(),
  nrow = 1
)
plot_simulation_expression(model, 22, "w")
plot_simulation_expression(model, 2, "w")
plot_simulation_expression(model, 3, "w")
plot_simulation_expression(model, 5)

plot_simulation_expression(model, 5, "w") + facet_wrap(~module_group, ncol = 3) + geom_vline(aes(xintercept = 0))

# reg_names <- sample(colnames(model$simulations$regulation), 10)
# df <- bind_cols(
#   model$simulations$meta,
#   as.data.frame(model$simulations$regulation[, reg_names]),
# ) %>% gather("regulation", "value", !!reg_names) %>% 
#   mutate(regulation = gsub("regulation_", "", regulation))
# ggplot(df) + 
#   geom_point(aes(sim_time, value, colour = paste0(from, "_", to)), size = .5) +
#   facet_wrap(~regulation, scales = "free_y") +
#   theme_bw()

model <- model %>%
  generate_experiment()


traj <- 
  model %>% 
  wrap_dataset()

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
g3 <- dynplot(traj, layout = layout_graph(traj)) +
  geom_cell_point(color = "grey80") +
  new_scale_fillcolour() +
  geom_trajectory_segments(aes(colour = milestone_percentages), size = 2) +
  geom_milestone_label(aes(fill = milestone_id)) +
  scale_milestones_fillcolour()

patchwork::wrap_plots(g1, g2, g3, nrow = 1)



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
  scale_milestones_fillcolour()

expr <- log2(as.matrix(model$experiment$xcounts) + 1)
regulators <- model$feature_info %>% filter(is_tf) %>% pull(feature_id)
# handle <- bred::qsub_submit_bred(
#   expr = expr,
#   samples = rownames(expr),
#   regulators = regulators,
#   targets = colnames(expr)
# )
# write_rds(lst(handle, model), "~/handle.rds")
# list2env(read_rds("~/handle.rds"), .GlobalEnv)
grn_pred <- bred::qsub_retrieve_bred(handle, force = TRUE, interaction_importance_filter = .1)
grn_pred <- grn_pred %>% filter(adj_importance > .01)

targets_finished <- unique(grn_pred$target)
gold_network <- model$feature_network %>% filter(to %in% targets_finished)
gold_regulation <- model$experiment$regulation[, paste0("regulation_", gold_network$from, "_", gold_network$to)]
regulation_index <- 
  grn_pred %>% select(regulator, target) %>% unique() %>% mutate(
  paste = paste0("regulation_", regulator, "_", target),
  index = match(paste, colnames(gold_regulation)))

grn_pred2 <- 
  grn_pred %>% 
  left_join(regulation_index %>% select(regulator, target, index), by = c("regulator", "target")) %>% 
  mutate(
    x = ifelse(is.na(index), 0, gold_regulation[cbind(sample, index)])
  ) %>% 
  arrange(desc(adj_importance))

# ONE
grn_pred2_one <- grn_pred2 %>% filter(sample == "cell1")
eval_one <- GENIE3::evaluate_ranking_direct(
  values = grn_pred2_one$adj_importance, 
  are_true = grn_pred2_one$x > 0,
  num_positive_interactions = sum(gold_regulation[1,] > 0),
  num_possible_interactions = length(regulators) * length(targets_finished)
)
eval_one$area_under
GENIE3::plot_curves(eval_one)

# ALL
grn_pred2_all <- grn_pred2 %>%
  group_by(regulator, target) %>% 
  summarise(adj_importance = sum(adj_importance) / nrow(gold_regulation)) %>% 
  ungroup() %>% 
  left_join(gold_network %>% transmute(from, to, x = strength), by = c("regulator"="from", "target"="to")) %>% 
  mutate(x = ifelse(is.na(x), 0, x)) %>% 
  arrange(desc(adj_importance))
eval_all <- GENIE3::evaluate_ranking_direct(
  values = grn_pred2_all$adj_importance, 
  are_true = grn_pred2_all$x > 0,
  num_positive_interactions = nrow(gold_network),
  num_possible_interactions = length(regulators) * length(targets_finished)
)
eval_all$area_under
GENIE3::plot_curves(eval_all)

# edge
# cell_grouping <- factor(dynwrap::group_onto_trajectory_edges(traj))
grn_mat <- grn_pred %>% mutate(j = as.integer(regulator) + length(regulators) * (as.integer(target) - 1)) %>% reshape2::acast(sample ~ j, value.var = "adj_importance", fill = 0)
grn_mat[grn_mat < 1e-2] <- 0
space <- dyndimred::dimred_landmark_mds(Matrix::Matrix(grn_mat, sparse = TRUE), ndim = 4, distance_method = "pearson")
cell_grouping <- kmeans(space, centers = 7)$cluster
SCORPIUS::draw_trajectory_plot(space, factor(cell_grouping))

cell_group_counts <- table(cell_grouping)
grn_pred2_edge <- grn_pred2 %>%
  mutate(
    edge = cell_grouping[sample],
    num_cells = cell_group_counts[edge]
  ) %>% 
  group_by(edge, regulator, target) %>% 
  summarise(
    adj_importance = sum(adj_importance / num_cells),
    x = sum(x / num_cells)
  ) %>% 
  ungroup() %>% 
  arrange(desc(adj_importance))

eval_edges <- map(names(cell_group_counts), function(grp) {
  ev <- GENIE3::evaluate_ranking_direct(
    values = grn_pred2_edge %>% filter(edge == grp) %>% .$adj_importance, 
    are_true = grn_pred2_edge %>% filter(edge == grp) %>% .$x > 0,
    num_positive_interactions = nrow(gold_network),
    num_possible_interactions = length(regulators) * length(targets_finished)
  )
  ev$area_under <- ev$area_under %>% mutate(name = grp)
  ev$metrics <- ev$metrics %>% mutate(name = grp)
  ev
})
eval_edge <- lst(
  area_under = map_df(eval_edges, "area_under"),
  metrics = map_df(eval_edges, "metrics")
)

bind_rows(eval_all$area_under %>% mutate(name = "all"), eval_edge$area_under)
GENIE3::plot_curves(list(metrics = bind_rows(eval_all$metrics %>% mutate(name = "all"), eval_edge$metrics)))
