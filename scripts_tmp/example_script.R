library(tidyverse)
library(dyngen)

set.seed(2)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 60,
    num_targets = 200,
    num_hks = 500,
    distance_metric = "pearson",
    backbone = backbone_bifurcating(),
    tf_network_params = tf_network_random(min_tfs_per_module = 3),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_custom(),
    gold_standard_params = gold_standard_default(time_per_edge = 2),
    simulation_params = simulation_default(total_time = 20, num_simulations = 8),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    num_cores = 8,
    download_cache_dir = "~/.cache/dyngen"
  ) %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_kinetics()
  
plot_backbone(model)
plot_feature_network(model, tfs_only = TRUE)
plot_feature_network(model)

model <- model %>% 
  generate_gold_standard()

plot_gold_simulations(model, mapping = aes(comp_1, comp_2)) + scale_colour_brewer(palette = "Dark2")
plot_gold_expression(model)

model <- model %>% 
  generate_cells() %>% 
  generate_experiment()

plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")
plot_simulations(model)
plot_simulation_expression(model)

traj <- 
  model %>% 
  wrap_dyngen_dataset()

g1 <- dynplot::plot_dimred(traj)
g2 <- dynplot::plot_default(traj)
patchwork::wrap_plots(g1, g2, nrow = 1)
dynplot::plot_heatmap(traj, features_oi = 100)






# TODO: also make plot for simulations