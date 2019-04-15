library(tidyverse)
library(dyngen)

set.seed(3)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 50,
    num_targets = 200,
    num_hks = 300,
    distance_metric = "pearson",
    backbone = backbone_bifurcating_converging(),
    tf_network_params = tf_network_random(min_tfs_per_module = 3),
    feature_network_params = feature_network_default(),
    kinetics_params = kinetics_default(),
    simulation_params = simulation_default(total_time = 10, num_simulations = 32),
    gold_standard_params = gold_standard_default(),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    num_cores = 8,
    download_cache_dir = "~/.cache/dyngen"
  ) %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_kinetics() %>% 
  generate_gold_standard() %>% 
  generate_cells() %>% 
  generate_experiment()

plot_backbone(model)
plot_feature_network(model)
plot_simulations(model)
plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")


g <- 
  ggplot(model$simulations$meta %>% filter(t >= 0, simulation_i <= 16) %>% mutate(edge = paste0(from, "_", to))) + 
  geom_path(aes(t, simulation_i + time, colour = edge, group = simulation_i), size = 1) +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = seq(1, 17))

ggsave(plot = g, filename = "~/gs.pdf", width = 20, height = 16)

traj <- wrap_dyngen_dataset(model)
dynplot::plot_dimred(traj)
dynplot::plot_default(traj)
dynplot::plot_heatmap(traj, features_oi = 100)
