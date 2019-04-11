library(tidyverse)
library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 1000,
    num_features = 200,
    pct_tfs = .3,
    pct_hks = .3,
    dist_metric = "cosine",
    modulenet = modulenet_bifurcating_converging(),
    tfgen_params = tfgen_random(min_tfs_per_module = 3),
    simulation_params = simulation_default(total_time = 10, num_simulations = 32),
    simulation_setup = simulation_setup_custom(),
    verbose = TRUE,
    num_cores = 8
  ) %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_simulation_setup() %>% 
  simulate_cells() %>% 
  simulate_goldstandard() %>% 
  simulate_experiment()


plot_module_network(model)
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
