library(tidyverse)
library(dyngen)

set.seed(1)
model <- 
  initialise_model(
    num_cells = 1000,
    num_tfs = 60,
    num_targets = 200,
    num_hks = 0,
    distance_metric = "pearson",
    backbone = backbone_bifurcating_cycle(),
    tf_network_params = tf_network_random(min_tfs_per_module = 1),
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

model <- model  %>% 
  generate_gold_standard()

plot_gold_simulations(model, mapping = aes(comp_1, comp_2)) + scale_colour_brewer(palette = "Dark2")

model <- model %>% 
  generate_cells() %>% 
  generate_experiment()

plot_gold_mappings(model) + scale_colour_brewer(palette = "Dark2")
plot_simulations(model)

traj <- 
  model %>% 
  wrap_dyngen_dataset()

g1 <- dynplot::plot_dimred(traj)
g2 <- dynplot::plot_default(traj)
patchwork::wrap_plots(g1, g2, nrow = 1)
dynplot::plot_heatmap(traj, features_oi = 100)




edge_levels <- 
  model$gold_standard$mod_changes %>% 
  mutate(edge = paste0(from_, "->", to_)) %>% 
  pull(edge)

molecules <- model$feature_info %>% filter(is_tf) %>% gather(mol, val, x, y) %>% pull(val)
df <- bind_cols(
  model$gold_standard$meta,
  as.data.frame(as.matrix(model$gold_standard$counts))[,molecules]
) %>% 
  gather(molecule, value, one_of(molecules)) %>% 
  mutate(edge = factor(paste0(from_, "->", to_), levels = edge_levels)) %>% 
  left_join(model$feature_info %>% select(x, y, module_id) %>% gather(type, molecule, x, y), by = "molecule") %>% 
  group_by(module_id, t, simulation_i, burn, from, to, from_, to_, time, edge, type) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()



ggplot(df) +
  geom_line(aes(t, value, colour = module_id, linetype = type), size = 2) +
  facet_wrap(~edge) +
  theme_bw()# +
  # scale_color_manual(values = model$backbone$module_info %>% select(module_id, color) %>% deframe())
  # scale_colour_brewer(palette = "Set3")

model$gold_standard$mod_changes


# TODO: also make plot for simulations