devtools::load_all(".")

model <- 
  initialise_model(
    modulenet = modulenet_bifurcating_converging(),
    platform = platform_simple(n_cells = 1000, n_features = 50, pct_main_features = 1),
    tfgen_params = tfgen_random(percentage_tfs = 0.2),
    verbose = TRUE,
    num_cores = 8
  ) %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_simulation_setup() %>% 
  simulate_cells()

plot_module_network(model)
plot_feature_network(model)



expr <- model$simulations %>% select(one_of(model$simulation_system$molecule_ids)) %>% as.matrix
expr <- expr[,colSums(expr) != 0]
sim_f <- model$simulations[rowSums(expr) != 0,]
expr <- expr[rowSums(expr) != 0,]
space <- SCORPIUS::reduce_dimensionality(expr, dist_fun = SCORPIUS::correlation_distance)
plot_df <- bind_cols(sim_f %>% select(1:2), as.data.frame(space))




ggplot(plot_df %>% filter(t >= 0)) +
  geom_path(aes(Comp1, Comp2, colour = t, group = simulation_i)) +
  viridis::scale_color_viridis() +
  theme_bw()



ggplot(plot_df %>% filter(t >= 6)) +
  geom_path(aes(Comp1, Comp2, colour = t, group = simulation_i)) +
  viridis::scale_color_viridis() +
  theme_bw()


library(gganimate)

# get feature info
feature_info <- 
  model$feature_info %>% 
  mutate(
    main_index = match(is_main, c(TRUE, FALSE, NA)),
    size = c(4, 1, 1)[main_index],
    label = ""
  ) %>% 
  filter(is_main)

# get feature network
feature_network <- 
  model$feature_network %>% 
  filter(from %in% feature_info$feature_id & to %in% feature_info$feature_id)

# add extra edges invisible between regulators from the same module
feature_network <- 
  bind_rows(
    feature_network,
    feature_info %>%
      filter(is_tf) %>% 
      select(module_id, feature_id) %>%
      group_by(module_id) %>%
      do({
        crossing(from = .$feature_id, to = .$feature_id) %>%
          mutate(effect = -2)
      }) %>% 
      ungroup() %>% 
      filter(from < to)
  )

# determine colors
feature_network <-
  feature_network %>% 
  mutate(
    color = effect_colour(effect)
  )

# construct graph
graph <- igraph::graph_from_data_frame(
  feature_network,
  vertices = feature_info
)

# layout
layout <- igraph::layout.fruchterman.reingold(graph) %>% as.data.frame()
colnames(layout) <- c("comp1", "comp2")
rownames(layout) <- feature_info$feature_id
layout <- layout %>% rownames_to_column("feature_id")

edges <- 
  feature_network %>% 
  left_join(layout %>% rename_all(function(x) paste0("from_", x)), by = c("from" = "from_feature_id")) %>% 
  left_join(layout %>% rename_all(function(x) paste0("to_", x)), by = c("to" = "to_feature_id"))

plot_df <- 
  feature_info %>% 
  left_join(layout, by = "feature_id") %>% 
  left_join(
    model$simulations %>% 
      # filter(simulation_i == 1) %>% 
      select(t, simulation_i, one_of(feature_info$x)) %>% 
      gather(x, value, -t, -simulation_i),
    by = "x"
  )

ggplot(plot_df %>% filter(t == sample(t, 1), simulation_i == 1)) +
  geom_segment(aes(x = from_comp1, xend = to_comp1, y = from_comp2, yend = to_comp2), edges, arrow = arrow(), colour = "darkgray") +
  geom_point(aes(comp1, comp2, colour = value, size = is_tf)) +
  scale_colour_distiller(palette = "RdBu") +
  theme_bw() +
  coord_equal() + 
  facet_wrap(~simulation_i)

g <-
  ggplot(plot_df %>% filter(simulation_i == 1)) +
  geom_segment(aes(x = from_comp1, xend = to_comp1, y = from_comp2, yend = to_comp2), edges, arrow = arrow(), colour = "darkgray") +
  geom_point(aes(comp1, comp2, colour = value, size = is_tf)) +
  scale_colour_distiller(palette = "RdBu") +
  theme_bw() +
  coord_equal() +
  facet_wrap(~simulation_i) + 
  labs(title = "time = {current_frame}") +
  transition_manual(t)
# animate(g, fig.width = 1000, fig.height = 1000)
options(gganimate.dev_args = list(width = 1000, height = 1000))
print(g)
anim_save(filename = "~/help.gif")
