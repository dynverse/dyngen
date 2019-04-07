devtools::load_all(".")

set.seed(1)
model <- 
  initialise_model(
    num_cells = 1000,
    num_features = 50,
    pct_tfs = 1,
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
  simulate_cells()


plot_module_network(model)
# plot_feature_network(model, tfs_only = TRUE)
plot_feature_network(model)
plot_simulations(model)

model <- 
  model %>% 
  generate_goldstandard()


plot_gold_simulations(model)
plot_gold_mappings(model)

# write_rds(model, "~/yay3.rds")
# write_rds(model, "~/yay2.rds")
# model <- read_rds("~/yay.rds")

# library(gganimate)
# 
# # get feature info
# feature_info <- 
#   model$feature_info %>% 
#   mutate(
#     main_index = match(is_main, c(TRUE, FALSE, NA)),
#     size = c(4, 1, 1)[main_index],
#     label = ""
#   ) %>% 
#   filter(is_main)
# 
# # get feature network
# feature_network <- 
#   model$feature_network %>% 
#   filter(from %in% feature_info$feature_id & to %in% feature_info$feature_id)
# 
# # add extra edges invisible between regulators from the same module
# feature_network <- 
#   bind_rows(
#     feature_network,
#     feature_info %>%
#       filter(is_tf) %>% 
#       select(module_id, feature_id) %>%
#       group_by(module_id) %>%
#       do({
#         crossing(from = .$feature_id, to = .$feature_id) %>%
#           mutate(effect = -2)
#       }) %>% 
#       ungroup() %>% 
#       filter(from < to)
#   )
# 
# # determine colors
# feature_network <-
#   feature_network %>% 
#   mutate(
#     color = effect_colour(effect)
#   )
# 
# # construct graph
# graph <- igraph::graph_from_data_frame(
#   feature_network,
#   vertices = feature_info
# )
# 
# # layout
# #layout <- igraph::layout.fruchterman.reingold(graph) %>% as.data.frame()
# layout <- igraph::layout_with_kk(graph) %>% as.data.frame()
# colnames(layout) <- c("comp1", "comp2")
# rownames(layout) <- feature_info$feature_id
# layout <- layout %>% rownames_to_column("feature_id")
# 
# edges <- 
#   feature_network %>% 
#   filter(effect != -2) %>% 
#   left_join(layout %>% rename_all(function(x) paste0("from_", x)), by = c("from" = "from_feature_id")) %>% 
#   left_join(layout %>% rename_all(function(x) paste0("to_", x)), by = c("to" = "to_feature_id"))
# 
# plot_df <- 
#   feature_info %>% 
#   left_join(layout, by = "feature_id") %>% 
#   left_join(
#     model$simulations %>% 
#       # filter(simulation_i == 1) %>% 
#       select(t, simulation_i, one_of(feature_info$x)) %>% 
#       gather(x, value, -t, -simulation_i),
#     by = "x"
#   ) %>% 
#   mutate(
#     colour_expr = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)[cut(value, breaks = 100)]
#   )
# 
# ggplot(plot_df %>% filter(t == sample(t, 1), simulation_i == 1)) +
#   geom_segment(aes(x = from_comp1, xend = to_comp1, y = from_comp2, yend = to_comp2, linetype = factor(effect, level = c(1, -1))), edges, arrow = arrow(), colour = "darkgray") +
#   geom_point(aes(comp1, comp2, colour = colour_expr), size = 5) +
#   geom_text(aes(comp1, comp2, colour = color, label = feature_id), nudge_y = .2) +
#   scale_colour_identity() +
#   theme_bw() +
#   coord_equal() + 
#   facet_wrap(~simulation_i)
# 
# g <-
#   ggplot(plot_df %>% filter(simulation_i == 1)) +
#   geom_segment(aes(x = from_comp1, xend = to_comp1, y = from_comp2, yend = to_comp2, linetype = factor(effect, level = c(1, -1))), edges, arrow = arrow(), colour = "darkgray") +
#   geom_point(aes(comp1, comp2, colour = colour_expr), size = 5) +
#   geom_text(aes(comp1, comp2, colour = color, label = feature_id), nudge_y = .2) +
#   scale_colour_identity() +
#   theme_bw() +
#   coord_equal() + 
#   facet_wrap(~simulation_i) +
#   labs(title = "time = {current_frame}") +
#   transition_manual(t)
# # animate(g, fig.width = 1000, fig.height = 1000)
# options(gganimate.dev_args = list(width = 1000, height = 1000))
# print(g)
# anim_save(filename = "~/help.gif")
