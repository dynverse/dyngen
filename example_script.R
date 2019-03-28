devtools::load_all(".")

model <- 
  initialise_model(
    modulenet = modulenet_consecutive_bifurcating(),
    platform = platform_simple(n_cells = 1000, n_features = 2000, pct_main_features = 0.2),
    verbose = TRUE
  ) %>% 
  generate_tf_network() %>% 
  generate_feature_network()

plot_module_network(model)
plot_feature_network(model)

