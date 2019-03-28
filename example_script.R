devtools::load_all(".")

model <- 
  initialise_model(
    modulenet = modulenet_consecutive_bifurcating(),
    platform = platform_simple(n_cells = 1000, n_features = 10000, pct_trajectory_features = 0.2),
    verbose = TRUE
  ) %>% 
  generate_tfnet() %>% 
  generate_feature_network()
