devtools::load_all(".")

model <- 
  generator_config(
    modulenet = modulenet_consecutive_bifurcating(),
    platform = platform_simple(n_cells = 1000, n_features = 10000),
    verbose = TRUE
  ) %>% 
  generate_tfnet() %>% 
  generate_targetnet()
