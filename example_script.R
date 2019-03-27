devtools::load_all(".")

model <- 
  generator_config(
    modulenet = modulenet_consecutive_bifurcating(),
    verbose = TRUE
  ) %>% 
  generate_tfnet()
