library(tidyverse)

out <- 
  initialise_model(
    backbone = backbone_bifurcating(),
    num_tfs = 40,
    num_targets = 20,
    num_hks = 0,
    verbose = FALSE,
    num_cells = 100,
    gold_standard_params = gold_standard_default(tau = .1, census_interval = 10),
    simulation_params = simulation_default(
      ssa_algorithm = ssa_etl(tau = .1),
      census_interval = 10,
      experiment_params = simulation_type_wild_type(num_simulations = 4)
    )
  ) %>% 
  generate_dataset()

example_model <- out$model

map_df(
  names(example_model),
  function(nm) {
    tibble(name = nm, size = (pryr::object_size(example_model[[nm]])))
  }
) %>% arrange(size)

usethis::use_data(example_model, compress = "xz", overwrite = TRUE)
