library(tidyverse)

set.seed(1)

backbone <- backbone_bifurcating()
total_time <- simtime_from_backbone(backbone)

out <- 
  initialise_model(
    backbone = backbone,
    num_tfs = nrow(backbone$module_info),
    num_targets = 10,
    num_hks = 10,
    verbose = FALSE,
    num_cells = 100,
    gold_standard_params = gold_standard_default(tau = .1, census_interval = 100),
    simulation_params = simulation_default(
      ssa_algorithm = ssa_etl(tau = .1),
      census_interval = 500,
      experiment_params = simulation_type_wild_type(num_simulations = 100),
      compute_cellwise_grn = TRUE,
      compute_rna_velocity = TRUE
    )
  ) %>% 
  generate_dataset(format = "none")

example_model <- out$model
example_model$num_cores <- NULL
example_model$download_cache_dir <- NULL

map_df(
  names(example_model),
  function(nm) {
    tibble(name = nm, size = (pryr::object_size(example_model[[nm]])))
  }
) %>% arrange(size)

usethis::use_data(example_model, compress = "xz", overwrite = TRUE)
