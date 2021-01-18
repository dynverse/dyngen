#' Combine multiple dyngen models
#'
#' Assume the given models have the exact same feature ids and ran up until the `generate_cells()` step.
#' In addition, the user is expected to run `generate_experiment()` on the combined models.
#'
#' @param models A named list of models. The names of the list will be used to 
#'   prefix the different cellular states in the combined model.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # create the config object
#' init_config <- initialise_model(
#'   backbone = backbone_linear(),
#'   num_cells = 500,
#'   num_targets = 250,
#'   num_hks = 50,
#'   gold_standard_params = gold_standard_default(census_interval = 1, tau = 0.05),
#'   simulation_params = simulation_default(
#'     # burn_time = 10,
#'     # total_time = 10,
#'     census_interval = 1,
#'     ssa_algorithm = ssa_etl(tau = 0.05),
#'     experiment_params = simulation_type_wild_type(num_simulations = 40)
#'   )
#' )
#'
#' # generate the genes and their kinetics
#' model_common <-
#'   init_config %>%
#'   generate_tf_network() %>%
#'   generate_feature_network()
#' plot_feature_network(model_common)
#'
#' # run the simulation once
#' model_a <- model_common %>%
#'   generate_kinetics() %>%
#'   generate_gold_standard() %>%
#'   generate_cells()
#'
#' # run the simulation once more
#' model_b <- model_common %>%
#'   generate_kinetics() %>%
#'   generate_gold_standard() %>%
#'   generate_cells()
#'
#' # combine models, do experiment afterwards
#' model_ab <- combine_models(list("left" = model_a, "right" = model_b) %>%
#'   generate_experiment()
#'
#' # show a dimensionality reduction
#' plot_simulations(model_ab)
#' plot_gold_mappings(model_ab, do_facet = FALSE)
#'
#' # create a dynwrap dataset
#' dataset <- wrap_dataset(model_ab)
#' }
combine_models <- function(models) {
  assert_that(
    is.list(models),
    length(models) >= 1,
    !is.null(names(models))
  )
  
  model_combined <- models[[1]]
  model_combined$gold_standard <- list()
  model_combined$simulations <- list()
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    model_prefix <- names(models)[[i]]
    prefix_fun <- function(x) paste0(model_prefix, "_", x)
    
    if (model$verbose) cat("Merging model ", i, "/", length(models), " ", model_prefix, "\n", sep = "")
    
    # combine gold standard
    model_combined$gold_standard$mod_changes <- bind_rows(
      model_combined$gold_standard$mod_changes,
      model$gold_standard$mod_changes %>% mutate_at(c("from", "to", "from_", "to_"), prefix_fun)
    )
    model_combined$gold_standard$meta <- bind_rows(
      model_combined$gold_standard$meta,
      model$gold_standard$meta %>% mutate_at(c("from", "to", "from_", "to_"), prefix_fun)
    )
    model_combined$gold_standard$counts <- rbind(
      model_combined$gold_standard$counts,
      model$gold_standard$counts
    )
    model_combined$gold_standard$network <- bind_rows(
      model_combined$gold_standard$network,
      model$gold_standard$network %>% mutate_at(c("from", "to"), prefix_fun)
    )
    
    # combine simulations
    simulation_i_offset <-
      if (!is.null(model_combined$simulations$meta)) {
        max(model_combined$simulations$meta$simulation_i)
      } else {
        0
      }
    model_combined$simulations$meta <- bind_rows(
      model_combined$simulations$meta,
      model$simulations$meta %>% mutate_at(c("from", "to"), prefix_fun) %>%
        mutate(simulation_i = simulation_i + simulation_i_offset)
    )
    model_combined$simulations$counts <- rbind(
      model_combined$simulations$counts,
      model$simulations$counts
    )
    # could also combine propensities, rna velocity, etc
  }
  
  # recalculate the dimred
  if (model$verbose) cat("Recomputing dimred\n")
  model_combined <- model_combined %>%
    dyngen:::calculate_dimred()
  
  # return output
  model_combined
}
