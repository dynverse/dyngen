#' Combine multiple dyngen models
#'
#' Assume the given models have the exact same feature ids and ran up until the `generate_cells()` step.
#' In addition, the user is expected to run `generate_experiment()` on the combined models.
#' 
#' See the [vignette on simulating batch effects](https://dyngen.dynverse.org/articles/advanced_topics/simulating_knockouts.html) on how to use this function.
#'
#' @param models A named list of models. The names of the list will be used to 
#'   prefix the different cellular states in the combined model.
#' @param duplicate_gold_standard Whether or not the gold standards of the models are 
#'   different and should be duplicated and prefixed.
#'
#' @export
#'
#' @examples
#' data("example_model")
#' model_ab <- combine_models(list("left" = example_model, "right" = example_model))
#'
#' \donttest{
#' # show a dimensionality reduction
#' plot_simulations(model_ab)
#' plot_gold_mappings(model_ab, do_facet = FALSE)
#' }
combine_models <- function(models, duplicate_gold_standard = TRUE) {
  assert_that(
    is.list(models),
    length(models) >= 1,
    !is.null(names(models))
  )
  
  model_combined <- models[[1]]
  if (duplicate_gold_standard) {
    model_combined$gold_standard <- list()
  }
  model_combined$simulations <- list()
  if ("experiment" %in% names(model_combined)) {
    model_combined$experiment <- list(
      feature_info = models[[1]]$feature_info
    )
  }
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    model_prefix <- names(models)[[i]]
    prefix_fun <-
      if (duplicate_gold_standard) {
        function(x) paste0(model_prefix, "_", x)
      } else {
        identity
      }
    
    if (model$verbose) cat("Merging model ", i, "/", length(models), " ", model_prefix, "\n", sep = "")
    
    # combine gold standard
    if (duplicate_gold_standard) {
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
    }
    
    # combine simulations
    simulation_i_offset <- nrow(model_combined$simulations$perturbed_parameters) %||% 0
    step_ix_offset <- nrow(model_combined$simulations$counts) %||% 0 # compute if need be for experiment
    
    model_combined$simulations$meta <- bind_rows(
      model_combined$simulations$meta,
      model$simulations$meta %>% 
        mutate_at(c("from", "to"), prefix_fun) %>%
        mutate(
          model = model_prefix,
          simulation_i = .data$simulation_i + simulation_i_offset
        )
    )
    for (name in c("counts", "rna_velocity", "cellwise_grn", "reaction_firings", "reaction_propensities", "perturbed_parameters")) {
      model_combined$simulations[[name]] <- rbind(
        model_combined$simulations[[name]],
        model$simulations[[name]]
      )
    }
    if (!is.null(model$simulations$kd_multiplier)) {
      model_combined$simulations$kd_multiplier <- bind_rows(
        model_combined$simulations$kd_multiplier,
        model$simulations$kd_multiplier %>% 
          mutate(simulation_i = .data$simulation_i + simulation_i_offset)
      )
    }
    
    # combine experiment
    if ("experiment" %in% names(model_combined)) {
      cell_id_offset <- nrow(model_combined$experiment$cell_info) %||% 0
      new_cell_ids <- paste0("cell", seq_len(nrow(model$experiment$cell_info)) + cell_id_offset)
      
      for (name in c("counts_premrna", "counts_mrna", "counts_protein", "cellwise_grn", "rna_velocity")) {
        new_mat <- model$experiment[[name]]
        if (!is.null(new_mat)) {
          rownames(new_mat) <- new_cell_ids
          model_combined$experiment[[name]] <- rbind(
            model_combined$experiment[[name]],
            new_mat
          )
        }
      }
      
      model_combined$experiment$cell_info <- bind_rows(
        model_combined$experiment$cell_info,
        model$experiment$cell_info %>% 
          mutate_at(c("from", "to"), prefix_fun) %>%
          mutate(
            model = model_prefix,
            cell_id = new_cell_ids,
            step_ix = .data$step_ix + step_ix_offset,
            simulation_i = .data$simulation_i + simulation_i_offset
          )
      )
    }
  }
  
  # recalculate the dimred
  if (model_combined$simulation_params$compute_dimred) {
    if (model_combined$verbose) cat("Recomputing dimred\n")
    model_combined <- calculate_dimred(model_combined)
  }
  
  # return output
  model_combined
}
