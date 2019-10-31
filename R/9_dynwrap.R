#' Wrap a simulation into a dynwrap object
#' 
#' The output of this object can be used with \pkg{dyno}.
#' 
#' @param model A dyngen output model for which the experiment has been emulated with [generate_experiment()].
#' @param store_grn Whether or not to also store GRN information.
#' @param store_dimred Whether or not to store the dimensionality reduction constructed on the true counts.
#' 
#' @export
#' @importFrom dynwrap wrap_expression add_trajectory add_dimred
wrap_dataset <- function(model, store_grn = FALSE, store_dimred = FALSE) {
  dataset <- wrap_expression(
    counts = model$experiment$xcounts,
    expression = as(log2(model$experiment$xcounts + 1), "dgCMatrix"),
    expression_unspliced = as(log2(model$experiment$wcounts + 1), "dgCMatrix"),
    expression_protein = as(log2(model$experiment$ycounts + 1), "dgCMatrix"),
    cell_info = model$experiment$cell_info %>% select(-from, -to, -time),
    feature_info = model$experiment$feature_info
  ) %>% 
    add_trajectory(
      milestone_network = model$gold_standard$network,
      progressions = model$experiment$cell_info %>% select(cell_id, from, to, percentage = time)
    )
  
  if (store_dimred) {
    dataset <- dataset %>% 
      add_dimred(
        dimred = model$simulations$dimred[model$experiment$cell_info$step_ix, ] %>% 
          magrittr::set_rownames(model$experiment$cell_info$cell_id),
        dimred_segment_points = model$gold_standard$dimred[!model$gold_standard$meta$burn,],
        dimred_segment_progressions = model$gold_standard$meta[!model$gold_standard$meta$burn,] %>% 
          select(from, to, percentage = time),
        pair_with_velocity = TRUE
      )
  }
  
  # add a few more values
  if (store_grn) {
    dataset$feature_network <- 
      model$feature_network
    dataset$feature_network_sc <- 
      model$experiment$regulation
  }
  
  dataset
}