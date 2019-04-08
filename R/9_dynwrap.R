#' @export
wrap_dyngen_dataset <- function(model) {
  traj <- 
    dynwrap::wrap_expression(
      counts = model$experiment$counts,
      expression = as(log2(model$experiment$counts + 1), "dgCMatrix"), # TODO: add normalisation
      cell_info = model$experiment$cell_info %>% select(-from, -to, -time),
      feature_info = model$experiment$feature_info
    ) %>% 
    dynwrap::add_trajectory(
      milestone_network = model$goldstandard$network,
      progressions = model$experiment$cell_info %>% select(cell_id, from, to, percentage = time)
    ) %>% 
    dynwrap::add_dimred(
      dimred = model$simulations$dimred[model$experiment$cell_info$step_ix, ] %>% magrittr::set_rownames(model$experiment$cell_info$cell_id),
      dimred_segment_points = model$goldstandard$dimred,
      dimred_segment_progressions = model$goldstandard$meta %>% select(from, to, percentage = time)
    )
}