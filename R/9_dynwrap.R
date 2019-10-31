#' Wrap a simulation into a dynwrap object
#' 
#' The output of this object can be used with \pkg{dyno}.
#' 
#' @param model A dyngen output model for which the experiment has been emulated with [generate_experiment()].
#' 
#' @export
#' @importFrom dynwrap wrap_expression add_trajectory add_dimred
wrap_dataset <- function(model) {
  dynwrap::wrap_expression(
    counts = model$experiment$xcounts,
    expression = as(log2(model$experiment$xcounts + 1), "dgCMatrix"),
    expression_unspliced = as(log2(model$experiment$wcounts + 1), "dgCMatrix"),
    expression_protein = as(log2(model$experiment$ycounts + 1), "dgCMatrix"),
    cell_info = model$experiment$cell_info %>% select(-from, -to, -time),
    feature_info = model$experiment$feature_info
  ) %>% 
    dynwrap::add_trajectory(
      milestone_network = model$gold_standard$network,
      progressions = model$experiment$cell_info %>% select(cell_id, from, to, percentage = time)
    ) %>% 
    dynwrap::add_dimred(
      dimred = model$simulations$dimred[model$experiment$cell_info$step_ix, ] %>% 
        magrittr::set_rownames(model$experiment$cell_info$cell_id),
      dimred_segment_points = model$gold_standard$dimred[!model$gold_standard$meta$burn,],
      dimred_segment_progressions = model$gold_standard$meta[!model$gold_standard$meta$burn,] %>% 
        select(from, to, percentage = time)
    )
}