#' @export
calculate_dimred <- function(model) {
  if (!model %has_name% "simulations") {
    stop("First run some simulations before calculating the dimred")
  }
  
  sim_counts <- model$simulations$counts
  sim_ix <- seq_len(nrow(sim_counts))
  
  has_gs <- model %has_name% "goldstandard" && model$goldstandard %has_name% "counts"
  if (has_gs) {
    gs_counts <- model$goldstandard$counts
    counts <- Matrix::rbind2(sim_counts, gs_counts)
  } else {
    counts <- sim_counts
  }
  
  
  dimred <- dyndimred::dimred_landmark_mds(counts, distance_metric = model$simulation_params$dimred_method)
  
  model$simulations$dimred <- dimred[sim_ix, , drop = FALSE]
  if (has_gs) {
    model$goldstandard$dimred <- dimred[-sim_ix, , drop = FALSE]
  }
  
  model
}