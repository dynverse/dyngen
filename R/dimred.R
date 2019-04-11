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
    gs_ix <- seq_len(nrow(counts))[-sim_ix]
    landmark_ix <- if (length(gs_ix) > 1000) sample(gs_ix, 1000) else gs_ix
  } else {
    counts <- sim_counts
    landmark_ix <- if (length(sim_ix) > 1000) sample(sim_ix, 1000) else sim_ix
  }
  
  counts <- counts[, grep("TF", colnames(counts)), drop = FALSE]
  
  dist_metrics <- dynutils::list_distance_metrics()
  dist_metric <- model$dist_metric
  assert_that(dist_metric %all_in% dist_metrics)
  
  dist_2lm <- dynutils::calculate_distance(counts[landmark_ix, , drop = FALSE], counts, metric = dist_metric)
  dist_lm <- dist_2lm[, landmark_ix, , drop = FALSE]
  dimred <- dyndimred:::.lmds_cmdscale(dist_lm, dist_2lm, ndim = 3, rescale = TRUE)
  attr(dimred, "landmark_space") <- NULL
  dimred <- dyndimred:::process_dimred(dimred)
  
  model$simulations$dimred <- dimred[sim_ix, , drop = FALSE]
  if (has_gs) {
    model$goldstandard$dimred <- dimred[-sim_ix, , drop = FALSE]
  }
  
  model
}