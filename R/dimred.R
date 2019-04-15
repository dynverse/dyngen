#' @export
calculate_dimred <- function(model) {
  if (!model %has_name% "simulations") {
    stop("First run some simulations before calculating the dimred")
  }
    
  sim_counts <- model$simulations$counts
  sim_ix <- seq_len(nrow(sim_counts))
  
  has_gs <- model %has_name% "gold_standard" && model$gold_standard %has_name% "counts"
  if (has_gs) {
    gs_counts <- model$gold_standard$counts
    counts <- rbind(sim_counts, gs_counts)
    gs_ix <- seq_len(nrow(counts))[-sim_ix]
    landmark_ix <- if (length(gs_ix) > 1000) sample(gs_ix, 1000) else gs_ix
  } else {
    counts <- sim_counts
    landmark_ix <- if (length(sim_ix) > 1000) sample(sim_ix, 1000) else sim_ix
  }
  
  counts <- counts[, grep("TF", colnames(counts)), drop = FALSE]
  
  dist_2lm <- as.matrix(dynutils::calculate_distance(
    x = counts[landmark_ix, , drop = FALSE], 
    y = counts, 
    metric = model$distance_metric
  ))
  dist_lm <- dist_2lm[, landmark_ix, drop = FALSE]
  dimred <- as.matrix(dyndimred:::.lmds_cmdscale(
    dist_lm, 
    dist_2lm, 
    ndim = 3, 
    rescale = TRUE
  ))
  dimnames(dimred) <- list(rownames(counts), paste0("comp_", seq_len(ncol(dimred))))
  
  model$simulations$dimred <- dimred[sim_ix, , drop = FALSE]
  if (has_gs) {
    model$gold_standard$dimred <- dimred[-sim_ix, , drop = FALSE]
  }
  
  model
}