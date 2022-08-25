#' @importFrom lmds lmds
calculate_dimred <- function(
  model, 
  num_landmarks = 1000,
  dimred_simulations = TRUE,
  dimred_gold = TRUE
) {
  # satisfy r cmd check
  is_tf <- NULL
  
  # set rcpp thread options to model$num_cores
  prev_num_cores <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = model$num_cores)
  
  # check whether the simulations have been run
  has_sim <- dimred_simulations && model %has_name% "simulations" && model$simulations %has_name% "counts"
  if (has_sim) {
    sim_counts <- model$simulations$counts
    sim_ix <- seq_len(nrow(sim_counts))
  } else {
    sim_counts <- NULL
    sim_ix <- numeric(0)
  }
  
  # check whether the gold standard has been run
  has_gs <- dimred_gold && model %has_name% "gold_standard" && model$gold_standard %has_name% "counts"
  if (has_gs) {
    gs_counts <- model$gold_standard$counts
    gs_ix <- length(sim_ix) + seq_len(nrow(gs_counts))
  } else {
    gs_counts <- NULL
    gs_ix <- numeric(0)
  }
  
  # throw error if neither have run
  if (!has_gs && !has_sim) {
    stop("First generate the gold standard or the simulations before calculating a dimred")
  }
  
  # use only tf data
  tf_info <- model$feature_info %>% filter(is_tf)
  
  if (!is.null(sim_counts)) {
    sim_wcounts <- sim_counts[, tf_info$mol_premrna, drop = FALSE]
    sim_xcounts <- sim_counts[, tf_info$mol_mrna, drop = FALSE]
    sim_ycounts <- sim_counts[, tf_info$mol_protein, drop = FALSE]
  } else {
    sim_wcounts <- sim_xcounts <- sim_ycounts <- NULL
  }
  
  if (!is.null(gs_counts)) {
    gs_wcounts <- gs_counts[, tf_info$mol_premrna, drop = FALSE]
    gs_xcounts <- gs_counts[, tf_info$mol_mrna, drop = FALSE]
    gs_ycounts <- gs_counts[, tf_info$mol_protein, drop = FALSE]
  } else {
    gs_wcounts <- gs_xcounts <- gs_ycounts <- NULL
  }
  
  # combine data and select landmarks
  # WARNING: do not change the order of this rbind without changing the sim_ix and gs_ix objects
  # (and vice versa)
  counts <- rbind(
    sim_xcounts,
    gs_xcounts
  )
  
  # log2 transformation
  counts@x <- log2(counts@x + 1)
  
  # normalise
  max_cols <- apply(counts, 2, quantile, .99)
  max_cols[max_cols == 0] <- 1
  counts <- sweep(counts, 2, max_cols, "/") %>%
    Matrix::Matrix(sparse = TRUE)

  # TODO: sample gold standard ix with dynwrap
  # waypoint selection method?
  
  # sample matching indices from w and x  
  # WARNING: do not change the order of this rbind without changing the sim_ix and gs_ix objects
  landmark_ix <- c(gs_ix, sim_ix)
  
  if (length(landmark_ix) > num_landmarks) {
    landmark_ix <- sample(landmark_ix, num_landmarks)
  }
  
  # calculate distances to lndmarks
  dist_2lm <- as.matrix(suppressWarnings(dynutils::calculate_distance(
    x = counts[landmark_ix, , drop = FALSE], 
    y = counts, 
    method = model$distance_metric
  )))
  
  # # calculate distances between landmarks
  # dist_lm <- dist_2lm[, landmark_ix, drop = FALSE]
  attr(dist_2lm, "landmark_ix") <- landmark_ix
  
  # calculate dimred
  dimred <- lmds::cmdscale_landmarks(dist_2lm, ndim = 3)
  dimnames(dimred) <- list(rownames(counts), paste0("comp_", seq_len(ncol(dimred))))
  
  # separate out sim dimred and gs dimred
  if (has_sim) {
    model$simulations$dimred <- dimred[sim_ix, , drop = FALSE]
  }
  if (has_gs) {
    model$gold_standard$dimred <- dimred[gs_ix, , drop = FALSE]
  }
  
  # restore previous setting
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = prev_num_cores)
  
  # return model
  model
}