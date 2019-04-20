#' @export
calculate_dimred <- function(model) {
  # check whether the simulations have been run
  has_sim <- model %has_name% "simulations" && model$simulations %has_name% "counts"
  if (has_sim) {
    sim_counts <- model$simulations$counts
    sim_ix <- seq_len(nrow(sim_counts))
  } else {
    sim_counts <- NULL
    sim_ix <- numeric(0)
  }
  
  # check whether the gold standard has been run
  has_gs <- model %has_name% "gold_standard" && model$gold_standard %has_name% "counts"
  if (has_gs) {
    gs_counts <- model$gold_standard$counts
    gs_ix <- seq(length(sim_ix) + 1, length(sim_ix) + nrow(gs_counts))
    
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
  tf_molecules <- tf_info %>% select(w, x, y) %>% gather(col, val) %>% pull(val)
  
  if (!is.null(sim_counts)) {
    sim_counts <- sim_counts[, tf_molecules, drop = FALSE]
  }
  
  if (!is.null(gs_counts)) {
    gs_counts <- gs_counts[, tf_molecules, drop = FALSE]
  }
  
  # combine data and select landmarks
  counts <- rbind(sim_counts, gs_counts)
  
  # normalise
  max_cols <- apply(counts, 2, quantile, .99)
  max_cols[max_cols == 0] <- 1
  counts <- sweep(counts, 2, max_cols, "/") %>% 
    Matrix::Matrix(sparse = TRUE)
  
  landmark_ix <- 
    if (length(gs_ix) > 0) {
      if (length(gs_ix) > 1000) {
        sample(gs_ix, 1000)
      } else {
        gs_ix
      }
    } else {
      if (length(sim_ix) > 1000) {
        sample(sim_ix, 1000)
      } else {
        sim_ix
      }
    }
  
  # calculate distances to lndmarks
  dist_2lm <- as.matrix(dynutils::calculate_distance(
    x = counts[landmark_ix, , drop = FALSE], 
    y = counts, 
    method = model$distance_metric
  ))
  
  # calculate distances between landmarks
  dist_lm <- dist_2lm[, landmark_ix, drop = FALSE]
  
  # calculate dimred
  dimred <- as.matrix(dyndimred:::.lmds_cmdscale(
    dist_lm, 
    dist_2lm, 
    ndim = 3, 
    rescale = TRUE
  ))
  dimnames(dimred) <- list(rownames(counts), paste0("comp_", seq_len(ncol(dimred))))
  
  # separate out sim dimred and gs dimred
  if (has_sim) {
    model$simulations$dimred <- dimred[sim_ix, , drop = FALSE]
  }
  if (has_gs) {
    model$gold_standard$dimred <- dimred[gs_ix, , drop = FALSE]
  }
  
  # return model
  model
}