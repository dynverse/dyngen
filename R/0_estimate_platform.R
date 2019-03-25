#' Estimate a platform from a dataset
#' 
#' Altenatively, [get_simple_platform()] returns a simple platform parameter configuration.
#'  
#' @param counts The counts with cells in rows and genes in columns.
#' @param grouping A named vector representing a grouping of the cells.
#' @param subsample The number of cells to subsample.
#' 
#' @rdname platforms
#' 
#' @export
get_platform_from_counts <- function(counts, grouping, subsample = 500) {
  requireNamespace("splatter")
  
  # add a try catch to the splatEstDropout function because it errors too frequently
  old_fun <- splatter:::splatEstDropout
  new_fun <- function(...) {
    tryCatch({
      old_fun(...)
    }, error = function(e) {
      warning("Could not estimate dropout parameters, defaulting to mid = 0.01 and shape = 1.")
      splatter::setParams(params, dropout.mid = 0.01, dropout.shape = 1)
    })
    
  }
  assignInNamespace("splatEstDropout", new_fun, asNamespace("splatter"))
  on.exit(assignInNamespace("splatEstDropout", old_fun, asNamespace("splatter")))
  
  # remove genes that are not sufficiently expressed
  min.pct <- 0.05
  counts <- counts[, apply(counts, 2, function(x) mean(x > 0) > min.pct), drop = FALSE]
  
  # sample the number of 
  if (!is.null(subsample)) {
    ix <- sample.int(nrow(counts), min(nrow(counts), subsample))
  } else {
    ix <- seq_len(nrow(counts))
  }
  
  # estimate splatter params
  estimate <- splatter::splatEstimate(t(counts[ix, , drop = FALSE]))
  attr(class(estimate), "package") <- NULL # make sure scater won't get loaded when the platform is loaded
      
  # determine how many features change between trajectory stages
  group_ids <- unique(dataset_raw$grouping)
      
  # differential expression using wilcox test
  diffexp <- map_df(group_ids, function(group_id) {
    inside <- dataset_raw$grouping == group_id
    outside <- dataset_raw$grouping != group_id
    
    expression_inside <- log2(counts[inside, ] + 1)
    expression_outside <- log2(counts[outside, ] + 1)
    
    result <- map_df(colnames(expression_inside), function(feature_id) {
      pvalue <- wilcox.test(expression_inside[, feature_id], expression_outside[, feature_id])$p.value
      log2fc <- mean(expression_inside[, feature_id]) - mean(expression_outside[, feature_id])
      
      tibble(
        pvalue = pvalue,
        log2fc = log2fc,
        feature_id = feature_id,
        group_id = group_id
      )
    })
  }) %>%
    mutate(qvalue = p.adjust(pvalue, "fdr"))
  
  qvalue_cutoff <- 0.05
  log2fc_cutoff <- 1
  diffexp_features <- diffexp %>% filter(
    qvalue < qvalue_cutoff,
    abs(log2fc) > log2fc_cutoff
  ) %>%
    pull(feature_id) %>%
    unique()
  
  trajectory_dependent_features <- length(diffexp_features) / ncol(dataset_raw$counts)
  
  # create platform object
  lst(
    estimate,
    trajectory_dependent_features,
    n_cells = nrow(counts),
    n_features = ncol(counts)
  )
}

#' @param n_cells The number of cells
#' @param n_features The number of features
#' @param trajectory_dependent_features The percentage of features that are being driven by the trajectory (or vice versa)
#' @param dropout_mean_rate The mean rate of dropouts
#' @param dropout_mean_shape The shape of dropouts
#' 
#' @rdname platforms
#' 
#' @export
get_simple_platform <- function(
  n_cells = 100L,
  n_features = 100L,
  trajectory_dependent_features = 0.1,
  dropout_rate = 0.01,
  dropout_shape = 1
) {
  list(
    platform_id = "simple",
    estimate = splatter::newSplatParams(mean.rate = dropout_rate, mean.shape = dropout_shape),
    n_cells = n_cells,
    n_features = n_features,
    trajectory_dependent_features = trajectory_dependent_features
  )
}