#' Generate a target network 
#' 
#' [generate_feature_network()] generates a network of target genes that are regulated
#' by the previously generated TFs, and also a separate network of housekeeping genes (HKs).
#' [feature_network_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the transcription network has been generated with [generate_tf_network()].
#' @param realnet The name of a gene regulatory network (GRN) in [realnets]. 
#'   Alternatively, a custom GRN can be used by passing a weighted sparse matrix.
#' @param damping A damping factor used for the page rank algorithm used to subsample the realnet.
#' @param target_resampling How many targets / HKs to sample from the realnet per iteration.
#' @param max_in_degree The maximum in-degree of a target / HK.
#' 
#' @export
generate_feature_network <- function(
  model
) {
  if (model$verbose) cat("Sampling feature network from real network\n")
  
  # verify that this function has all the data it needs
  assert_that(
    model %has_names% c("feature_info", "feature_network"),
    msg = "Execute generate_tf_network() before executing generate_feature_network()"
  )
  
  # process realnet
  realnet <- model$feature_network_params$realnet
  
  if (is.character(realnet)) {
    # download realnet (~2MB)
    realnet <- .feature_network_fetch_realnet(model)
  }
  
  # check realnet
  assert_that(
    dynutils::is_sparse(realnet),
    !is.null(rownames(realnet)),
    !is.null(colnames(realnet))
  )
  assert_that(
    nrow(model$feature_info) <= nrow(realnet),
    msg = paste0("Number of regulators in realnet (", nrow(realnet), ") is not large enough (>= ", nrow(model$feature_info), ")")
  )
  
  sample_tfs_per <- min(model$feature_network_params$target_resampling, model$numbers$num_targets)
  sample_hks_per <- min(model$feature_network_params$target_resampling, model$numbers$num_hks)
  ncol_at_least <- model$numbers$num_tfs + max(sample_tfs_per, sample_hks_per)
  assert_that(
    ncol(realnet) >= ncol_at_least,
    msg = paste0(
      "Number of genes in realnet (", ncol(realnet), ") is not large enough (>= ", ncol_at_least, "). ",
      "Choose a larger realnet or make use of the `feature_network_params$target_resampling` parameter."
    )
  )
  
  # sample target network from realnet
  if (sample_tfs_per > 0) {
    num_target_start <- seq(1, model$numbers$num_targets, by = sample_tfs_per)
    num_target_stop <- c((num_target_start - 1) %>% tail(-1), model$numbers$num_targets)
    downstreams <- 
      pmap(lst(start = num_target_start, stop = num_target_stop), function(start, stop) {
        .feature_network_sample_downstream(model, realnet, num_targets = stop - start + 1, target_index_offset = start - 1)
      })
    target_info <- map_df(downstreams, "target_info")
    target_network <- map_df(downstreams, "target_network")
  } else {
    target_info <- NULL
    target_network <- NULL
  }
  
  # sample house keeping
  if (sample_hks_per > 0) {
    num_hk_start <- seq(1, model$numbers$num_hks, by = sample_hks_per)
    num_hk_stop <- c((num_hk_start - 1) %>% tail(-1), model$numbers$num_hks)
    housekeepings <- 
      pmap(lst(start = num_hk_start, stop = num_hk_stop), function(start, stop) {
        .feature_network_sample_housekeeping(model, realnet, num_hks = stop - start + 1, hk_index_offset = start - 1)
      })
    hk_info <- map_df(housekeepings, "hk_info")
    hk_network <- map_df(housekeepings, "hk_network")
  } else {
    hk_info <- NULL
    hk_network <- NULL
  }
  
  # return output
  model$feature_info <- 
    bind_rows(
      model$feature_info,
      target_info,
      hk_info
    )
  
  model$feature_network <- 
    bind_rows(
      model$feature_network,
      target_network,
      hk_network
    )
  
  model
}

#' @export
#' @rdname generate_feature_network
feature_network_default <- function(
  realnet = sample(realnets$name, 1),
  damping = 0.05,
  target_resampling = Inf,
  max_in_degree = 5
) {
  lst(
    realnet,
    damping,
    target_resampling,
    max_in_degree
  )
}

.feature_network_fetch_realnet <- function(model) {
  realnet <- model$feature_network_params$realnet
  
  data(realnets, package = "dyngen", envir = environment())
  assert_that(realnet %all_in% realnets$name)
  
  realnet_url <- realnets$url[[match(realnet, realnets$name)]]
  
  .download_cacheable_file(
    url = realnet_url, 
    cache_dir = model$download_cache_dir, 
    verbose = model$verbose
  )
}
#' @importFrom Matrix summary
#' @importFrom igraph graph_from_data_frame page_rank E
.feature_network_sample_downstream <- function(
  model,
  realnet,
  num_targets = model$numbers$num_targets, 
  target_index_offset = 0
) {
  requireNamespace("igraph")
  
  # determine desired number of targets for each tf
  tf_info <- 
    model$feature_info %>% 
    mutate(
      num_targets = .generate_partitions(
        num_elements = num_targets,
        num_groups = n(), 
        min_elements_per_group = 0
      )
    )
  tf_names <- tf_info$feature_id
  
  # map tfs to rownames of realnet randomly
  tf_mapper <- set_names(tf_names, sample(rownames(realnet), length(tf_names)))
  rownames(realnet) <- ifelse(rownames(realnet) %in% names(tf_mapper), tf_mapper[rownames(realnet)], rownames(realnet))
  colnames(realnet) <- ifelse(colnames(realnet) %in% names(tf_mapper), tf_mapper[colnames(realnet)], colnames(realnet))
  
  # convert to igraph
  gr <- 
    realnet %>% 
    Matrix::summary() %>% 
    as.data.frame() %>% 
    group_by(j) %>% 
    sample_n(min(n(), sample(1:model$feature_network_params$max_in_degree, 1)), weight = x) %>% 
    ungroup() %>% 
    transmute(
      i = rownames(realnet)[i], 
      j = colnames(realnet)[j],
      weight = x
    ) %>%
    igraph::graph_from_data_frame(vertices = colnames(realnet))
  
  personalized <- rep(0, ncol(realnet)) %>% set_names(colnames(realnet))
  personalized[tf_info$feature_id] <- tf_info$num_targets
  page_rank <- igraph::page_rank(
    gr,
    personalized = personalized, 
    directed = TRUE,
    weights = igraph::E(gr)$weight,
    damping = model$feature_network_params$damping
  )
  features_sel <- 
    enframe(page_rank$vector, "feature_id", "score") %>% 
    mutate(score = score + stats::runif(n(), 0, 1e-15)) %>% 
    filter(!feature_id %in% tf_names) %>% 
    sample_n(num_targets, weight = score) %>% 
    pull(feature_id) %>% 
    unique()
  
  # get induced subgraph
  subgr <-
    gr %>% 
    igraph::induced_subgraph(c(features_sel, tf_names))
  target_regnet <- 
    subgr %>% 
    igraph::as_data_frame() %>% 
    as_tibble() %>% 
    select(from, to) %>% 
    # remove connections to tfs (to avoid ruining the given module network)
    filter(!to %in% tf_names)
  
  # some targets may be missing
  missing_targets <- setdiff(features_sel, target_regnet$to)
  if (length(missing_targets) > 0) {
    deg <- igraph::degree(subgr)
    new_regs <- sample(names(deg), length(missing_targets), prob = deg, replace = TRUE)
    target_regnet <- 
      bind_rows(
        target_regnet, 
        tibble(from = new_regs, to = missing_targets)
      )
  }
  
  # rename non-tf features
  target_names <- unique(c(target_regnet$from, target_regnet$to)) %>% setdiff(tf_names)
  target_mapper <- 
    set_names(
      paste0("Target", seq_along(target_names) + target_index_offset),
      target_names
    )
  
  target_regnet <- 
    target_regnet %>% 
    mutate_all(~ ifelse(. %in% names(target_mapper), target_mapper[.], .))

  # create target info
  target_info <- 
    tibble(
      feature_id = target_mapper,
      is_tf = FALSE,
      is_hk = FALSE,
      burn = TRUE # extra genes should be available during burn in
    )
  
  lst(
    target_network = target_regnet,
    target_info = target_info
  )
}

#' @importFrom Matrix summary
#' @importFrom igraph graph_from_data_frame page_rank E
.feature_network_sample_housekeeping <- function(
  model,
  realnet,
  num_hks = model$numbers$num_hks, 
  hk_index_offset = 0
) {
  requireNamespace("igraph")
  
  # convert to igraph
  gr <- 
    realnet %>% 
    Matrix::summary() %>% 
    as.data.frame() %>% 
    group_by(j) %>% 
    sample_n(min(n(), sample(1:model$feature_network_params$max_in_degree, 1)), weight = x) %>% 
    ungroup() %>% 
    transmute(
      i = rownames(realnet)[i], 
      j = colnames(realnet)[j],
      weight = x
    ) %>%
    igraph::graph_from_data_frame(vertices = colnames(realnet))
  
  hk_names <- gr %>% 
    igraph::bfs(sample.int(ncol(realnet), 1), neimode = "all") %>% 
    .[["order"]] %>% 
    .[seq_len(num_hks)] %>% 
    names()
  
  hk_regnet <- 
    gr %>% 
    igraph::induced_subgraph(hk_names) %>% 
    igraph::as_data_frame() %>% 
    as_tibble() %>% 
    select(from, to) %>% 
    filter(from != to)
  
  # rename hk features
  hk_mapper <- 
    set_names(
      paste0("HK", seq_along(hk_names) + hk_index_offset),
      hk_names
    )
  
  hk_regnet <- 
    hk_regnet %>% 
    mutate(from = hk_mapper[from], to = hk_mapper[to])
  
  # create target info
  hk_info <- 
    tibble(
      feature_id = hk_mapper,
      is_tf = FALSE,
      is_hk = TRUE,
      burn = TRUE # extra genes should be available during burn in
    )
  
  lst(
    hk_network = hk_regnet,
    hk_info = hk_info
  )
}
