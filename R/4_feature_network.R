#' @export
feature_network_default <- function(
  realnet_name = sample(realnets$name, 1),
  min_targets_per_tf = 0L,
  damping = 0.05
) {
  lst(
    type = "realnet_sampler",
    realnet_name,
    min_targets_per_tf,
    damping
  )
}

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
  
  # determine number of targets for each tf
  tf_info <- 
    model$feature_info %>% 
    mutate(
      num_targets = .generate_partitions(
        num_elements = model$numbers$num_targets,
        num_groups = model$numbers$num_tfs, 
        min_elements_per_group = model$feature_network_params$min_targets_per_tf
      )
    )
  tf_names <- tf_info$feature_id
  
  # download realnet (~2MB)
  realnet <- .feature_network_fetch_realnet(model)
  
  # verify that realnet is large enough
  assert_that(
    nrow(model$feature_info) <= nrow(realnet),
    msg = paste0("Number of regulators in realnet (", nrow(realnet), ") is not large enough; should be >= ", nrow(model$feature_info))
  )
  
  # map tfs to rownames of realnet randomly
  tf_mapper <- set_names(tf_names, sample(rownames(realnet), length(tf_names)))
  rownames(realnet) <- ifelse(rownames(realnet) %in% names(tf_mapper), tf_mapper[rownames(realnet)], rownames(realnet))
  colnames(realnet) <- ifelse(colnames(realnet) %in% names(tf_mapper), tf_mapper[colnames(realnet)], colnames(realnet))
  
  # sample target network from realnet
  downstream <- .feature_network_sample_downstream(model, tf_info, realnet)
  housekeeping <- .feature_network_sample_housekeeping(model, realnet)
  
  # return output
  model$feature_info <- 
    bind_rows(
      model$feature_info,
      downstream$target_info,
      housekeeping$hk_info
    )
  
  model$feature_network <- 
    bind_rows(
      model$feature_network,
      downstream$target_network,
      housekeeping$hk_network
    )
  
  model
}

.feature_network_fetch_realnet <- function(model) {
  realnet_name <- model$feature_network_params$realnet_name
  
  data(realnets, package = "dyngen", envir = environment())
  assert_that(realnet_name %all_in% realnets$name)
  
  realnet_url <- realnets$url[[match(realnet_name, realnets$name)]]
  
  .download_cacheable_file(realnet_url, model)
}
#' @importFrom Matrix summary
#' @importFrom igraph graph_from_data_frame page_rank E
.feature_network_sample_downstream <- function(model, tf_info, realnet) {
  requireNamespace("igraph")
  
  tf_names <- tf_info$feature_id
  
  # convert to igraph
  gr <- 
    realnet %>% 
    Matrix::summary() %>% 
    as.data.frame() %>% 
    transmute(
      i = rownames(realnet)[i], 
      j = colnames(realnet)[j],
      weight = x
    ) %>%
    igraph::graph_from_data_frame(vertices = colnames(realnet))
  
  # derive targets per tf
  target_regnet <- map2_df(
    tf_info$feature_id,
    tf_info$num_targets,
    function(tf_name, num_targets) {
      if (num_targets > 0) {
        personalized <- rep(0, ncol(realnet)) %>% set_names(colnames(realnet))
        personalized[[tf_name]] <- 1
        
        # calculate page rank from regulator
        page_rank <- igraph::page_rank(
          gr,
          personalized = personalized, 
          directed = TRUE,
          weights = igraph::E(gr)$weight,
          damping = model$feature_network_params$damping
        )
        
        # select top targets
        features_sel <- 
          enframe(page_rank$vector, "feature_id", "score") %>% 
          mutate(score = score + stats::runif(n(), 0, 1e-15)) %>% 
          sample_n(num_targets, weight = score) %>% 
          pull(feature_id) %>% 
          c(tf_name)
        
        # get induced subgraph
        gr %>% 
          igraph::induced_subgraph(features_sel) %>% 
          igraph::ego(10, tf_name, mode = "out") %>% 
          first() %>% 
          names() %>% 
          igraph::induced_subgraph(gr, .) %>% 
          igraph::as_data_frame() %>% 
          as_tibble() %>% 
          select(from, to)
      } else {
        tibble(from = character(0), to = character(0))
      }
    }
  )
  
  # filter duplicates and connections between tfs
  target_regnet <- 
    target_regnet %>% 
    # can contain duplicates
    unique() %>% 
    # remove connections to tfs (to avoid ruining the given module network)
    filter(!to %in% tf_names)
  
  # rename non-tf features
  if (nrow(target_regnet) > 0) {
    target_names <- unique(c(target_regnet$from, target_regnet$to)) %>% setdiff(tf_names)
    target_mapper <- 
      set_names(
        paste0("Target", seq_along(target_names)),
        target_names
      )
    
    target_regnet <- 
      target_regnet %>% 
      mutate_all(~ ifelse(. %in% names(target_mapper), target_mapper[.], .))
  } else {
    target_mapper <- character(0)
  }
  
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
.feature_network_sample_housekeeping <- function(model, realnet) {
  requireNamespace("igraph")
  
  num_hks <- model$numbers$num_hks
  
  # convert to igraph
  gr <- 
    realnet %>% 
    Matrix::summary() %>% 
    as.data.frame() %>% 
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
    filter(from != to) %>% 
    group_by(to) %>% 
    slice(1:4) %>% 
    ungroup()
  
  # rename hk features
  if (nrow(hk_regnet) > 0) {
    hk_mapper <- 
      set_names(
        paste0("HK", seq_along(hk_names)),
        hk_names
      )
    
    hk_regnet <- 
      hk_regnet %>% 
      mutate(from = hk_mapper[from], to = hk_mapper[to])
  } else {
    hk_mapper <- character(0)
  }
  
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
