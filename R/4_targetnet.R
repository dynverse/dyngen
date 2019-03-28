#' @export
networkgen_realnet_sampler <- function(
  realnet_name = realnets$name,
  min_targets_per_tf = 5L,
  damping = 0.05
) {
  realnet_name <- match.arg(realnet_name)
  
  data(realnets, package = "dyngen")
  assert_that(realnet_name %all_in% realnets$name)
  
  realnet_url <- realnets$url[[match(realnet_name, realnets$name)]]
  
  lst(
    type = "realnet_sampler",
    realnet_name,
    realnet_url,
    min_targets_per_tf,
    damping
  )
}

#' @export
generate_feature_network <- function(
  model
) {
  if (model$verbose) cat("Sampling feature network from real network\n")
  
  # determine number of targets for each tf
  tf_info <- model$tf_info %>% 
    mutate(
      num_targets = .generate_partitions(
        num_elements = max(model$feature_numbers$num_tfs, model$feature_numbers$num_targets),
        num_groups = model$feature_numbers$num_tfs, 
        min_elements_per_group = model$networkgen_params$min_targets_per_tf
      )
    )
  tf_names <- tf_info$feature_id
  
  # download realnet (~2MB)
  realnet <- .networkgen_fetch_realnet(model)
  
  # map tfs to rownames of realnet randomly
  tf_mapper <- set_names(tf_names, sample(rownames(realnet), length(tf_names)))
  rownames(realnet) <- ifelse(rownames(realnet) %in% names(tf_mapper), tf_mapper[rownames(realnet)], rownames(realnet))
  colnames(realnet) <- ifelse(colnames(realnet) %in% names(tf_mapper), tf_mapper[colnames(realnet)], colnames(realnet))
  
  # sample target network from realnet
  out <- .networkgen_sample_from_realnet(model, tf_info, realnet)

  # return output
  model$feature_info <- 
    bind_rows(
      model$tf_info,
      out$target_info
    )
  
  model$feature_network <- 
    bind_rows(
      model$tf_network,
      out$target_network
    )
  
  model
}

.networkgen_fetch_realnet <- function(model) {
  tmpfile <- tempfile()
  on.exit(file.remove(tmpfile))
  download.file(model$networkgen_params$realnet_url, destfile = tmpfile, quiet = !model$verbose)
  realnet <- read_rds(tmpfile)
  
  assert_that(
    nrow(model$tf_info) <= nrow(realnet),
    msg = paste0("Number of regulators in realnet (", nrow(realnet), ") is not large enough; should be >= ", nrow(model$tf_info))
  )
  
  realnet
}

.networkgen_sample_from_realnet <- function(model, tf_info, realnet) {
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
      personalized <- rep(0, ncol(realnet)) %>% set_names(colnames(realnet))
      personalized[[tf_name]] <- 1
      
      # calculate page rank from regulator
      page_rank <- igraph::page_rank(
        gr,
        personalized = personalized, 
        directed = TRUE,
        weights = igraph::E(gr)$weight,
        damping = model$networkgen_params$damping
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
        select(regulator = from, target = to)
    }
  )
  
  # filter duplicates and connections between tfs
  target_regnet <- 
    target_regnet %>% 
    # can contain duplicates
    unique() %>% 
    # remove connections between tfs (to avoid ruining the given module network)
    filter(!(regulator %in% tf_names && target %in% tf_names))
  
  # rename non-tf features
  tmp_names <- unique(c(target_regnet$regulator, target_regnet$target))
  other_mapper <- setdiff(tmp_names, tf_names) %>% 
  {set_names(paste0("Target", seq_along(.)), .)}
  
  target_regnet <- 
    target_regnet %>% 
    mutate_all(~ ifelse(. %in% names(other_mapper), other_mapper[.], .))
  
  # create target info
  target_info <- 
    tibble(
      feature_id = other_mapper,
      is_tf = FALSE,
      is_trajectory_driven = TRUE,
      a0 = NA, # a0 is decided later based on regulation
      burn = FALSE # extra genes should be available during burn in ; TODO: should this not be the other way around then? 
    )
  
  lst(
    target_network = target_regnet,
    target_info = target_info
  )
}