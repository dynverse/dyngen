#' Convert simulation output to different formats.
#' 
#' For use with other packages compatible with dyno / anndata.
#' 
#' @param model A dyngen output model for which the experiment has been emulated with [generate_experiment()].
#' @param store_cellwise_grn Whether or not to also store cellwise GRN information.
#' @param store_dimred Whether or not to store the dimensionality reduction constructed on the true counts.
#' @param store_rna_velocity WHether or not to store the log propensity ratios.
#' 
#' @return A dyno/anndata object
#' 
#' @export
#' @rdname convert
#' 
#' @examples
#' model <- initialise_model(
#'   backbone = backbone_bifurcating()
#' )
#' \dontshow{
#' # actually use a smaller example 
#' # to reduce execution time during
#' # testing of the examples
#' model <- initialise_model(
#'   backbone = model$backbone,
#'   num_cells = 5,
#'   num_targets = 0,
#'   num_hks = 0,
#'   gold_standard_params = gold_standard_default(census_interval = 1, tau = 0.1),
#'   simulation_params = simulation_default(
#'     burn_time = 10,
#'     total_time = 10,
#'     census_interval = 1,
#'     ssa_algorithm = ssa_etl(tau = 0.1),
#'     experiment_params = simulation_type_wild_type(num_simulations = 1)
#'   )
#' )
#' }
#' \donttest{
#' model <- model %>%
#'   generate_tf_network() %>%
#'   generate_feature_network() %>%
#'   generate_kinetics() %>%
#'   generate_gold_standard() %>%
#'   generate_cells() %>%
#'   generate_experiment()
#'   
#' dataset <- as_dyno(model)
#' 
#' # dynplot::plot_dimred(dataset)
#' }
as_dyno <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  requireNamespace("dynwrap")
  
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  dataset <- dynwrap::wrap_expression(
    id = model$id,
    counts = counts,
    counts_spliced = model$experiment$counts_mrna,
    counts_unspliced = model$experiment$counts_premrna,
    counts_protein = model$experiment$counts_protein,
    expression = as(log2(counts + 1), "dgCMatrix"),
    cell_info = model$experiment$cell_info %>% select(-.data$from, -.data$to, -.data$time),
    feature_info = model$experiment$feature_info
  ) %>% 
    dynwrap::add_trajectory(
      milestone_network = model$gold_standard$network,
      progressions = model$experiment$cell_info %>%
        select(.data$cell_id, .data$from, .data$to, percentage = .data$time)
    )
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    
    dataset <- dataset %>% 
      dynwrap::add_dimred(
        dimred = dimred,
        dimred_segment_points = model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE],
        dimred_segment_progressions = model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE] %>% 
          select(.data$from, .data$to, percentage = .data$time)
      )
  }
  
  # add a few more values
  if (store_cellwise_grn) {
    regulatory_network <- model$feature_network %>% 
      select(regulator = .data$from, target = .data$to, .data$strength, .data$effect)
    regulation_sc <- model$experiment$cellwise_grn
    
    regulators <- unique(regulatory_network$regulator)
    targets <- colnames(dataset$counts)
    regulatory_network_sc <- 
      Matrix::summary(regulation_sc) %>% 
      transmute(
        cell_id = factor(dataset$cell_ids[.data$i], levels = dataset$cell_ids),
        regulator = factor(regulatory_network$regulator[.data$j], levels = regulators),
        target = factor(regulatory_network$target[.data$j], levels = targets),
        strength = .data$x
      ) %>% 
      as_tibble()
    
    
    dataset <- dataset %>% 
      dynwrap::add_regulatory_network(
        regulatory_network = regulatory_network,
        regulatory_network_sc = regulatory_network_sc,
        regulators = regulators,
        targets = targets
      )
  }
  
  if (store_rna_velocity) {
    dataset <- dataset %>% extend_with(
      "dynwrap::with_rna_velocity",
      rna_velocity = model$experiment$rna_velocity
    )
  }
  
  dataset
}

#' @rdname convert
#' @export
wrap_dataset <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  .Deprecated("as_dyno")
  as_dyno(
    model = model, 
    store_dimred = store_dimred,
    store_cellwise_grn = store_cellwise_grn,
    store_rna_velocity = store_rna_velocity
  )
}


#' @rdname convert
#' @importFrom tibble column_to_rownames
#' @export
as_anndata <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  requireNamespace("anndata")
  
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  
  ad <- anndata::AnnData(
    X = counts,
    layers = list(
      counts_spliced = model$experiment$counts_mrna,
      counts_unspliced = model$experiment$counts_premrna,
      counts_protein = model$experiment$counts_protein
    ),
    obs = model$experiment$cell_info %>%
      select(-.data$from, -.data$to, -.data$time) %>% 
      as.data.frame() %>% 
      column_to_rownames("cell_id"),
    var = model$experiment$feature_info %>% 
      as.data.frame() %>% 
      column_to_rownames("feature_id"),
    uns = list(
      traj_milestone_network = model$gold_standard$network,
      traj_progressions = model$experiment$cell_info %>%
        select(.data$cell_id, .data$from, .data$to, percentage = .data$time) %>% 
        as.data.frame()
    )
  )
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    ad$obsm$dimred <- dimred
    
    ad$uns$traj_dimred_segments <- bind_cols(
      model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE],
      as.data.frame(model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE])
    )
  }
  
  if (store_cellwise_grn) {
    ad$uns$regulatory_network <- model$feature_network %>% 
      select(regulator = .data$from, target = .data$to, .data$strength, .data$effect)
    ad$obsm$regulatory_network_sc <- model$experiment$cellwise_grn
    ad$uns$regulatory_network_regulators <- unique(model$feature_network$from)
    ad$uns$regulatory_network_targets <- colnames(counts)
  }
  
  if (store_rna_velocity) {
    ad$obsm$rna_velocity <- model$experiment$rna_velocity
  }
  
  ad
}