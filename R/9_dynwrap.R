#' Wrap a simulation into a dynwrap object
#' 
#' The output of this object can be used with \pkg{dyno}.
#' 
#' @param model A dyngen output model for which the experiment has been emulated with [generate_experiment()].
#' @param store_cellwise_grn Whether or not to also store cellwise GRN information.
#' @param store_dimred Whether or not to store the dimensionality reduction constructed on the true counts.
#' @param store_rna_velocity WHether or not to store the log propensity ratios.
#' 
#' @export
#' @importFrom dynwrap wrap_expression add_trajectory add_dimred add_regulatory_network
wrap_dataset <- function(
  model,
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_dimred = !is.null(model$experiment$dimred),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  dataset <- wrap_expression(
    id = model$id,
    counts = counts,
    counts_spliced = model$experiment$counts_mrna,
    counts_unspliced = model$experiment$counts_premrna,
    counts_protein = model$experiment$counts_protein,
    expression = as(log2(counts + 1), "dgCMatrix"),
    cell_info = model$experiment$cell_info %>% select(-from, -to, -time),
    feature_info = model$experiment$feature_info
  ) %>% 
    add_trajectory(
      milestone_network = model$gold_standard$network,
      progressions = model$experiment$cell_info %>% select(cell_id, from, to, percentage = time)
    )
  
  if (store_dimred) {
    dataset <- dataset %>% 
      add_dimred(
        dimred = model$simulations$dimred[model$experiment$cell_info$step_ix, ] %>% 
          magrittr::set_rownames(model$experiment$cell_info$cell_id),
        dimred_segment_points = model$gold_standard$dimred[!model$gold_standard$meta$burn,],
        dimred_segment_progressions = model$gold_standard$meta[!model$gold_standard$meta$burn,] %>% 
          select(from, to, percentage = time)
      )
  }
  
  # add a few more values
  if (store_cellwise_grn) {
    regulatory_network <- model$feature_network %>% 
      select(regulator = from, target = to, strength, effect)
    regulation_sc <- model$experiment$cellwise_grn
    
    regulators <- unique(regulatory_network$regulator)
    targets <- colnames(dataset$counts)
    regulatory_network_sc <- 
      Matrix::summary(regulation_sc) %>% 
      transmute(
        cell_id = factor(dataset$cell_ids[i], levels = dataset$cell_ids),
        regulator = factor(regulatory_network$regulator[j], levels = regulators),
        target = factor(regulatory_network$target[j], levels = targets),
        strength = x
      ) %>% 
      as_tibble()
    
    
    dataset <- dataset %>% 
      add_regulatory_network(
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