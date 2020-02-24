#' Wrap a simulation into a dynwrap object
#' 
#' The output of this object can be used with \pkg{dyno}.
#' 
#' @param model A dyngen output model for which the experiment has been emulated with [generate_experiment()].
#' @param store_grn Whether or not to also store GRN information.
#' @param store_dimred Whether or not to store the dimensionality reduction constructed on the true counts.
#' 
#' @export
#' @importFrom dynwrap wrap_expression add_trajectory add_dimred
wrap_dataset <- function(model, store_grn = FALSE, store_dimred = FALSE) {
  dataset <- wrap_expression(
    id = model$id,
    counts = model$experiment$counts_mrna,
    expression = as(log2(model$experiment$counts_mrna + 1), "dgCMatrix"),
    expression_unspliced = as(log2(model$experiment$counts_premrna + 1), "dgCMatrix"),
    expression_protein = as(log2(model$experiment$counts_protein + 1), "dgCMatrix"),
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
  if (store_grn) {
    regulatory_network <- model$feature_network %>% 
      select(regulator = from, target = to, strength, effect)
    regulation_sc <- model$experiment$regulation
    
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
  
  dataset
}

#' Add a GRN to a dynwrap object
#'
#' @inheritParams dynwrap::common_param
#' @param regulatory_network A data frame consisting of three columns: `"regulator"`, `"target"`, `"strength"`.
#' @param regulatory_network_sc A data frame consisting of four columns: `"cell_id"`, `"regulator"`, `"target"`, `"strength"`.
#' @param regulators The feature ids of the regulators.
#' @param targets The feature ids of the targets.
#' @param ... Extra arguments to be saved in the model.
#'
#' @export
add_regulatory_network <- function(dataset, regulatory_network, regulatory_network_sc = NULL, regulators = NULL, targets = NULL, ...) {
  # check regulatory network
  assert_that(
    is.data.frame(regulatory_network),
    regulatory_network %has_names% c("regulator", "target", "strength"),
    is.character(regulatory_network$regulator) || is.factor(regulatory_network$regulator),
    is.character(regulatory_network$target) || is.factor(regulatory_network$target),
    is.numeric(regulatory_network$strength),
    !is.null(regulators),
    !is.null(targets),
    all(regulatory_network$regulator %in% regulators),
    all(regulatory_network$target %in% targets)
  )
  
  if (!is.factor(regulatory_network$regulator)) {
    regulatory_network$regulator <- factor(regulatory_network$regulator, regulators)
  }
  if (!is.factor(regulatory_network$target)) {
    regulatory_network$target <- factor(regulatory_network$target, targets)
  }
  
  # check sc regulatory network
  cell_ids <- dataset$cell_ids
  
  assert_that(
    is.data.frame(regulatory_network_sc),
    regulatory_network_sc %has_names% c("cell_id", "regulator", "target", "strength"),
    is.character(regulatory_network_sc$cell_id) || is.factor(regulatory_network_sc$cell_id),
    is.character(regulatory_network_sc$regulator) || is.factor(regulatory_network_sc$regulator),
    is.character(regulatory_network_sc$target) || is.factor(regulatory_network_sc$target),
    is.numeric(regulatory_network_sc$strength),
    !is.null(dataset$cell_ids),
    all(regulatory_network_sc$cell_id %in% dataset$cell_ids),
    all(regulatory_network_sc$regulator %in% regulators),
    all(regulatory_network_sc$target %in% targets)
  )
  
  if (!is.factor(regulatory_network_sc$cell_id)) {
    regulatory_network_sc$cell_id <- factor(regulatory_network_sc$cell_id, cell_ids)
  }
  if (!is.factor(regulatory_network_sc$regulator)) {
    regulatory_network_sc$regulator <- factor(regulatory_network_sc$regulator, regulators)
  }
  if (!is.factor(regulatory_network_sc$target)) {
    regulatory_network_sc$target <- factor(regulatory_network_sc$target, targets)
  }
  
  dataset <- dataset %>% extend_with(
    "dynwrap::with_regulatory_network",
    regulatory_network = regulatory_network,
    regulatory_network_sc = regulatory_network_sc,
    regulators = regulators,
    targets = targets,
    ...
  )
  
}
