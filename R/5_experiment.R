#' Extract cells from a simulation, either simulating a snapshot experiment or a synchronized experiment
#' 
#' @param simulation List containing the expression of simulation(s)
#' @param model The model of the simulation
#' @param samplesettings List containing the information of how the cells should be sampled
#' 
#' @export
take_experiment_cells <- function(simulation, model, samplesettings = list(type = "snapshot", ncells = 500)) {
  # sample
  if (samplesettings$type == "snapshot"){
    
  } else if (samplesettingss$type == "synchronized") {
  }
  
  # get cellinfo
  rownames(expression) <- paste0("C", seq_len(nrow(expression)))
  cellinfo <- cellinfo %>% mutate(cell_id = rownames(expression))
  
  # filter
  experiment <- filter_expression(expression, cellinfo, model$geneinfo)
  
  experiment
}

sample_snapshot <- function(simulation, gs, ncells = 500) {
  sample_ids <- gs$progressions %>% 
    filter(!burn) %>% 
    group_by(step_id) %>% 
    summarise() %>% 
    sample_n(ncells) %>% 
    pull(step_id)
  expression <- simulation$expression[sample_ids, ]
  
  lst(expression, cellinfo = tibble(step_id = sample_ids))
}
snapshot_sampler <- function(ncells=10) {function(simulation, gs) {sample_synchronized(simulation, gs, ncells=ncells)}}

sample_synchronized <- function(simulation, gs, ntimepoints = 10, timepoints = seq(0, max(simulation$stepinfo$simulationtime), length.out=ntimepoints), ncells_per_timepoint = 12) {
  ncells_per_timepoint <- min(ncells_per_timepoint, length(unique(simulation$stepinfo$simulation_id)))
  
  non_burn_step_ids <- gs$progressions %>% 
    filter(!burn) %>% 
    pull(step_id) %>% 
    unique()
  
  stepinfo <- simulation$stepinfo %>% filter(step_id %in% non_burn_step_ids)
  
  sample_step_info <- map_dfr(timepoints, function(timepoint) {
    stepinfo %>% group_by(simulation_id) %>% 
      summarise(step_id = step_id[which.min(abs(simulationtime - timepoint))]) %>% 
      mutate(timepoint = timepoint) %>% 
      sample_n(ncells_per_timepoint)
  })
  
  lst(
    expression = simulation$expression[sample_step_info$step_id, ],
    cellinfo = sample_step_info
  )
}
ssynchronized_sampler <- function(ntimepoints=10) {function(simulation, gs) sample_synchronized(simulation, gs, ntimepoints=ntimepoints)}

#' Checks the expression for certain properties
#' 
#' @param expression Expression matrix
#' @export
check_expression <- function(expression) {
  checks <- list(
    contains_na = any(is.na(expression)),
    contains_zero_cells = any(apply(expression, 1, max) == 0),
    contains_zero_genes = any(apply(expression, 2, max) == 0),
    contains_nonchanging_cells = any(is.na(apply(expression, 1, sd))),
    contains_nonchanging_genes = any(is.na(apply(expression, 2, sd)))
  )
  
  checks
}

#' Filter expression
#' 
#' @param expression Expression matrix
#' @param cellinfo Cell info dataframe
#' @param geneinfo Gene info dataframe
#' 
#' @export
filter_expression <- function(expression, cellinfo, geneinfo) {
  remove_cells <- (apply(expression, 1, max) == 0) | is.na(apply(expression, 1, sd))
  
  expression <- expression[!remove_cells, ]
  cellinfo <- cellinfo %>% slice(match(rownames(expression), cell_id))
  lst(expression, cellinfo, geneinfo)
}



#' Get housekeeping reference means
#' 
#' @param counts Expression matrix containing counts
#' @export
get_housekeeping_reference_means <- function(counts) colMeans(counts)

#' Add housekeeping genes
#' 
#' @param expression The original expression data.
#' @param geneinfo The original gene info
#' @param housekeeping_reference_means The mean expression of a set of genes in the reference dataset
#' @param n_housekeeping_genes The number of genes to add
#' @param overallaverage Overall average expression of the original dataset, this keeps the overall average expression the same even with housekeeping genes
#' @param gene_id_generator Function to generate gene_ids
#' 
#' @export
#' @importFrom utils data
#' @importFrom magrittr set_colnames
add_housekeeping_poisson <- function(
  expression, 
  geneinfo, 
  housekeeping_reference_means, 
  n_housekeeping_genes=200, 
  overallaverage = mean(expression),
  gene_id_generator = function(n) {paste0("GH", seq_len(n))}
) {
  if(is.null(housekeeping_reference_means)) stop("Reference means required!!")
  
  meanpoissons <- overallaverage * housekeeping_reference_means/mean(housekeeping_reference_means)
  
  additional_expression <- purrr::map(sample(meanpoissons, n_housekeeping_genes), ~rpois(nrow(expression), .)) %>%
    invoke(cbind, .) %>% 
    magrittr::set_colnames(gene_id_generator(n_housekeeping_genes))
  
  geneinfo <- dplyr::bind_rows(
    geneinfo %>% dplyr::mutate(housekeeping=FALSE), 
    tibble(gene=colnames(additional_expression), housekeeping=TRUE)
  )
  
  list(expression=cbind(expression, additional_expression), geneinfo=geneinfo)
}