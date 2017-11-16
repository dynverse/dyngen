# Experiment =====================================================
#' Run an experiment from a simulation
#' @param simulation The simulation
#' @param gs The gold standard
#' @param sampler Function telling how the cells should be sampled
#' @param platform Platform
run_experiment <- function(
  simulation, 
  gs,
  sampler,
  platform
) {
  # mutate <- dplyr::mutate
  # filter <- dplyr::filter
  
  # first sample the cells from the sample, using the given number of cells from the platform
  n_cells <- platform$n_cells
  sampled <- sampler(simulation, gs, n_cells)
  expression_simulated <- sampled$expression
  rownames(expression_simulated) <- paste0("C", seq_len(nrow(expression_simulated)))
  cellinfo <- sampled$cellinfo %>% mutate(cell_id = rownames(expression_simulated))
  
  n_genes_simulated <- ncol(expression_simulated)
  
  # generate housekeeping expression
  # number of genes housekeeping depends on the fraction in the reference
  n_genes_housekeeping <- round((n_genes_simulated / platform$pct_changing) * (1 - platform$pct_changing))
  n_genes <- n_genes_simulated + n_genes_housekeeping
  
  # we now extract the splatter estimates
  estimate <- platform$estimate
  estimate@nGenes <- n_genes_housekeeping
  
  estimate@nCells <- n_cells;estimate@groupCells <- n_cells
  
  housekeeping_simulation <- splatter::splatSimulateSingle(estimate)
  
  # we only use the earlier steps in the splatSImulateSingle, but then you need to dig deep into splat::: 's ...
  # we now combine the genemeans from splatter with the simulated expression values
  # then use the libsizes from splatter to estimate the "true" expression from each cell, which will then be used to estimate the tre counts
  
  # see splatter:::splatSimSingleCellMeans
  exp.lib.sizes <- Biobase::pData(housekeeping_simulation)$ExpLibSize
  cell.means.gene <- rep(Biobase::fData(housekeeping_simulation)$GeneMean, n_cells) %>% matrix(ncol=n_cells)
  cell.means.gene <- rbind(cell.means.gene, t(expression_simulated / mean(expression_simulated) * mean(cell.means.gene)))
  cell.props.gene <- t(t(cell.means.gene)/colSums(cell.means.gene))
  expression <- t(t(cell.props.gene) * exp.lib.sizes)
  geneinfo <- tibble(gene_id = ifelse(rownames(expression) == "", paste0("H", seq_len(nrow(expression))), rownames(expression)), housekeeping = rownames(expression) == "")
  rownames(expression) <- geneinfo$gene_id
  
  # see splatter:::splatSimTrueCounts
  true_counts <- matrix(rpois(n_genes * n_cells, lambda = expression), nrow = n_genes, ncol = n_cells)
  dimnames(true_counts) <- dimnames(expression)
  
  # finally, if present, dropouts will be simulated
  # see splatter:::splatSimDropout
  logistic <- function (x, x0, k) {1/(1 + exp(-k * (x - x0)))}
  if (estimate@dropout.present | TRUE) {
    drop.prob <- sapply(seq_len(n_cells), function(idx) {
      eta <- log(expression[, idx])
      return(logistic(eta, x0 = estimate@dropout.mid, k = estimate@dropout.shape))
    })
    keep <- matrix(rbinom(n_cells * n_genes, 1, 1 - drop.prob), 
                   nrow = n_genes, ncol = n_cells)
    
    counts <- true_counts
    counts[!keep] <- 0
  } else {
    counts <- true_counts
  }
  dimnames(counts) <- dimnames(true_counts)
  
  experiment <- lst(
    cellinfo,
    expression_simulated,
    expression = t(expression),
    true_counts = t(true_counts),
    counts = t(counts),
    geneinfo
  )
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

#' Filter counts
#' 
#' @param experiment Experiment list, containing expression, cellinfo and geneinfo
#' 
#' @export
filter_experiment <- function(experiment) {
  remove_cells <- (apply(experiment$expression, 1, max) == 0) | is.na(apply(experiment$expression, 1, sd))
  
  experiment$expression <- experiment$expression[!remove_cells, ]
  experiment$cellinfo <- experiment$cellinfo %>% slice(match(rownames(experiment$expression), cell_id))
  experiment
}


#' 
#' #' Get housekeeping reference means
#' #' 
#' #' @param counts Expression matrix containing counts
#' #' @export
#' get_housekeeping_reference_means <- function(counts) colMeans(counts)
#' 
#' #' Add housekeeping genes
#' #' 
#' #' @param expression The original expression data.
#' #' @param geneinfo The original gene info
#' #' @param housekeeping_reference_means The mean expression of a set of genes in the reference dataset
#' #' @param n_housekeeping_genes The number of genes to add
#' #' @param overallaverage Overall average expression of the original dataset, this keeps the overall average expression the same even with housekeeping genes
#' #' @param gene_id_generator Function to generate gene_ids
#' #' 
#' #' @export
#' #' @importFrom utils data
#' #' @importFrom magrittr set_colnames
#' add_housekeeping_poisson <- function(
#'   expression, 
#'   geneinfo, 
#'   housekeeping_reference_means, 
#'   n_housekeeping_genes=200, 
#'   overallaverage = mean(expression),
#'   gene_id_generator = function(n) {paste0("GH", seq_len(n))}
#' ) {
#'   if(is.null(housekeeping_reference_means)) stop("Reference means required!!")
#'   
#'   meanpoissons <- overallaverage * housekeeping_reference_means/mean(housekeeping_reference_means)
#'   
#'   additional_expression <- purrr::map(sample(meanpoissons, n_housekeeping_genes), ~rpois(nrow(expression), .)) %>%
#'     invoke(cbind, .) %>% 
#'     magrittr::set_colnames(gene_id_generator(n_housekeeping_genes))
#'   
#'   geneinfo <- dplyr::bind_rows(
#'     geneinfo %>% dplyr::mutate(housekeeping=FALSE), 
#'     tibble(gene=colnames(additional_expression), housekeeping=TRUE)
#'   )
#'   
#'   list(expression=cbind(expression, additional_expression), geneinfo=geneinfo)
#' }