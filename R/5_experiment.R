#' Extract cells from a simulation, either simulating a snapshot experiment or a synchronized experiment
#' @param simulation List containing the expression of simulation(s)
#' @param takesettings List containing the information of how the cells should be taken
#' @export
take_experiment_cells <- function(simulation, model, takesettings = list(type="snapshot", ncells=500)) {
  # sample
  if (takesettings$type == "snapshot"){
    
    sample_ids <- sample(seq_len(nrow(simulation$expression)), takesettings$ncells)
    expression <- simulation$expression[sample_ids, ]
    cellinfo <- simulation$stepinfo[sample_ids, ]
    
  } else if (takesettings$type == "synchronized") {
    
    totaltime <- max(simulation$stepinfo$simulationtime)
    timepoints <- seq(0, totaltime, length.out=takesettings$ntimepoints)
    
    sample_steps <- map(timepoints, function(timepoint) {
      simulation$stepinfo %>% group_by(simulation_id) %>% 
        summarise(step_id = step_id[which.min(abs(simulationtime - timepoint))]) %>% 
        mutate(timepoint=timepoint)
    }) %>% bind_rows()
    
    expression <- simulation$expression[sample_steps$step_id, ]
    cellinfo <- sample_steps %>% left_join(simulation$stepinfo, by="step_id")
    
  }
  # get cellinfo
  rownames(expression) <- paste0("C", seq_len(nrow(expression)))
  cellinfo <- cellinfo %>% mutate(cell_id=rownames(expression))
  
  # get geneinfo
  geneinfo <- model$geneinfo
  
  # filter
  experiment <- filter_expression(expression, cellinfo, geneinfo)
  
  experiment
}

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
#' @param ngenes The number of genes to add
#' @param overallaverage Overall average expression of the original dataset, this keeps the overall average expression the same even with housekeeping genes
#' @param meanpoissons The mean expression of a set of genes in the reference dataset
#' 
#' @export
#' @importFrom utils data
add_housekeeping_poisson <- function(expression, geneinfo, housekeeping_reference_means, n_housekeeping_genes=200, overallaverage = mean(expression)) {
  if(is.null(housekeeping_reference_means)) stop("Reference means required!!")
  
  meanpoissons <- overallaverage * housekeeping_reference_means/mean(housekeeping_reference_means)
  
  additional_expression <- purrr::map(sample(meanpoissons, n_housekeeping_genes), ~rpois(nrow(expression), .)) %>%
    invoke(cbind, .) %>% 
    magrittr::set_colnames(seq_len(n_housekeeping_genes)+ncol(expression))
  
  geneinfo <- dplyr::bind_rows(
    geneinfo %>% dplyr::mutate(housekeeping=F), 
    tibble(gene=colnames(additional_expression) %>% as.numeric(), housekeeping=T)
  )
  
  list(expression=cbind(expression, additional_expression), geneinfo=geneinfo)
}
