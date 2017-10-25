#' Extract cells from a simulation, either simulating a snapshot experiment or a synchronized experiment
#' @param simulation List containing the expression of simulation(s)
#' @param takesettings List containing the information of how the cells should be taken
#' @export
take_experiment_cells <- function(simulation, takesettings = list(type="snapshot", ncells=500)) {
  if(takesettings$type == "snapshot"){
    sample_ids <- sample(seq_len(nrow(simulation$expression)), takesettings$ncells)
    experiment <- list(
      expression=simulation$expression[sample_ids, ],
      cellinfo=simulation$stepinfo[sample_ids, ]
    )
  } else if(takesettings$type == "synchronized") {
    totaltime <- max(simulation$stepinfo$simulationtime)
    timepoints <- seq(0, totaltime, length.out=takesettings$ntimepoints)
    
    sample_steps <- map(timepoints, function(timepoint) {
      simulation$stepinfo %>% group_by(simulation_id) %>% 
        summarise(step_id = step_id[which.min(abs(simulationtime - timepoint))]) %>% 
        mutate(timepoint=timepoint)
    }) %>% bind_rows()
    
    experiment <- list(
      expression=simulation$expression[sample_steps$step_id, ],
      cellinfo=sample_steps %>% left_join(simulation$stepinfo, by="step_id")
    )
  }
  rownames(experiment$expression) <- paste0("C", seq_len(nrow(experiment$expression)))
  cellinfo <- experiment$cellinfo %>% mutate(cell_id=rownames(experiment$expression))
  
  tibble::lst(cellinfo, expression)
}

#' Add housekeeping genes
#' 
#' @param expression The original expression data.
#' @param geneinfo The original gene info
#' @param ngenes The number of genes to add
#' @param overallaverage TODO: Zouter/wouters help?
#' 
#' @export
#' @importFrom utils data
add_housekeeping_poisson <- function(expression, geneinfo, ngenes=200, overallaverage = mean(expression)) {
  utils::data(ginhoux, envir = environment())
  reference_expression <- (2^ginhoux$expression)-1
  meanpoissons <- colMeans(reference_expression) %>% {./mean(reference_expression)*overallaverage}
  
  additional_expression <- purrr::map(sample(meanpoissons, ngenes), ~rpois(nrow(expression), .)) %>% 
    invoke(cbind, .) %>% 
    magrittr::set_colnames(seq_len(ngenes)+ncol(expression))
  
  geneinfo <- dplyr::bind_rows(geneinfo %>% dplyr::mutate(housekeeping=F), tibble(gene=colnames(additional_expression) %>% as.numeric(), housekeeping=T))
  
  list(expression=cbind(expression, additional_expression), geneinfo=geneinfo)
}
