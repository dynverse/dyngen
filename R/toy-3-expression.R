#' Generate toy expression based on a milestone_network and cell progression information
#' @export
generate_expression <- function(milestone_network, progressions, ngenes=100, noise_std=0.05) {
  nedges <- nrow(milestone_network)
  milestone_expressions <- list()
  milestone_network <- milestone_network %>% mutate(splinefuns=map(seq_len(n()), ~NULL))
  
  for (edge_id in seq_len(nedges)) {
    edge <- dyneval:::extract_row_to_list(milestone_network, edge_id)
    
    # check whether the starting and ending milestones have already been visited, otherwise the start and end are random
    start <- if (edge$from %in% names(milestone_expressions)) milestone_expressions[[edge$from]] else runif(ngenes)
    end <- if (edge$to %in% names(milestone_expressions)) milestone_expressions[[edge$to]] else runif(ngenes)
    
    milestone_expressions[[edge$from]] <- start
    milestone_expressions[[edge$to]] <- end
    
    xs <- map(seq_len(ngenes), ~c(0, sort(runif(3)), 1))
    ys <- pmap(list(x=xs, start=start, end=end), function(x, start, end) c(start, start, runif(length(x) - 4), end, end))
    
    milestone_network$splinefuns[edge_id] <- map2(xs, ys, function(x, y) {
      approxfun(x, y)
    }) %>% list()
  }
  
  filtered_progression <- progressions %>% # a cell can only be in one edge (maximum in tents)
    group_by(cell_id) %>%
    arrange(-percentage) %>%
    #mutate(percentage = sum(percentage)) %>%
    filter(row_number() == 1)
  
  # extract expression for each edge
  expression <- filtered_progression %>%
    group_by(from, to) %>%
    summarise(percentages = list(percentage), cell_ids=list(cell_id)) %>%
    left_join(milestone_network, by=c("from", "to")) %>%
    rowwise() %>%
    do(
      expression=map(.$splinefuns, function(f) f(.$percentage)) %>% invoke(rbind, .),
      cell_id=.$cell_id
    ) %>% {
      magrittr::set_colnames(invoke(cbind, .$expression), unlist(.$cell_id))
    } %>% t
  
  expression <- expression + rnorm(length(expression), 0, noise_std)
  
  expression <- expression[unique(progressions$cell_id), ]
  
  expression
}


# library(tidyverse)
# 
# milestone_network = dyngen::generate_toy_milestone_network("bifurcating")
# progressions = dyngen::random_progressions_tented(milestone_network, 500)
# expression <- generate_expression(milestone_network, progressions, 100)
# 
# cell_order <- order_cells(milestone_network, progressions)
# expression[cell_order$order, ] %>% pheatmap::pheatmap(cluster_rows=F, gaps_row = cell_order$edge_id %>% diff() %>% {which(. != 0)})
# 
# source("../dyneval/R_modular_TI/dimred_wrappers.R")
# space <- dimred_ica(expression, ndim=2)
# ggplot(space %>% as.data.frame()) + geom_point(aes(Comp1, Comp2))
# 
# library(rgl)
# plot3d(space[, 1], space[, 2], space[, 3])
