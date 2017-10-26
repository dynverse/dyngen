#' Plot a modulenet
#' 
#' @param model The model
#' 
#' @importFrom igraph layout.graphopt graph_from_data_frame plot.igraph
#' @export
plot_modulenet <- function(model) {
  graph <- igraph::graph_from_data_frame(model$modulenet, vertices = model$modulenodes)
  
  modulenames <- unique(model$geneinfo$module_id)
  colors <- rainbow(length(modulenames))
  igraph::V(graph)$color <- colors
  igraph::E(graph)$color <- c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(1,-1, 0)))]
  
  layout <- igraph::layout.graphopt(graph, charge=0.01, niter=10000)
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    vertex.size = 20, 
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1, 
    vertex.label.color = "black"
  )
}


#' Plot a gene network
#' @param model The model
#' @param colorby By what to color
#' @param main_only Whether to only draw the main network
#' @param label Whether to label genes
#' @export
plot_net <- function(model, colorby=c("module", "main"), main_only=TRUE, label=FALSE) {
  colorby <- match.arg(colorby)
  
  geneinfo  <- model$geneinfo
  if(main_only) geneinfo <- geneinfo %>% filter(main)
  net <- model$net %>% filter((from %in% geneinfo$gene_id) & (to %in% geneinfo$gene_id))
  
  graph <- igraph::graph_from_data_frame(net %>% select(from, to), vertices=geneinfo$gene_id)
  
  # layout
  set.seed(1) # to get same layout
  layout <- igraph::layout.fruchterman.reingold(graph)
  main_filter <- as.numeric(factor(geneinfo$main, levels = c(FALSE, TRUE, NA), exclude=NULL))
  
  # change vertex/edge colors and sizes
  igraph::V(graph)$size <- c(1, 4, 1)[main_filter]
  
  if (!label) {
    igraph::V(graph)$label <- rep("", length(igraph::V(graph)))
  }
  
  if (colorby == "main") {
    igraph::V(graph)$color <- c("black", "white", "grey")[main_filter]
  } else if (colorby == "module") {
    modulenames <- model$modulenodes$module_id
    colors <- rainbow(length(modulenames)) %>% set_names(modulenames)
    igraph::V(graph)$color <- colors[geneinfo$module_id]
  }
  igraph::E(graph)$color <- c("#d63737", "#3793d6", "#7cd637")[as.numeric(factor(net$effect, levels = c(1,-1, 0)))]
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1
  )
}

#' Plot the overlap in targets in a heatmap
#' @param model The model
#' @importFrom pheatmap pheatmap
plot_net_overlaps <- function(model) {
  jaccard <- function(x, y) {length(intersect(x, y))/length(union(x,y))}
  pheatmap::pheatmap(sapply(model$geneinfo$gene_id, function(i) sapply(model$geneinfo$gene_id, function(j) jaccard(model$net$from[model$net$to==i], model$net$from[model$net$to==j]))))
}





#' Plot the simulations
#' @importFrom grDevices rainbow
#' @importFrom stats sd
plot_simulations = function(simulations, samplingrate=0.1) {
  requireNamespace("rgl")
  for (i in seq_len(length(simulations))) {
    sample = sample(c(T, F), size=nrow(simulations[[i]]$expression), T, c(samplingrate, 1-samplingrate))
    
    sample[simulations[[i]]$expression[sample,] %>% apply(1, mean) %>% {. == 0}] = FALSE
    
    simulations[[i]]$subcellinfo = simulations[[i]]$cellinfo[sample,]
    simulations[[i]]$subexpression = simulations[[i]]$expression[sample,]
    
    simulations[[i]]$expression[sample,] %>% apply(1, stats::sd) %>% {. == 0} %>% sum %>% print
    simulations[[i]]$subexpression %>% {apply(., 1, stats::sd) == 0} %>% sum %>% print
  }
  
  overallexpression = map(simulations, "subexpression") %>% do.call(rbind, .)
  overallexpression = overallexpression %>% set_rownames(1:nrow(overallexpression))
  overallcellinfo = assign_progression(overallexpression, reference)
  overallcellinfo$simulation_id = map(seq_len(length(simulations)), ~rep(., nrow(simulations[[.]]$expression))) %>% unlist %>% factor()
  overallcellinfo$observed = T
  
  space =  lmds(overallexpression, 3) %>% as.data.frame() %>% bind_cols(overallcellinfo)
  space %>% filter(observed) %>% ggplot() + geom_path(aes(Comp1, Comp2, group=simulation_id, color=progression)) + geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3)
  space %>% filter(observed) %>% ggplot() + geom_path(aes(Comp1, Comp2, group=simulation_id, color=progression)) + geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3) + viridis::scale_color_viridis(option="A")
  space %>% filter(observed) %>% ggplot() + geom_path(aes(V1, V2, group=simulation_id, color=state_id)) + geom_path(aes(V1, V2), data=space %>% filter(!observed), size=3)
  
  for (i in unique(as.numeric(space$simulation_id))) 
    rgl::lines3d(space %>% filter(simulation_id == i), col = grDevices::rainbow(length(unique(space$simulation_id)))[[i]])
  
}