#' Plotting
#' 
#' @param model The model to plot
#' 
#' @importFrom igraph layout.graphopt graph_from_data_frame plot.igraph
plot_modulenet <- function(model) {
  graph <- igraph::graph_from_data_frame(model$modulenet, vertices = model$modulenodes)
  
  modulenames <- unique(model$geneinfo$module)
  colors <- rainbow(length(modulenames))
  V(graph)$color <- colors[match(names(V(graph)), modulenames)]
  
  layout <- igraph::layout.graphopt(graph, charge=0.01, niter=10000)
  #layout <- ForceAtlas2::layout.forceatlas2(graph, iterations=3000, plotstep=100)
  #png(file.path(imagefolder, "net_consecutive_bifurcating.png"), pointsize = 30, width=1000, height=1000)
  #igraph::plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = 20, edge.arrow.size=0.5, edge.loop.angle=0.1, vertex.color=c("#222222", "#662222")[model$modulenodes$a0+1], vertex.label.color="white")
  
  igraph::plot.igraph(
    graph,
    edge.color = c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(-1,1, 0)))], 
    layout = layout,
    vertex.size = 20, 
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1, 
    vertex.color = V(graph)$color, 
    vertex.label.color = "black"
  )
  #dev.off()
}

plot_net <- function(model) {
  graph <- igraph::graph_from_data_frame(model$net)
  layout <- igraph::layout_with_fr(graph)
  ldtf_filter <- as.numeric(factor(model$geneinfo$ldtf, levels = c(F, T)))
  igraph::plot.igraph(
    graph,
    edge.color = c("#d63737", "#3793d6", "#7cd637")[as.numeric(factor(model$net$effect, levels = c(-1,1, 0)))],
    layout = layout,
    vertex.size = c(1,5)[ldtf_filter], 
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1,
    vertex.color = c("black", "white")[ldtf_filter]
  )
}

#' @importFrom pheatmap pheatmap
plot_net_overlaps <- function(model) {
  jaccard <- function(x, y) {length(intersect(x, y))/length(union(x,y))}
  pheatmap::pheatmap(sapply(model$geneinfo$gene, function(i) sapply(model$geneinfo$gene, function(j) jaccard(model$net$from[model$net$to==i], model$net$from[model$net$to==j]))))
}

#' @importFrom igraph graph_from_data_frame V E plot.igraph
plot_net <- function(model) {
  goi <- unique(c(model$net$from, model$net$to))
  subnet <- model$net %>% filter(from %in% goi) %>% filter(to %in% goi)
  graph <- igraph::graph_from_data_frame(subnet, vertices=model$geneinfo$gene)
  
  modulenames <- model$modulenodes$module
  
  colors <- rainbow(length(modulenames))
  V(graph)$color <- colors[match(model$geneinfo$module[match(names(V(graph)), model$geneinfo$gene)], modulenames)]
  #V(graph)$label <- model$geneinfo$gene
  V(graph)$label <- ""
  E(graph)$color <- subnet$effect %>% factor(levels=c(1, -1)) %>% as.numeric() %>% c("blue", "red")[.]
  igraph::plot.igraph(graph, vertex.size=6, edge.arrow.size=0.5)
}

#' @importFrom igraph graph_from_data_frame V E plot.igraph
plot_net_tfs <- function(model) {
  goi <- unique(model$net$from)
  subnet <- model$net %>% filter(from %in% goi) %>% filter(to %in% goi)
  graph <- igraph::graph_from_data_frame(subnet)
  
  modulenames <- model$modulenodes$module
  
  colors <- rainbow(length(modulenames))
  V(graph)$color <- colors[match(model$geneinfo$module[match(names(V(graph)), model$geneinfo$gene)], modulenames)]
  E(graph)$color <- subnet$effect %>% factor(levels=c(1, -1)) %>% as.numeric() %>% c("blue", "red")[.]
  igraph::plot.igraph(graph)
}

