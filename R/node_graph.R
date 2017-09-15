#' @importFrom stats setNames
make_node_graph <- function(linegraph_dataframe, prefix="S") {
  # maak een lijst van alle edges
  edges <- c(linegraph_dataframe$from, linegraph_dataframe$to) %>% unique %>% sort
  
  # plak 2 nodes aan elke edge; voor een node i heten deze i_start en i_end
  edges_str <- sort(c(paste0(edges, "_start"), paste0(edges, "_end")))
  
  # laad de linegraph in in igraph, en bekijk welke nodes equivalent zijn aan elkaar
  gr <- linegraph_dataframe %>% 
    transmute(from_str = paste0(from, "_end"), to_str = paste0(to, "_start")) %>% 
    igraph::graph_from_data_frame(vertices = edges_str)
  
  # if you want: plot(gr)
  edge_to_node <- igraph::components(gr)$membership
  
  # geef elke node een unieke naam
  edge_to_node <- stats::setNames(paste0(prefix, edge_to_node), names(edge_to_node))
  
  # vertaal elke edge naar zijn bijhorende from en to node
  edge_df <- bind_rows(lapply(edges, function(i) {
    data_frame(from = edge_to_node[paste0(i, "_start")], to = edge_to_node[paste0(i, "_end")], edge_id = i)
  }))
  
  edge_df
}