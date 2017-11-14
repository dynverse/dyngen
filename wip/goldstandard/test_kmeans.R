smooth_expression <- function(expression) {
  expression %>% 
    zoo::rollmean(50, c("extend", "extend", "extend")) %>%
    magrittr::set_rownames(rownames(expression))
}

preprocess_simulation_for_gs <- function(simulation, model) {
  print("smoothing...")
  
  simulation$expression_smooth <- simulation$expression %>% as.data.frame() %>% split(simulation$stepinfo$simulation_id) %>% map(smooth_expression) %>% do.call(rbind, .)
  dimnames(simulation$expression_smooth) <- dimnames(simulation$expression)
  
  print("normalizing...")
  simulation$expression_normalized <- t(dynutils::scale_quantile(t(simulation$expression_smooth), outlier_cutoff=0.05))
  
  print("calculating module expression...")
  simulation$expression_modules <- simulation$expression_normalized %>% t %>% as.data.frame() %>% split(model$geneinfo$module_id) %>% map(~apply(., 2, mean)) %>% do.call(rbind, .) %>% t
  
  simulation
}

# preprocessing
simulation <- preprocess_simulation_for_gs(simulation, model)

expression <- simulation$expression_modules
stepinfo <- simulation$stepinfo

# subsample
# samplexpression <- expression[sample(nrow(expression), min(4000, nrow(expression))), ]
samplexpression <- expression[((stepinfo$step)%%2) == 1, ]
samplestepinfo <- stepinfo %>% slice(match(rownames(samplexpression), step_id))

# kmeans
centers <- samplexpression %>% kmeans(200, iter.max = 100)
samplestepinfo$center_id <- centers$cluster

samplestepinfo %>% ggplot(aes(step, center_id, color=simulation_id, group=simulation_id)) + geom_line()

center_distances <- centers$centers %>% dist()

# # do not allow flipping (due to noise)
# samplestepinfo$run <- rle(samplestepinfo$center_id) %>% {rep(seq_along(.$lengths), .$lengths)}
# samplestepinfo <- samplestepinfo %>% group_by(simulation_id, run) %>% 
#   nest() %>% 
#   group_by(simulation_id) %>% 
#   mutate(
#     center_id = map_dbl(data, ~.$center_id[[1]]),
#     next_center_id = lead(center_id, 1),
#     prev_center_id = lag(center_id, 1)
#   ) %>% 
#   mutate(center_id = ifelse(next_center_id == prev_center_id & !is.na(next_center_id == prev_center_id), next_center_id, center_id)) %>% 
#   select(-next_center_id, -prev_center_id) %>% 
#   mutate(data = map(data, ~select(., -center_id))) %>% 
#   unnest(data) %>% 
#   ungroup()
# 
# samplestepinfo %>% ggplot(aes(step, center_id, color=factor(simulation_id), group=simulation_id)) + geom_line() + facet_wrap(~simulation_id)

# extract center network from stepinfo
extract_center_network <- function(stepinfo, center_column_name = "center_id") {
  stepinfo$center_id <- stepinfo[[center_column_name]]
  stepinfo$run <- rle(stepinfo$center_id) %>% {rep(seq_along(.$lengths), .$lengths)}
  stepinfo_runs <- stepinfo %>% group_by(simulation_id, run) %>% 
    summarise(
      center_id = center_id[[1]],
      length = n()
    )
  
  center_network <- stepinfo_runs %>% group_by(simulation_id) %>% 
    mutate(from = center_id, to = lead(center_id, 1)) %>% 
    drop_na() %>% 
    select(from, to, simulation_id) %>% 
    ungroup() %>% 
    group_by(from, to) %>% 
    summarise(simulation_ids = list(unique(simulation_id)), n_simulations = length(unique(simulation_id)), n_runs = n()) %>% 
    ungroup()
  
  center_network
}

center_network <- extract_center_network(samplestepinfo)

if(!igraph::is.connected(center_network %>% igraph::graph_from_data_frame())) { stop("not connected!")}

center_network$width <- center_network$n_simulations / max(center_network$n_simulations) * 10
center_graph <- center_network %>% igraph::graph_from_data_frame()
center_graph %>% plot


cluster_ids <- unique(c(center_network$from, center_network$to))

samplestepinfo

##

calculate_passthrough_similarity <- function(graph) {
  adj <- as_adjacency_matrix(graph) %>% as.matrix()
  adj <- adj / apply(adj, 1, sum)
  newadj <- adj
  adjs <- list()
  for (i in 1:20) {
    newadj <- newadj / apply(newadj, 1, sum)
    newadj <- newadj %*% adj
    adjs <- c(adjs, list(newadj))
  }
  adjs <- map2(seq_along(adjs), adjs, function(step, adj) {adj %>% reshape2::melt(varnames=c("from", "to"), value.name="prob") %>% mutate(step=step)}) %>% bind_rows()
  probs <- cbind(adjs$prob %>% split(adjs$from) %>% do.call(rbind, .), adjs$prob %>% split(adjs$to) %>% do.call(rbind, .))
  #cor(probs) %>% pheatmap::pheatmap()
  cor(t(probs)) %>% magrittr::set_colnames(colnames(adj)) %>% magrittr::set_rownames(rownames(adj))
}

samplestepinfo$new_center_id <- samplestepinfo$center_id

n <- 20
settings <- tibble(center_cutoff = seq(0.1, 0.5, length.out = n), passthrough_similarity_cutoff = 1)

center_network <- extract_center_network(samplestepinfo, "new_center_id")
center_graph <- center_network %>% graph_from_data_frame()

center_graph %>% plot()

setting_id <- 1

for (setting_id in seq_len(nrow(settings))) {
  dynutils::extract_row_to_list(settings, setting_id) %>% list2env(.GlobalEnv)
  
  # calculate time distances
  cluster_ids <- unique(c(center_network$from, center_network$to))
  
  # calculate vertex similarity
  passthrough_similarity <- calculate_passthrough_similarity(center_graph)[as.character(cluster_ids), as.character(cluster_ids)]
  
  # pheatmap::pheatmap(diffsinks)
  
  time_distances <- center_network %>%
    mutate(from=factor(from, levels=cluster_ids), to=factor(to, levels=cluster_ids)) %>%
    reshape2::acast(from~to, fill=0, value.var="n_runs", drop=F)
  
  # calculate expression distances
  center_distances <- limma::avearrays(t(samplexpression), samplestepinfo$new_center_id) %>% t() %>% dist() %>% as.matrix() %>% .[as.character(cluster_ids), as.character(cluster_ids)]
  
  combined <- (center_distances < center_cutoff) & (passthrough_similarity > passthrough_similarity_cutoff)
  # combined %>% graph_from_adjacency_matrix() %>% plot()
  components <- combined %>% igraph::graph_from_adjacency_matrix() %>% igraph::components() %>% .$membership
  
  center_graph <- igraph::contract(center_graph, components, toString) %>% igraph::simplify()
  
  center_id_map <- names(V(center_graph)) %>% tibble(old_center_id = ., new_center_id = seq_along(.)) %>% mutate(old_center_id=strsplit(old_center_id, ", ")) %>% unnest(old_center_id) %>% mutate_all(funs(as.integer))
  
  samplestepinfo <- samplestepinfo %>% select(-matches("old_center_id")) %>% rename(old_center_id = new_center_id) %>% left_join(center_id_map, by="old_center_id")
  
  center_network <- extract_center_network(samplestepinfo, "new_center_id")
  center_graph <- center_network %>% graph_from_data_frame()

  center_graph %>% plot()
}







center_distances <- limma::avearrays(t(samplexpression), samplestepinfo$new_center_id) %>% t() %>% dist() %>% as.matrix()


samplestepinfo$run <- rle(samplestepinfo$new_center_id) %>% {rep(seq_along(.$lengths), .$lengths)}
samplestepinfo %>% group_by(samplestepinfo)







score_simplicity <- function(graph) graph %>% igraph::as.undirected() %>% degree() %>% mean()
score_simplicity <- function(graph) sum(igraph::degree(graph %>% igraph::as.undirected()) > 2)
#score_simplicity <- function(graph) length(get.elementary.circuits(graph %>% igraph::as.undirected() %>% igraph::simplify()) %>% keep(~length(.) > 3))

changed <- TRUE
distance_cutoff <- 2
# simplicity_cutoff <- 5
simplicity_cutoff <- 2.5
simplicity_cutoff <- 4
distance_cutoff_step <- 0.2

while(score_simplicity(center_graph) > simplicity_cutoff) {
  distances <- center_graph %>% igraph::distances()
  distance_distances <- dist(distances) %>% as.matrix()
  
  components <- (distance_distances <= distance_cutoff) %>% igraph::graph_from_adjacency_matrix() %>% igraph::components() %>% .$membership
  
  center_graph <- igraph::contract(center_graph, components, vertex.attr.comb=toString) %>% igraph::simplify()
  center_graph <- igraph::set.vertex.attribute(center_graph, "name", value=seq_along(V(center_graph))) # rename them, necessary for elementary circuits, center_ids are kept in the center_ids attribute
  
  distance_cutoff <- distance_cutoff + distance_cutoff_step
  
  center_graph %>% plot()
  
  print(score_simplicity(center_graph))
}
center_graph %>% plot()

# map old center ids to new ones
center_id_map <- V(center_graph)$center_ids %>% tibble(center_id = ., new_center_id = seq_along(.)) %>% mutate(center_id=strsplit(center_id, ", ")) %>% unnest(center_id) %>% mutate_all(funs(as.integer))

# recalculate edge direction
center_network <- center_graph %>% igraph::as.undirected() %>% igraph::as_data_frame()
samplestepinfo <- samplestepinfo %>% left_join(center_id_map, by="center_id")

extract_center_network <- function(stepinfo, center_column_name = "center_id") {
  samplestepinfo$center_id <- stepinfo[[center_column_name]]
  samplestepinfo$run <- rle(samplestepinfo$center_id) %>% {rep(seq_along(.$lengths), .$lengths)}
  samplestepinfo_runs <- samplestepinfo %>% group_by(simulation_id, run) %>% 
    summarise(
      center_id = center_id[[1]],
      length = n()
    )
  
  center_network <- samplestepinfo_runs %>% group_by(simulation_id) %>% 
    mutate(from = center_id, to = lead(center_id, 1)) %>% 
    drop_na() %>% 
    select(from, to, simulation_id) %>% 
    ungroup() %>% 
    group_by(from, to) %>% 
    summarise(simulation_ids = list(unique(simulation_id)), n_simulations = length(unique(simulation_id)), n_runs = n()) %>% 
    ungroup()
  
  center_network
}

center_network <- extract_center_network(samplestepinfo, "new_center_id")

# look at each edge, in every direction, and decide which direction to choose or whether to collapse
fraction_cutoff <- 5
center_network <- center_network %>% 
  mutate(edge_id = map2_chr(from, to, ~paste0(sort(c(.x, .y)), collapse="_"))) %>%  # unique, sorted, edge_id
  arrange(edge_id) %>% 
  group_by(edge_id) %>% 
  mutate(
    n = n(),
    fraction = max(n_runs)/min(n_runs)
  ) %>% 
  arrange(n, fraction) %>% 
  filter(
    (row_number() == which.max(n_runs)) |  # either maximal
    (n == 2 & fraction < (row_number() == which.max(n_runs))) # or indicisive
  ) %>% 
  mutate(n = n()) %>% 
  ungroup()

# collapse
components <- center_network %>% 
  filter(n == 2) %>% 
  select(from, to) %>% 
  igraph::graph_from_data_frame(vertices=unique(c(center_network$from, center_network$to))) %>% 
  igraph::components() %>% 
  .$membership

center_graph <- center_network %>% igraph::graph_from_data_frame()
V(center_graph)$center_ids <- names(V(center_graph))
center_graph <- igraph::contract(center_graph, components, vertex.attr.comb=toString) %>% igraph::simplify()
center_id_map <- V(center_graph)$center_ids %>% tibble(center_id = ., new_center_id = seq_along(.)) %>% mutate(center_id=strsplit(center_id, ", ")) %>% unnest(center_id) %>% mutate_all(funs(as.integer))


center_graph %>% plot()




















samplestepinfo$state_id <- factor(center_grouping[samplestepinfo$center_id])
samplestepinfo$center_id <- factor(samplestepinfo$center_id)
source("../dynmodular/dimred_wrappers.R")
spaces <- map(c(dimred_pca, dimred_ica), ~samplexpression %>% {. + runif(length(.), 0, 0.01)} %>% .(ndim=2) %>% as.data.frame() %>% bind_cols(samplestepinfo))
map(spaces, function(space) {
  map(c("center_id", "state_id"), function(colorize) {
    space$color <- space[[colorize]]
    ggplot(space, aes(Comp1, Comp2)) + geom_point(aes_string(color=colorize))
  })
}) %>% unlist(recursive = F) %>% cowplot::plot_grid(plotlist=.)
