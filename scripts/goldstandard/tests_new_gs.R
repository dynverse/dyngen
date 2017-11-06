preprocess_simulation_for_gs <- function(simulation, model) {
  simulation$expression_normalized <- t(dynutils::scale_quantile(t(simulation$expression), outlier_cutoff=0.05))
  simulation$expression_modules <- simulation$expression_normalized %>% t %>% as.data.frame() %>% split(model$geneinfo$module_id) %>% map(~apply(., 2, mean)) %>% do.call(rbind, .) %>% t
  
  simulation
}


simulation <- preprocess_simulation_for_gs(simulation, model)


expression <- simulation$expression_modules
simulation_ids <- unique(simulation$stepinfo$simulation_id)

simulation_id_oi <- 5

verbose <- FALSE

# calculate for every step the minimal distance to other simulations
distances <- map(simulation_ids, function(simulation_id_oi) {
  stepinfo_oi <- simulation$stepinfo %>% filter(simulation_id == simulation_id_oi)
  expression_oi <- expression[stepinfo_oi$step_id, ]
  
  distances <- map(setdiff(simulation_ids, simulation_id_oi), function(simulation_id_ref) {
    stepinfo_or <- simulation$stepinfo %>% filter(simulation_id == simulation_id_ref)
    expression_or <- expression[stepinfo_or$step_id, ]
    
    FNN::knnx.dist(expression_or, expression_oi, 1)[, 1] %>% tibble(distance = ., cell_id = rownames(expression_oi), simulation_id_ref = simulation_id_ref, step_id=stepinfo_oi$step_id, simulation_id = simulation_id_oi)
  }) %>% bind_rows()
  
  if (verbose) {
    distances_matrix <- distances %>% mutate(step_id = factor(step_id, levels=unique(step_id))) %>% reshape2::acast(step_id~simulation_id_ref, value.var="distance")
    distances_matrix %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  }
  
  distances
}) %>% bind_rows()

# calculate whether two simulations are connected
distances <- distances %>% 
  group_by(simulation_id, simulation_id_ref) %>% 
  mutate(connected = zoo::rollapply(distance, 100, function(x) all(x < 0.4), align="left", partial=FALSE, fill=NA))

connected <- distances %>% reshape2::acast(cell_id~simulation_id_ref, value.var="connected", mean, fill=1)
connected[is.na(connected)] <- 1

connected2 <- connected[sample(nrow(connected), 10000), ]

connected2 %>% pheatmap::pheatmap()

cell_distances <- SCORPIUS::euclidean_distance(connected2)
clust <- fastcluster::hclust(as.dist(cell_distances))
plot(clust)
labels <- cutree(clust, h=0)

label_counts <- labels %>% table %>% sort(TRUE)
clusters <- map2(names(label_counts), label_counts, function(cluster_id, count) {
  tibble(connected = list(connected2[which(labels == cluster_id)[[1]], ]), ncells = count, cluster_id=as.numeric(cluster_id))
}) %>% bind_rows()

clusters$selected <- FALSE
for (i in seq_len(nrow(clusters))) {
  cluster <- clusters %>% dynutils::extract_row_to_list(i)
  
  connected_cluster <- cluster$connected
  connected_similarities <- map_dbl(clusters %>% filter(selected) %>% pull(connected), function(connected_ref_cluster) {
    (connected_cluster - connected_ref_cluster) %>% abs %>% sum
  })
  
  if(length(connected_similarities) == 0 || min(connected_similarities) > 2) {
    if(cluster$ncells > 100) {
      clusters[i, "selected"] <- TRUE
    }
  }
}
clusters

selected_cluster_ids <- clusters %>% filter(selected) %>% pull(cluster_id)
labels_selected <- ifelse(labels %in% selected_cluster_ids, labels, NA)
stepinfo <- simulation$stepinfo %>% slice(match(names(labels), step_id)) %>% mutate(cluster_id=labels_selected) %>% arrange(simulation_id, step)


cluster_links <- stepinfo %>% group_by(simulation_id) %>% drop_na() %>% summarise(cluster_ids=list(cluster_id)) %>% mutate(cluster_ids = map(cluster_ids, function(cluster_ids) cluster_ids[c(1, which(diff(cluster_ids) != 0)+1)])) %>% pull(cluster_ids)

cluster_network <- map(cluster_links, function(cluster_links) {
  if (length(cluster_links) > 1) {
    tibble(from=cluster_links[1:length(cluster_links)-1], to=cluster_links[2:length(cluster_links)])
  }
}) %>% bind_rows() %>% group_by(from, to) %>% summarise(n=n())

cluster_network %>% igraph::graph_from_data_frame() %>% plot













simulation$stepinfo %>% slice(match(names(labels), step_id)) %>% mutate(cluster_id=labels_selectedsimulation$stepinfo %>% slice(match(names(labels), step_id)) %>% mutate(cluster_id=labels) %>% arrange(simulation_id, step) %>% Viewexpression <- simulation$expression
simulation_ids <- unique(simulation$stepinfo$simulation_id)

simulation_id_oi <- 5

verbose <- TRUE

step_infos <- map(simulation_ids, function(simulation_id_oi) {
  stepinfo_oi <- simulation$stepinfo %>% filter(simulation_id == simulation_id_oi)
  expression_oi <- expression[stepinfo_oi$step_id, ]
  
  distances <- map(setdiff(simulation_ids, simulation_id_oi), function(simulation_id_ref) {
    stepinfo_or <- simulation$stepinfo %>% filter(simulation_id == simulation_id_ref)
    expression_or <- expression[stepinfo_or$step_id, ]
    
    FNN::knnx.dist(expression_or, expression_oi, 1)[, 1] %>% tibble(distance = ., cell_id = rownames(expression_oi), simulation_id_ref = simulation_id_ref, step_id=stepinfo_oi$step_id)
  }) %>% bind_rows()
  window_smoothing <- 50
  distances_smoothened <- distances %>% mutate(step_id = factor(step_id, levels=unique(step_id))) %>% reshape2::acast(step_id~simulation_id_ref, value.var="distance") %>% zoo::rollapply(window_smoothing, mean, fill=NA, align="left", partial=TRUE)
  
  if (verbose) distances_smoothened %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  
  window_connections <- 500
  connections <- distances_smoothened %>% zoo::rollapply(window_connections, function(x) all(x > 0.5), fill=NA, align="left", partial=TRUE) %>% apply(1, function(x) c(1, 2)[as.numeric(x) + 1])
  if (verbose) connections %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F, breaks=c(0, 1))
  
  window_change <- 300
  stable_regions <- connections %>% t() %>% zoo::rollapply(window_change, function(x) all(x == x[[1]]), fill=NA, align="left") %>% t() %>% apply(2, all)
  
  stepinfo_oi %>% mutate(
    stable_region = stable_regions,
    piece_id = c(1, cumsum(diff(stable_region) == 1) + 1)
  )
}) %>% bind_rows()

step_infos %>% bind_rows() %>% ggplot() + geom_line(aes(simulationtime, piece_id)) + facet_wrap(~simulation_id)

step_infos %>% filter(stable_region) %>% {split(., .$piece_id)}




library(GNG)
gng <- GNG::gng(simulation$expression, verbose = TRUE, max_nodes=100)
gng$edges %>% igraph::graph_from_data_frame() %>% plot


gng$edges %>% igraph::graph_from_data_frame()

# simulation$molecules[simulation$stepinfo$simulation_id == 1, ] %>% diff(100) %>% zoo::rollmean(100) %>% apply(1, mean) %>% zoo::rollapply(100, sd) %>% plot
# 
# 
# 
# convergence <- simulation$molecules[simulation$stepinfo$simulation_id == 1, ] %>% zoo::rollapply(200, sd, 1) %>% apply(1, max)
# 
# convergence %>% plot()
# 
# convergence
# 
# 
# simulation$molecules[simulation$stepinfo$simulation_id == 1, ] %>% pheatmap::pheatmap(cluster_rows=F, cluster_cols=F)
