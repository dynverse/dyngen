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
samplexpression <- expression[((stepinfo$step)%%5) == 1, ]
samplestepinfo <- stepinfo %>% slice(match(rownames(samplexpression), step_id))
samplexpression_split <- samplexpression %>% as.data.frame() %>% split(samplestepinfo$simulation_id) %>% map(as.matrix)

# cluster simulations, extract prototypes
library(dtwclust)
clust <- dtwclust::tsclust(samplexpression_split, type="hierarchical")
plot(clust)

labels <- cutree(clust, h=500)
prototypes <- map(samplexpression_split %>% split(labels), function(samplexpression_split_clust) {
  prototype <- DBA(samplexpression_split_clust, centroid = samplexpression_split_clust[[samplexpression_split_clust %>% map_int(nrow) %>% which.max()]])
  prototype
})

walk(prototypes, ~pheatmap::pheatmap(., cluster_cols=F, cluster_rows=F))
prototypes <- map(prototypes, smooth_expression)
walk(prototypes, ~pheatmap::pheatmap(., cluster_cols=F, cluster_rows=F))

expression_prototypes <- do.call(rbind, prototypes)
stepinfo_prototypes <- map2(
  prototypes, 
  seq_along(prototypes), 
  function(prototype, prototype_id) 
    tibble(prototype_id, step_id = rownames(prototype), step = seq_len(nrow(prototype)))
) %>% bind_rows() %>% mutate(prototype_id = factor(prototype_id))

# calculate membership
module_bytes <- 2**(1:ncol(expression_prototypes))
membership <- (expression_prototypes > 0.4) %>% as.numeric() %>% matrix(nrow=nrow(expression_prototypes))
stepinfo_prototypes$membership <- (expression_prototypes > 0.5) %>% apply(1, function(x) sum(module_bytes[as.logical(x)]))

stepinfo_prototypes %>% ggplot(aes(step, membership, color=prototype_id, group=prototype_id)) + geom_line()

# smooth membership
# window <- 100
# for (i in 1:5) {
#   stepinfo_prototypes$membership <- stepinfo_prototypes %>% group_by(prototype_id) %>%
#     mutate(membership = zoo::rollapply(membership, window, function(x) {table(x, exclude=NULL) %>% {names(.)[which.max(.)]}}, partial=TRUE)) %>%
#     pull(membership) %>%
#     as.integer()
# }

stepinfo_prototypes %>% ggplot(aes(step, membership, color=prototype_id, group=prototype_id)) + geom_line()

# calculate binary network
stepinfo_prototypes$run <- rle(stepinfo_prototypes$membership) %>% {rep(seq_along(.$lengths), .$lengths)}
stepinfo_prototypes_runs <- stepinfo_prototypes %>% group_by(prototype_id, run) %>% 
  summarise(
    membership = membership[[1]],
    length = n()
  )

# construct the membership network,...
# by first grouping on prototype and creating a network
membership_network <- stepinfo_prototypes_runs %>% group_by(prototype_id) %>% 
  mutate(from = membership, to = lead(membership, 1)) %>% 
  drop_na() %>% 
  select(from, to, prototype_id) %>% 
  ungroup() %>% 
  group_by(from, to) %>% 
  summarise(prototypes = list(prototype_id)) %>% 
  ungroup()

membership_network %>% igraph::graph_from_data_frame() %>% plot












# find closeness of prototypes
min_link_distance <- 0.5

stepinfo_prototypes <- tibble(step_id = map(prototypes, rownames) %>% unlist(), prototype_id = map2(seq_along(prototypes), prototypes, ~rep(.x, nrow(.y))) %>% unlist(), step = map(prototypes, ~seq_len(nrow(.))) %>% unlist())

prototype_distances <- prototypes %>% do.call(rbind, .) %>% dist %>% as.matrix() %>% reshape2::melt(varnames=c("step_id_from", "step_id_to"), value.name="distance")
prototype_distances <- prototype_distances %>% left_join(stepinfo_prototypes %>% setNames(paste0(names(.), "_from"))) %>% left_join(stepinfo_prototypes %>% setNames(paste0(names(.), "_to")))

prototype_closeness <- prototype_distances %>% group_by(step_id_from, prototype_id_to) %>% summarise(min_distance = min(distance)) %>% ungroup()
distance_connected_cutoff <- 0.5 # depends on noise levels
prototype_closeness$connected <- prototype_closeness$min_distance < distance_connected_cutoff

connected <- prototype_closeness %>% reshape2::acast(step_id_from~prototype_id_to, value.var="connected", mean, fill=1)
connected[is.na(connected)] <- 1 # same simulation -> always connected duhh

connected <- connected[stepinfo_prototypes$step_id, unique(stepinfo_prototypes$prototype_id)]
binary <- 2**(0:(ncol(connected)-1))
stepinfo_prototypes <- stepinfo_prototypes %>% mutate(cluster_id = apply(connected,1, function(x) sum(binary[as.logical(x)]))) 

stepinfo_prototypes$cluster_id %>% plot

# smooth cluster assignment using majority voting
window <- 100
stepinfo_prototypes$cluster_id <- stepinfo_prototypes %>% group_by(prototype_id) %>%
  mutate(cluster_id = zoo::rollapply(cluster_id, window, function(x) {table(x, exclude=NULL) %>% {names(.)[which.max(.)]}}, partial=TRUE)) %>%
  pull(cluster_id)
stepinfo_prototypes$cluster_id <- as.numeric(stepinfo_prototypes$cluster_id)

stepinfo_prototypes %>% ggplot() + geom_line(aes(step, as.numeric(cluster_id), color=factor(prototype_id))) + facet_wrap(~prototype_id)














distances <- map(simulation_ids, function(simulation_id_oi) {
  stepinfo_oi <- simulation$stepinfo %>% filter(simulation_id == simulation_id_oi)
  expression_oi <- expression[stepinfo_oi$step_id, ]
  
  # for each simulation, calculate closest distances to the simulation of interest
  distances <- map(setdiff(simulation_ids, simulation_id_oi), function(simulation_id_ref) {
    stepinfo_or <- simulation$stepinfo %>% filter(simulation_id == simulation_id_ref)
    expression_or <- expression[stepinfo_or$step_id, ]
    
    FNN::knnx.dist(expression_or, expression_oi, 1)[, 1] %>% tibble(distance = ., cell_id = rownames(expression_oi), simulation_id_ref = simulation_id_ref, step_id=stepinfo_oi$step_id, step=stepinfo_oi$step, simulation_id = simulation_id_oi)
  }) %>% bind_rows()
  
  if (verbose) {
    distances_matrix <- distances %>% mutate(step_id = factor(step_id, levels=unique(step_id))) %>% reshape2::acast(step_id~simulation_id_ref, value.var="distance")
    distances_matrix %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  }
  
  distances
}) %>% bind_rows()



distances <- map(simulation_ids, function(simulation_id_oi) {
  stepinfo_oi <- simulation$stepinfo %>% filter(simulation_id == simulation_id_oi)
  expression_oi <- expression[stepinfo_oi$step_id, ]
  
  # for each simulation, calculate closest distances to the simulation of interest
  distances <- map(setdiff(simulation_ids, simulation_id_oi), function(simulation_id_ref) {
    stepinfo_or <- simulation$stepinfo %>% filter(simulation_id == simulation_id_ref)
    expression_or <- expression[stepinfo_or$step_id, ]
    
    FNN::knnx.dist(expression_or, expression_oi, 1)[, 1] %>% tibble(distance = ., cell_id = rownames(expression_oi), simulation_id_ref = simulation_id_ref, step_id=stepinfo_oi$step_id, step=stepinfo_oi$step, simulation_id = simulation_id_oi)
  }) %>% bind_rows()
  
  if (verbose) {
    distances_matrix <- distances %>% mutate(step_id = factor(step_id, levels=unique(step_id))) %>% reshape2::acast(step_id~simulation_id_ref, value.var="distance")
    distances_matrix %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  }
  
  distances
}) %>% bind_rows()













distances <- map(simulation_ids, function(simulation_id_oi) {
  stepinfo_oi <- simulation$stepinfo %>% filter(simulation_id == simulation_id_oi)
  expression_oi <- expression[stepinfo_oi$step_id, ]
  
  # for each simulation, calculate closest distances to the simulation of interest
  distances <- map(setdiff(simulation_ids, simulation_id_oi), function(simulation_id_ref) {
    stepinfo_or <- simulation$stepinfo %>% filter(simulation_id == simulation_id_ref)
    expression_or <- expression[stepinfo_or$step_id, ]
    
    FNN::knnx.dist(expression_or, expression_oi, 1)[, 1] %>% tibble(distance = ., cell_id = rownames(expression_oi), simulation_id_ref = simulation_id_ref, step_id=stepinfo_oi$step_id, step=stepinfo_oi$step, simulation_id = simulation_id_oi)
  }) %>% bind_rows()
  
  if (verbose) {
    distances_matrix <- distances %>% mutate(step_id = factor(step_id, levels=unique(step_id))) %>% reshape2::acast(step_id~simulation_id_ref, value.var="distance")
    distances_matrix %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  }
  
  distances
}) %>% bind_rows()





















changed <- TRUE
for (cutoff in c(0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6)) {
  adj <- center_network %>% reshape2::acast(from~to, value.var="from", fill=-1) %>% {.!=-1}
  
  A <- adj %*% t(adj)
  D <- (!adj) %*% t(!adj)
  n <- nrow(adj)
  
  jaccards <- (A+D) / (A+D+2 * (n - (A+D)))
  
  # (jaccards > 0.6) %>% igraph::graph_from_adjacency_matrix() %>% plot
  center_map <- (jaccards > cutoff) %>% igraph::graph_from_adjacency_matrix() %>% igraph::components() %>% .$membership
  
  center_network$from <- center_map[center_network$from]
  center_network$to <- center_map[center_network$to]
  
  if(nrow(center_network) == nrow(center_network %>% distinct(from, to))) {
    changed <- FALSE
  }
  
  print(length(unique(center_map)))
  
  center_network <- center_network %>% distinct(from, to)
  
  center_network %>% igraph::graph_from_data_frame() %>% plot()
}
























module_bytes <- 2**(1:ncol(samplexpression))
binary <- (samplexpression > 0.5) %>% as.numeric() %>% matrix(nrow=nrow(samplexpression))
samplestepinfo$binary <- (samplexpression > 0.5) %>% apply(1, function(x) sum(module_bytes[as.logical(x)]))

library(dtwclust)

binary_split <- binary %>% as.data.frame() %>% split(samplestepinfo$simulation_id) %>% map(as.matrix)
clust <- binary_split %>% tsclust(type = "hierarchical")

plot(clust)

# calculate membership
min_expression <- 0.5
module_bytes <- 2**(1:ncol(samplexpression))
membership <- (samplexpression > min_expression) %>% as.numeric() %>% matrix(nrow=nrow(samplexpression))
samplestepinfo$membership <- (samplexpression > min_expression) %>% apply(1, function(x) sum(module_bytes[as.logical(x)]))

samplestepinfo %>% ggplot(aes(step, membership, color=simulation_id, group=simulation_id)) + geom_line()

# do not allow flipping (due to noise)
samplestepinfo$run <- rle(samplestepinfo$membership) %>% {rep(seq_along(.$lengths), .$lengths)}
samplestepinfo <- samplestepinfo %>% group_by(simulation_id, run) %>% 
  nest() %>% 
  group_by(simulation_id) %>% 
  mutate(
    membership = map_dbl(data, ~.$membership[[1]]),
    next_membership = lead(membership, 1),
    prev_membership = lag(membership, 1)
  ) %>% 
  mutate(membership = ifelse(next_membership == prev_membership & !is.na(next_membership == prev_membership), next_membership, membership)) %>% 
  select(-next_membership, -prev_membership) %>% 
  mutate(data = map(data, ~select(., -membership))) %>% 
  unnest(data) %>% 
  ungroup()


# calculate binary network
samplestepinfo$run <- rle(samplestepinfo$membership) %>% {rep(seq_along(.$lengths), .$lengths)}
samplestepinfo_runs <- samplestepinfo %>% group_by(simulation_id, run) %>% 
  summarise(
    membership = membership[[1]],
    length = n()
  )

# construct the membership network,...
# by first grouping on prototype and creating a network
membership_network <- samplestepinfo_runs %>% group_by(simulation_id) %>% 
  mutate(from = membership, to = lead(membership, 1)) %>% 
  drop_na() %>% 
  select(from, to, simulation_id) %>% 
  ungroup() %>% 
  group_by(from, to) %>% 
  summarise(simulation_ids = list(unique(simulation_id)), n_simulations = length(unique(simulation_id)), n_runs = n()) %>% 
  ungroup()

membership_network$width <- membership_network$n_simulations / max(membership_network$n_simulations) * 10

membership_network %>% igraph::graph_from_data_frame() %>% plot




membership_network %>% igraph::graph_from_data_frame() %>% 
  igraph::distances(mode="out") %>% 
  pheatmap::pheatmap()

membership_network %>% igraph::graph_from_data_frame() %>% 
  igraph::distances(mode="out") %>% 
  t() %>% 
  dist() %>% 
  hclust() %>% 
  plot

membership_network %>% group_by(from, to) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  arrange(n)
