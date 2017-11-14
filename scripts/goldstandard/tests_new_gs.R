pdf("hi.pdf")

plot_simulation_modules(simulation, model)

preprocess_simulation_for_gs <- function(simulation, model) {
  print("smoothing...")
  smooth_expression <- function(expression) {
    expression %>% 
      zoo::rollmean(20, c("extend", "extend", "extend")) %>%
      magrittr::set_rownames(rownames(expression))
  }
  
  # simulation$expression_smooth <- simulation$expression %>% as.data.frame() %>% split(simulation$stepinfo$simulation_id) %>% map(smooth_expression) %>% do.call(rbind, .)
  # dimnames(simulation$expression_smooth) <- dimnames(simulation$expression)
  
  print("normalizing...")
  simulation$expression_normalized <- t(dynutils::scale_quantile(t(simulation$expression), outlier_cutoff=0.05))
  
  print("calculating module expression...")
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

# calculate whether two simulations are connected
distance_connected_cutoff <- 0.5 # depends on noise levels
distances <- distances %>% 
  group_by(simulation_id, simulation_id_ref) %>% 
  mutate(connected = distance < distance_connected_cutoff) %>% 
  ungroup()
  #mutate(connected = zoo::rollapply(distance, 100, function(x) all(x < distance_connected_cutoff), align="left", partial=FALSE, fill=NA)) %>% 
  

distances %>% ggplot() + geom_histogram(aes(distance)) + geom_vline(xintercept = distance_connected_cutoff)
distances %>% sample_n(50000) %>% ggplot() + geom_line(aes(step, distance, color=factor(simulation_id_ref))) + facet_wrap(~simulation_id) + geom_hline(yintercept=distance_connected_cutoff)

connected <- distances %>% reshape2::acast(cell_id~simulation_id_ref, value.var="connected", mean, fill=1)
connected[is.na(connected)] <- 1 # same simulation -> always connected duhh

# connected2 <- connected[sample(nrow(connected), 5000), ]
# connected2 %>% pheatmap::pheatmap()

# cluster based on connectedness
stepinfo <- simulation$stepinfo
connected <- connected[stepinfo$step_id, unique(stepinfo$simulation_id)]
binary <- 2**(1:ncol(connected))
stepinfo <- stepinfo %>% mutate(cluster_id = apply(connected,1, function(x) sum(binary[as.logical(x)]))) 

# look at the consecutive length of the same connected cluster in a simulation (a "run")
stepinfo$run <- rle(stepinfo$cluster_id) %>% {rep(seq_along(.$lengths), .$lengths)}
run_info <- stepinfo %>% group_by(simulation_id, cluster_id, run) %>% summarise(run_length=n()) %>% ungroup()# %>% ggplot() + geom_point(aes(cluster_id, run_length, group=cluster_id, color=factor(simulation_id)))
run_info %>% ggplot() + geom_point(aes(cluster_id, log10(run_length), group=cluster_id, color=factor(simulation_id)))

# select cluster ids based on their run length
selected_cluster_ids <- run_info %>% 
  group_by(cluster_id) %>% 
  summarise(top_run_length = max(run_length)) %>% 
  mutate(nsimulations = binaryLogic::as.binary(cluster_id, size=length(binary)) %>% map_int(sum)) %>% 
  filter(top_run_length > 100, nsimulations > length(unique(run_info$simulation_id)) * 0.1) %>% 
  pull(cluster_id)

# now assign some parts as fuzzy
stepinfo$cluster_id[!(stepinfo$cluster_id %in% selected_cluster_ids)] <- NA
cluster_order <- c(NA, unique(na.omit(stepinfo$cluster_id)))
stepinfo$cluster_id <- match(stepinfo$cluster_id, cluster_order) - 1

stepinfo %>% ggplot() + geom_line(aes(step, as.numeric(cluster_id), color=factor(simulation_id))) + facet_wrap(~simulation_id)

# now reassign fuzzy if start and end are the same
stepinfo <- stepinfo %>% 
  group_by(simulation_id) %>% 
  mutate(cluster_id_na = ifelse(cluster_id == 0, NA, cluster_id)) %>% 
  mutate(
    cluster_id_notna_prev = zoo::na.locf(cluster_id_na, fromLast = FALSE, na.rm=FALSE),
    cluster_id_notna_next = zoo::na.locf(cluster_id_na, fromLast = TRUE, na.rm=FALSE)
  ) %>% 
  mutate(cluster_id = ifelse(cluster_id_notna_prev==cluster_id_notna_next, cluster_id_notna_next, cluster_id)) %>% 
  select(-cluster_id_na, -cluster_id_notna_prev, -cluster_id_notna_next) %>% 
  ungroup()

stepinfo %>% ggplot() + geom_line(aes(step, as.numeric(cluster_id), color=factor(simulation_id))) + facet_wrap(~simulation_id)

# use majority rule to reassign cluster_id, smoothing this out
window <- 100
stepinfo$cluster_id <- stepinfo %>% group_by(simulation_id) %>%
  mutate(cluster_id = zoo::rollapply(cluster_id, window, function(x) {table(x, exclude=NULL) %>% {names(.)[which.max(.)]}}, partial=TRUE)) %>%
  pull(cluster_id)
stepinfo$cluster_id <- as.numeric(stepinfo$cluster_id)

stepinfo %>% ggplot() + geom_line(aes(step, as.numeric(cluster_id), color=factor(simulation_id))) + facet_wrap(~simulation_id)

# add fuzzy to previous
stepinfo <- stepinfo %>% 
  group_by(simulation_id) %>% 
  mutate(cluster_id_na = ifelse(cluster_id == 0, NA, cluster_id),
         cluster_id= zoo::na.locf(cluster_id_na, fromLast = FALSE)) %>% 
  select(-cluster_id_na) %>% 
  ungroup()

stepinfo %>% ggplot() + geom_line(aes(step, as.numeric(cluster_id), color=factor(simulation_id))) + facet_wrap(~simulation_id)

# divide in consecutive pieces
stepinfo$run <- rle(stepinfo$cluster_id) %>% {rep(seq_along(.$lengths), .$lengths)}
stepinfo$piece_id <- group_indices(stepinfo %>% group_by(simulation_id, run))

pieces <- stepinfo %>% group_by(piece_id, simulation_id, cluster_id) %>% 
  summarise(expression = list(expression[step_id, ])) %>% 
  ungroup()
pieces <- pieces %>% 
  group_by(simulation_id) %>% 
  mutate(next_cluster_id = lead(cluster_id, 1)) %>% 
  ungroup()

# split pieces if self-similar (oscillatory)
pieces <- map(split(pieces, pieces$piece_id), function(piece) {
  subsample <- 10
  window <- round(500/subsample)
  piece_expression <- piece$expression[[1]]
  piece_expression <- piece_expression[seq(1, nrow(piece_expression), subsample), ]
  dist <- dist(piece_expression) %>% as.matrix()
  
  restarts <- c()
  prevstart <- 1
  for (i in (window*2):nrow(dist)) {
    if (i > prevstart + window*2) {
      local_dist <- dist[i, prevstart:(i-window * 2)]
      prev_dist <- dist[i, (i-window * 2):(i-window)]
      
      if(all(prev_dist > distance_connected_cutoff) & any(local_dist < distance_connected_cutoff)) {
        print(i)
        restarts <- c(restarts, i)
        prevstart <- i
      }
    }
  }
  
  if (piece$piece_id == 1) {
    dist %>% pheatmap::pheatmap(cluster_rows=F, cluster_cols=F, gaps_col=restarts)
  }
  
  subpieces_membership <- cumsum(seq_len(nrow(dist)) %in% restarts)
  subpieces_expression <- piece$expression[[1]] %>% as.data.frame() %>% split(subpieces_membership) %>% map(as.matrix)
  
  #dist[nrow(dist)/1, ] %>% zoo::rollapply(window, function(x) all(x > distance_connected_cutoff)) %>% plot
  #dist %>% pheatmap::pheatmap(cluster_rows=F, cluster_cols=F)
  
  subpieces <- map(subpieces_expression, ~mutate(piece, expression=list(.))) %>% 
    bind_rows() %>% 
    mutate(next_cluster_id = lead(cluster_id, 1))
  
  if(nrow(subpieces) > 1) {
    # remove first and last cycle
    #subpieces <- subpieces[c(-1, -nrow(subpieces)),]
  }
  subpieces
}) %>% bind_rows()
pieces$piece_id <- seq_len(nrow(pieces))

# look at each cluster and divide it in substates
library(dtwclust)
skip <- 20
median_cutoff <- 100
pieces <- map(split(pieces, pieces$cluster_id), function(pieces) {
  if (nrow(pieces) <= 2) {
    pieces$cluster_substate_id <- 1
  } else {
    ts <- pieces$expression
    ts <- ts %>% map(~.[seq(1, nrow(.), skip), , drop=F])
    
    distances <- map(ts, function(t1) {
      map_dbl(ts, function(t2) {
        min(dist(t1, t2))
      })
    }) %>% unlist() %>% matrix(nrow=length(ts))
    
    pieces$cluster_substate_id <- distances %>% {. < distance_connected_cutoff} %>% igraph::graph_from_adjacency_matrix() %>% igraph::components() %>% .$membership
    
    # using dtw clust, problems with incomplete cycles...
    # clust <- tsclust(ts, type="hierarchical", control=hierarchical_control(method="average"))
    # plot(clust)
    # labels <- cutree(clust, h = quantile(clust$height %>% keep(~. != 0), 0.05) * median_cutoff)
    # pieces$cluster_substate_id <- labels
  }
  pieces
}) %>% bind_rows() %>% arrange(piece_id) # again sort on piece_id to keep order

# define states and state network
pieces$state_id <- pieces %>% group_by(cluster_id, cluster_substate_id) %>% group_indices()
pieces <- pieces %>% 
  group_by(simulation_id) %>% 
  mutate(next_state_id = lead(state_id, 1))

state_network <- pieces %>% 
  select(state_id, next_state_id) %>% 
  rename(from = state_id, to = next_state_id) %>% 
  drop_na() %>% 
  group_by(from, to) %>% 
  summarise(n_transitions = n()) %>% 
  ungroup() %>% 
  mutate(edge_id = seq_along(from))

state_network %>% mutate(width = n_transitions/sum(n_transitions) * 10) %>% igraph::graph_from_data_frame(vertices=pieces$state_id %>% unique()) %>% plot

# add piece_id and state_id to stepinfo
stepinfo <- stepinfo %>% 
  select(-piece_id) %>% 
  left_join(
    pieces %>% 
      ungroup() %>% 
      mutate(step_id = map(expression, rownames)) %>% 
      select(piece_id, step_id, state_id) %>% 
      unnest(step_id), by="step_id"
  ) # add piece_id to stepinfo

# now extract time for every state taking into account the magnitude of changes in expression
# stepinfo <- map(pieces %>% split(pieces$state_id), function(subpieces) {
#   substepinfo <- subpieces %>% mutate(step_id=map(expression, rownames)) %>% select(-expression) %>% unnest(step_id) %>% ungroup() %>% left_join(stepinfo %>% select(step, simulationtime, step_id), by="step_id")
#   subexpression <- subpieces$expression %>% do.call(rbind, .)
#   
#   substepinfo$diff <- subexpression %>% 
#     as.data.frame() %>% 
#     split(substepinfo$piece_id) %>% 
#     map(~c(0, rowMeans(abs(diff(as.matrix(.)))))) %>% 
#     map(~./sum(.)) %>% 
#     unlist()
#   
#   substepinfo <- substepinfo %>% group_by(piece_id) %>% 
#     mutate(time = cumsum(diff))
#   
#   plot(substepinfo$simulationtime, substepinfo$time)
#   
#   substepinfo
# }) %>% bind_rows() %>% ungroup()

# load dimred for plotting
source("../dynmodular/dimred_wrappers.R")

# fix startCircle of princurve package
startCircle <- function (x) 
{
  d <- dim(x)
  n <- d[1]
  p <- d[2]
  xbar <- apply(x, 2, "mean")
  ray <- sqrt((scale(x, xbar, FALSE)^2) %*% rep(1, p))
  radius <- mean(ray)
  s <- cbind(radius * sin((pi * (1:101))/50), radius * cos((pi * 
                                                              (1:101))/50))
  if (p > 2) 
    s <- cbind(s, matrix(0, 101, p - 2))
  princurve::get.lam(x, s)
}
assignInNamespace("startCircle", startCircle, "princurve")


timeinfo <- map(pieces %>% split(pieces$state_id), function(subpieces) {
  substepinfo <- subpieces %>% mutate(step_id=map(expression, rownames)) %>% select(-expression) %>% unnest(step_id) %>% ungroup() %>% left_join(stepinfo %>% select(step, simulationtime, step_id), by="step_id")
  subexpression <- subpieces$expression %>% do.call(rbind, .)
  space <- dimred_pca(subexpression, ndim=2)
  
  samplexpression <- subexpression[sample(nrow(subexpression), min(4000, nrow(subexpression))), ]
  samplestepinfo <- substepinfo %>% slice(match(rownames(samplexpression), step_id))
  samplespace <- space[rownames(samplexpression),]
  
  # fit <- princurve::principal.curve(samplespace, plot.true = TRUE, trace = FALSE, smoother="periodic.lowess")
  # fit2 <- princurve::get.lam(space, fit$s, fit$tag)
  # times <- fit2$lambda
  
  # check if circular
  if(state_network %>% filter(from == subpieces$state_id[[1]], to == subpieces$state_id[[2]]) %>% nrow()) {
    smoother = "periodic.lowess"
  } else {
    smoother = "lowess"
  }
  
  init <- substepinfo$expression[substepinfo$expression %>% map(nrow) %>% which.max]
  
  fit <- princurve::principal.curve(samplexpression, start = init, plot.true = FALSE, trace = FALSE, smoother=smoother)
  fit2 <- princurve::get.lam(subexpression, fit$s, fit$tag)
  times <- fit2$lambda
  
  simulationtime_correlation <- substepinfo %>% mutate(time = times) %>% group_by(piece_id) %>% summarise(cor = cor(simulationtime, time)) %>% {mean(.$cor)}
  
  if(abs(simulationtime_correlation) < 0.5) {
    warning(glue::glue("low simulationtime correlation for state {subpieces$state_id[[1]]}: {simulationtime_correlation}"))
  }
  if(simulationtime_correlation < 0) {
    times <- max(times) - times
  }
  
  print(
    ggplot(samplespace %>% as.data.frame() %>% mutate(time=fit$lambda)) + geom_point(aes(Comp1, Comp2, color=time))
  )
  print(
    ggplot(space %>% as.data.frame() %>% mutate(time=times)) + geom_point(aes(Comp1, Comp2, color=time))
  )
  
  tibble(step_id = rownames(subexpression), time=times)
}) %>% bind_rows()

stepinfo <- stepinfo %>% select(-matches("^time$")) %>% left_join(timeinfo, by="step_id")

# final check & plotting

source("../dynmodular/dimred_wrappers.R")
samplestepinfo <- stepinfo %>% sample_n(2000) %>% arrange(step)
samplexpression <- simulation$expression[samplestepinfo$step_id, ]

spaces <- map(c(dimred_mds, dimred_pca), ~samplexpression %>% {. + runif(length(.), 0, 0.01)} %>% .(ndim=2) %>% as.data.frame() %>% bind_cols(samplestepinfo) %>% mutate(state_id = factor(state_id)))
map(spaces, function(space) {
  map(c("simulationtime", "state_id", "time"), function(colorize) {
    space$color <- space[[colorize]]
    ggplot(space, aes(Comp1, Comp2)) + geom_point(aes_string(color=colorize))
  })
}) %>% unlist(recursive = F) %>% cowplot::plot_grid(plotlist=.)

dev.off()
dev.off()
