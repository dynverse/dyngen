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
  simulation$expression_normalized <- t(dynutils::scale_quantile(t(simulation$expression), outlier_cutoff=0.05))
  
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
samplexpression <- expression[((stepinfo$step)%%50) == 1, ]
samplestepinfo <- stepinfo %>% slice(match(rownames(samplexpression), step_id))

samplexpression




ncells <- 100
nsimulations <- 3
samplestepinfo <- tibble(step = seq_len(ncells * nsimulations)) %>% mutate(simulationtime = rep(seq(0, 4, length.out = ncells), nsimulations), simulation_id = rep(seq_len(nsimulations), each=ncells))
samplexpression <- tibble("4" = samplestepinfo$simulationtime/4, "3"=`4`**(samplestepinfo$simulation_id/2)) %>% as.matrix() %>% magrittr::set_rownames(samplestepinfo$step)


# samplexpression <- c(1, 0, 0, 1, 0, 0.5) %>% matrix(nrow=3)
# samplestepinfo <- tibble(cell_id = seq_len(nrow(samplexpression))) %>% mutate(simulationtime = runif(n()))
# rownames(samplexpression) <- samplestepinfo$cell_id

source("../dynmodular/dimred_wrappers.R")

dimred <- function(expression) dimred_pca(expression, 2)
dimred <- function(expression) expression[, c("4", "3")] %>% magrittr::set_colnames(c("Comp1", "Comp2"))

plot_space <- function(space) {
  ggplot(space %>% as.data.frame() %>% bind_cols(samplestepinfo), aes(Comp1, Comp2)) + geom_path(aes(group=simulation_id)) + geom_point(aes(color=simulationtime))+ theme(legend.position="none")
}
space <- samplexpression %>% dimred()
plot_space(space)

# cell_dist_kernel <- function(x, original, rate=4) (dexp(x, rate)/rate - exp(-rate * original)) * 0.5
# cell_dist_kernel <- function(x, original, rate = 2) dnorm(original - x) * rate
# cell_dist_kernel <- function(x, original, rate=2) (dexp(x, rate=rate)) / rate
cell_dist_kernel <- function(x, original, rate=2) 

curexpression <- samplexpression
originalexpression <- curexpression
speed <- 1
settings <- tibble()
spaces <- list()
original_cell_dist <- curexpression %>% dist() %>% as.matrix()
# time_dist_kernel <- function(x, original, rate=4) (dexp(x, rate)/rate - exp(-rate * original)) * 0.1
time_dist_kernel <- function(x, original, rate=0.5) 1 - (dexp(abs(original - x), rate=rate)) / rate
time_dist_kernel <- function(x, original, rate=20) (dexp(abs(original - x), rate=rate)) / rate
time_dist_kernel <- function(x, original, sd=0.5) (dnorm(original - x, sd=sd) / dnorm(0, sd=sd))
# xs <- seq(0, 1, 0.01);map2_dbl(xs, 0.5, time_dist_kernel) %>% plot(xs, .)
time_dist_mask <- (as.matrix(dist(seq_along(samplestepinfo$step), "maximum")) == 1) & (as.matrix(dist(samplestepinfo$simulation_id)) == 0)
# time_dist_mask[lower.tri(time_dist_mask)] <- 0
# matrix(as.numeric(time_dist_mask), nrow=nrow(cell_dist)) %>% pheatmap::pheatmap(cluster_rows=F, cluster_cols=F)

unit_direction <- function(diffs_cell) {
  diffs_cell <- diffs_cell / sqrt(rowSums(t(apply(diffs_cell, 1, function(x) x **2))))
  diffs_cell[is.na(diffs_cell)] <- 0
  diffs_cell
}

for (i in 1:200) {
  cell_dist <- curexpression %>% dist() %>% as.matrix()
  time_weighing <- time_dist_kernel(original_cell_dist, cell_dist)
  time_weighing[!time_dist_mask] <- NA
  change_weighing <- time_weighing %>% rowMeans(na.rm=TRUE)
  
  cell_weighing <- cell_dist_kernel(original_cell_dist, original_cell_dist)
  
  diffs <- map(rownames(curexpression), function(cell_id) t(t(curexpression) - curexpression[cell_id, ]))
  names(diffs) <- rownames(curexpression)
  
  change_func <- function(cell_id) (diffs[[cell_id]] * cell_weighing[cell_id, ]) %>% colMeans
  change_diffs <- map(rownames(curexpression), change_func) %>% do.call(rbind, .)
  
  # time_weighing <- time_dist_kernel(cell_dist, original_cell_dist)
  # time_weighing[!time_dist_mask] <- 0
  # change_func <- function(cell_id) {
  #   weighing <- time_weighing[cell_id,]
  #   weighing <- weighing[weighing > 0]
  #   if (length(weighing) == 0) {
  #     weighing <- time_weighing[cell_id,1:2]
  #   }
  #   (unit_direction(diffs[[cell_id]][names(weighing), ,drop=F]) * weighing) %>% colMeans()
  # }
  # change_time <- map(rownames(curexpression), change_func) %>% do.call(rbind, .)
  # change_time %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  
  # change <- (change_diffs - change_time) * speed
  weights <- matrix(rep(change_weighing, ncol(change_diffs)), ncol=ncol(change_diffs))
  
  change <- (change_diffs * speed)
  
  curexpression <- curexpression + change
  
  print(i)
  
  if ((i %% 10) == 1) {
    space <- curexpression %>% dimred()
    spaces <- c(spaces, list(plot_space(space)))
  }
}
spaces %>% cowplot::plot_grid(plotlist = .)

curexpression %>% dist() %>% as.matrix() %>% hist()












space <- curexpression %>% t() %>% scale() %>% t %>% dimred_mds(2)
ggplot(space %>% as.data.frame() %>% bind_cols(samplestepinfo)) + geom_point(aes(Comp1, Comp2, color=simulationtime))



samplexpression %>% pheatmap::pheatmap(scale="column")
samplexpression %>% pheatmap::pheatmap()
curexpression %>% pheatmap::pheatmap(scale="column")
curexpression %>% pheatmap::pheatmap()

curexpression[, 5] %>% plot(samplestepinfo$simulationtime, .)
