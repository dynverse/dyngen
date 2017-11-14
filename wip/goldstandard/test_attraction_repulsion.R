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
samplexpression <- expression[((stepinfo$step)%%30) == 1, ]
samplestepinfo <- stepinfo %>% slice(match(rownames(samplexpression), step_id))

samplexpression


##

ncells <- 100
nsimulations <- 3
samplestepinfo <- tibble(step = seq_len(ncells * nsimulations)) %>% mutate(simulationtime = rep(seq(0, 4, length.out = ncells), nsimulations), simulation_id = rep(seq_len(nsimulations), each=ncells))
samplexpression <- tibble("4" = samplestepinfo$simulationtime/4, "3"=`4`**(samplestepinfo$simulation_id/2)) %>% as.matrix() %>% magrittr::set_rownames(samplestepinfo$step)


# samplexpression <- c(1, 0, 0, 1, 0, 0.5) %>% matrix(nrow=3)
# samplestepinfo <- tibble(cell_id = seq_len(nrow(samplexpression))) %>% mutate(simulationtime = runif(n()))
# rownames(samplexpression) <- samplestepinfo$cell_id

source("../dynmodular/dimred_wrappers.R")

# dimred <- function(expression) dimred_ica(expression, 2)
dimred <- function(expression) expression[, c("4", "3")] %>% magrittr::set_colnames(c("Comp1", "Comp2"))

plot_space <- function(space) {
  ggplot(space %>% as.data.frame() %>% bind_cols(samplestepinfo), aes(Comp1, Comp2)) + geom_path(aes(group=simulation_id)) + geom_point(aes(color=simulationtime))+ theme(legend.position="none")
}
space <- samplexpression %>% dimred()
plot_space(space)

# cell_dist_kernel <- function(x, original, rate=4) (dexp(x, rate)/rate - exp(-rate * original)) * 0.5
# cell_dist_kernel <- function(x, original, rate = 2) dnorm(original - x) * rate
cell_dist_kernel <- function(x, original, rate=2) (dexp(x, rate=rate)) / rate
# cell_dist_kernel <- function(x, original, rate=2) dnorm(original - x) * rate
  
curexpression <- samplexpression
originalexpression <- curexpression
speed <- 1
settings <- tibble()
spaces <- list()
original_cell_dist <- curexpression %>% dist() %>% as.matrix()

unit_direction <- function(diffs_cell) {
  diffs_cell <- diffs_cell / sqrt(rowSums(t(apply(diffs_cell, 1, function(x) x **2))))
  diffs_cell[is.na(diffs_cell)] <- 0
  diffs_cell
}

for (i in 1:20) {
  cell_dist <- curexpression %>% dist() %>% as.matrix()
  
  diffs <- map(rownames(curexpression), function(cell_id) t(t(curexpression) - curexpression[cell_id, ]))
  names(diffs) <- rownames(curexpression)
  
  cell_weighing <- cell_dist_kernel(original_cell_dist, cell_dist)
  
  change_func <- function(cell_id) (diffs[[cell_id]] * cell_weighing[cell_id, ]) %>% colMeans
  change_diffs <- map(rownames(curexpression), change_func) %>% do.call(rbind, .)
  
  densities <- apply(cell_dist, 1, function(y) quantile(y, 0.9))
  densities <- densities/max(densities)
  
  change <- (change_diffs * speed) * densities **10
  
  curexpression <- curexpression + change
  
  print(i)
  
  if ((i %% 2) == 1) {
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
