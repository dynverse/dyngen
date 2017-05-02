library(tidyverse)

.datasets_location <- "results"

datasetsinfo <- readRDS("results/datasets.rds")
dataset <- load_dataset(datasetsinfo$id[[1]], contents = contents_dataset(experiment=contents_experiment(simulations=TRUE)))

gs <- dataset$gs
experiment <- dataset$experiment
model <- dataset$model
experiment$expression_modules <- dyngen:::get_module_counts(experiment$expression, model$modulemembership)

task <- dyneval::wrap_ti_task_data(
  ti_type = dataset$model$modulenetname, 
  name = dataset$info$id,
  expression = dataset$experiment$expression,
  state_names = gs$milestone_names, 
  state_net = gs$milestone_net, 
  state_percentages = gs$milestone_percentages
)

#pheatmap::pheatmap(SCORPIUS::quant.scale(experiment$expression_modules, 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, annotation_col=gs$cellinfo %>% select(piecestateid) %>% mutate(piecestateid=factor(piecestateid)) %>% as.data.frame %>% set_rownames(gs$cellinfo$cell))

gs$reference$expression %>% dyngen:::get_module_counts(., model$modulemembership) %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=T, annotation_col=gs$reference$cellinfo, gaps_col = which(diff(as.numeric(gs$reference$cellinfo$piecestateid)) != 0))

dyneval::plotLearner.ti.default(task)
