library(tidyverse)

datasetsinfo <- readRDS("results/datasets.rds")
dataset <- load_dataset(datasetsinfo$id[[1]])

gs <- dataset$gs
experiment <- dataset$experiment
model <- dataset$model
experiment$expression_modules <- get_module_counts(experiment$expression, model$modulemembership)
pheatmap::pheatmap(SCORPIUS::quant.scale(experiment$expression_modules, 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, annotation_col=gs$cellinfo %>% select(piecestateid) %>% mutate(piecestateid=factor(piecestateid)) %>% as.data.frame %>% set_rownames(gs$cellinfo$cell))

gs$reference$expression %>% get_module_counts(., model$modulemembership) %>% t %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=T, annotation_col=gs$reference$cellinfo, gaps_col = which(diff(as.numeric(gs$reference$cellinfo$piecestateid)) != 0))


task <- dyneval::wrap_ti_task_data("blablabluh", "blub", expression = dataset$expression, state_names = dataset$gs$milestone_names, state_net = dataset$gs$milestone_net, state_percentages = dataset$gs$milestone_percentages)
dyneval::plotLearner.ti.default(task)