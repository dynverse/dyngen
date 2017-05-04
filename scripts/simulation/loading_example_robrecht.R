library(tidyverse)
library(cowplot)

.datasets_location <- "results"

datasetsinfo <- readRDS("results/datasets.rds")
dataset <- load_dataset(datasetsinfo$id[[1]], contents = contents_dataset(experiment=contents_experiment(simulations=TRUE)))

gs <- dataset$gs
experiment <- dataset$experiment
model <- dataset$model
experiment$expression_modules <- dyngen:::get_module_counts(experiment$expression, model$modulemembership)

task <- with(dataset, dyneval::wrap_ti_task_data(
  ti_type = model$modulenetname,
  name = info$id,
  expression = log2(counts+1),
  state_names = gs$milestone_names,
  state_net = gs$milestone_net,
  state_percentages = gs$milestone_percentages_notent %>% slice(match(rownames(counts), id))
))

pred_output1 <- dyneval::trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)

task_emdist <- dyneval::compute_emlike_dist(task)

pred_emdist1 <- dyneval::compute_emlike_dist(pred_output1)

corank1 <- dyneval::compute_coranking(task_emdist, pred_emdist1)

cowplot::plot_grid(
  dyneval::plotLearner.ti.default(task),
  dyneval::plotLearner.ti.default(pred_output1)
)

corank1$summary