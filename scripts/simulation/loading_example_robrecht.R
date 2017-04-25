library(tidyverse)

datasetsinfo = readRDS("results/datasets.rds")
dataset = load_dataset(datasetsinfo$id[[20]])



task = dyneval::wrap_ti_task_data("blablabluh", "blub", expression = dataset$expression, state_names = dataset$gs$milestone_names, state_net = dataset$gs$milestone_net, state_percentages = dataset$gs$milestone_percentages)
dyneval::plotLearner.ti.default(task)
