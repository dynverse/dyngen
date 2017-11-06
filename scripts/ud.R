library(tidyverse)
library(dynutils)
#library(dyngen)

params <- simple_params
params$model$modulenet_name <- "trifurcating"
params$simulation$local <- 8

# model
model <- invoke(generate_model_from_modulenet, params$model)
plot_net(model, label=TRUE, main_only = FALSE)
plot_modulenet(model)

# simulation
simulation <- invoke(simulate_multiple, params$simulation, model$system)
plot_simulation_modules(simulation, model)

# experiment
experiment <- invoke(run_experiment, params$experiment, simulation, model)

module_counts <- get_module_counts(experiment$expression, experiment$geneinfo)

source("../dynmodular/dimred_wrappers.R")
space <- dimred_ica(experiment$expression, ndim = 2)
ggplot(space %>% as.data.frame() %>% bind_cols(experiment$cellinfo)) + geom_point(aes(Comp1, Comp2, color=simulationtime)) + viridis::scale_color_viridis(option="A")

pheatmap::pheatmap(t(dynutils::scale_quantile(experiment$expression[order(experiment$cellinfo$simulationtime), ])), cluster_cols=F, cluster_rows=T)
pheatmap::pheatmap(t(dynutils::scale_quantile(module_counts[order(experiment$cellinfo$simulationtime), ])), cluster_cols=F, cluster_rows=F)


expression_df <- experiment$expression %>% 
  reshape2::melt(varnames=c("cell_id", "gene_id"), value.name="expression") %>% 
  left_join(experiment$geneinfo %>% rename(celltype_id = cell_id), by="gene_id") %>% 
  left_join(experiment$cellinfo, by="cell_id") %>% 
  filter(!is.na(module_id)) %>% 
  filter(main)
expression_df %>% ggplot() + geom_smooth(aes(simulationtime, expression, group=gene_id, color=factor(module_id))) + facet_wrap(~simulation_id)
