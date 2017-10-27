library(tidyverse)
library(dyngen)

params = list(
  model = list(
    # modulenet
    modulenet_name = "linear",
    # treeseed = 1,
    
    # network between tfs
    ngenes_per_module= function(n) sample(1, n, replace=TRUE), 
    edge_retainment = function(n) max(c(round(n/2), 1)),
    
    # extra targets
    target_adder_name = "realnet",
    realnet_name = "regulatorycircuits",
    damping = 0.05,
    ntargets_sampler = function() {sample(0, 1)}
  ),
  simulation = list(
    totaltime = 10,
    burntime = 0,
    local=TRUE
  ),
  experiment = list(
    # experiment setting
    samplesettings = list(type = "snapshot", ncells = 500),
    add_housekeeping = TRUE,
    n_housekeeping_genes = 100
  )
)# %>% list2env(.GlobalEnv)
params %>% walk(~list2env(., .GlobalEnv))

model <- invoke(generate_model_from_modulenet, params$model)
plot_net(model, label=TRUE, main_only = FALSE)
plot_modulenet(model)

simulation <- invoke(simulate_multiple, params$simulation, model$system)

housekeeping_reference_dataset <- readRDS("../dynalysis/analysis/data/datasets/real/cortical_interneuron_differentiation_frazer.rds")
params$experiment$housekeeping_reference_means <- get_housekeeping_reference_means(housekeeping_reference_dataset$counts)

experiment <- invoke(run_experiment, params$experiment, simulation, model)

module_counts <- get_module_counts(experiment$expression, experiment$geneinfo)

source("../dynmodular/dimred_wrappers.R")
space <- dimred_ica(experiment$expression, ndim = 2)
ggplot(space %>% as.data.frame() %>% bind_cols(experiment$cellinfo)) + geom_point(aes(Comp1, Comp2, color=simulationtime))

pheatmap::pheatmap(t(dynutils::scale_quantile(experiment$expression[order(experiment$cellinfo$simulationtime), ])), cluster_cols=F, cluster_rows=T)
pheatmap::pheatmap(t(dynutils::scale_quantile(module_counts[order(experiment$cellinfo$simulationtime), ])), cluster_cols=F, cluster_rows=F)
