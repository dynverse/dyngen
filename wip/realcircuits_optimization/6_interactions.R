base_model <- load_circuit("cell_cycle_arabidopsis")

plot_net(base_model)

base_model <- generate_formulae(base_model$net, base_model$geneinfo, base_model$celltypes) %>% c(base_model)


name2gene <- base_model$geneinfo$gene %>% set_names(base_model$geneinfo$name)
gene2name <- base_model$geneinfo$name %>% set_names(base_model$geneinfo$gene)

system <- base_model

degradations <- read_tsv("data/genecircuits/cell_cycle_arabidopsis/degradations.tsv") %>% mutate(from=name2gene[from], to=name2gene[to]) %>% filter(!is.na(special))
for (i in seq_len(nrow(degradations))) {
  degradation <- degradations[i, ]
  to <- degradation$to
  from <- degradation$from
  strength <- ifelse(is.null(degradation$strength), 1, degradation$strength)
  
  formula_name <- paste0("protein_degradation_", to)
  degradation_K <- add.variable(fvar(paste0("q_", to, "_", from)), type="degradation", gene=to, system=system, strength=strength, randomize=degradation$randomize)
  system$formulae[[formula_name]] <- system$formulae[[formula_name]] + fvar(paste0("y_", to)) * (degradation_K * fvar(paste0("y_", from))) / (2 + fvar(paste0("y_", from)))
}

associations <- read_tsv("data/genecircuits/cell_cycle_arabidopsis/associations.tsv") %>% mutate(from=name2gene[from], to=name2gene[to])
associations$assoc <- map2_chr(associations$from, associations$to, ~system$geneinfo %>% filter(name==paste0(gene2name[.x], "_", gene2name[.y])) %>% .$gene)


for (i in seq_len(nrow(associations))) {
  association <- associations[i, ]
  to <- association$to
  from <- association$from
  assocgene <- association$assoc
  strength <- ifelse(is.null(association$strength), 1, association$strength)
  
  print(from)
  print(to)
  print(assocgene)
  
  association_K <- add.variable(fvar(paste0("assoc_", to, "_", from)), type="association", from=from, to=to, system=system, strength=strength, randomize=association$randomize)
  add.formula(fvar(paste0("y_", from)) * fvar(paste0("y_", to)) * association_K, set_names(c(-1, -1, 1), c(paste0("y_", from), paste0("y_", to), paste0("y_", assocgene))), name=paste0("association_", assocgene), system=system)
  
  deassociation_K <- add.variable(fvar(paste0("deassoc_", to, "_", from)), type="deassociation", from=from, to=to, system=system)
  add.formula(fvar(paste0("y_", assocgene)) * deassociation_K, set_names(c(1, 1, -1), c(paste0("y_", from), paste0("y_", to), paste0("y_", assocgene))), paste0("deassociation_", assocgene), system=system)
}

system$formulae.strings <- map_chr(system$formulae, function(fl) fl@string)

base_model <- system
base_model <- generate_kinetics(base_model$vargroups, base_model$variables, base_model$nus.changes) %>% c(base_model)
base_model$burngenes <- determine_burngenes(base_model)

experiment <- run_cycle(base_model)
fitness_cycle(experiment)

#experiment$molecules %>% pheatmap(cluster_rows=F)

Goi = base_model$geneinfo$gene[match(c("RBR", "E2Fb", "RBR_E2Fb", "CYCB", "CDKA_CYCB", "SCF", "CYCD", "CDKA_CYCD", "CDKA"), base_model$geneinfo$name)]
experiment$molecules[, paste0("y_", Goi)] %>% pheatmap(cluster_cols=F, cluster_rows=F, labels_col=paste(Goi, gene2name[Goi]))

Goi = base_model$geneinfo$gene
experiment$molecules[, paste0("y_", Goi)] %>% pheatmap(cluster_cols=F, cluster_rows=F, labels_col=paste(Goi, gene2name[Goi]))

##

randomize <- function() {
  randomizer_functions <- list(k = function(k) 10^runif(1, -2, 2), q = function(q) exp(runif(1, log(50), log(1000))), assoc = function(assoc) runif(1, 10, 100), deassoc = function(assoc) runif(1, 0.1, 2))
  to_randomize <- keep(base_model$variables, function(.) (!is.null(.$randomize)) && .$randomize) %>% map_chr("name")# %>% sample(2)
  #notto_randomize <- keep(base_model$variables, function(.) (is.null(.$randomize)) || .$randomize) %>% map_chr("name")
  new_params = map_dbl(to_randomize, function(paramname) {
    paramtype = gsub("([A-Za-z0-9]*)_.*", "\\1", paramname)
    randomizer_functions[[paramtype]]()
  })
  base_model$params[names(new_params)] = new_params
  base_model
}

models = map(1:100, ~randomize())
#saveRDS(models, "models.rds")
#models = readRDS("models.rds")


scores = mclapply(models, function(model) {
  library(tidyverse)
  library(magrittr)
  experiment <- run_cycle(model, totaltime=50, burntime=5)
  fitness_cycle(experiment)
}, mc.cores = 8)
handle = qsub_lapply(models, function(model) {
  library(tidyverse)
  library(magrittr)
  experiment <- run_cycle(model, totaltime=50, burntime=5)
  fitness_cycle(experiment)
}, override_qsub_config(wait=F))
saveRDS(handle, file="handle.rds")
handle = readRDS("handle.rds")
scores = PRISM::qsub_retrieve(handle)
plot(unlist(scores))

scores <- unlist(scores)
scores[is.na(scores)] <- 0

model <-models[[which.max(scores)]]
experiment <- run_cycle(model)
fitness_cycle(experiment)

model$info = generate_model_info()

model %>% saveRDS("~/tmp/cycle_model.rds")














model1 = readRDS("~/tmp/cycle_model.rds")
model2 = add_targets_shared(model1$net, model1$geneinfo, tfs=model1$geneinfo$gene, add_to_existing_net = FALSE)
model2 = generate_formulae(model2$net, model2$geneinfo) %>% c(model2)
model2 = generate_kinetics(model2$vargroups, model2$variables, model2$nus.changes) %>% c(model2)
model = merge_models(model1, model2)
model$info = generate_model_info()

experiment = run_cycle(model, totaltime = 100)





model3 = generate_model("linear", genestartid = max(model2$geneinfo$gene))
model = merge_models(model1, model2)
model = merge_models(model, model3)

experiment = run_cycle(model, totaltime = 20)
experiment = run_experiment(model, totaltime = 20, burntime=5)

filter_expression = function(expression, ncells=min(nrow(expression), 1000), Goi=colnames(expression)) {
  if(!is.null(ncells)) {
    expression = expression[sample(seq_len(nrow(expression)), ncells), ]
  }
  expression[, colnames(expression) %in% as.character(Goi)]
}

source("scripts/evaluation/methods/dimred.R")
expression = filter_expression(experiment$expression, NULL, model2$geneinfo$gene) %>% {log2(.+1)}
space = ica(expression, ndim = 2)

plotdata = bind_cols(space %>% as.data.frame, experiment$cellinfo %>% slice(match(rownames(space), cell)))
ggplot(plotdata) + geom_point(aes(Comp2, Comp1, color=simulationtime)) + viridis::scale_color_viridis(option="A")

expression = filter_expression(experiment$expression)
space1 = ica(filter_expression(experiment$expression, NULL, model3$geneinfo$gene), ndim = 1)
space2 = ica(filter_expression(experiment$expression, NULL, model2$geneinfo$gene), ndim = 2)


space = cbind(space1, space2)

library(rgl)
rgl::plot3d(space[, 1], space[, 2], space[, 3])









graph <- base_model$net %>% select(from, to) %>% mutate(type="regulation") %>% bind_rows(associations %>% mutate(type="association")) %>% bind_rows(degradations %>% mutate(type="degradation")) %>% graph_from_data_frame()
plot(graph, edge.color=RColorBrewer::brewer.pal(9, "Set1")[factor(edge.attributes(graph)$type)])




add_targets_shared(model$net, ldtfs = unique(model$net$from))
