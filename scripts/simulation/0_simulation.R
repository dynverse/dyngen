scriptfolder = "./scripts/simulation/"
walk(c("1_network_generation.R", "2_formulae.R", "3_kinetics.R", "4_simulate.R", "5_goldstandard.R"), ~source(file.path(scriptfolder, .)))

model = load_modulenet("consecutive_bifurcating")
class(model) = "dyngen::model"
model = modulenet_to_genenet(model$modulenet, model$modulenodes) %>% c(model)
model$geneinfo = list_genes(model$geneinfo, model$modulemembership, model$net, model$modulenodes)

model = generate_formulae(model$net, model$geneinfo, model$celltypes) %>% c(model)
model = generate_kinetics(model$vargroups, model$variables, model$nus.changes) %>% c(model)

model$burngenes = determine_burngenes(model)

totaltime = 8
burntime = 2
#ncells = 500

#newtime = estimate_convergence(8, verbose=T, totaltime=15) %>% max() %>% {.+1}

dataset = simulate_multiple_cells(model, burntime, totaltime, ncells=1000)
dataset = simulate_one_cell(model, burntime, totaltime)
dataset = simulate_multiple_cells_split(model, burntime, totaltime, nsimulations = 40, ncellspersimulation = 20)

dataset$model = model
dataset$counts = simulate_scrnaseq(dataset$expression)

pheatmap(SCORPIUS::quant.scale(dataset$counts) %>% t, cluster_cols=F)

dataset$expression_modules = get_module_counts(dataset$expression, dataset$model$modulemembership)

pheatmap(SCORPIUS::quant.scale(dataset$expression_modules) %>% t, cluster_cols=F)
dataset$expression_modules %>% reshape2::melt(varnames=c("cell", "module"), value.name="expression") %>% ggplot() + geom_histogram(aes(expression)) + facet_wrap(~module) + geom_vline(xintercept = 2)

dataset$cellinfo = get_states_and_progression(dataset$expression_modules) %>% left_join(dataset$cellinfo, "cell")

# determine cell states & cell orderings within states
pheatmap(SCORPIUS::quant.scale(dataset$expression_modules, 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, annotation_col=dataset$cellinfo %>% select(state) %>% mutate(state=factor(state)) %>% as.data.frame %>% set_rownames(dataset$cellinfo$cell))
ggplot(dataset$cellinfo) + geom_area(aes(simulationtime, group=state, fill=factor(state)), position="fill", stat="bin", bins=20)
