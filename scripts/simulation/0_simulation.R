modulenetname="linear"; totaltime=5; burntime=2; nsimulations = 40; ncells = 500
run_experiment = function(modulenetname, totaltime, burntime=2, nsimulations = 40, ncells = 500, qsub.conf=NULL) {
  scriptfolder = "./scripts/simulation/"
  walk(c("1_network_generation.R", "2_formulae.R", "3_kinetics.R", "4_simulate.R", "5_scrnaseq.R", "6_goldstandard.R"), ~source(file.path(scriptfolder, .)))
  
  model = load_modulenet(modulenetname)
  class(model) = "dyngen::model"
  model = modulenet_to_genenet(model$modulenet, model$modulenodes) %>% c(model)
  model$geneinfo = list_genes(model$geneinfo, model$modulemembership, model$net, model$modulenodes)
  
  model = generate_formulae(model$net, model$geneinfo, model$celltypes) %>% c(model)
  model = generate_kinetics(model$vargroups, model$variables, model$nus.changes) %>% c(model)
  
  model$burngenes = determine_burngenes(model)
  
  #newtime = estimate_convergence(8, verbose=T, totaltime=15) %>% max() %>% {.+1}
  
  #experiment = simulate_multiple_cells(model, burntime, totaltime, ncells=1000)
  #experiment = simulate_one_cell(model, burntime, totaltime)
  experiment = simulate_multiple_cells_split(model, burntime, totaltime, nsimulations = nsimulations, ncellspersimulation = ceiling(ncells/nsimulations), qsub.conf=qsub.conf)
  class(experiment) = "dyngen::experiment"
  
  experiment$model = model
  
  # determine cell states & cell orderings within states
  experiment$expression_modules = get_module_counts(experiment$expression, experiment$model$modulemembership)
  
  pheatmap(SCORPIUS::quant.scale(experiment$expression_modules) %>% t, cluster_cols=F)
  experiment$expression_modules %>% reshape2::melt(varnames=c("cell", "module"), value.name="expression") %>% ggplot() + geom_histogram(aes(expression)) + facet_wrap(~module) + geom_vline(xintercept = 2)
  
  experiment$cellinfo = get_states_and_progression(experiment$expression_modules, experiment$model$modulenodes, ) %>% left_join(experiment$cellinfo, "cell")
  
  pheatmap(SCORPIUS::quant.scale(experiment$expression_modules, 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, annotation_col=experiment$cellinfo %>% select(state) %>% mutate(state=factor(state)) %>% as.data.frame %>% set_rownames(experiment$cellinfo$cell))
  ggplot(experiment$cellinfo) + geom_area(aes(simulationtime, group=state, fill=factor(state)), position="fill", stat="bin", bins=20)
  
  experiment
}


platforms = read_tsv("data/platforms.tsv")

run_scrnaseq = function(experiment, platform) {
  scriptfolder = "./scripts/simulation/"
  walk(c("5_scrnaseq.R"), ~source(file.path(scriptfolder, .)))
  
  dataset = experiment
  class(dataset) = "dyngen::dataset"
  
  dataset$counts = simulate_scrnaseq(dataset$expression, platform)
  dataset$platform = platform
  
  pheatmap(SCORPIUS::quant.scale(dataset$counts) %>% t, cluster_cols=F)
  
  dataset
}

generate_goldstandard = function(dataset) {
  
}