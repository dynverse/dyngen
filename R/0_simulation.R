#' @import tidyverse
generate_model = function(modulenetname, verbose=F) {
  model = load_modulenet(modulenetname)
  class(model) = "dyngen::model"
  
  model = modulenet_to_genenet(model$modulenet, model$modulenodes) %>% c(model)
  #mnet_wanted = model$modulenet
  #model = extract_network_from_modulenet(real_modulenet, mnet_wanted) %>% c(model)
  
  model$geneinfo = list_genes(model$geneinfo, model$modulemembership, model$net, model$modulenodes)
  if(verbose) plot_net_tfs(model)
  
  model = generate_formulae(model$net, model$geneinfo, model$celltypes) %>% c(model)
  model = generate_kinetics(model$vargroups, model$variables, model$nus.changes) %>% c(model)
  
  model$burngenes = determine_burngenes(model)
  
  model$info = list(date=date(), version=1, id=dambiutils:::random_time_string())
  
  model$ti = list(type=modulenetname, generator="modulenet_barabasi", )
  
  model
}


#' @import tidyverse
run_experiment = function(model, totaltime, burntime=2, nsimulations = 40, ncells = 500) {
  #newtime = estimate_convergence(8, verbose=T, totaltime=15) %>% max() %>% {.+1}
  
  #experiment = simulate_multiple_cells(model, burntime, totaltime, ncells=1000)
  #experiment = simulate_one_cell(model, burntime, totaltime, ncells=NULL, ssa.algorithm=ssa.em(noise_strength=2))
  experiment = simulate_multiple_cells_split(model, burntime, totaltime, nsimulations = nsimulations, ncellspersimulation = ceiling(ncells/nsimulations), ssa.algorithm=ssa.em(noise_strength=2), local=F)
  class(experiment) = "dyngen::experiment"
  
  experiment$model = model
  
  # get module counts
  experiment$expression_modules = get_module_counts(experiment$expression, experiment$model$modulemembership)
  
  pheatmap::pheatmap(SCORPIUS::quant.scale(experiment$expression_modules) %>% t, cluster_cols=F)
  pheatmap::pheatmap(SCORPIUS::quant.scale(experiment$expression) %>% t, cluster_cols=F, labels_row=experiment$model$geneinfo$name)
  experiment$expression_modules %>% reshape2::melt(varnames=c("cell", "module"), value.name="expression") %>% ggplot() + geom_histogram(aes(expression)) + facet_wrap(~module) + geom_vline(xintercept = 2)
  
  experiment$info = list(date=date(), version=1, id =dambiutils:::random_time_string(), modelid=model$info$id)
  
  experiment
}

#' @import tidyverse
run_scrnaseq = function(experiment, platform) {
  dataset = list()
  class(dataset) = "dyngen::dataset"
  dataset$counts = simulate_scrnaseq(experiment$expression, platform)
  dataset$platform = platform
  
  #pheatmap::pheatmap(SCORPIUS::quant.scale(dataset$counts) %>% t, cluster_cols=F)
  
  dataset$info = list(experimentid=experiment$info$id, platformname=platform$platformname, id=dambiutils:::random_time_string(), version=1)
  
  dataset
}

#' @import tidyverse
extract_goldstandard = function(experiment, verbose=F, seed=get.seed()) {
  if(is.null(experiment$simulations)) stop("requires multiple individual simulations to align to gold standard")
  gs = list()
  gs$piecenet = get_piecenet(experiment$model$statenet)
  gs$piecestates = get_piecestates(gs$piecenet)
  gs$piecestatenet = get_piecestatenet(gs$piecenet)
  
  if(verbose) print(">> smoothing")
  experiment$simulations = smoothe_simulations(experiment$simulations, experiment$model) # add smoothened expression values
  
  if(verbose) print(">> dividing")
  gs$pieces = divide_simulations(experiment$simulations, gs$piecestates, experiment$model)
  gs$pieces = gs$pieces %>% mutate(ncells=map_int(cells, length)) %>% group_by(piecestateid) %>% mutate(outlierness = abs((ncells - mean(ncells))/sd(ncells))) %>% filter(is.na(outlierness) | outlierness < 2) %>% ungroup # remove long pieces (probably skipped a state)
  
  if(verbose) print(">> generating reference")
  gs$reference = get_reference_expression(gs$pieces, experiment$model)
  
  if(verbose) print(">> assigning")
  gs$cellinfo = assign_progression(experiment$expression, gs$reference)
  
  pheatmap::pheatmap(SCORPIUS::quant.scale(experiment$expression_modules, 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, annotation_col=gs$cellinfo %>% select(piecestateid) %>% mutate(piecestateid=factor(piecestateid)) %>% as.data.frame %>% set_rownames(gs$cellinfo$cell))
  ggplot(gs$cellinfo) + geom_area(aes(simulationtime, group=piecestateid, fill=factor(piecestateid)), position="fill", stat="bin", bins=10)
  
  # if(verbose) print(">> calculating cell distances")
  # statenodes = gs$piecenet %>% rename(state=piece) %>% left_join(gs$reference$cellinfo %>% group_by(piecestateid) %>% summarise(maxprogression=max(progression)) %>% mutate(state=as.integer(piecestateid)), by="state")
  # snet = get_branchconnected_statenet(get_line_graph(gs$piecenet), statenodes)
  # gs$celldistances = get_cell_distances(gs$cellinfo %>% mutate(state=piecestateid), snet)
  
  if(verbose) print(">> extracting milestones")
  gs <- c(get_milestones(experiment, gs), gs)
  
  gs$info = list(experimentid=experiment$info$id, version=1, id=dambiutils:::random_time_string())
  
  gs
}