#' Model initializations
#' @export
generate_model_from_modulenet <- function(model) {
  if(!(c("modulenodes", "modulenet", "celltypes", "states") %in% names(model))) stop("Invalid model")
  
  # convert module network to gene network
  model = modulenet_to_genenet(model$modulenet, model$modulenodes, 0) %>% c(model)
  
  # list information of the genes, such as their module relationships
  model$geneinfo = list_genes(model$geneinfo, model$modulemembership, model$net, model$modulenodes)
  
  # generate thermodynamics formulae & kinetics
  model = generate_formulae(model$net, model$geneinfo, model$celltypes) %>% c(model)
  model = generate_kinetics(model$vargroups, model$variables, model$nus.changes) %>% c(model)
  
  # determine which genes will be burn in genes
}



generate_model = function(modulenetname=NULL, treeseed=NULL, genestart_id=0, verbose=F) {
  if(!is.null(modulenetname)) {
    model = load_modulenet(modulenetname)
  } else if(!is.null(treeseed)) {
    statenet = generate_random_tree(treeseed)
    model = from_states_to_modulenet(statenet)
  } else {
    stop("Can't generate model if you don't tell me how!")
  }
  
  model = modulenet_to_genenet(model$modulenet, model$modulenodes, genestart_id) %>% c(model)
  #mnet_wanted = model$modulenet
  #model = extract_network_from_modulenet(real_modulenet, mnet_wanted) %>% c(model)
  
  model$geneinfo = list_genes(model$geneinfo, model$modulemembership, model$net, model$modulenodes)
  if(verbose) plot_net_tfs(model)
  
  model = generate_formulae(model$net, model$geneinfo, model$celltypes) %>% c(model)
  model = generate_kinetics(model$vargroups, model$variables, model$nus.changes) %>% c(model)
  
  model$burngenes = determine_burngenes(model)
  
  model$info = list(date=date())
  model$ti = list(type=modulenetname, generator="modulenet_barabasi", generatorinfo = list())
  
  model
}

run_simulations <- function(model, totaltime, burntime = 2, nsimulations = 40, local = F) {
  simulations <- simulate_multiple(model, burntime, totaltime, nsimulations, local, ssa.algorithm=fastgssa::ssa.em(noise_strength=2))
  lst(simulations, model)
}

run_experiment <- function(simulation, takesettings, add_housekeeping=TRUE) {
  experiment = take_experiment_cells(simulation$simulations, takesettings)
  
  additional_data = add_housekeeping_poisson(experiment$expression, simulation$model$geneinfo)
  experiment$expression = additional_data$expression
  experiment$geneinfo = additional_data$geneinfo
  
  experiment
}

#' @importFrom pheatmap pheatmap
#' @importFrom pheatmap pheatmap
plot_experiment <- function(experiment) {
  experiment$expression_modules = get_module_counts(experiment$expression, experiment$model$modulemembership)
  
  pheatmap::pheatmap(dynutils::scale_quantile(experiment$expression_modules) %>% t, cluster_cols=F)
  pheatmap::pheatmap(dynutils::scale_quantile(experiment$expression) %>% t, cluster_cols=F, labels_row=experiment$model$geneinfo$name)
  experiment$expression_modules %>% reshape2::melt(varnames=c("cell", "module"), value.name="expression") %>% ggplot() + geom_histogram(aes(expression)) + facet_wrap(~module) + geom_vline(xintercept = 2)
}

run_scrnaseq = function(experiment, platform) {
  dataset = list()
  class(dataset) = "dyngen::dataset"
  dataset$counts = simulate_scrnaseq(experiment$expression, platform)
  dataset$platform = platform
  
  #pheatmap::pheatmap(dynutils::scale_quantile(dataset$counts) %>% t, cluster_cols=F)
  
  dataset
}

#' @importFrom pheatmap pheatmap
#' @importFrom dynutils random_time_string
#' @importFrom stats sd
extract_goldstandard = function(experiment, verbose=F, seed=get.seed()) {
  if(is.null(experiment$simulations)) stop("requires multiple individual simulations to align to gold standard")
  gs = list()
  gs$piecenet = get_piecenet(experiment$model$statenet)
  gs$states = get_states(gs$piecenet)
  gs$statenet = get_statenet(gs$piecenet)
  
  if(verbose) print(">> smoothing")
  experiment$simulations = smoothe_simulations(experiment$simulations, experiment$model) # add smoothened expression values
  
  if(verbose) print(">> dividing")
  gs$pieces = divide_simulations(experiment$simulations, gs$states, experiment$model)
  gs$pieces = gs$pieces %>%
    mutate(ncells=map_int(cells, length)) %>% 
    group_by(state_id) %>%
    mutate(outlierness = abs((ncells - mean(ncells))/stats::sd(ncells))) %>% 
    filter(is.na(outlierness) | outlierness < 2) %>%
    ungroup # remove long pieces (probably skipped a state)
  
  if(verbose) print(">> generating reference")
  gs$reference = get_reference_expression(gs$pieces, experiment$model)
  
  if(verbose) print(">> assigning")
  gs$cellinfo = assign_progression(experiment$expression, gs$reference)
  
 pheatmap::pheatmap(dynutils::scale_quantile(experiment$expression_modules, 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, annotation_col=gs$cellinfo %>% select(state_id) %>% mutate(state_id=factor(state_id)) %>% as.data.frame %>% set_rownames(gs$cellinfo$cell))
  ggplot(gs$cellinfo %>% left_join(experiment$cellinfo, by="cell")) + geom_area(aes(simulationtime, group=state_id, fill=factor(state_id)), position="fill", stat="bin", bins=10)
  
  # if(verbose) print(">> calculating cell distances")
  # statenodes = gs$piecenet %>% rename(state=piece) %>% left_join(gs$reference$cellinfo %>% group_by(stateid) %>% summarise(maxprogression=max(progression)) %>% mutate(state=as.integer(state_id)), by="state")
  # snet = get_branchconnected_statenet(get_line_graph(gs$piecenet), statenodes)
  # gs$celldistances = get_cell_distances(gs$cellinfo %>% mutate(state=state_id), snet)
  
  if(verbose) print(">> extracting milestones")
  #gs <- c(get_milestones(experiment, gs), gs)
  
  gs$info = list(experiment_id=experiment$info$id, id=paste0(.version, "/", dambiutils::random_time_string()))
  
  gs
}