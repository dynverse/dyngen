#' Model generation from modulenet
#' Default parameters are defined in individual functions and are inherited in `zzz_param_inheritance.R`
#' 
#' @param model the module model -- is a list which must contain modulenodes, modulenet, celltypes and states
#' 
#' @inheritParams add_targets_realnet
#' @inheritParams modulenet_to_genenet
#' @export
generate_model_from_modulenet <- function(
  # module net ---------------
  modulenet_name, 
  
  # params for modulenet_name == "tree"
  treeseed,
  decay,
  
  # module net to net ------------
  ngenes_per_module,
  edge_retainment,
  ntargets_sampler,
  
  # add targets -----------------------
  target_adder_name, 
  
  # params for target_adder_name == "realnet_name"
  realnet_name,
  damping
) {
  # load modulenet
  if (modulenet_name == "tree") {
    stagenet = generate_random_tree(treeseed)
    model = from_stages_to_modulenet(stagenet)
  } else {
    model <- load_modulenet(modulenet_name)
  }
  
  # convert module network to gene network between modules
  model <- modulenet_to_genenet(
    model$modulenet, 
    model$modulenodes, 
    ngenes_per_module=ngenes_per_module, 
    edge_retainment=edge_retainment
  ) %>% c(model)
  
  # add some targets
  if (target_adder_name == "realnet") {
    model <- add_targets_realnet(
      model$net, model$geneinfo, 
      realnet_name = realnet_name, 
      damping=damping,
      ntargets_sampler = ntargets_sampler
    ) %>% dynutils::merge_lists(model, .)
  } else {
    print("not supported")
  }
  
  # randomize parameters 
  model$net <- randomize_network_parameters(model$net)
  
  # generate thermodynamics formulae & kinetics
  model <- generate_formulae(model$net, model$geneinfo, model$cells) %>% c(model)
  model <- generate_kinetics(model$vargroups, model$variables, model$nus_df, model$formulae) %>% c(model)
  model$burn_variables <- determine_burn_variables(model)
  
  model
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