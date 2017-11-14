#' Model generation from modulenet
#' 
#' Default parameters are defined in individual functions and are inherited in `zzz_param_inheritance.R`
#' 
#' @param target_adder_name The method for adding targets. Only "realnet" is currently supported.
#' @inheritParams load_modulenet
#' @inheritParams generate_system
#' @inheritParams add_targets_realnet
#' @inheritParams modulenet_to_genenet 
#' @inheritParams generate_random_tree
#' 
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
  damping,
  
  # params for system
  samplers
) {
  # load modulenet
  if (modulenet_name == "tree") {
    stagenet <- generate_random_tree(treeseed)
    model <- from_stages_to_modulenet(stagenet)
  } else {
    model <- load_modulenet(modulenet_name)
  }
  
  # convert module network to gene network between modules
  model <- modulenet_to_genenet(
    model$modulenet, 
    model$modulenodes, 
    ngenes_per_module = ngenes_per_module, 
    edge_retainment = edge_retainment
  ) %>% c(model)
  
  # add some targets
  if (target_adder_name == "realnet") {
    model <- add_targets_realnet(
      model$net, model$geneinfo, 
      realnet_name = realnet_name, 
      damping = damping,
      ntargets_sampler = ntargets_sampler
    ) %>% dynutils::merge_lists(model, .)
  } else {
    print("not supported")
  }
  
  # randomize parameters 
  model$net <- randomize_network_parameters(model$net)
  
  # generate thermodynamics formulae & kinetics
  model$system <- generate_system(model$net, model$geneinfo, model$cells, samplers)
  
  model
}

run_experiment <- function(simulation, model, samplesettings, n_housekeeping_genes, housekeeping_reference_means, add_housekeeping=TRUE) {
  experiment <- take_experiment_cells(simulation, model, samplesettings)
  
  if (add_housekeeping) {
    additional_data <- add_housekeeping_poisson(experiment$expression, model$geneinfo, housekeeping_reference_means, n_housekeeping_genes)
    experiment$expression <- additional_data$expression
    experiment$geneinfo <- additional_data$geneinfo
  }
  
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