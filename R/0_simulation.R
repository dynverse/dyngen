# Model generation ==========================================================================
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
  # module net --
  modulenet_name, 
  
  # params for modulenet_name == "tree"
  treeseed,
  decay,
  
  # module net to net --
  ngenes_per_module,
  edge_retainment,
  ntargets_sampler,
  
  # add targets --
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

# Simulation =================================================================
#' Simulate multiple cells
#' 
#' @param system The system to simulate
#' @param burntime The burn-in time before sampling
#' @param totaltime The total simulation time
#' @param nsimulations The number of simulations
#' @param local Whether or not to use \code{\link[PRISM]{qsub_lapply}}
#' @param ssa_algorithm Which GSSA algorithm to use
#' 
#' @importFrom PRISM qsub_lapply
#' @importFrom pbapply pblapply
#' @export
simulate_multiple <- function(system, burntime, totaltime, nsimulations = 16, local=FALSE, ssa_algorithm = fastgssa::ssa.em(noise_strength=4)) {
  force(system) # force the evaluation of the system argument, as the qsub environment will be empty except for existing function arguments
  if(!local) {
    multilapply = function(x, fun) {PRISM::qsub_lapply(x, fun, qsub_environment = list2env(list()))}
  } else {
    multilapply = function(x, fun) {pbapply::pblapply(x, fun, cl = getOption("ncores"))}
  }
  
  seeds <- sample.int(nsimulations * 10, nsimulations, replace = F)
  simulations = multilapply(seq_len(nsimulations), function(i) {
    set.seed(seeds[[i]]) # set seed, to avoid the same seeds in multiple cells (the case when eg. using pblapply)
    
    cell = simulate_cell(system, totaltime=totaltime, burntime=burntime, ssa_algorithm=ssa_algorithm)
    
    rownames(cell$molecules) = paste0(i, "_", seq_len(nrow(cell$molecules)))
    
    expression = cell$molecules[,stringr::str_detect(colnames(cell$molecules), "x_")]
    colnames(expression) = gsub("x_(.*)", "\\1", colnames(expression))
    stepinfo = tibble(step_id = rownames(cell$molecules), step=seq_along(cell$times), simulationtime=cell$times, simulation_id=i)
    
    lst(molecules=cell$molecules, stepinfo=stepinfo, expression=expression)
  })
  
  molecules <- map(simulations, "molecules") %>% do.call(rbind, .)
  expression <- map(simulations, "expression") %>% do.call(rbind, .)
  stepinfo <- map(simulations, "stepinfo") %>% do.call(bind_rows, .)
  
  lst(molecules, expression, stepinfo)
}








# Gold standard =============================================================
extract_goldstandard <- function(simulation, model, reference_length, max_path_length, smooth_window, verbose=TRUE, preprocess=TRUE) {
  if(preprocess) {
    if (verbose) print("Preprocessing")
    simulation <- preprocess_simulation_for_gs(simulation, model, smooth_window=smooth_window)
  }
  
  expression <- simulation$expression_modules
  stepinfo <- simulation$stepinfo
  
  start_milestones <- unique(model$edge_operations %>% filter(start) %>% pull(from))
  milestone_network <- model$edge_operations %>% mutate(edge_id = seq_len(n())) %>% select(from, to, edge_id, burn)
  
  if (verbose) print("Extracting milestone paths")
  paths <- get_milestone_paths(milestone_network, start_milestones, max_path_length)
  
  if (verbose) print("Processing operations")
  operations <- process_operations(model$edge_operations, model$modulenodes$module_id)
  
  path_operations <- extract_path_operations(operations, paths)
  
  if (verbose) print("Extracting references")
  references <- extract_references(path_operations, milestone_network, reference_length=reference_length)
  
  if (verbose) print("Mapping simulations onto reference")
  simulation_expressions <- expression %>% as.data.frame() %>% split(stepinfo$simulation_id)
  times <- map_to_reference(simulation_expressions, references)
  
  if (verbose) print("Postprocessing")
  gs <- list()
  gs$progressions <- stepinfo %>% left_join(times, by="step_id") %>% left_join(milestone_network, by=c("edge_id"))
  gs$milestone_network <- milestone_network
  gs$references <- references
  gs$expression_modules <- simulation$expression_modules
  
  gs
  
  # when things go wrong:
  # pheatmap::pheatmap(references[[8]]$reference_expression, cluster_rows=F, cluster_cols=F, show_rownames = FALSE)
  # pheatmap::pheatmap(simulation_expressions[[1]], cluster_rows=F, cluster_cols=F, show_rownames = FALSE)
  # 
}

check_goldstandard <- function(gs) {
  gs$progressions$edge_id <- factor(gs$progressions$edge_id, levels=gs$milestone_network$edge_id)
  edge_counts <- table(gs$progressions$edge_id)
  if (any(edge_counts == 0)) {warning("Some edges not represented!")}
  
  lst(edge_counts)
}

#' @importFrom cowplot plot_grid
plot_goldstandard <- function(simulation, model, gs) {
  if(is.null(simulation$expression_modules)) simulation <- preprocess_simulation_for_gs(simulation, model)
  
  # Dimensionality reduction
  source("../dynmodular/dimred_wrappers.R")
  # samplexpression <- gs$expression_modules[gs$progressions$step_id[(gs$progressions$step)%%10 == 0], ]
  samplexpression <- gs$expression_modules[sample(gs$progressions$step_id, 1000), ]
  colnames(samplexpression) <- paste0("M", colnames(samplexpression))
  sampleprogressions <- gs$progressions %>% slice(match(rownames(samplexpression), step_id)) %>% arrange(simulation_id, step)
  samplexpression <- samplexpression[sampleprogressions$step_id, ]
  spaces <- map(c(dimred_pca, dimred_ica, dimred_lmds), ~samplexpression %>% {. + runif(length(.), 0, 0.01)} %>% .(ndim=2) %>% as.data.frame() %>% bind_cols(sampleprogressions) %>% bind_cols(samplexpression %>% as.data.frame()))
  
  # Extract milestone grouping
  milestone_ids <- unique(c(gs$milestone_network$from, gs$milestone_network$to))
  grouping <- dynutils::convert_progressions_to_milestone_percentages(
    space$step_id, 
    milestone_ids, 
    gs$milestone_network, 
    gs$progressions %>% rename(cell_id = step_id) %>% slice(match(sampleprogressions$step_id, cell_id))
  ) %>% 
    group_by(cell_id) %>% filter(row_number() == which.max(percentage)) %>% rename(step_id = cell_id) %>% 
    ungroup()
  
  generate_milestone_space <- function(space, grouping, milestone_network) {
    space <- space %>% left_join(grouping, by="step_id")
    space_milestone <- space %>% group_by(milestone_id) %>% summarise(Comp1=mean(Comp1), Comp2=mean(Comp2))
    space_milestone_network <- milestone_network %>% 
      left_join(space_milestone %>% rename_all(~paste0(., "_from")) %>% rename(from = milestone_id_from), by="from") %>% 
      left_join(space_milestone %>% rename_all(~paste0(., "_to")) %>% rename(to = milestone_id_to), by="to")
    
    lst(space_milestone, space_milestone_network, space)
  }
  
  # Plot gold standard
  map(spaces, function(space) {
    space <- space %>% mutate_at(., vars(simulation_id, edge_id), funs(factor))
    
    map(c("percentage", "simulationtime", "simulation_id", "edge_id", "milestones"), function(colorize) {
      space$color <- space[[colorize]]
      if(colorize %in% c("percentage", "simulationtime")) {
        color_scale <- viridis::scale_color_viridis(option="A")
      } else if(colorize %in% c("simulation_id", "edge_id")) {
        color_scale <- scale_color_discrete()
      }
      
      if (colorize != "milestones") {
        ggplot(space, aes(Comp1, Comp2)) + geom_path(aes_string(color=colorize, group="simulation_id")) + geom_point(aes_string(color=colorize)) + color_scale + theme(legend.position = "none") + ggtitle(colorize)
      } else {
        milestone_space <- generate_milestone_space(space, grouping, gs$milestone_network)
        
        ggplot() + 
          geom_point(aes(Comp1, Comp2, color=factor(milestone_id, levels=milestone_ids)), data=milestone_space$space, alpha=0.1) +
          ggnetwork::geom_edges(aes(Comp1_from, Comp2_from, xend=Comp1_to, yend = Comp2_to), milestone_space$space_milestone_network, color = "grey50") +
          ggnetwork::geom_edges(aes(Comp1_from+(Comp1_to - Comp1_from)/2.1, Comp2_from+(Comp2_to - Comp2_from)/2.1, xend=Comp1_from+(Comp1_to - Comp1_from)/1.9, yend = Comp2_from+(Comp2_to - Comp2_from)/1.9), milestone_space$space_milestone_network, color = "grey50", arrow=arrow(type="closed", length=unit(0.1, "inches"))) +
          ggnetwork::geom_nodelabel(aes(Comp1, Comp2, color = factor(milestone_id, levels=milestone_ids), label=milestone_id), milestone_space$space_milestone, size=5) +
          ggnetwork::theme_blank() + 
          theme(legend.position = "none") + 
          ggtitle(colorize)
      }
    })
  }) %>% unlist(recursive = F) %>% cowplot::plot_grid(plotlist=., nrow=length(spaces)) %>% print()
  
  # Plot heatmap
  sampleprogressions_ordered <- sampleprogressions %>% arrange(edge_id, percentage)
  samplexpression_ordered <- samplexpression[sampleprogressions_ordered$step_id, ]
  
  samplexpression_ordered %>% t() %>% pheatmap::pheatmap(
    cluster_cols=F, 
    cluster_rows=F, 
    gaps_col = which(diff(sampleprogressions_ordered$edge_id) != 0),
    annotation_col = sampleprogressions_ordered %>% 
      select(edge_id, simulation_id, percentage) %>% 
      mutate(edge_id = factor(edge_id), simulation_id = simulation_id) %>% 
      as.data.frame() %>% 
      magrittr::set_rownames(sampleprogressions_ordered$step_id)
  )
  
  # Plot modules
  space <- spaces[[1]]
  space <- space %>% mutate_at(., vars(simulation_id, edge_id), funs(factor))
  map(colnames(samplexpression), function(colorize) {
    space$color <- space[[colorize]]
    color_scale <- viridis::scale_color_viridis(option="A")
    ggplot(space, aes(Comp1, Comp2)) + geom_path(aes_string(color=colorize, group="simulation_id")) + geom_point(aes_string(color=colorize)) + color_scale + theme(legend.position = "none") + ggtitle(colorize)
  }) %>% cowplot::plot_grid(plotlist=.) %>% print()
  
  # reference and simulation expression changes
  reference_data <- imap(gs$references, function(reference, i) {
    reference$reference_expression %>% reshape2::melt(varnames=c("step", "module"), value.name="expression") %>% mutate(reference_id = i)
  }) %>% bind_rows()
  
  reference_data %>% 
    {ggplot(.) + 
      geom_raster(aes(step, factor(module), fill=expression, group=module)) + 
      facet_wrap(~reference_id) + 
      scale_fill_distiller(palette="RdBu")} %>% 
    print()
  
  # 
  samplexpression <- simulation$expression_modules[simulation$stepinfo$step_id[(simulation$stepinfo$step)%%10 == 1], ]
  sampleprogressions <- simulation$stepinfo %>% slice(match(rownames(samplexpression), step_id)) %>% arrange(simulation_id, step)
  samplexpression <- samplexpression[sampleprogressions$step_id, ]
  
  samplexpression %>% 
    reshape2::melt(varnames=c("step_id", "module"), value.name="expression") %>% 
    left_join(sampleprogressions, by="step_id") %>% 
    ggplot() + 
      geom_raster(aes(step, factor(module), fill=expression, group=module)) + 
      facet_wrap(~simulation_id) + 
      scale_fill_distiller(palette="RdBu")
}






# Experiment =====================================================
#' Run an experiment from a simulation
#' @param simulation The simulation
#' @param model The model
#' @param samplesettings How the cells should be sampled
#' @param n_housekeeping_genes Number of housekeeping genes to add
#' @param housekeeping_reference_means Reference means for housekeeping genes
run_experiment <- function(
  simulation, 
  model, 
  samplesettings, 
  n_housekeeping_genes, 
  housekeeping_reference_means
) {
  experiment <- take_experiment_cells(simulation, model, samplesettings)
  
  if (n_housekeeping_genes > 0) {
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
