## Plotting models
#' Model plotting
#' @param model The model
#' @param colorby By what to color
#' @param main_only Whether to only draw the main network
#' @param label Whether to label genes
#' @export
plot_model <- function(model) {
  list(
    plot_net(model, label=FALSE, main_only = FALSE),
    plot_modulenet(model),
    plot_net_overlaps(model)
  )
}

#' @rdname plot_model
#' 
#' @importFrom igraph layout.graphopt graph_from_data_frame plot.igraph
#' @importFrom grDevices rainbow
#' @export
plot_modulenet <- function(model) {
  graph <- igraph::graph_from_data_frame(model$modulenet, vertices = model$modulenodes)
  
  modulenames <- unique(model$geneinfo$module_id)
  colors <- rainbow(length(modulenames))
  igraph::V(graph)$color <- colors
  igraph::E(graph)$color <- c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(1,-1, 0)))]
  
  layout <- igraph::layout.graphopt(graph, charge=0.01, niter=10000)
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    vertex.size = 20, 
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1, 
    vertex.label.color = "black"
  )
}

#' @rdname plot_model
#' @importFrom grDevices rainbow rgb
#' @export
plot_net <- function(model, colorby=c("module", "main"), main_only=TRUE, label=FALSE) {
  colorby <- match.arg(colorby)
  
  geneinfo  <- model$geneinfo
  if(main_only) geneinfo <- geneinfo %>% filter(main)
  net <- model$net %>% filter((from %in% geneinfo$gene_id) & (to %in% geneinfo$gene_id))
  
  # add extra edges invisible between regulators from the same module
  net <- bind_rows(
    net, 
    geneinfo %>% group_by(module_id) %>% filter(n() > 1) %>% {split(., .$module_id)} %>% map(~as.data.frame(t(combn(.$gene_id, 2)), stringsAsFactors=FALSE)) %>% bind_rows() %>% bind_rows(tibble(V1=character(), V2=character())) %>% mutate(effect = -2) %>% rename(from = V1, to=V2)
  )
  
  
  graph <- igraph::graph_from_data_frame(net %>% select(from, to), vertices=geneinfo$gene_id)
  
  # layout
  set.seed(1) # to get same layout
  layout <- igraph::layout.fruchterman.reingold(graph)
  main_filter <- as.numeric(factor(geneinfo$main, levels = c(FALSE, TRUE, NA), exclude=NULL))
  
  # change vertex/edge colors and sizes
  igraph::V(graph)$size <- c(1, 4, 1)[main_filter]
  
  if (!label) {
    igraph::V(graph)$label <- rep("", length(igraph::V(graph)))
  }
  
  if (colorby == "main") {
    igraph::V(graph)$color <- c("black", "white", "grey")[main_filter]
  } else if (colorby == "module") {
    modulenames <- model$modulenodes$module_id
    colors <- rainbow(length(modulenames)) %>% set_names(modulenames)
    igraph::V(graph)$color <- colors[geneinfo$module_id]
  }
  igraph::E(graph)$color <- c("#d63737", "#3793d6", "#7cd637", grDevices::rgb(0, 0, 0, alpha=0))[as.numeric(factor(net$effect, levels = c(1,-1, 0, -2)))]
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1
  )
}

#' @rdname plot_model
#' @importFrom pheatmap pheatmap
#' @export
plot_net_overlaps <- function(model) {
  jaccard <- function(x, y) {length(intersect(x, y))/length(union(x,y))}
  pheatmap::pheatmap(sapply(model$geneinfo$gene_id, function(i) sapply(model$geneinfo$gene_id, function(j) jaccard(model$net$from[model$net$to==i], model$net$from[model$net$to==j]))))
}


## Plotting simulation
#' Simulation plotting
#' @param simulation The simulation
#' @param spaces Spaces dataframe
#' @param subsample Subsampled expression data
#' @export
plot_simulation <- function(simulation) {
  list(
    plot_simulation_space_individual(simulation),
    plot_simulation_space_time(simulation),
    plot_simulation_space_modules(simulation),
    plot_simulation_simulations_lines(simulation),
    plot_simulation_modules_lines(simulation),
    plot_simulation_modules_heatmap(simulation)
  )
}

subsample_simulation <- function(simulation) {
  if(is.null(simulation$expression_modules)) simulation <- preprocess_simulation_for_gs(simulation, model)
  
  sample_n <- round(nrow(simulation$stepinfo)/2000)
  samplestepinfo <- simulation$stepinfo %>% group_by(simulation_id) %>% filter((step %% sample_n) == 0) %>% ungroup()
  samplexpression <- simulation$expression[samplestepinfo$step_id, ]
  samplexpression_modules <- simulation$expression_modules[samplestepinfo$step_id, ]
  
  samplexpression <- samplexpression + runif(length(samplexpression), 0, 0.01)
  samplexpression_modules <- samplexpression_modules + runif(length(samplexpression_modules), 0, 0.01)
  
  lst(samplexpression, samplexpression_modules, samplestepinfo)
}

dimred_pca = function(x, ndim=3) {
  space = prcomp(t(x))$rotation[,seq_len(ndim)]
  process_dimred(space)
}
dimred_mds = function(x, ndim=3) {
  requireNamespace("SCORPIUS")
  space = SCORPIUS::reduce_dimensionality(SCORPIUS::correlation_distance(x),ndim = ndim)
  process_dimred(space)
}
dimred_ica = function(x, ndim=3) {
  requireNamespace("fastICA")
  space = fastICA::fastICA(t(scale(t(x))), ndim)$S
  process_dimred(space)
}
dimred_tsne = function(x, ndim=3) {
  requireNamespace("Rtsne")
  space = Rtsne::Rtsne(as.dist(SCORPIUS::correlation_distance(x)), dims = ndim, is_distance = TRUE,  perplexity = 20)$Y
  rownames(space) = rownames(x)
  process_dimred(space)
}
process_dimred = function(space) {
  space = as.matrix(space)
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}
dimred <- function(x, dimred_name="pca", ndim=2) {
  get(paste0("dimred_", dimred_name))(x, ndim=ndim)
}

dimred_simulation <- function(simulation, subsample = subsample_simulation(simulation), dimred_names = c("pca", "ica", "mds"), expression_names=c("samplexpression", "samplexpression_modules")) {
  map(
    dimred_names, 
    function(dimred_name) {
      map(
        expression_names, 
        function(expression_name) {
          expression <- subsample[[expression_name]]
          space <- expression %>% 
            dimred(dimred_name = dimred_name, ndim=2) %>% 
            as.data.frame() %>% 
            bind_cols(subsample$samplestepinfo)
          
          tibble(dimred_name = dimred_name, expression_name = expression_name, space=list(space))
        }
      ) %>% bind_rows()
    }
  ) %>% bind_rows()
}

#' @rdname plot_simulation
#' @export
plot_simulation_space_individual <- function(simulation, spaces=dimred_simulation(simulation)) {
  spaces %>% 
    unnest(space) %>% 
    ggplot(aes(Comp1, Comp2, group=simulation_id, color=factor(simulation_id))) + 
      geom_path() + 
      facet_wrap(dimred_name~expression_name, scales = "free") +
      theme_void() +
      theme(legend.position="none")
}

#' @rdname plot_simulation
#' @export
plot_simulation_space_time <- function(simulation, spaces=dimred_simulation(simulation)) {
  spaces %>% 
    unnest(space) %>% 
    ggplot(aes(Comp1, Comp2, group=simulation_id, color=simulationtime)) + 
      geom_path() + 
      facet_wrap(dimred_name~expression_name, scales = "free") +
      viridis::scale_color_viridis() +
      theme_void() +
      theme(legend.position="none")
}

#' @rdname plot_simulation
#' @export
plot_simulation_space_modules <- function(
  simulation, 
  subsample = subsample_simulation(simulation), 
  spaces=dimred_simulation(simulation, subsample, "ica", "samplexpression_modules")
) {
  module_ids <- colnames(simulation$expression_modules)
  samplexpression_modules_df <- subsample$samplexpression_modules %>% 
    reshape2::melt(varnames=c("step_id", "module_id"), value.name="expression") %>% 
    mutate_if(is.factor, as.character)
  
  spaces$space <- map(spaces$space, ~bind_cols(., ))
  
  spaces %>% 
    unnest(space) %>% 
    left_join(samplexpression_modules_df, by="step_id") %>% 
    ggplot(aes(Comp1, Comp2, group=simulation_id, color=expression)) + 
    geom_path() + 
    facet_wrap(dimred_name~module_id, scales = "free") +
    viridis::scale_color_viridis() +
    theme_void() +
    theme(legend.position="none")
}

#' @rdname plot_simulation
#' @export
plot_simulation_simulations_lines <- function(simulation) {
  subsample <- subsample_simulation(simulation)
  
  expression_df <- subsample$samplexpression_modules %>% 
    reshape2::melt(varnames=c("step_id", "module_id"), value.name="expression") %>% 
    mutate_if(is.factor, as.character) %>% 
    left_join(subsample$samplestepinfo, by="step_id")
  
  expression_df %>% 
    ggplot() + 
      geom_line(aes(simulationtime, expression, group=module_id, color=factor(module_id))) + 
      facet_wrap(~simulation_id) +
      theme_void()
}

#' @rdname plot_simulation
#' @export
plot_simulation_modules_lines <- function(simulation) {
  subsample <- subsample_simulation(simulation)
  
  expression_df <- subsample$samplexpression_modules %>% 
    reshape2::melt(varnames=c("step_id", "module_id"), value.name="expression") %>% 
    mutate_if(is.factor, as.character) %>% 
    left_join(subsample$samplestepinfo, by="step_id")
  
  expression_df %>% 
    ggplot() + 
    geom_line(aes(simulationtime, expression, group=simulation_id, color=factor(simulation_id))) + 
    facet_wrap(~module_id) +
    cowplot::theme_cowplot() +
    theme(legend.position="none", panel.grid.major.x=element_line(color="grey"))
}

#' @rdname plot_simulation
#' @export
plot_simulation_modules_heatmap <- function(simulation) {
  subsample <- subsample_simulation(simulation)
  
  subsample$samplexpression_modules %>% 
    reshape2::melt(varnames=c("step_id", "module_id"), value.name="expression") %>% 
    mutate_if(is.factor, as.character) %>% 
    left_join(subsample$samplestepinfo, by="step_id") %>% 
    ggplot() + 
        geom_raster(aes(step, factor(module_id), fill=expression, group=module_id)) + 
        facet_wrap(~simulation_id) + 
        scale_fill_distiller(palette="RdBu")
}

#' @rdname plot_simulation
#' @export
plot_simulation_3D <- function(simulation) {
  subsample <- subsample_simulation(simulation)
  
  space <- dimred_mds(subsample$samplexpression) %>% bind_cols(subsample$stepinfo)
  for (i in unique(as.numeric(space$simulation_id))) {
    rgl::lines3d(space %>% filter(simulation_id == i), col = grDevices::rainbow(length(unique(space$simulation_id)))[[i]])
  }
}


#' Gold standard plotting
#' @param simulation The simulation
#' @param gs The gold standard
#' @param spaces Spaces dataframe
#' @export
plot_goldstandard <- function(simulation, gs) {
  list(
    plot_goldstandard_edges(simulation, gs),
    plot_goldstandard_burn(simulation, gs),
    plot_goldstandard_network(simulation, gs),
    plot_goldstandard_heatmap(simulation, gs),
    plot_goldstandard_references(gs)
  )
}

dimred_goldstandard <- function(simulation, gs) {
  subsample <- subsample_simulation(simulation)
  spaces <- dimred_simulation(simulation, subsample = subsample, expression_names=c("samplexpression_modules"))
  step_ids <- subsample$samplestepinfo$step_id
  
  sampleprogressions <- gs$progressions %>% 
    slice(match(step_ids, step_id))
  
  spaces$space <- map(spaces$space, left_join, sampleprogressions)
  
  spaces
}

#' @rdname plot_goldstandard
#' @export
plot_goldstandard_edges <- function(simulation, gs, spaces = dimred_goldstandard(simulation, gs)) {
  spaces %>% 
    unnest(space) %>% 
    ggplot(aes(Comp1, Comp2, group=simulation_id, color=factor(edge_id))) + 
    geom_path() + 
    facet_wrap(dimred_name~expression_name, scales = "free") +
    theme_void()
}

#' @rdname plot_goldstandard
#' @export
plot_goldstandard_burn <- function(simulation, gs, spaces = dimred_goldstandard(simulation, gs)) {
  spaces %>% 
    unnest(space) %>% 
    ggplot(aes(Comp1, Comp2, group=simulation_id, color=!burn)) + 
    geom_path() + 
    facet_wrap(dimred_name~expression_name, scales = "free") +
    theme_void()
}

#' @rdname plot_goldstandard
#' @export
plot_goldstandard_network <- function(simulation, gs, spaces = dimred_goldstandard(simulation, gs)) {
  space <- spaces$space[[1]]
  
  milestone_ids <- unique(c(gs$milestone_network$from, gs$milestone_network$to))
  grouping <- convert_progressions_to_milestone_percentages(
    space$step_id, 
    milestone_ids, 
    gs$milestone_network, 
    gs$progressions %>% rename(cell_id = step_id) %>% slice(match(space$step_id, cell_id))
  ) %>% 
    group_by(cell_id) %>% filter(row_number() == which.max(percentage)) %>% rename(step_id = cell_id) %>% 
    ungroup()
  
  generate_milestone_space <- function(space, grouping, milestone_network) {
    space <- space %>% left_join(grouping, by="step_id")
    space_milestone <- space %>% group_by(milestone_id) %>% summarise(Comp1=mean(Comp1), Comp2=mean(Comp2))
    space_milestone_network <- milestone_network %>% 
      left_join(space_milestone %>% rename_all(~paste0(., "_from")) %>% rename(from = milestone_id_from), by="from") %>% 
      left_join(space_milestone %>% rename_all(~paste0(., "_to")) %>% rename(to = milestone_id_to), by="to")
    
    tibble(space_milestone = list(space_milestone), space_milestone_network = list(space_milestone_network))
  }
  
  spaces_milestone <- map(spaces$space, generate_milestone_space, grouping, gs$milestone_network) %>% bind_rows() %>% bind_cols(spaces)
  
  pmap(as.list(spaces_milestone), function(space, space_milestone, space_milestone_network, ...) {
    space <- space %>% mutate_at(., vars(simulation_id, edge_id), funs(factor))
    
    map(c("percentage", "simulationtime", "simulation_id", "edge_id", "milestones"), function(colorize) {
      space$color <- space[[colorize]]
      if(colorize %in% c("percentage", "simulationtime")) {
        color_scale <- viridis::scale_color_viridis(option="A")
      } else if(colorize %in% c("simulation_id", "edge_id")) {
        color_scale <- scale_color_discrete()
      }
      
      if (colorize != "milestones") {
        ggplot(space, aes(Comp1, Comp2)) + 
          geom_path(aes_string(color=colorize, group="simulation_id")) + 
          geom_point(aes_string(color=colorize)) + 
          color_scale + 
          theme_void() +
          theme(legend.position = "none") + 
          ggtitle(colorize)
      } else {
        ggplot() + 
          geom_point(aes(Comp1, Comp2, color=factor(milestone_id, levels=milestone_ids)), data=space_milestone, alpha=0.1) +
          ggnetwork::geom_edges(aes(Comp1_from, Comp2_from, xend=Comp1_to, yend = Comp2_to), space_milestone_network, color = "grey50") +
          ggnetwork::geom_edges(aes(Comp1_from+(Comp1_to - Comp1_from)/2.1, Comp2_from+(Comp2_to - Comp2_from)/2.1, xend=Comp1_from+(Comp1_to - Comp1_from)/1.9, yend = Comp2_from+(Comp2_to - Comp2_from)/1.9), space_milestone_network, color = "grey50", arrow=arrow(type="closed", length=unit(0.1, "inches"))) +
          ggnetwork::geom_nodelabel(aes(Comp1, Comp2, color = factor(milestone_id, levels=milestone_ids), label=milestone_id), space_milestone, size=5) +
          ggnetwork::theme_blank() + 
          theme(legend.position = "none") + 
          ggtitle(colorize)
      }
    })
  }) %>% unlist(recursive = F) %>% cowplot::plot_grid(plotlist=., nrow=length(spaces)) %>% print()
}

#' @rdname plot_goldstandard
#' @export
plot_goldstandard_heatmap <- function(simulation, gs) {
  subsample <- subsample_simulation(simulation)
  step_ids <- rownames(subsample$samplexpression)
  sampleprogressions <- gs$progressions %>% 
    slice(match(step_ids, step_id))
  
  sampleprogressions_ordered <- sampleprogressions %>% arrange(edge_id, percentage)
  samplexpression_ordered <- subsample$samplexpression[sampleprogressions_ordered$step_id, ]
  samplexpression_modules_ordered <- subsample$samplexpression_modules[sampleprogressions_ordered$step_id, ]
  
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
  
  samplexpression_modules_ordered %>% t() %>% pheatmap::pheatmap(
    cluster_cols=F, 
    cluster_rows=F, 
    gaps_col = which(diff(sampleprogressions_ordered$edge_id) != 0),
    annotation_col = sampleprogressions_ordered %>% 
      select(edge_id, simulation_id, percentage) %>% 
      mutate(edge_id = factor(edge_id), simulation_id = simulation_id) %>% 
      as.data.frame() %>% 
      magrittr::set_rownames(sampleprogressions_ordered$step_id)
  )
}

#' @rdname plot_goldstandard
#' @export
plot_goldstandard_references <- function(gs) {
  reference_data <- imap(gs$references, function(reference, i) {
    reference$reference_expression %>% reshape2::melt(varnames=c("step", "module"), value.name="expression") %>% mutate(reference_id = i)
  }) %>% bind_rows()
  
  reference_data %>% 
    ggplot() + 
      geom_raster(aes(step, factor(module), fill=expression, group=module)) + 
      facet_wrap(~reference_id) + 
      scale_fill_distiller(palette="RdBu") +
      cowplot::theme_cowplot()
}




#' Plotting functions for experiment
#'
#' @param experiment The experiment
#' @export
plot_experiment <- function(experiment) {
  plots <- map(c("expression_simulated", "expression", "true_counts", "counts"), function(expression_name) {
    expr <- experiment[[expression_name]]
    
    if(ncol(expr) > 500) {
      expr <- expr[, sample(colnames(expr), 500)]
    }    
    if(nrow(expr) > 500) {
      expr <- expr[sample(rownames(expr), 500), ]
    }
    
    space <- dimred_ica(expr, ndim=2)
    space %>% as.data.frame() %>% ggplot() + geom_point(aes(Comp1, Comp2)) + ggtitle(expression_name)
  })
  cowplot::plot_grid(plotlist = plots) %>% print()
}

#' Plotting functions for normalisation
#'
#' @param experiment The experiment
#' @param normalisation The normalisation
#' @export
plot_normalisation <- function(experiment, normalisation) {
  plots <- map(c("experiment$expression_simulated", "experiment$expression", "experiment$true_counts", "experiment$counts", "normalisation$count", "normalisation$expression"), function(expression_name) {
    expr <- eval(parse(text=expression_name))
    
    if(ncol(expr) > 500) {
      expr <- expr[, sample(colnames(expr), 500)]
    }    
    if(nrow(expr) > 500) {
      expr <- expr[sample(rownames(expr), 500), ]
    }
    
    space <- dimred_ica(expr, ndim=2)
    space %>% as.data.frame() %>% ggplot() + geom_point(aes(Comp1, Comp2)) + ggtitle(expression_name)
  })
  cowplot::plot_grid(plotlist = plots) %>% print()
}