#' Plot a modulenet
#' 
#' @param model The model
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


#' Plot a gene network
#' @param model The model
#' @param colorby By what to color
#' @param main_only Whether to only draw the main network
#' @param label Whether to label genes
#' @importFrom grDevices rainbow rgb
#' @export
plot_net <- function(model, colorby=c("module", "main"), main_only=TRUE, label=FALSE) {
  colorby <- match.arg(colorby)
  
  geneinfo  <- model$geneinfo
  if(main_only) geneinfo <- geneinfo %>% filter(main)
  net <- model$net %>% filter((from %in% geneinfo$gene_id) & (to %in% geneinfo$gene_id))
  
  net <- bind_rows(
    net, 
    geneinfo %>% group_by(module_id) %>% filter(n() > 1) %>% {split(., .$module_id)} %>% map(~as.data.frame(t(combn(.$gene_id, 2)))) %>% bind_rows() %>% mutate(effect = -2) %>% rename(from = V1, to=V2)
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

#' Plot the overlap in targets in a heatmap
#' @param model The model
#' @importFrom pheatmap pheatmap
plot_net_overlaps <- function(model) {
  jaccard <- function(x, y) {length(intersect(x, y))/length(union(x,y))}
  pheatmap::pheatmap(sapply(model$geneinfo$gene_id, function(i) sapply(model$geneinfo$gene_id, function(j) jaccard(model$net$from[model$net$to==i], model$net$from[model$net$to==j]))))
}




plot_simulation_modules <- function(simulation, model) {
  expression_df <- simulation$expression %>% 
    reshape2::melt(varnames=c("step_id", "gene_id"), value.name="expression") %>% 
    left_join(model$geneinfo %>% rename(celltype_id = cell_id), by="gene_id") %>% 
    left_join(simulation$stepinfo, by="step_id") %>% 
    filter((step %% 10) == 1) %>% 
    filter(!is.na(module_id)) %>% 
    filter(main)
  expression_df %>% ggplot() + geom_smooth(aes(simulationtime, expression, group=gene_id, color=factor(module_id))) + facet_wrap(~simulation_id)
}


#' 
#' #' Plot the simulations
#' #' 
#' #' @param simulations A list of simulations
#' #' @param samplingrate A sampling rate for sampling which cells to plot
#' #' 
#' #' @importFrom grDevices rainbow
#' #' @importFrom stats sd
#' #' @importFrom magrittr set_rownames
#' plot_simulations <- function(simulations, samplingrate=0.1) {
#'   requireNamespace("rgl")
#'   for (i in seq_len(length(simulations))) {
#'     sample <- sample(c(T, F), size=nrow(simulations[[i]]$expression), T, c(samplingrate, 1-samplingrate))
#'     
#'     ix <- simulations[[i]]$expression[sample,] %>% apply(1, mean) %>% {. == 0}
#'     sample[ix] <- FALSE
#'     
#'     simulations[[i]]$subcellinfo <- simulations[[i]]$cellinfo[sample,]
#'     simulations[[i]]$subexpression <- simulations[[i]]$expression[sample,]
#'     
#'     # simulations[[i]]$expression[sample,] %>% apply(1, stats::sd) %>% {. == 0} %>% sum %>% print
#'     # simulations[[i]]$subexpression %>% {apply(., 1, stats::sd) == 0} %>% sum %>% print
#'   }
#'   
#'   overallexpression <- map(simulations, "subexpression") %>% do.call(rbind, .)
#'   overallexpression <- overallexpression %>% magrittr::set_rownames(1:nrow(overallexpression))
#'   overallcellinfo <- assign_progression(overallexpression, reference)
#'   overallcellinfo$simulation_id <- map(seq_len(length(simulations)), ~rep(., nrow(simulations[[.]]$expression))) %>% unlist %>% factor()
#'   overallcellinfo$observed <- T
#'   
#'   space <- lmds(overallexpression, 3) %>% as.data.frame() %>% bind_cols(overallcellinfo)
#'   
#'   ggplot() + 
#'     geom_path(aes(Comp1, Comp2, group=simulation_id, color=progression), space %>% filter(observed)) +
#'     geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3)
#'   ggplot() + 
#'     geom_path(aes(Comp1, Comp2, group=simulation_id, color=progression), space %>% filter(observed)) + 
#'     geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3) + 
#'     viridis::scale_color_viridis(option="A")
#'   ggplot() + 
#'     geom_path(aes(V1, V2, group=simulation_id, color=state_id), space %>% filter(observed)) + 
#'     geom_path(aes(V1, V2), data=space %>% filter(!observed), size=3)
#'   
#'   for (i in unique(as.numeric(space$simulation_id))) {
#'     rgl::lines3d(space %>% filter(simulation_id == i), col = grDevices::rainbow(length(unique(space$simulation_id)))[[i]])
#'   }
#'   
#' }




#' Plot the gold standard
#' @importFrom cowplot plot_grid
#' @param simulation The simulation
#' @param model The model
#' @param gs The gold standard
plot_goldstandard <- function(simulation, model, gs) {
  if(is.null(simulation$expression_modules)) simulation <- preprocess_simulation_for_gs(simulation, model)
  
  # Dimensionality reduction
  source(paste0(find.package("dyngen"), "/ext_data/dimred_wrappers.R"))
  # samplexpression <- gs$expression_modules[gs$progressions$step_id[(gs$progressions$step)%%10 == 0], ]
  samplexpression <- gs$expression_modules[sample(gs$progressions$step_id, 1000), ]
  colnames(samplexpression) <- paste0("M", colnames(samplexpression))
  sampleprogressions <- gs$progressions %>% slice(match(rownames(samplexpression), step_id)) %>% arrange(simulation_id, step)
  samplexpression <- samplexpression[sampleprogressions$step_id, ]
  spaces <- map(c(dimred_pca, dimred_ica, dimred_mds), ~samplexpression %>% {. + runif(length(.), 0, 0.01)} %>% .(ndim=2) %>% as.data.frame() %>% bind_cols(sampleprogressions) %>% bind_cols(samplexpression %>% as.data.frame()))
  
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
  
  # reference expression changes
  reference_data <- imap(gs$references, function(reference, i) {
    reference$reference_expression %>% reshape2::melt(varnames=c("step", "module"), value.name="expression") %>% mutate(reference_id = i)
  }) %>% bind_rows()
  
  reference_data %>% 
  {ggplot(.) + 
      geom_raster(aes(step, factor(module), fill=expression, group=module)) + 
      facet_wrap(~reference_id) + 
      scale_fill_distiller(palette="RdBu")} %>% 
    print()
  
  # module expression changes
  samplexpression <- simulation$expression_modules[simulation$stepinfo$step_id[(simulation$stepinfo$step)%%10 == 1], ]
  sampleprogressions <- simulation$stepinfo %>% slice(match(rownames(samplexpression), step_id)) %>% arrange(simulation_id, step)
  samplexpression <- samplexpression[sampleprogressions$step_id, ]
  
  samplexpression %>% 
    reshape2::melt(varnames=c("step_id", "module"), value.name="expression") %>% 
    left_join(sampleprogressions, by="step_id") %>% 
    {ggplot(.) + 
    geom_raster(aes(step, factor(module), fill=expression, group=module)) + 
    facet_wrap(~simulation_id) + 
    scale_fill_distiller(palette="RdBu")} %>% 
    print()
}

plot_experiment <- function(experiment) {
  plots <- map(c("expression_simulated", "expression", "true_counts", "counts"), function(expression_name) {
    expr <- experiment[[expression_name]]
    source(paste0(find.package("dyngen"), "/ext_data/dimred_wrappers.R"))
    
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

plot_normalisation <- function(experiment, normalisation) {
  source(paste0(find.package("dyngen"), "/ext_data/dimred_wrappers.R"))
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
