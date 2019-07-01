#' @importFrom tidygraph tbl_graph activate
#' @importFrom ggraph circle ggraph geom_edge_loop geom_edge_fan geom_node_circle geom_node_text geom_node_point theme_graph scale_edge_width_continuous
#' @importFrom ggplot2 ggplot scale_fill_manual coord_equal scale_colour_manual scale_size_manual coord_equal labs geom_path theme_bw aes facet_wrap geom_line
#' @importFrom viridis scale_color_viridis
#' @importFrom grid arrow unit
NULL

effect_colour <- function(effect) {
  index <- match(effect, c(1, -1, -2, NA, 0))
  c(
    "#3793d6",     # activating => blue
    "#d63737",     # repressing => red
    "#00000000",   # not real => invisible
    "darkgray",    # not decided yet
    "#7cd637"      # ??? :P => green
  )[index]
}

#' Visualise the backbone of a model
#' 
#' @param model A dyngen initial model created with [initialise_model()].
#' 
#' @importFrom igraph layout.graphopt V E
#' @export
plot_backbone <- function(model) {
  node_legend <- model$backbone$module_info %>% select(module_id, color) %>% deframe()
  
  nodes <- model$backbone$module_info %>% rename(name = module_id)
  edges <- model$backbone$module_network %>% arrange(from == to)
  
  gr <- tbl_graph(nodes = nodes, edges = edges)
  layout <- 
    gr %>% 
    igraph::layout.graphopt(charge = .01, niter = 10000) %>% 
    dynutils::scale_minmax() %>% 
    magrittr::set_rownames(nodes$name) %>% 
    magrittr::set_colnames(c("x", "y")) %>% 
    as.data.frame() 
  
  r <- .03
  cap <- circle(8, "mm")
  str <- .2
  arrow <- grid::arrow(type = "closed", angle = ifelse(igraph::E(gr)$effect == 1, 30, 89), length = grid::unit(3, "mm"))
  
  ggraph(gr, layout = "manual", node.positions = layout) +
    geom_edge_loop(aes(width = strength, strength = str), arrow = arrow, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(width = strength), arrow = arrow, start_cap = cap, end_cap = cap) +
    geom_node_circle(aes(r = r, fill = name)) +
    geom_node_text(aes(label = name)) +
    theme_graph(base_family = 'Helvetica') +
    scale_fill_manual(values = node_legend) +
    scale_edge_width_continuous(trans = "log10", range = c(.5, 3)) +
    coord_equal()
}

#' Visualise the feature network of a model
#' 
#' @param model A dyngen intermediary model for which the feature network has been generated with [generate_feature_network()].
#' @param show_tfs Whether or not to show the transcription factors.
#' @param show_targets Whether or not to show the targets.
#' @param show_hks Whether or not to show the housekeeping genes.
#' 
#' @importFrom igraph layout_with_fr V E
#' @export
plot_feature_network <- function(
  model,
  show_tfs = TRUE,
  show_targets = TRUE,
  show_hks = FALSE
) {
  # get feature info
  feature_info <- 
    model$feature_info %>%
    mutate(
      color_by = case_when(
        !is.na(module_id) ~ module_id,
        is_tf ~ "TF",
        !is_hk ~ "Target",
        is_hk ~ "HK"
      )
    )
  
  color_legend <- c(
    "TF" = "black",
    "Target" = "darkgray",
    "HK" = "lightgray",
    model$backbone$module_info %>% select(module_id, color) %>% deframe()
  )
  
  # remove unwanted features and convert NA into "NA"
  feature_info <- 
    feature_info %>% 
    filter(
      show_tfs & is_tf |
      show_targets & !is_tf & !is_hk |
      show_hks & is_hk
    ) %>% 
    mutate(
      color_by = ifelse(is.na(color_by), "NA", color_by)
    )
  
  # filter feature network
  feature_network <- 
    model$feature_network %>% 
    filter(from %in% feature_info$feature_id & to %in% feature_info$feature_id) %>% 
    arrange(from == to)
  
  # add extra edges invisible between regulators from the same module
  feature_network <- 
    bind_rows(
      feature_network,
      feature_info %>%
        filter(is_tf) %>% 
        select(module_id, feature_id) %>%
        group_by(module_id) %>%
        do({
          crossing(from = .$feature_id, to = .$feature_id) %>%
            mutate(effect = -2)
        }) %>% 
        ungroup() %>% 
        filter(from < to)
    )
  
  gr <- tbl_graph(nodes = feature_info, edges = feature_network)
  layout <- igraph::layout_with_fr(gr) %>% 
    dynutils::scale_minmax() %>%
    magrittr::set_rownames(feature_info$feature_id) %>%
    magrittr::set_colnames(c("x", "y")) %>%
    as.data.frame()
  gr <- gr %>% activate(edges) %>% filter(is.na(effect) | effect != -2)
  
  cap <- circle(2.5, "mm")
  arrow <- grid::arrow(
    type = "closed",
    length = grid::unit(3, "mm"),
    angle = case_when(
      igraph::E(gr)$effect == 1 ~ 30,
      igraph::E(gr)$effect == -1 ~ 89,
      TRUE ~ 0
    )
  )
  
  ggraph(gr, layout = "manual", node.positions = layout) +
    geom_edge_fan(
      arrow = arrow, 
      start_cap = cap, 
      end_cap = cap
    ) +
    geom_node_point(aes(colour = color_by, size = as.character(is_tf))) +
    theme_graph(base_family = 'Helvetica') +
    scale_colour_manual(values = color_legend) +
    scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3)) +
    coord_equal() +
    labs(size = "is TF", color = "Module group")
}

#' Visualise the simulations using the dimred
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_cells()].
#' @param mapping Which components to plot.
#' 
#' @export
plot_simulations <- function(model, mapping = aes(comp_1, comp_2)) {
  plot_df <- 
    bind_cols(
      model$simulations$meta,
      model$simulations$dimred %>% as.data.frame
    )
  
  ggplot(plot_df %>% filter(sim_time >= 0), mapping) +
    geom_path(aes(colour = sim_time, group = simulation_i)) +
    viridis::scale_color_viridis() +
    theme_bw()
}

#' Visualise the simulations using the dimred
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_cells()].
#' @param detailed Whether or not to colour according to each separate sub-edge in the gold standard.
#' @param mapping Which components to plot.
#' 
#' @export
plot_gold_simulations <- function(model, detailed = FALSE, mapping = aes(comp_1, comp_2)) {
  plot_df <- 
    bind_cols(
      model$gold_standard$meta,
      model$gold_standard$dimred %>% as.data.frame
    ) %>% filter(!burn)
  
  if (!detailed && model %has_names% "simulations" && model$simulations %has_names% "dimred") {
    plot_df <- plot_df %>% 
      bind_rows(
        bind_cols(
          model$simulations$meta,
          model$simulations$dimred %>% as.data.frame
        ) %>% filter(sim_time >= 0)
      )
  }
  
  if (detailed) {
    plot_df <- plot_df %>% mutate(edge = paste0(from_, "_", to_))
  } else {
    plot_df <- plot_df %>% mutate(edge = paste0(from, "_", to))
  }
  
  ggplot(mapping = mapping) +
    geom_path(aes(group = simulation_i), plot_df %>% filter(simulation_i > 0), colour = "darkgray") +
    geom_path(aes(colour = edge, group = paste0(from_, "_", to_)), plot_df %>% filter(simulation_i == 0), size = 2) +
    theme_bw()
}

#' Visualise the mapping of the simulations to the gold standard
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_cells()].
#' @param selected_simulations Which simulation indices to visualise.
#' @param do_facet Whether or not to facet according to simulation index.
#' @param mapping Which components to plot.
#' 
#' @export
plot_gold_mappings <- function(model, selected_simulations = NULL, do_facet = TRUE, mapping = aes(comp_1, comp_2)) {
  plot_df <- 
    bind_rows(
      bind_cols(
        model$simulations$meta,
        model$simulations$dimred %>% as.data.frame
      ) %>% filter(sim_time >= 0),
      bind_cols(
        model$gold_standard$meta,
        model$gold_standard$dimred %>% as.data.frame
      ) %>% filter(!burn)
    )
  
  if (!is.null(selected_simulations)) {
    plot_df <- plot_df %>% filter(simulation_i %in% selected_simulations)
  }
  
  plot_df <- plot_df %>% mutate(edge = paste0(from, "->", to))
  
  g <- ggplot(mapping = mapping) +
    geom_path(aes(colour = edge), plot_df %>% filter(simulation_i == 0)) +
    geom_path(aes(colour = edge, group = simulation_i), plot_df %>% filter(simulation_i != 0)) +
    theme_bw() 
  
  if (do_facet) {
    g <- g +
    facet_wrap(~ simulation_i)
  }
  
  g
}


#' Visualise the expression of the gold standard over simulation time
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_gold_standard()].
#' @param what Which molecule types to visualise.
#' 
#' @export
plot_gold_expression <- function(model, what = c("w", "x", "y")) {
  edge_levels <- 
    model$gold_standard$mod_changes %>% 
    mutate(edge = paste0(from_, "->", to_)) %>% 
    pull(edge)
  
  molecules <- model$feature_info %>% filter(is_tf) %>% gather(mol, val, w, x, y) %>% pull(val)
  df <- bind_cols(
    model$gold_standard$meta,
    as.data.frame(as.matrix(model$gold_standard$counts))[,molecules]
  ) %>% 
    gather(molecule, value, one_of(molecules)) %>% 
    mutate(edge = factor(paste0(from_, "->", to_), levels = edge_levels)) %>% 
    left_join(model$feature_info %>% select(w, x, y, module_id) %>% gather(type, molecule, w, x, y), by = "molecule") %>% 
    group_by(module_id, sim_time, simulation_i, burn, from, to, from_, to_, time, edge, type) %>% 
    summarise(value = mean(value)) %>% 
    ungroup() %>% 
    filter(type %in% what)
  
  ggplot(df) +
    geom_line(aes(sim_time, value, colour = module_id, linetype = type, size = type)) +
    scale_size_manual(values = c(w = .5, x = 1, y = .5)) +
    facet_wrap(~edge) +
    theme_bw()
}


#' Visualise the expression of the simulations over simulation time
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_cells()].
#' @param simulation_i Which simulation to visualise.
#' @param what Which molecule types to visualise.
#' 
#' @export
plot_simulation_expression <- function(model, simulation_i = 1, what = c("w", "x", "y")) {
  edge_levels <-
    model$gold_standard$network %>%
    mutate(edge = paste0(from, "->", to)) %>%
    pull(edge)
  
  molecules <- model$feature_info %>% filter(is_tf) %>% gather(mol, val, w, x, y) %>% pull(val)
  df <- bind_cols(
    model$simulations$meta,
    as.data.frame(as.matrix(model$simulations$counts)[,molecules])
  ) %>% 
    filter(simulation_i == !!simulation_i) %>% 
    gather(molecule, value, one_of(molecules)) %>% 
    mutate(edge = factor(paste0(from, "->", to), levels = edge_levels)) %>%
    left_join(model$feature_info %>% select(w, x, y, module_id) %>% gather(type, molecule, w, x, y), by = "molecule") %>% 
    group_by(module_id, sim_time, simulation_i, from, to, time, edge, type) %>% 
    summarise(value = mean(value)) %>% 
    ungroup() %>% 
    mutate(module_group = gsub("[0-9]*$", "", module_id)) %>% 
    filter(type %in% what)
  
  ggplot(df) +
    geom_line(aes(sim_time, value, linetype = type, colour = module_id, size = type)) +
    scale_size_manual(values = c(w = .5, x = 1, y = .5)) +
    facet_wrap(~module_group, ncol = 1) +
    theme_bw()
}
