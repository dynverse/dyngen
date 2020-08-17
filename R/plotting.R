#' @importFrom tidygraph tbl_graph activate
#' @importFrom ggplot2 ggsave labs aes coord_equal scale_colour_manual scale_size_manual facet_wrap geom_line geom_text ggplot theme_bw geom_path geom_point geom_step
#' @importFrom ggraph circle geom_edge_fan geom_edge_loop geom_node_circle geom_node_text ggraph scale_edge_width_continuous theme_graph geom_node_label geom_node_point
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

#' Visualise the backbone state network of a model
#' 
#' @param model A dyngen initial model created with [initialise_model()].
#' @param detailed Whether or not to also plot the substates of transitions.
#' 
#' @export
#' 
#' @examples
#' data("example_model")
#' plot_backbone_statenet(example_model)
plot_backbone_statenet <- function(model, detailed = FALSE) {
  # satisfy r cmd check 
  from_ <- to_ <- from <- to <- mod_diff <- name <- time <- module_progression <- 
    from_cap <- to_cap <- main <- NULL
  
  edges <- model$backbone$expression_patterns
  
  large_cap <- 4
  small_cap <- 1
  if (detailed) {
    edges <- .generate_gold_standard_mod_changes(edges) %>% 
      rename(from = from_, to = to_, from_ = from, to_ = to, module_progression = mod_diff) %>% 
      mutate(
        from_cap = ifelse(from == from_, large_cap, small_cap),
        to_cap = ifelse(to == to_, large_cap, small_cap)
      )
    nodes <- tibble(
      name = unique(c(edges$from, edges$to)),
      main = name %in% c(edges$from_, edges$to_)
    )
  } else {
    nodes <- tibble(
      name = unique(c(edges$from, edges$to)),
      main = TRUE
    )
    edges <- edges %>% mutate(
      from_cap = large_cap,
      to_cap = large_cap
    )
  }
  
  gr <- tbl_graph(edges = edges %>% rename(weight = time), nodes = nodes)
  
  r <- .05
  arrow <- grid::arrow(type = "closed", length = grid::unit(3, "mm"))
  
  ggraph(gr, layout = "igraph", algorithm = "kk") +
    geom_edge_fan(
      aes(
        label = module_progression,
        start_cap = circle(from_cap, "mm"),
        end_cap = circle(to_cap, "mm")
      ), 
      arrow = arrow, 
      colour = "gray"
    ) +
    geom_node_point(data = function(df) df %>% filter(!main)) +
    geom_node_label(aes(label = name), function(df) df %>% filter(main)) +
    theme_graph(base_family = 'Helvetica') +
    coord_equal()
}

#' Visualise the backbone of a model
#' 
#' @param model A dyngen initial model created with [initialise_model()].
#' 
#' @importFrom igraph layout.graphopt V E
#' @export
#' 
#' @examples
#' data("example_model")
#' plot_backbone_modulenet(example_model)
plot_backbone_modulenet <- function(model) {
  # satisfy r cmd check
  module_id <- color <- from <- to <- strength <- effect <- name <- NULL
  
  node_legend <- model$backbone$module_info %>% select(module_id, color) %>% deframe()
  
  nodes <- model$backbone$module_info %>% rename(name = module_id)
  edges <- model$backbone$module_network %>% arrange(from == to)
  
  gr <- tbl_graph(nodes = nodes, edges = edges)
  layout <- 
    gr %>% 
    igraph::layout.graphopt(charge = .01, niter = 10000) %>% 
    dynutils::scale_minmax() %>% 
    as.data.frame() 
  rownames(layout) <- nodes$name
  colnames(layout) <- c("x", "y")
  
  r <- .03
  cap <- circle(4, "mm")
  str <- .2
  arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(3, "mm"))
  arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(3, "mm"))
  
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
    geom_edge_loop(aes(width = strength, strength = str, filter = effect >= 0), arrow = arrow_up, start_cap = cap, end_cap = cap) +
    geom_edge_loop(aes(width = strength, strength = str, filter = effect < 0), arrow = arrow_down, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(width = strength, filter = effect >= 0), arrow = arrow_up, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(width = strength, filter = effect < 0), arrow = arrow_down, start_cap = cap, end_cap = cap) +
    geom_node_circle(aes(r = r, colour = name), fill = "white") +
    geom_node_text(aes(label = name)) +
    theme_graph(base_family = 'Helvetica') +
    scale_colour_manual(values = node_legend) +
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
#' 
#' @examples
#' data("example_model")
#' plot_feature_network(example_model)
plot_feature_network <- function(
  model,
  show_tfs = TRUE,
  show_targets = TRUE,
  show_hks = FALSE
) {
  # satisfy r cmd check
  module_id <- color <- is_tf <- is_hk <- color_by <- from <- to <- 
    feature_id <- `.` <- edges <- effect <- NULL
  
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
    as.data.frame()
  rownames(layout) <- feature_info$feature_id
  colnames(layout) <- c("x", "y")
  
  gr <- gr %>% activate(edges) %>% filter(is.na(effect) | effect != -2)
  
  cap <- circle(2.5, "mm")
  str <- .2
  
  arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(3, "mm"))
  arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(3, "mm"))
  
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
    geom_edge_loop(aes(strength = str, filter = !is.na(effect) & effect >= 0 & from == to), arrow = arrow_up, start_cap = cap, end_cap = cap) +
    geom_edge_loop(aes(strength = str, filter = !is.na(effect) & effect < 0 & from == to), arrow = arrow_down, start_cap = cap, end_cap = cap) +
    geom_edge_loop(aes(strength = str, filter = is.na(effect))) +
    geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0 & from != to), arrow = arrow_up, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(filter = !is.na(effect) & effect < 0), arrow = arrow_down, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(filter = is.na(effect))) +
    geom_node_point(aes(colour = color_by, size = as.character(is_tf))) +
    theme_graph(base_family = "Helvetica") +
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
#' 
#' @examples
#' data("example_model")
#' plot_simulations(example_model)
plot_simulations <- function(model, mapping = aes(comp_1, comp_2)) {
  # satisfy r cmd check
  comp_1 <- comp_2 <- sim_time <- simulation_i <- NULL
  
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
#' @param highlight Which simulation to highlight. If highlight == 0 then the gold simulation will be highlighted.
#' 
#' @export
#' 
#' @examples
#' data("example_model")
#' plot_gold_simulations(example_model)
plot_gold_simulations <- function(model, detailed = FALSE, mapping = aes(comp_1, comp_2), highlight = 0) {
  # satisfy r cmd check
  comp_1 <- comp_2 <- burn <- sim_time <- from_ <- to_ <- from <- to <- simulation_i <- edge <- NULL
  
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
    geom_path(aes(group = simulation_i), plot_df %>% filter(simulation_i != highlight), colour = "darkgray") +
    geom_path(aes(colour = edge, group = paste0(from_, "_", to_)), plot_df %>% filter(simulation_i == highlight), size = 2) +
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
#' 
#' @examples
#' data("example_model")
#' plot_gold_mappings(example_model)
plot_gold_mappings <- function(model, selected_simulations = NULL, do_facet = TRUE, mapping = aes(comp_1, comp_2)) {
  # satisfy r cmd check
  comp_1 <- comp_2 <- sim_time <- burn <- simulation_i <- from <- to <- edge <- NULL
  
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
#' @param label_changing Whether or not to add a label next to changing molecules.
#' 
#' @export
#' 
#' @examples
#' data("example_model")
#' plot_gold_expression(example_model, what = "mol_mrna", label_changing = FALSE)
plot_gold_expression <- function(
  model, 
  what = c("mol_premrna", "mol_mrna", "mol_protein"),
  label_changing = TRUE
) {
  # satisfy r cmd check
  from_ <- to_ <- edge <- is_tf <- mol <- val <- molecule <- value <- module_id <- 
    type <- sim_time <- simulation_i <- burn <- from <- to <- time <- color <- NULL
  
  assert_that(what %all_in% c("mol_premrna", "mol_mrna", "mol_protein"))
  
  edge_levels <- 
    model$gold_standard$mod_changes %>% 
    mutate(edge = paste0(from_, "->", to_)) %>% 
    pull(edge)
  
  molecules <- model$feature_info %>% filter(is_tf) %>% gather(mol, val, !!!what) %>% pull(val)
  df <- bind_cols(
    model$gold_standard$meta,
    as.data.frame(as.matrix(model$gold_standard$counts))[,molecules]
  ) %>% 
    gather(molecule, value, one_of(molecules)) %>% 
    mutate(edge = factor(paste0(from_, "->", to_), levels = edge_levels)) %>% 
    left_join(model$feature_info %>% select(module_id, !!!what) %>% gather(type, molecule, !!!what), by = "molecule") %>% 
    group_by(module_id, sim_time, simulation_i, burn, from, to, from_, to_, time, edge, type) %>% 
    summarise(value = mean(value)) %>% 
    ungroup() %>% 
    filter(type %in% what)
  
  g <- ggplot(df, aes(sim_time, value, colour = module_id)) +
    geom_line(aes(linetype = type, size = type)) +
    scale_size_manual(values = c(mol_premrna = .5, mol_mrna = 1, mol_protein = .5)) +
    scale_colour_manual(values = model$backbone$module_info %>% select(module_id, color) %>% deframe) +
    facet_wrap(~edge) +
    theme_bw()
  
  if (label_changing) {
    g <- g +
      geom_text(
        aes(label = paste0(module_id, "_", type)), 
        df %>% group_by(edge, module_id, type) %>% filter(sim_time == max(sim_time) & any(diff(value) > 0.01)) %>% ungroup,
        hjust = 1, 
        vjust = 0,
        nudge_y = .15
      )
  }
  
  g
}


#' Visualise the expression of the simulations over simulation time
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_cells()].
#' @param simulation_i Which simulation to visualise.
#' @param what Which molecule types to visualise.
#' @param facet What to facet on.
#' @param label_nonzero Plot labels for non-zero molecules.
#' 
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats approx
#' 
#' @export
#' 
#' @examples
#' data("example_model")
#' plot_simulation_expression(example_model)
plot_simulation_expression <- function(
  model, 
  simulation_i = 1:4,
  what = c("mol_premrna", "mol_mrna", "mol_protein"),
  facet = c("simulation", "module_group", "module_id", "none"),
  label_nonzero = FALSE
) {
  # satisfy r check
  is_tf <- mol <- val <- molecule <- value <- module_id <- type <- sim_time <- 
    color <- module_group <- x <- y <- `.` <- NULL
  
  assert_that(what %all_in% c("mol_premrna", "mol_mrna", "mol_protein"))
  facet <- match.arg(facet)
  
  molecules <- model$feature_info %>% filter(is_tf) %>% gather(mol, val, !!!what) %>% pull(val)
  df <- bind_cols(
    model$simulations$meta,
    as.data.frame(as.matrix(model$simulations$counts)[,molecules])
  ) %>% 
    filter(simulation_i %in% !!simulation_i) %>% 
    gather(molecule, value, one_of(molecules)) %>% 
    left_join(
      model$feature_info %>%
        select(!!!what, module_id) %>% 
        gather(type, molecule, !!!what) %>% 
        mutate(type = factor(type, levels = what)), 
      by = "molecule"
    ) %>% 
    group_by(module_id, sim_time, simulation_i, type) %>% 
    summarise(value = mean(value)) %>% 
    ungroup() %>% 
    mutate(module_group = gsub("[0-9]*$", "", module_id)) %>% 
    filter(type %in% what)
  
  g <- ggplot(df, aes(sim_time, value)) +
    geom_step(aes(linetype = type, size = type, colour = module_id)) +
    scale_size_manual(values = c(mol_premrna = .5, mol_mrna = 1, mol_protein = .5)) +
    scale_colour_manual(values = model$backbone$module_info %>% select(module_id, color) %>% deframe) +
    theme_bw()
  
  if (label_nonzero) {
    pts <- seq(0, max(model$simulations$meta$sim_time), by = 5)
    df_labels <- 
      df %>% 
      group_by(module_id, type, module_group) %>% do({
        df2 <- .
        approx(x = df2$sim_time, y = df2$value, xout = pts) %>%
          as_tibble() %>% 
          rename(sim_time = x, value = y)
      }) %>% 
      ungroup() %>% 
      filter(value > 0)
    g <- g +
      geom_point(data = df_labels) +
      ggrepel::geom_text_repel(
        aes(label = paste0(module_id, "_", type)), 
        df_labels
      )
  }
  
  if (facet == "simulation") {
    g <- g + facet_wrap(~simulation_i, ncol = 1)
  } else if (facet == "module_group") {
    g <- g + facet_wrap(~module_group, ncol = 1)
  } else if (facet == "module_id") {
    g <- g + facet_wrap(~module_id, ncol = 1)
  }
  
  g
}
