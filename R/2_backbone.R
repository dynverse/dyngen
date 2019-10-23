#' Backbone of the simulation model
#' 
#' A module is a group of genes which, to some extent, shows the same
#' expression behaviour. Several modules are connected together such that
#' one or more genes from one module will regulate the expression of
#' another module. By creating chains of modules, a dynamic behaviour in gene 
#' regulation can be created.
#' 
#' @param module_info A tibble containing meta information on the modules themselves.
#' 
#' * module_id (character): the name of the module
#' * basal (numeric): basal expression level of genes in this module, must be between \[0, 1\]
#' * burn (logical): whether or not outgoing edges of this module will 
#'   be active during the burn in phase
#' * independence (numeric): the independence factor between regulators of this module, must be between \[0, 1\]
#'   
#' @param module_network A tibble describing which modules regulate which other modules.
#' 
#' * from (character): the regulating module
#' * to (character): the target module
#' * effect (integer): `1L` if the regulating module upregulates 
#'   the target module, `-1L` if it downregulates
#' * strength (numeric): the strength of the interaction
#' * cooperativity (numeric): cooperativity factor, larger 1 if positive cooperativity,
#'   between 0 and 1 for negative cooperativity
#'   
#' @param expression_patterns A tibble describing the expected expression pattern
#' changes when a cell is simulated by dyngen. Each row represents one transition between
#' two cell states. 
#' 
#' * from (character): name of a cell state
#' * to (character): name of a cell state
#' * module_progression (character): differences in module expression between the two states. 
#'   Example: `"+4,-1|+9|-4"` means the expression of module 4 will go up at the same time as module 1 goes down; 
#'   afterwards module 9 expression will go up, and afterwards module 4 expression will go down again.
#' * start (logical): Whether or not this from cell state is the start of the trajectory
#' * burn (logical): Whether these cell states are part of the burn in phase. Cells will 
#'   not get sampled from these cell states.
#' * time (numeric): The duration of an transition.
#' 
#' @importFrom grDevices rainbow
#' @export
#' 
#' @seealso [list_backbones()] for a list of all backbone methods.
backbone <- function(
  module_info,
  module_network,
  expression_patterns
) {
  assert_that(
    is.data.frame(module_info),
    module_info %has_names% c("module_id", "basal", "burn", "independence"),
    is.character(module_info$module_id),
    is.numeric(module_info$basal),
    all(0 <= module_info$basal & module_info$basal <= 1),
    is.logical(module_info$burn),
    is.numeric(module_info$independence),
    all(0 <= module_info$independence & module_info$independence <= 1),
    
    is.data.frame(module_network),
    module_network %has_names% c("from", "to", "effect", "strength", "cooperativity"),
    is.character(module_network$from),
    is.character(module_network$to),
    is.integer(module_network$effect),
    all(module_network$effect == 1L | module_network$effect == -1L),
    is.numeric(module_network$strength),
    is.numeric(module_network$cooperativity),
    module_network$from %all_in% module_info$module_id,
    module_network$to %all_in% module_info$module_id,
    
    is.data.frame(expression_patterns),
    expression_patterns %has_names% c("from", "to", "module_progression", "start", "burn", "time"),
    is.character(expression_patterns$from),
    is.character(expression_patterns$to),
    is.character(expression_patterns$module_progression),
    is.logical(expression_patterns$start),
    is.logical(expression_patterns$burn),
    is.numeric(expression_patterns$time)
  )
  
  tmp_modules <-
    expression_patterns$module_progression %>%
    strsplit("[,|]") %>%
    unlist() %>%
    unique() %>% 
    gsub("[+-]", "", .)
  assert_that(tmp_modules %all_in% module_info$module_id)
  
  if (! module_info %has_name% "color") {
    if (all(grepl("^[a-zA-Z]*[0-9]*$", module_info$module_id))) {
      # module names are divided up into groups
      module_info <-
        module_info %>% 
        mutate(
          group = gsub("[0-9]*$", "", module_id)
        )
      group_names <- unique(module_info$group)
      group_colours <- grDevices::rainbow(length(group_names))
      
      module_info <- module_info %>% 
        group_by(group) %>% 
        mutate(
          group_colour = group_colours[match(group, group_names)],
          color = colour_brighten(group_colour, rev(seq(1, .4, length.out = n())))
        ) %>% 
        ungroup() %>% 
        select(-group, -group_colour)
    } else {
      module_info <-
        module_info %>% 
        mutate(color = grDevices::rainbow(n()))
    }
  }
  
  # add burn modules to burn transition
  extra_modules <- module_info %>% filter(basal > 0) %>% pull(module_id)
  if (length(extra_modules) > 0) {
    expression_patterns$module_progression[[1]] <- 
      c(
        strsplit(expression_patterns$module_progression[[1]], ",")[[1]],
        paste0("+", extra_modules)
      ) %>% 
      unique() %>% 
      paste(collapse = ",")
  }
  
  # make sure time is not an integer
  expression_patterns$time <- as.numeric(expression_patterns$time)
  
  lst(
    module_info, 
    module_network, 
    expression_patterns
  ) %>% 
    add_class("dyngen::backbone")
}


#' List of all predefined backbone models
#' 
#' A module is a group of genes which, to some extent, shows the same
#' expression behaviour. Several modules are connected together such that
#' one or more genes from one module will regulate the expression of
#' another module. By creating chains of modules, a dynamic behaviour in gene 
#' regulation can be created.
#' 
#' @export
#' @rdname backbone_models
#' 
#' @seealso [backbone()] for more information on the data structures that define the backbone.
list_backbones <- function() {
  list(
    bifurcating = backbone_bifurcating,
    bifurcating_converging = backbone_bifurcating_converging,
    bifurcating_cycle = backbone_bifurcating_cycle,
    bifurcating_loop = backbone_bifurcating_loop,
    binary_tree = backbone_binary_tree,
    branching = backbone_branching,
    consecutive_bifurcating = backbone_consecutive_bifurcating,
    converging = backbone_converging,
    cycle = backbone_cycle,
    disconnected = backbone_disconnected,
    linear = backbone_linear,
    trifurcating = backbone_trifurcating
  )
}