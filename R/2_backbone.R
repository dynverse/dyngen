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
#' * a0 (numeric): basal expression level of genes in this module
#' * burn (logical): whether or not outgoing edges of this module will 
#'   be active during the burn in phase
#'   
#' @param module_network A tibble describing which modules regulate which other modules.
#' 
#' * from (character): the regulating module
#' * to (character): the target module
#' * effect (integer): `1` if the regulating module upregulates 
#'   the target module, `-1` if it downregulates
#' * strength (numeric): the strength of the interaction
#' * cooperativity (numeric): cooperativity factor, larger 1 if positive cooperativity,
#'   between 0 and 1 for negative cooperativity
#' @param expression_patterns A tibble describing the expected expression pattern
#' changes when a cell is simulated by dyngen.
#' 
#' * from (character): name of a cell state
#' * to (character): name of a cell state
#' * module_progression (character): differences in module expression between the two states. 
#'   Example: `"+4,-1|+9|-4"` means the expression of module 4 will go up at the same time as module 1 goes down; 
#'   afterwards module 9 expression will go up, and afterwards module 4 expression will go down again.
#' * start (logical): Whether or not this from cell state is the start of the trajectory
#' * burn (logical): Whether these cell states are part of the burn in phase. Cells will 
#'   not get sampled from these cell states.
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
    module_network$from %all_in% module_info$module_id,
    module_network$to %all_in% module_info$module_id
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
  expression_patterns$module_progression[[1]] <- 
    c(
      strsplit(expression_patterns$module_progression[[1]], ",")[[1]],
      paste0("+", module_info %>% filter(a0 > 0) %>% pull(module_id))
    ) %>% 
    unique() %>% 
    paste(collapse = ",")
  
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