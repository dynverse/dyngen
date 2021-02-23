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
#' * hill (numeric): hill coefficient, larger than 1 for positive cooperativity,
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
#' @return A dyngen backbone.
#' 
#' @seealso [dyngen] on how to run a dyngen simulation
#' 
#' @examples
#' library(tibble)
#' backbone <- backbone(
#'   module_info = tribble(
#'     ~module_id, ~basal, ~burn, ~independence,
#'     "M1",       1,      TRUE,  1,
#'     "M2",       0,      FALSE, 1,
#'     "M3",       0,      FALSE, 1
#'   ),
#'   module_network = tribble(
#'     ~from, ~to,  ~effect, ~strength, ~hill,
#'     "M1",  "M2", 1L,      1,         2,
#'     "M2",  "M3", 1L,      1,         2
#'   ), 
#'   expression_patterns = tribble(
#'     ~from, ~to,  ~module_progression, ~start, ~burn, ~time,
#'     "s0",  "s1", "+M1",               TRUE,   TRUE,  30,
#'     "s1",  "s2", "+M2,+M3",           FALSE,  FALSE, 80
#'   )
#' )
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
    module_network %has_names% c("from", "to", "effect", "strength", "hill"),
    is.character(module_network$from),
    is.character(module_network$to),
    is.integer(module_network$effect),
    all(module_network$effect == 1L | module_network$effect == -1L),
    is.numeric(module_network$strength),
    is.numeric(module_network$hill),
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
    unique()
  assert_that(
    gsub("[+-]", "", tmp_modules) %all_in% module_info$module_id,
    msg = "Not all modules listed under `expression_Patterns$module_progression` are in `module_info$module_id`."
  )
  
  if (! module_info %has_name% "color") {
    if (all(grepl("^[a-zA-Z]*[0-9]*$", module_info$module_id))) {
      # module names are divided up into groups
      module_info <-
        module_info %>% 
        mutate(
          group = gsub("[0-9]*$", "", .data$module_id)
        )
      group_names <- unique(module_info$group)
      group_colours <- grDevices::rainbow(length(group_names))
      
      module_info <- module_info %>% 
        group_by(.data$group) %>% 
        mutate(
          group_colour = group_colours[match(.data$group, group_names)],
          color = colour_brighten(.data$group_colour, rev(seq(1, .4, length.out = n())))
        ) %>% 
        ungroup() %>% 
        select(-.data$group, -.data$group_colour)
    } else {
      module_info <-
        module_info %>% 
        mutate(color = grDevices::rainbow(n()))
    }
  }
  
  # add burn modules to burn transition
  extra_modules <- module_info %>% filter(.data$basal > 0) %>% pull(.data$module_id)
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
#' @return A list of all the available backbone generators.
#' 
#' @seealso [dyngen] on how to run a dyngen simulation
#' 
#' @examples
#' names(list_backbones())
#' 
#' bb <- backbone_bifurcating()
#' bb <- backbone_bifurcating_converging()
#' bb <- backbone_bifurcating_cycle()
#' bb <- backbone_bifurcating_loop()
#' bb <- backbone_binary_tree()
#' bb <- backbone_branching()
#' bb <- backbone_consecutive_bifurcating()
#' bb <- backbone_converging()
#' bb <- backbone_cycle()
#' bb <- backbone_cycle_simple()
#' bb <- backbone_disconnected()
#' bb <- backbone_linear()
#' bb <- backbone_linear_simple()
#' bb <- backbone_trifurcating()
#' 
#' model <- initialise_model(
#'   backbone = bb
#' )
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
    cycle_simple = backbone_cycle_simple,
    disconnected = backbone_disconnected,
    linear = backbone_linear,
    linear_simple = backbone_linear_simple,
    trifurcating = backbone_trifurcating
  )
}