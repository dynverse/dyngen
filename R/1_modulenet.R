#' Module network models for creating a gene regulatory network
#' 
#' A module is a group of genes which, to some extent, shows the same
#' expression behaviour. Several modules are connected together such that
#' one or more genes from one module will regulate the expression of
#' another module. By creating chains of modules, a dynamic behaviour in gene 
#' regulation can be created.
#' 
#' A modulenet model contains three tibbles, namely `module_info`,
#' `module_network`, and `expression_patterns`.
#' 
#' `module_info` contains meta information on the modules themselves.
#' 
#' * module_id (character): the name of the module
#' * a0 (numeric): basal expression level of genes in this module
#' * burn (logical): whether or not outgoing edges of this module will 
#'   be active during the burn in phase
#'     
#' `module_network` describes which modules regulate which other modules.
#' 
#' * from (character): the regulating module
#' * to (character): the target module
#' * effect (integer): `1` if the regulating module upregulates 
#'   the target module, `-1` if it downregulates
#' * strength (numeric): the strength of the interaction
#' * cooperativity (numeric): cooperativity factor, larger 1 if positive cooperativity,
#'   between 0 and 1 for negative cooperativity
#' 
#' `expression_patterns` describes the expected expression pattern
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
#' @export
modulenet <- function(
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
    module_info <-
      module_info %>% 
      mutate(color = grDevices::rainbow(n()))
  }
  
  lst(
    module_info, 
    module_network, 
    expression_patterns
  ) %>% 
    add_class("dyngen::modulenet")
}

#' @export
#' @rdname modulenet
list_modulenet_generators <- function() {
  list(
    bifurcating = modulenet_bifurcating,
    bifurcating_converging = modulenet_bifurcating_converging,
    bifurcating_cycle = modulenet_bifurcating_cycle,
    bifurcating_loop = modulenet_bifurcating_loop,
    binary_tree = modulenet_binary_tree,
    consecutive_bifurcating = modulenet_consecutive_bifurcating,
    converging = modulenet_converging,
    cycle = modulenet_cycle,
    linear = modulenet_linear,
    linear_long = modulenet_linear_long,
    trifurcating = modulenet_trifurcating
  )
}

#' @export
#' @rdname modulenet
modulenet_bifurcating <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 1, TRUE,
    "M8", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M8", 1, 1, 2,
    "M8", "M2", 1, 0.2, 2,
    "M2", "M3", 1, 1, 2,
    "M2", "M4", 1, 1, 2,
    "M3", "M3", 1, 10, 2,
    "M4", "M4", 1, 10, 2,
    "M4", "M3", -1, 10, 2,
    "M3", "M4", -1, 10, 2,
    "M4", "M1", -1, 10, 2,
    "M3", "M1", -1, 10, 2,
    "M3", "M6", 1, 1, 2,
    "M3", "M7", -1, 10, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1,+M7", TRUE, TRUE,
    "S1", "S2", "+M8|+M2", FALSE, FALSE,
    "S2", "S3", "+M3|-M1|-M8,+M6|-M7,-M2", FALSE, FALSE,
    "S2", "S4", "+M4|-M1|-M8|-M2", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_bifurcating_converging <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 0, FALSE,
    "M8", 0, FALSE,
    "M9", 0, FALSE,
    "M10", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M4", 1, 10, 2,
    "M3", "M5", 1, 1, 2,
    "M5", "M5", 1, 10, 2,
    "M4", "M5", -1, 100, 2,
    "M5", "M4", -1, 100, 2,
    "M4", "M6", 1, 0.25, 5,
    "M5", "M7", 1, 0.25, 5,
    "M6", "M8", 1, 1, 2,
    "M7", "M8", 1, 1, 2,
    "M8", "M8", 1, 10, 2,
    "M8", "M4", -1, 10, 2,
    "M8", "M5", -1, 10, 2,
    "M8", "M9", 1, 1, 2,
    "M9", "M1", -1, 10, 2,
    "M9", "M10", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2|+M3", FALSE, FALSE,
    "S2", "S3", "+M5|+M7|+M8", FALSE, FALSE,
    "S2", "S4", "+M4|+M6|+M8", FALSE, FALSE,
    "S3", "S5", "+M9,-M5|-M7,+M10", FALSE, FALSE,
    "S4", "S5", "+M9,-M4|-M6,+M10", FALSE, FALSE,
    "S5", "S6", "-M1,-M8|-M9,-M2|-M10,-M3", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_bifurcating_cycle <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 0, FALSE,
    "M8", 0, FALSE,
    "M9", 0, FALSE,
    "M10", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M4", 1, 10, 2,
    "M3", "M5", 1, 1, 2,
    "M5", "M5", 1, 10, 2,
    "M4", "M5", -1, 100, 2,
    "M5", "M4", -1, 100, 2,
    "M4", "M6", 1, 0.1, 5,
    "M5", "M7", 1, 0.1, 5,
    "M6", "M8", 1, 1, 2,
    "M7", "M8", 1, 1, 2,
    "M8", "M8", 1, 10, 2,
    "M8", "M4", -1, 10, 2,
    "M8", "M5", -1, 10, 2,
    "M8", "M9", 1, 1, 2,
    "M9", "M1", -1, 10, 2,
    "M9", "M10", 1, 1, 2,
    "M10", "M8", -1, 100, 2,
    "M10", "M4", -1, 100, 2,
    "M10", "M5", -1, 100, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2|+M3", FALSE, FALSE,
    "S2", "S3", "+M4|+M6|+M8", FALSE, FALSE,
    "S2", "S4", "+M5|+M7|+M8", FALSE, FALSE,
    "S3", "S5", "-M4,+M9,-M6|+M10,-M1", FALSE, FALSE,
    "S4", "S5", "-M5,+M9,-M7|+M10,-M1", FALSE, FALSE,
    "S5", "S1", "-M2,-M8|-M3,-M9|-M10|+M1", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_bifurcating_loop <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 0, FALSE,
    "M8", 0, FALSE,
    "M9", 0, FALSE,
    "M10", 0, FALSE,
    "M11", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M4", 1, 10, 2,
    "M3", "M5", 1, 1, 2,
    "M5", "M5", 1, 10, 2,
    "M4", "M5", -1, 100, 2,
    "M5", "M4", -1, 100, 2,
    "M4", "M6", 1, 0.1, 5,
    "M6", "M8", 1, 1, 2,
    "M7", "M8", 1, 1, 2,
    "M8", "M8", 1, 10, 2,
    "M8", "M4", -1, 10, 2,
    "M8", "M9", 1, 1, 2,
    "M9", "M1", -1, 10, 2,
    "M9", "M10", 1, 1, 2,
    "M10", "M8", -1, 100, 2,
    "M10", "M4", -1, 100, 2,
    "M5", "M11", 1, 10, 2,
    "M5", "M8", -1, 100, 2,
    "M10", "M5", -1, 100, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2|+M3", FALSE, FALSE,
    "S2", "S4", "+M4|+M6|+M8|+M9,-M4", FALSE, FALSE,
    "S4", "S1", "-M1,+M10,-M6|-M2,-M8|-M3,-M9|-M10|+M1", FALSE, FALSE,
    "S2", "S3", "+M5|+M11", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_binary_tree <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1A", 1, TRUE,
    "M1B", 0, FALSE,
    "M2A", 0, FALSE,
    "M2B", 0, FALSE,
    "M2C", 0, FALSE,
    "M3A", 0, FALSE,
    "M3B", 0, FALSE,
    "M3C", 0, FALSE,
    "M4A", 0, FALSE,
    "M5A", 0, FALSE,
    "M4B", 0, FALSE,
    "M5B", 0, FALSE,
    "M5C", 0, FALSE,
    "M6A", 0, FALSE,
    "M6B", 0, FALSE,
    "M7A", 0, FALSE,
    "M7B", 0, FALSE,
    "M8A", 0, FALSE,
    "M8B", 0, FALSE,
    "M8C", 0, FALSE,
    "M9A", 0, FALSE,
    "M9B", 0, FALSE,
    "M10A", 0, FALSE,
    "M10B", 0, FALSE,
    "M10C", 0, FALSE,
    "M10D", 0, FALSE,
    "M10E", 0, FALSE,
    "M11A", 0, FALSE,
    "M11B", 0, FALSE,
    "M11C", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1A", "M1B", 1, 0.1, 2,
    "M1B", "M2A", 1, 1, 2,
    "M2A", "M2A", 1, 10, 2,
    "M1B", "M3A", 1, 1, 2,
    "M3A", "M3A", 1, 10, 2,
    "M3A", "M2A", -1, 10, 2,
    "M2A", "M3A", -1, 10, 2,
    "M2A", "M2B", 1, 1, 2,
    "M2B", "M2C", 1, 1, 2,
    "M2C", "M4A", 1, 1, 2,
    "M2C", "M5A", 1, 1, 2,
    "M4A", "M5A", -1, 10, 2,
    "M5A", "M4A", -1, 10, 2,
    "M4A", "M4B", 1, 1, 2,
    "M5A", "M5B", 1, 1, 2,
    "M5B", "M5C", 1, 1, 2,
    "M5C", "M6A", 1, 1, 2,
    "M5C", "M7A", 1, 1, 2,
    "M6A", "M7A", -1, 10, 2,
    "M7A", "M6A", -1, 10, 2,
    "M6A", "M6B", 1, 1, 2,
    "M7A", "M7B", 1, 1, 2,
    "M3A", "M3B", 1, 1, 2,
    "M3B", "M3C", 1, 1, 2,
    "M3C", "M8A", 1, 1, 2,
    "M3C", "M9A", 1, 1, 2,
    "M8A", "M9A", -1, 10, 2,
    "M9A", "M8A", -1, 10, 2,
    "M8A", "M8B", 1, 1, 2,
    "M8B", "M8C", 1, 1, 2,
    "M9A", "M9B", 1, 1, 2,
    "M4B", "M10A", 1, 1, 2,
    "M4B", "M11A", 1, 1, 2,
    "M10A", "M11A", -1, 10, 2,
    "M11A", "M10A", -1, 10, 2,
    "M11A", "M11B", 1, 1, 2,
    "M11B", "M11C", 1, 1, 2,
    "M10A", "M10B", 1, 1, 2,
    "M10B", "M10C", 1, 1, 2,
    "M10C", "M10D", 1, 1, 2,
    "M10D", "M10E", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1A", TRUE, TRUE,
    "S1", "S2", "+M1B", FALSE, FALSE,
    "S2", "S3", "+M2A|+M2B|+M2C", FALSE, FALSE,
    "S2", "S4", "+M3A|+M3B|+M3C", FALSE, FALSE,
    "S3", "S5", "+M4A|+M4B", FALSE, FALSE,
    "S3", "S6", "+M5A|+M5B|+M5C", FALSE, FALSE,
    "S4", "S7", "+M8A|+M8B|+M8C", FALSE, FALSE,
    "S4", "S8", "+M9A|+M9B", FALSE, FALSE,
    "S6", "S9", "+M6A|+M6B", FALSE, FALSE,
    "S6", "S10", "+M7A|+M7B", FALSE, FALSE,
    "S5", "S11", "+M10A|+M10B|+M10C|+M10D|+M10E", FALSE, FALSE,
    "S5", "S12", "+M11A|+M11B|+M11C", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_consecutive_bifurcating <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 0, FALSE,
    "M8", 0, FALSE,
    "M9", 0, FALSE,
    "M10", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 0.1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M3", 1, 10, 2,
    "M2", "M4", 1, 1, 2,
    "M4", "M4", 1, 10, 2,
    "M4", "M3", -1, 10, 2,
    "M3", "M4", -1, 10, 2,
    "M3", "M5", 1, 0.1, 2,
    "M5", "M6", 1, 1, 2,
    "M5", "M7", 1, 1, 2,
    "M6", "M6", 1, 10, 2,
    "M7", "M7", 1, 10, 2,
    "M6", "M7", -1, 10, 2,
    "M7", "M6", -1, 10, 2,
    "M4", "M7", -1, 10, 2,
    "M4", "M6", -1, 10, 2,
    "M7", "M8", 1, 1, 2,
    "M6", "M9", 1, 1, 2,
    "M4", "M10", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2", FALSE, FALSE,
    "S2", "S3", "+M3|+M5", FALSE, FALSE,
    "S2", "S4", "+M4|+M10", FALSE, FALSE,
    "S3", "S5", "+M7|+M8", FALSE, FALSE,
    "S3", "S6", "+M6|+M9", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_converging <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, TRUE,
    "M3", 0, TRUE,
    "M4", 0, TRUE,
    "M5", 0, TRUE,
    "M6", 0, TRUE,
    "M7", 0, TRUE,
    "M8", 0, FALSE,
    "M9", 0, FALSE,
    "M10", 0, FALSE,
    "M11", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M4", 1, 10, 2,
    "M3", "M5", 1, 1, 2,
    "M5", "M5", 1, 10, 2,
    "M4", "M5", -1, 100, 2,
    "M5", "M4", -1, 100, 2,
    "M4", "M6", 1, 0.25, 5,
    "M5", "M7", 1, 0.25, 5,
    "M6", "M8", 1, 1, 2,
    "M7", "M8", 1, 1, 2,
    "M8", "M8", 1, 10, 2,
    "M8", "M4", -1, 10, 2,
    "M8", "M5", -1, 10, 2,
    "M8", "M9", 1, 1, 2,
    "M9", "M1", -1, 10, 2,
    "M9", "M10", 1, 1, 2,
    "M10", "M11", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2|+M3", FALSE, TRUE,
    "S2", "S3", "+M5|+M7", FALSE, TRUE,
    "S2", "S4", "+M4|+M6", FALSE, TRUE,
    "S3", "S5", "+M8|+M9,-M5|-M7,+M10", FALSE, FALSE,
    "S4", "S5", "+M8|+M9,-M4|-M6,+M10", FALSE, FALSE,
    "S5", "S6", "+M11,-M1,-M8|-M9,-M2|-M10,-M3", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_cycle <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 1, FALSE,
    "M3", 1, TRUE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, TRUE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", -1, 100, 2,
    "M2", "M3", -1, 100, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M5", 1, 1, 2,
    "M5", "M1", -1, 100, 2,
    "M1", "M6", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M3,+M1|+M6", TRUE, TRUE,
    "S1", "S2", "-M1,+M4|+M5,-M6", FALSE, FALSE,
    "S2", "S3", "-M3,+M2|-M4|-M5", FALSE, FALSE,
    "S3", "S1", "+M1,-M2|+M6|+M3", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_linear <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M5", 1, 1, 2,
    "M5", "M6", 1, 1, 2,
    "M6", "M7", 1, 1, 2,
    "M7", "M7", 1, 5, 2,
    "M7", "M1", -1, 5, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2|+M3|+M4|+M5|+M6|+M7|-M1|-M2|-M3|-M4|-M5|-M6", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_linear_long <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 0, FALSE,
    "M8", 0, FALSE,
    "M9", 0, FALSE,
    "M10", 0, FALSE,
    "M11", 0, FALSE,
    "M12", 0, FALSE,
    "M13", 0, FALSE,
    "M14", 0, FALSE,
    "M15", 0, FALSE,
    "M16", 0, FALSE,
    "M17", 0, FALSE,
    "M18", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M5", 1, 1, 2,
    "M5", "M6", 1, 1, 2,
    "M6", "M7", 1, 1, 2,
    "M7", "M8", 1, 1, 2,
    "M8", "M9", 1, 1, 2,
    "M9", "M10", 1, 1, 2,
    "M10", "M11", 1, 1, 2,
    "M11", "M12", 1, 1, 2,
    "M12", "M13", 1, 1, 2,
    "M13", "M14", 1, 1, 2,
    "M14", "M15", 1, 1, 2,
    "M15", "M16", 1, 1, 2,
    "M16", "M17", 1, 1, 2,
    "M17", "M18", 1, 1, 2,
    "M18", "M18", 1, 5, 2,
    "M18", "M1", -1, 5, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1", TRUE, TRUE,
    "S1", "S2", "+M2|+M3|+M4|+M5|+M6|+M7|+M8|+M9|+M10|+M11|+M12|+M13|+M14|+M15|+M16|+M17|+M18|-M1|-M2|-M3|-M4|-M5|-M6|-M7|-M8|-M9|-M10|-M11|-M12|-M13|-M14|-M15|-M16|-M17", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname modulenet
modulenet_trifurcating <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "M1", 1, TRUE,
    "M2", 0, FALSE,
    "M3", 0, FALSE,
    "M4", 0, FALSE,
    "M5", 0, FALSE,
    "M6", 0, FALSE,
    "M7", 1, TRUE,
    "M8", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M8", 1, 1, 2,
    "M8", "M2", 1, 1, 2,
    "M2", "M3", 1, 1, 2,
    "M2", "M4", 1, 1, 2,
    "M2", "M5", 1, 1, 2,
    "M3", "M3", 1, 10, 2,
    "M4", "M4", 1, 10, 2,
    "M5", "M5", 1, 10, 2,
    "M4", "M3", -1, 10, 2,
    "M3", "M4", -1, 10, 2,
    "M4", "M5", -1, 10, 2,
    "M3", "M5", -1, 10, 2,
    "M5", "M3", -1, 10, 2,
    "M5", "M4", -1, 10, 2,
    "M4", "M1", -1, 10, 2,
    "M3", "M1", -1, 10, 2,
    "M5", "M1", -1, 10, 2,
    "M5", "M6", 1, 1, 2,
    "M3", "M7", -1, 10, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "S0", "S1", "+M1,+M7", TRUE, TRUE,
    "S1", "S2", "+M8|+M2", FALSE, FALSE,
    "S2", "S3", "+M3|-M7|-M1|-M8|-M2", FALSE, FALSE,
    "S2", "S4", "+M4|-M1|-M8|-M2", FALSE, FALSE,
    "S2", "S5", "+M5|+M6|-M1|-M8|-M2", FALSE, FALSE
  )
  
  modulenet(module_info, module_network, expression_patterns)
}

