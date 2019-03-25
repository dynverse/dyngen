#' Module network models for starting the simulation with
#' 
#' A module is a group of genes which, to some extent, shows the same
#' expression behaviour. Several modules are connected together such that
#' one or more genes from one module will regulate the expression of
#' another module. By creating chains of modules, a dynamic behaviour in gene 
#' regulation can be created.
module_network <- function(
  nodes,
  edges,
  operations
) {
  lst(nodes, edges, operations)
}

#' @export
#' @rdname module_network
module_network_bifurcating <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "6", 0, FALSE,
    "7", 1, TRUE,
    "8", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "8", 1, 1, 2,
    "8", "2", 1, 0.2, 2,
    "2", "3", 1, 1, 2,
    "2", "4", 1, 1, 2,
    "3", "3", 1, 10, 2,
    "4", "4", 1, 10, 2,
    "4", "3", -1, 10, 2,
    "3", "4", -1, 10, 2,
    "4", "1", -1, 10, 2,
    "3", "1", -1, 10, 2,
    "3", "6", 1, 1, 2,
    "3", "7", -1, 10, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1,+7", TRUE, TRUE,
    "1", "2", "+8|+2", FALSE, FALSE,
    "2", "3", "+3|-1|-8,+6|-7,-2", FALSE, FALSE,
    "2", "4", "+4|-1|-8|-2", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_bifurcating_converging <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 0, FALSE,
    "8", 0, FALSE,
    "9", 0, FALSE,
    "10", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "3", "4", 1, 1, 2,
    "4", "4", 1, 10, 2,
    "3", "5", 1, 1, 2,
    "5", "5", 1, 10, 2,
    "4", "5", -1, 100, 2,
    "5", "4", -1, 100, 2,
    "4", "6", 1, 0.25, 5,
    "5", "7", 1, 0.25, 5,
    "6", "8", 1, 1, 2,
    "7", "8", 1, 1, 2,
    "8", "8", 1, 10, 2,
    "8", "4", -1, 10, 2,
    "8", "5", -1, 10, 2,
    "8", "9", 1, 1, 2,
    "9", "1", -1, 10, 2,
    "9", "10", 1, 1, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2|+3", FALSE, FALSE,
    "2", "3", "+5|+7|+8", FALSE, FALSE,
    "2", "4", "+4|+6|+8", FALSE, FALSE,
    "3", "5", "+9,-5|-7,+10", FALSE, FALSE,
    "4", "5", "+9,-4|-6,+10", FALSE, FALSE,
    "5", "6", "-1,-8|-9,-2|-10,-3", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_bifurcating_cycle <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 0, FALSE,
    "8", 0, FALSE,
    "9", 0, FALSE,
    "10", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "3", "4", 1, 1, 2,
    "4", "4", 1, 10, 2,
    "3", "5", 1, 1, 2,
    "5", "5", 1, 10, 2,
    "4", "5", -1, 100, 2,
    "5", "4", -1, 100, 2,
    "4", "6", 1, 0.1, 5,
    "5", "7", 1, 0.1, 5,
    "6", "8", 1, 1, 2,
    "7", "8", 1, 1, 2,
    "8", "8", 1, 10, 2,
    "8", "4", -1, 10, 2,
    "8", "5", -1, 10, 2,
    "8", "9", 1, 1, 2,
    "9", "1", -1, 10, 2,
    "9", "10", 1, 1, 2,
    "10", "8", -1, 100, 2,
    "10", "4", -1, 100, 2,
    "10", "5", -1, 100, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2|+3", FALSE, FALSE,
    "2", "3", "+4|+6|+8", FALSE, FALSE,
    "2", "4", "+5|+7|+8", FALSE, FALSE,
    "3", "5", "-4,+9,-6|+10,-1", FALSE, FALSE,
    "4", "5", "-5,+9,-7|+10,-1", FALSE, FALSE,
    "5", "1", "-2,-8|-3,-9|-10|+1", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_bifurcating_loop <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 0, FALSE,
    "8", 0, FALSE,
    "9", 0, FALSE,
    "10", 0, FALSE,
    "11", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "3", "4", 1, 1, 2,
    "4", "4", 1, 10, 2,
    "3", "5", 1, 1, 2,
    "5", "5", 1, 10, 2,
    "4", "5", -1, 100, 2,
    "5", "4", -1, 100, 2,
    "4", "6", 1, 0.1, 5,
    "6", "8", 1, 1, 2,
    "7", "8", 1, 1, 2,
    "8", "8", 1, 10, 2,
    "8", "4", -1, 10, 2,
    "8", "9", 1, 1, 2,
    "9", "1", -1, 10, 2,
    "9", "10", 1, 1, 2,
    "10", "8", -1, 100, 2,
    "10", "4", -1, 100, 2,
    "5", "11", 1, 10, 2,
    "5", "8", -1, 100, 2,
    "10", "5", -1, 100, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2|+3", FALSE, FALSE,
    "2", "4", "+4|+6|+8|+9,-4", FALSE, FALSE,
    "4", "1", "-1,+10,-6|-2,-8|-3,-9|-10|+1", FALSE, FALSE,
    "2", "3", "+5|+11", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_binary_tree <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1A", 1, TRUE,
    "1B", 0, FALSE,
    "2A", 0, FALSE,
    "2B", 0, FALSE,
    "2C", 0, FALSE,
    "3A", 0, FALSE,
    "3B", 0, FALSE,
    "3C", 0, FALSE,
    "4A", 0, FALSE,
    "5A", 0, FALSE,
    "4B", 0, FALSE,
    "5B", 0, FALSE,
    "5C", 0, FALSE,
    "6A", 0, FALSE,
    "6B", 0, FALSE,
    "7A", 0, FALSE,
    "7B", 0, FALSE,
    "8A", 0, FALSE,
    "8B", 0, FALSE,
    "8C", 0, FALSE,
    "9A", 0, FALSE,
    "9B", 0, FALSE,
    "10A", 0, FALSE,
    "10B", 0, FALSE,
    "10C", 0, FALSE,
    "10D", 0, FALSE,
    "10E", 0, FALSE,
    "11A", 0, FALSE,
    "11B", 0, FALSE,
    "11C", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1A", "1B", 1, 0.1, 2,
    "1B", "2A", 1, 1, 2,
    "2A", "2A", 1, 10, 2,
    "1B", "3A", 1, 1, 2,
    "3A", "3A", 1, 10, 2,
    "3A", "2A", -1, 10, 2,
    "2A", "3A", -1, 10, 2,
    "2A", "2B", 1, 1, 2,
    "2B", "2C", 1, 1, 2,
    "2C", "4A", 1, 1, 2,
    "2C", "5A", 1, 1, 2,
    "4A", "5A", -1, 10, 2,
    "5A", "4A", -1, 10, 2,
    "4A", "4B", 1, 1, 2,
    "5A", "5B", 1, 1, 2,
    "5B", "5C", 1, 1, 2,
    "5C", "6A", 1, 1, 2,
    "5C", "7A", 1, 1, 2,
    "6A", "7A", -1, 10, 2,
    "7A", "6A", -1, 10, 2,
    "6A", "6B", 1, 1, 2,
    "7A", "7B", 1, 1, 2,
    "3A", "3B", 1, 1, 2,
    "3B", "3C", 1, 1, 2,
    "3C", "8A", 1, 1, 2,
    "3C", "9A", 1, 1, 2,
    "8A", "9A", -1, 10, 2,
    "9A", "8A", -1, 10, 2,
    "8A", "8B", 1, 1, 2,
    "8B", "8C", 1, 1, 2,
    "9A", "9B", 1, 1, 2,
    "4B", "10A", 1, 1, 2,
    "4B", "11A", 1, 1, 2,
    "10A", "11A", -1, 10, 2,
    "11A", "10A", -1, 10, 2,
    "11A", "11B", 1, 1, 2,
    "11B", "11C", 1, 1, 2,
    "10A", "10B", 1, 1, 2,
    "10B", "10C", 1, 1, 2,
    "10C", "10D", 1, 1, 2,
    "10D", "10E", 1, 1, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1A", TRUE, TRUE,
    "1", "2", "+1B", FALSE, FALSE,
    "2", "3", "+2A|+2B|+2C", FALSE, FALSE,
    "2", "4", "+3A|+3B|+3C", FALSE, FALSE,
    "3", "5", "+4A|+4B", FALSE, FALSE,
    "3", "6", "+5A|+5B|+5C", FALSE, FALSE,
    "4", "7", "+8A|+8B|+8C", FALSE, FALSE,
    "4", "8", "+9A|+9B", FALSE, FALSE,
    "6", "9", "+6A|+6B", FALSE, FALSE,
    "6", "10", "+7A|+7B", FALSE, FALSE,
    "5", "11", "+10A|+10B|+10C|+10D|+10E", FALSE, FALSE,
    "5", "12", "+11A|+11B|+11C", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_consecutive_bifurcating <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 0, FALSE,
    "8", 0, FALSE,
    "9", 0, FALSE,
    "10", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 0.1, 2,
    "2", "3", 1, 1, 2,
    "3", "3", 1, 10, 2,
    "2", "4", 1, 1, 2,
    "4", "4", 1, 10, 2,
    "4", "3", -1, 10, 2,
    "3", "4", -1, 10, 2,
    "3", "5", 1, 0.1, 2,
    "5", "6", 1, 1, 2,
    "5", "7", 1, 1, 2,
    "6", "6", 1, 10, 2,
    "7", "7", 1, 10, 2,
    "6", "7", -1, 10, 2,
    "7", "6", -1, 10, 2,
    "4", "7", -1, 10, 2,
    "4", "6", -1, 10, 2,
    "7", "8", 1, 1, 2,
    "6", "9", 1, 1, 2,
    "4", "10", 1, 1, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2", FALSE, FALSE,
    "2", "3", "+3|+5", FALSE, FALSE,
    "2", "4", "+4|+10", FALSE, FALSE,
    "3", "5", "+7|+8", FALSE, FALSE,
    "3", "6", "+6|+9", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_converging <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, TRUE,
    "3", 0, TRUE,
    "4", 0, TRUE,
    "5", 0, TRUE,
    "6", 0, TRUE,
    "7", 0, TRUE,
    "8", 0, FALSE,
    "9", 0, FALSE,
    "10", 0, FALSE,
    "11", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "3", "4", 1, 1, 2,
    "4", "4", 1, 10, 2,
    "3", "5", 1, 1, 2,
    "5", "5", 1, 10, 2,
    "4", "5", -1, 100, 2,
    "5", "4", -1, 100, 2,
    "4", "6", 1, 0.25, 5,
    "5", "7", 1, 0.25, 5,
    "6", "8", 1, 1, 2,
    "7", "8", 1, 1, 2,
    "8", "8", 1, 10, 2,
    "8", "4", -1, 10, 2,
    "8", "5", -1, 10, 2,
    "8", "9", 1, 1, 2,
    "9", "1", -1, 10, 2,
    "9", "10", 1, 1, 2,
    "10", "11", 1, 1, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2|+3", FALSE, TRUE,
    "2", "3", "+5|+7", FALSE, TRUE,
    "2", "4", "+4|+6", FALSE, TRUE,
    "3", "5", "+8|+9,-5|-7,+10", FALSE, FALSE,
    "4", "5", "+8|+9,-4|-6,+10", FALSE, FALSE,
    "5", "6", "+11,-1,-8|-9,-2|-10,-3", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_cycle <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 1, FALSE,
    "3", 1, TRUE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, TRUE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", -1, 100, 2,
    "2", "3", -1, 100, 2,
    "3", "4", 1, 1, 2,
    "4", "5", 1, 1, 2,
    "5", "1", -1, 100, 2,
    "1", "6", 1, 1, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+3,+1|+6", TRUE, TRUE,
    "1", "2", "-1,+4|+5,-6", FALSE, FALSE,
    "2", "3", "-3,+2|-4|-5", FALSE, FALSE,
    "3", "1", "+1,-2|+6|+3", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_linear <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "3", "4", 1, 1, 2,
    "4", "5", 1, 1, 2,
    "5", "6", 1, 1, 2,
    "6", "7", 1, 1, 2,
    "7", "7", 1, 5, 2,
    "7", "1", -1, 5, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2|+3|+4|+5|+6|+7|-1|-2|-3|-4|-5|-6", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_linear_long <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 0, FALSE,
    "8", 0, FALSE,
    "9", 0, FALSE,
    "10", 0, FALSE,
    "11", 0, FALSE,
    "12", 0, FALSE,
    "13", 0, FALSE,
    "14", 0, FALSE,
    "15", 0, FALSE,
    "16", 0, FALSE,
    "17", 0, FALSE,
    "18", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "3", "4", 1, 1, 2,
    "4", "5", 1, 1, 2,
    "5", "6", 1, 1, 2,
    "6", "7", 1, 1, 2,
    "7", "8", 1, 1, 2,
    "8", "9", 1, 1, 2,
    "9", "10", 1, 1, 2,
    "10", "11", 1, 1, 2,
    "11", "12", 1, 1, 2,
    "12", "13", 1, 1, 2,
    "13", "14", 1, 1, 2,
    "14", "15", 1, 1, 2,
    "15", "16", 1, 1, 2,
    "16", "17", 1, 1, 2,
    "17", "18", 1, 1, 2,
    "18", "18", 1, 5, 2,
    "18", "1", -1, 5, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1", TRUE, TRUE,
    "1", "2", "+2|+3|+4|+5|+6|+7|+8|+9|+10|+11|+12|+13|+14|+15|+16|+17|+18|-1|-2|-3|-4|-5|-6|-7|-8|-9|-10|-11|-12|-13|-14|-15|-16|-17", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


#' @export
#' @rdname module_network
module_network_trifurcating <- function() {
  nodes <- tribble(
    ~module_id, ~a0, ~burn,
    "1", 1, TRUE,
    "2", 0, FALSE,
    "3", 0, FALSE,
    "4", 0, FALSE,
    "5", 0, FALSE,
    "6", 0, FALSE,
    "7", 1, TRUE,
    "8", 0, FALSE
  )

  edges <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "1", "8", 1, 1, 2,
    "8", "2", 1, 1, 2,
    "2", "3", 1, 1, 2,
    "2", "4", 1, 1, 2,
    "2", "5", 1, 1, 2,
    "3", "3", 1, 10, 2,
    "4", "4", 1, 10, 2,
    "5", "5", 1, 10, 2,
    "4", "3", -1, 10, 2,
    "3", "4", -1, 10, 2,
    "4", "5", -1, 10, 2,
    "3", "5", -1, 10, 2,
    "5", "3", -1, 10, 2,
    "5", "4", -1, 10, 2,
    "4", "1", -1, 10, 2,
    "3", "1", -1, 10, 2,
    "5", "1", -1, 10, 2,
    "5", "6", 1, 1, 2,
    "3", "7", -1, 10, 2
  )

  operations <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn,
    "0", "1", "+1,+7", TRUE, TRUE,
    "1", "2", "+8|+2", FALSE, FALSE,
    "2", "3", "+3|-7|-1|-8|-2", FALSE, FALSE,
    "2", "4", "+4|-1|-8|-2", FALSE, FALSE,
    "2", "5", "+5|+6|-1|-8|-2", FALSE, FALSE
  )
  
  module_network(nodes, edges, operations)
}


