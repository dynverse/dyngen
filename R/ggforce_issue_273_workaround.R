# these functions are workarounds for issue 273: https://github.com/thomasp85/ggforce/issues/273

#' @importFrom ggplot2 geom_blank
#' @importFrom rlang eval_tidy new_data_mask
geom_edge_loop_workaround <- function(edges, mapping, ...) {
  if (is.null(edges) || !any(eval_tidy(mapping$filter, new_data_mask(list2env(edges))))) {
    geom_blank()
  } else {
    geom_edge_loop(mapping, ...)
  }
}

geom_edge_fan_workaround <- function(edges, mapping, ...) {
  if (is.null(edges) || !any(eval_tidy(mapping$filter, new_data_mask(list2env(edges))))) {
    geom_blank()
  } else {
    geom_edge_fan(mapping, ...)
  }
}