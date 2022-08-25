# these functions are workarounds for issue 273: https://github.com/thomasp85/ggforce/issues/273

#' @importFrom ggplot2 geom_blank
#' @importFrom rlang eval_tidy new_data_mask
geom_edge_loop_workaround <- function(edges, mapping, ...) {
  filter_result <- eval_tidy(mapping$filter, new_data_mask(list2env(edges)))
  if (any(filter_result)) {
    geom_edge_loop(mapping, ...)
  } else {
    geom_blank()
  }
}

geom_edge_fan_workaround <- function(edges, mapping, ...) {
  filter_result <- eval_tidy(mapping$filter, new_data_mask(list2env(edges)))
  if (any(filter_result)) {
    geom_edge_fan(mapping, ...)
  } else {
    geom_blank()
  }
}