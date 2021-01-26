.add_timing <- function(model, group, task) {
  entry <- tibble(group, task, time = Sys.time())
  model$timings <- bind_rows(model$timings, entry)
  model
}
#' Return the timings of each of the dyngen steps
#' @param model A dyngen object
#' 
#' @export
#' 
#' @examples 
#' model <- 
#'   initialise_model(backbone = backbone_linear()) %>% 
#'   generate_tf_network()
#' 
#' timings <- get_timings(model)
get_timings <- function(model) {
  model$timings %>%
    group_by(.data$group) %>% 
    summarise(
      task = .data$task[-length(.data$task)],
      time_elapsed = as.numeric(difftime(
        .data$time[-1],
        .data$time[-length(.data$time)],
        units = "secs"
      )),
      .groups = "drop"
    ) 
}
