.add_timing <- function(model, group, task) {
  entry <- tibble(group, task, time = Sys.time())
  model$timings <- bind_rows(model$timings, entry)
  model
}
#' Return the timings of each of the dyngen steps
#' 
#' @param model A dyngen object
#' 
#' @return A data frame with columns `"group"`, `"task"`, `"time_elapsed"`.
#' 
#' @export
#' 
#' @examples 
#' data("example_model")
#' timings <- get_timings(example_model)
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
