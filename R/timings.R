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
    do({
      dat <- .
      if (nrow(dat) > 1) {
        tibble(
          group = dat$group[[1]],
          task = dat$task[-length(dat$task)],
          time_elapsed = as.numeric(difftime(
            dat$time[-1],
            dat$time[-length(dat$time)],
            units = "secs"
          ))
        )
      } else {
        NULL
      }
    }) %>% 
    ungroup()
}
