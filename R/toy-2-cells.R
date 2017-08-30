#' Assign random progressions without tenting
#' @export
random_progressions <- function(milestone_network, ncells=100) {
  cell_ids <- paste0("C", seq_len(ncells))
  tibble(cell_id = cell_ids) %>%
    bind_cols(milestone_network[sample(seq_len(nrow(milestone_network)), length(cell_ids), replace=TRUE, prob=milestone_network$length), ]) %>%
    mutate(percentage = map_dbl(length, ~runif(1, 0, .))) %>%
    select(cell_id, from, to, percentage)
}

#' Assign random progressions with tenting
#' @export
random_progressions_tented <- function(milestone_network, ncells=100) {
  from_probabilities <- milestone_network %>% group_by(from) %>% summarise(prob=sqrt(sum(length^2)))
  cell_ids <- paste0("C", seq_len(ncells))
  
  tibble(cell_id = cell_ids) %>% 
    mutate(from=sample(from_probabilities$from, length(cell_ids), replace=TRUE, prob=from_probabilities$prob)) %>% 
    left_join(milestone_network, by="from") %>% 
    group_by(cell_id) %>% 
    mutate(
      from_percentage = runif(n()), # first calculate the from percentage
      percentage_relative = runif(n()), # use this from percentage to extract the to percentages
      percentage = (1-from_percentage)*(percentage_relative/sum(percentage_relative))
    ) %>% 
    ungroup() %>% 
    select(cell_id, from, to, percentage)
}

# progressions <- random_progressions_tented(milestone_network, 1000)
# task <- dyneval::wrap_ti_prediction("lol", "lol", unique(progressions$cell_id), unique(c(milestone_network$from, milestone_network$to)), milestone_network, progressions=progressions)
# dynplot::plot_strip_connections(task, task)
# dyneval::plot_default(task)
