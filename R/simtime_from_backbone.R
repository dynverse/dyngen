#' Determine simulation time from backbone
#' 
#' @param backbone A valid dyngen backbone object
#' @param burn Whether or not to compute the simtime for only the burn phase
#' 
#' @export
simtime_from_backbone <- function(backbone, burn = FALSE) {
  exp_pat <- backbone$expression_patterns
  if (burn) exp_pat <- exp_pat %>% filter(burn)
  sim_time_sum <- exp_pat %>% filter(start) %>% pull(from) %>% set_names(rep(0, length(.)), .)
  for (i in seq_len(nrow(exp_pat))) {
    sim_time_sum[[exp_pat$to[[i]]]] <- sim_time_sum[[exp_pat$from[[i]]]] + exp_pat$time[[i]]
  }
  total_time <- max(sim_time_sum * 1.2)
  total_time
}