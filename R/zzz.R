#' @importFrom Rcpp loadModule
#' @examples
#' library(Rcpp)
#' ssa <- new(dyngen:::SSA_EM, .01, 2)
#' ssa$simulate(initial_state, params, as.matrix(nu), final_time, max_duration, stop_on_neg_state, stop_on_neg_propensity, verbose, time_next_console)
loadModule("SSA_EM", TRUE)


