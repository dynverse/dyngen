#' @export
list_experiment_samplers <- function() {
  lst(
    snapshot = experiment_sampler_snapshot,
    synchronised = experiment_sampler_synchronised
  )
}

#' @export
experiment_sampler_snapshot <- function(
  weight_bw = 0.1
) {
  sampler_fun <- function(simulation, gold_standard, num_cells) {
    # do sampling
  }
  
  lst(
    sampler_type = "snapshot",
    sampler_fun,
    weight_bw
  )
}

#' @export
experiment_sampler_synchronised <- function(
  num_timepoints = 10
) {
  sampler_fun <- function(simulation, gold_standard, num_cells) {
    # do sampling
  }
  
  lst(
    sampler_type = "synchronised",
    sampler_fun,
    num_timepoints
  )
}