#' @export
simulation_default <- function(
  burn_time = 2,
  total_time = 10,
  num_simulations = 16,
  ssa_algorithm = fastgssa::ssa.em(noise_strength = 4)
) {
  lst(
    burn_time,
    total_time,
    num_simulations,
    ssa_algorithm
  )
}