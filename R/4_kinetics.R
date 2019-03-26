kinetics_samplers <- function(
  sample_r = function(n) runif(n, 10, 200), 
  sample_d = function(n) runif(n, 2, 8), 
  sample_p = function(n) runif(n, 2, 8), 
  sample_q = function(n) runif(n, 1, 5),
  calculate_a0 = function(effects) {
    if (all(effects != 1)) {
      1
    } else if (all(effects != -1)) {
      0.0001
    } else {
      0.5
    }
  },
  calculate_a = function(configuration_id, effects) {
    bound <- get_binding_configuration(configuration_id, length(effects)) # TODO: ??
    if(any(effects[bound] == -1)) {
      0
    } else {
      1
    }
  },
  sample_strength = function(n) runif(n, 1, 20),
  calculate_k = function(max_protein, strength) {
    max_protein / 2 / strength
  },
  sample_cooperativity = function(n) runif(n, 1, 4)
) {
  lst(
    sample_r,
    sample_d,
    sample_p,
    sample_q,
    calculate_a0,
    calculate_a,
    sample_strength,
    calculate_k,
    sample_cooperativity
  )
}