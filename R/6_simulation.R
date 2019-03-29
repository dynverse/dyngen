#' @export
simulation_default <- function(
  burn_time = 2,
  total_time = 20,
  num_simulations = 32,
  seeds = seq(1, num_simulations),
  ssa_algorithm = fastgssa::ssa.em(noise_strength = 4)
) {
  assert_that(length(seeds) == num_simulations)
  
  lst(
    burn_time,
    total_time,
    num_simulations,
    seeds,
    ssa_algorithm
  )
}

simulate_cells <- function(
  model
) {
  sim_params <- model$simulation_params
  sim_system <- model$simulation_system
  
  model$simulations <- 
    bind_rows(pbapply::pblapply(
      seq_len(sim_params$num_simulations),
      cl = model$num_cores,
      function(i) {
        set.seed(sim_params$seeds[[i]])
        
        propensity_funs <-
          sim_system$formulae %>% 
          select(formula_id, formula) %>% 
          deframe
        
        if (sim_params$burn_time > 0) {
          nus_burn <- sim_system$nus
          nus_burn[setdiff(rownames(nus_burn), sim_system$burn_variables),] <- 0
          
          # burn in
          out <- fastgssa::ssa(
            initial.state = sim_system$initial_state, 
            propensity.funs = propensity_funs,
            nu = nus_burn %>% as.matrix,
            final.time = sim_params$burn_time, 
            parms = sim_system$parameters,
            method = sim_params$ssa_algorithm,
            recalculate.all = TRUE, 
            stop.on.negstate = FALSE,
            stop.on.propensity = FALSE,
            verbose = FALSE
          )
          
          burn_out <- .simulate_cells_process_ssa(out)
          
          new_initial_state <- 
            burn_out %>% 
            slice(n()) %>% 
            select(one_of(sim_system$molecule_ids)) %>% 
            .[1, , drop = TRUE] %>% 
            unlist()
        } else {
          burn_out <- NULL
          new_initial_state <- system$initial_state
        }
        
        # actual simulation
        out <- fastgssa::ssa(
          initial.state = new_initial_state,
          propensity.funs = propensity_funs,
          nu = sim_system$nus %>% as.matrix,
          final.time = sim_params$total_time,
          parms = sim_system$parameters,
          method = sim_params$ssa_algorithm,
          recalculate.all = TRUE, 
          stop.on.negstate = FALSE,
          stop.on.propensity = FALSE,
          verbose = FALSE
        )
        
        # add both burnin as normal simulation together
        sim_out <- .simulate_cells_process_ssa(out)
        
        simulation <- 
          bind_rows(
            burn_out %>% mutate(t = t - max(t)), 
            sim_out
          ) %>% 
          mutate(simulation_i = i) %>% 
          select(simulation_i, t, everything())
        
        simulation
      }
    ))
  
  model
}

.simulate_cells_process_ssa <- function(out) {
  final_time <- out$args$final.time
  
  out$timeseries %>%
    filter(t <= final_time) %>%
    mutate(t = ifelse(row_number() == n(), final_time, t))
}
