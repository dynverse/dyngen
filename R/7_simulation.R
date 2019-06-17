#' @export
simulation_default <- function(
  burn_time = 2,
  total_time = 20,
  census_interval = .01,
  num_simulations = 32,
  seeds = sample.int(10 * num_simulations, num_simulations),
  ssa_algorithm = ssa_direct()
) {
  assert_that(length(seeds) == num_simulations)
  
  lst(
    burn_time,
    total_time,
    census_interval,
    num_simulations,
    seeds,
    ssa_algorithm
  )
}

#' @export
generate_cells <- function(
  model
) {
  sim_params <- model$simulation_params
  
  if (model$verbose) cat("Precompiling propensity functions for simulations\n")
  propensity_funs <- .generate_cells_precompile_propensity_funs(model, env = environment())
  
  # simulate cells one by one
  if (model$verbose) cat("Running ", sim_params$num_simulations, " simulations\n", sep = "")
  simulations <- 
    pbapply::pblapply(
      X = seq_len(sim_params$num_simulations),
      cl = model$num_cores,
      FUN = function(i) .generate_cells_simulate_cell(i, model = model, propensity_funs = propensity_funs)
    )
  
  # split up simulation data
  model$simulations <- lst(
    meta = map_df(simulations, "meta"),
    counts = do.call(rbind, map(simulations, "counts"))
  )
  
  # predict state
  if (model$verbose) cat("Mapping simulations to gold standard\n", sep = "")
  model$simulations$meta <- .generate_cells_predict_state(model)
  
  # perform dimred
  if (model$verbose) cat("Performing dimred\n", sep = "")
  model <- model %>% calculate_dimred()
  
  # return
  model
}


#' @importFrom fastgssa compile_propensity_functions
.generate_cells_precompile_propensity_funs <- function(model, env = parent.frame()) {
  # fetch paraneters and settings
  sim_system <- model$simulation_system
  
  # filter propensity functions
  propensity_funs <-
    sim_system$formulae %>% 
    select(formula_id, formula) %>% 
    deframe
  
  # compile prop funs
  comp_funs <- fastgssa::compile_propensity_functions(
    propensity_funs = propensity_funs,
    state = sim_system$initial_state,
    params = sim_system$parameters,
    env = env
  )
  
  comp_funs
}

.generate_cells_simulate_cell <- function(simulation_i, model, propensity_funs) {
  sim_params <- model$simulation_params
  sim_system <- model$simulation_system
  
  set.seed(sim_params$seeds[[simulation_i]])
  
  if (sim_params$burn_time > 0) {
    nus_burn <- sim_system$nus
    nus_burn[setdiff(rownames(nus_burn), sim_system$burn_variables),] <- 0
    
    # burn in
    out <- fastgssa::ssa(
      initial_state = sim_system$initial_state,
      propensity_funs = propensity_funs,
      nu = nus_burn,
      final_time = sim_params$burn_time,
      census_interval = sim_params$census_interval,
      params = sim_system$parameters,
      method = sim_params$ssa_algorithm,
      stop_on_neg_state = FALSE,
      verbose = FALSE
    )
    
    burn_meta <- 
      out$output %>%
      transmute(time = ifelse(n() == row_number(), sim_params$burn_time, time)) %>%
      mutate(time = time - max(time))
    burn_counts <- do.call(rbind, out$output$state) %>% Matrix::Matrix(sparse = TRUE)
    
    new_initial_state <-
      out$output$state %>% last()
  } else {
    burn_meta <- NULL
    burn_counts <- NULL
    new_initial_state <- system$initial_state
  }
  
  # actual simulation
  out <- fastgssa::ssa(
    initial_state = new_initial_state, 
    propensity_funs = propensity_funs,
    nu = sim_system$nus,
    final_time = sim_params$total_time, 
    census_interval = sim_params$census_interval,
    params = sim_system$parameters,
    method = sim_params$ssa_algorithm,
    stop_on_neg_state = FALSE,
    verbose = FALSE
  )
  
  meta <- out$output %>%
    transmute(time = ifelse(n() == row_number(), sim_params$total_time, time))
  counts <- do.call(rbind, out$output$state) %>%
    Matrix::Matrix(sparse = TRUE)
  
  if (!is.null(burn_meta)) {
    meta <- bind_rows(burn_meta, meta) %>% mutate(simulation_i)
    counts <- rbind(burn_counts, counts)
  }
  
  lst(meta, counts)
}

.generate_cells_predict_state <- function(model) {
  # fetch gold standard data
  gs_meta <- model$gold_standard$meta
  gold_ix <- !gs_meta$burn
  gs_meta <- gs_meta[gold_ix, , drop = FALSE]
  gs_counts <- model$gold_standard$counts[gold_ix, , drop = FALSE]
  gs_dimred <- model$gold_standard$dimred[gold_ix, , drop = FALSE]
  
  # fetch simulation data
  # (gold standard counts only contains TFs, so filter those)
  sim_meta <- model$simulations$meta
  sim_counts <- model$simulations$counts[, colnames(gs_counts), drop = FALSE]
  
  # calculate 1NN -> a full distance matrix could be avoided
  if (nrow(sim_counts) > 10000) {
    start <- seq(1, nrow(sim_counts), by = 10000)
    stop <- pmin(start + 9999, nrow(sim_counts))
    best_matches <- pmap(lst(start, stop), function(start, stop) {
      dis <- dynutils::calculate_distance(gs_counts, sim_counts[start:stop, , drop = FALSE], method = model$distance_metric)
      best_match <- apply(dis, 2, which.min)
    })
    best_match <- unlist(best_matches)
  } else {
    dis <- dynutils::calculate_distance(gs_counts, sim_counts, method = model$distance_metric)
    best_match <- apply(dis, 2, which.min)
  }
  
  # add predictions to sim_meta
  sim_meta <- 
    bind_cols(
      sim_meta %>% select(simulation_i, sim_time = time),
      gs_meta[best_match, , drop = FALSE] %>% select(from, to, time)
    ) %>% 
    group_by(from, to) %>% 
    mutate(time = dynutils::scale_minmax(time)) %>% 
    ungroup()
  
  # check if all branches are present
  sim_edges <- sim_meta %>% group_by(from, to) %>% summarise(n = n())
  network <- full_join(model$gold_standard$network, sim_edges, by = c("from", "to"))
  
  if (any(is.na(network$n))) {
    warning("Simulation does not contain all gold standard edges. This simulation likely suffers from bad kinetics; choose a different seed and rerun.")
  }
  
  # return output
  sim_meta
}
