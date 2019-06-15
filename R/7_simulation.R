#' @export
simulation_default <- function(
  burn_time = 2,
  total_time = 20,
  num_simulations = 32,
  seeds = sample.int(10 * num_simulations, num_simulations),
  ssa_algorithm = ssa_em(noise_strength = 10)
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

#' @export
generate_cells <- function(
  model
) {
  sim_params <- model$simulation_params
  
  if (model$verbose) cat("Running ", sim_params$num_simulations, " simulations\n", sep = "")
  
  # simulate cells one by one
  simulations <- 
    pbapply::pblapply(
      X = seq_len(sim_params$num_simulations),
      FUN = function(i) .generate_cells_simulate_cell(i, model = model)
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

.generate_cells_simulate_cell <- function(simulation_i, model) {
  sim_params <- model$simulation_params
  sim_system <- model$simulation_system
  
  set.seed(sim_params$seeds[[simulation_i]])
  
  propensity_funs <-
    sim_system$formulae %>% 
    select(formula_id, formula) %>% 
    deframe
  
  if (sim_params$burn_time > 0) {
    nus_burn <- sim_system$nus
    nus_burn[setdiff(rownames(nus_burn), sim_system$burn_variables),] <- 0
    
    # burn in
    out <- fastgssa(
      initial_state = sim_system$initial_state,
      propensity_funs = propensity_funs,
      nu = nus_burn,
      final_time = sim_params$burn_time,
      params = sim_system$parameters,
      method = sim_params$ssa_algorithm,
      stop_on_neg_state = FALSE,
      stop_on_neg_propensity = FALSE,
      verbose = FALSE
    )
    
    burn_meta <- 
      out$timeseries %>%
      transmute(time = ifelse(n() == row_number(), sim_params$burn_time, time)) %>%
      mutate(time = time - max(time))
    burn_counts <- do.call(rbind, out$timeseries$state) %>% Matrix::Matrix(sparse = TRUE)
    
    new_initial_state <-
      out$timeseries$state %>% last()
  } else {
    burn_meta <- NULL
    burn_counts <- NULL
    new_initial_state <- system$initial_state
  }
  
  # actual simulation
  out <- fastgssa(
    initial_state = new_initial_state, 
    propensity_funs = propensity_funs,
    nu = sim_system$nus,
    final_time = sim_params$total_time, 
    params = sim_system$parameters,
    method = sim_params$ssa_algorithm,
    stop_on_neg_state = FALSE,
    stop_on_neg_propensity = FALSE,
    verbose = FALSE
  )
  
  meta <- out$timeseries %>%
    transmute(time = ifelse(n() == row_number(), sim_params$total_time, time))
  counts <- do.call(rbind, out$timeseries$state) %>%
    Matrix::Matrix(sparse = TRUE)
  
  if (!is.null(burn_meta)) {
    meta <- bind_rows(burn_meta, meta) %>% mutate(simulation_i)
    counts <- rbind(burn_counts, counts)
  }
  
  lst(meta, counts)
}

.generate_cells_process_ssa <- function(out) {
  final_time <- out$args$final_time
  
  out$timeseries %>%
    filter(t <= final_time) %>%
    mutate(t = ifelse(row_number() == n(), final_time, t))
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
  dis <- dynutils::calculate_distance(gs_counts, sim_counts, method = model$distance_metric)
  best_match <- apply(dis, 2, which.min)
  
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
