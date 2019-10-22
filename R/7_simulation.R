#' Simulate the cells
#' 
#' [generate_cells()] runs simulations in order to determine the gold standard
#' of the simulations.
#' [simulation_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the gold standard been generated with [generate_gold_standard()].
#' @param burn_time The burn in time of the system, used to determine an initial state vector.
#' @param total_time The total simulation time of the system.
#' @param ssa_algorithm Which SSA algorithm to use for simulating the cells with [GillespieSSA2::ssa()]
#' @param census_interval A granularity parameter for the outputted simulation.
#' @param num_simulations The number of simulations to run.
#' @param seeds A set of seeds for each of the simulations.
#' @param store_grn Whether or not to store the GRN activation values.
#' @param store_reaction_firings Whether or not to store the number of reaction firings.
#' @param store_reaction_propensities Whether or not to store the propensity values of the reactions.
#' 
#' @importFrom GillespieSSA2 ssa
#' @export
generate_cells <- function(model) {
  if (model$verbose) cat("Precompiling reactions for simulations\n")
  reactions <- .generate_cells_precompile_reactions(model)
  
  # simulate cells one by one
  if (model$verbose) cat("Running ", model$simulation_params$num_simulations, " simulations\n", sep = "")
  simulations <- 
    pbapply::pblapply(
      X = seq_len(model$simulation_params$num_simulations),
      cl = model$num_cores,
      FUN = .generate_cells_simulate_cell,
      model = model,
      reactions = reactions
    )
  
  # split up simulation data
  model$simulations <- lst(
    meta = map_df(simulations, "meta"),
    counts = do.call(rbind, map(simulations, "counts")),
    regulation = do.call(rbind, map(simulations, "regulation")),
    reaction_firings = do.call(rbind, map(simulations, "reaction_firings")),
    reaction_propensities = do.call(rbind, map(simulations, "reaction_propensities"))
  )
  
  # predict state
  if (model$verbose) cat("Mapping simulations to gold standard\n", sep = "")
  if (!is.null(model[["gold_standard"]])) {
    model$simulations$meta <- .generate_cells_predict_state(model)
  } else {
    model$simulations$meta <- model$simulations$meta %>% rename(sim_time = time)
  }
  
  # perform dimred
  if (model$verbose) cat("Performing dimred\n", sep = "")
  model <- model %>% calculate_dimred()
  
  # return
  model
}

#' @export
#' @rdname generate_cells
#' @importFrom GillespieSSA2 ssa_etl
simulation_default <- function(
  burn_time = 2,
  total_time = 10,
  ssa_algorithm = ssa_etl(tau = .001),
  census_interval = .01,
  num_simulations = 32,
  seeds = sample.int(10 * num_simulations, num_simulations),
  store_grn = TRUE,
  store_reaction_firings = FALSE,
  store_reaction_propensities = FALSE
) {
  assert_that(length(seeds) == num_simulations)
  
  lst(
    burn_time,
    total_time,
    ssa_algorithm,
    census_interval,
    num_simulations,
    seeds,
    store_grn,
    store_reaction_firings,
    store_reaction_propensities
  )
}

#' @importFrom GillespieSSA2 compile_reactions
.generate_cells_precompile_reactions <- function(model) {
  # fetch paraneters and settings
  sim_system <- model$simulation_system
  reactions <- sim_system$reactions 
  
  buffer_ids <- unique(unlist(map(reactions, "buffer_ids")))
  
  # compile prop funs
  comp_funs <- GillespieSSA2::compile_reactions(
    reactions = reactions,
    buffer_ids = buffer_ids,
    state_ids = if (is.matrix(sim_system$initial_state)) colnames(sim_system$initial_state) else names(sim_system$initial_state),
    params = sim_system$parameters,
    hardcode_params = FALSE,
    fun_by = 1000L
  )
  
  comp_funs
}

.generate_cells_simulate_cell <- function(simulation_i, model, reactions, verbose = FALSE) {
  sim_params <- model$simulation_params
  sim_system <- model$simulation_system
  
  set.seed(sim_params$seeds[[simulation_i]])
  
  # get initial state
  initial_state <- sim_system$initial_state
  if (is.matrix(initial_state)) {
    initial_state <- initial_state[simulation_i, ]
  }
  
  if (sim_params$burn_time > 0) {
    burn_reaction_firings <- reactions
    rem <- setdiff(sim_system$molecule_ids, sim_system$burn_variables)
    if (length(rem) > 0) {
      burn_reaction_firings$state_change[match(rem, sim_system$molecule_ids), ] <- 0
    }
    
    # burn in
    out <- GillespieSSA2::ssa(
      initial_state = initial_state,
      reactions = burn_reaction_firings,
      final_time = sim_params$burn_time,
      census_interval = sim_params$census_interval,
      params = sim_system$parameters,
      method = sim_params$ssa_algorithm,
      stop_on_neg_state = FALSE,
      verbose = verbose,
      log_buffer = sim_params$store_grn,
      log_firings = sim_params$store_reaction_firings,
      log_propensity = sim_params$store_reaction_propensities
    )
    
    burn_meta <- 
      tibble(
        time = c(head(out$time, -1), sim_params$burn_time)
      ) %>%
      mutate(time = time - max(time))
    burn_counts <- out$state %>% Matrix::Matrix(sparse = TRUE)
    burn_regulation <- if (sim_params$store_grn) out$buffer %>% Matrix::Matrix(sparse = TRUE) else NULL
    burn_reaction_firings <- if (sim_params$store_reaction_firings) out$firings %>% Matrix::Matrix(sparse = TRUE) else NULL
    burn_reaction_propensities <- if (sim_params$store_reaction_propensities) out$propensity %>% Matrix::Matrix(sparse = TRUE) else NULL
    
    new_initial_state <- out$state[nrow(out$state), ]
  } else {
    burn_meta <- NULL
    burn_counts <- NULL
    new_initial_state <- initial_state
    burn_regulation <- NULL
    burn_reaction_firings <- NULL
  }
  
  # actual simulation
  out <- GillespieSSA2::ssa(
    initial_state = new_initial_state, 
    reactions = reactions,
    final_time = sim_params$total_time, 
    census_interval = sim_params$census_interval,
    params = sim_system$parameters,
    method = sim_params$ssa_algorithm,
    stop_on_neg_state = FALSE,
    verbose = verbose,
    log_buffer = sim_params$store_grn,
    log_firings = sim_params$store_reaction_firings,
    log_propensity = sim_params$store_reaction_propensities
  )
  
  meta <- 
    tibble(
      time = c(head(out$time, -1), sim_params$total_time)
    )
  counts <- out$state %>% Matrix::Matrix(sparse = TRUE)
  regulation <- if (sim_params$store_grn) out$buffer %>% Matrix::Matrix(sparse = TRUE) else NULL
  reaction_firings <- if (sim_params$store_reaction_firings) out$firings %>% Matrix::Matrix(sparse = TRUE) else NULL
  reaction_propensities <- if (sim_params$store_reaction_propensities) out$propensity %>% Matrix::Matrix(sparse = TRUE) else NULL
  
  if (!is.null(burn_meta)) {
    meta <- bind_rows(burn_meta, meta)
    counts <- rbind(burn_counts, counts)
    regulation <- rbind(burn_regulation, regulation)
    reaction_firings <- rbind(burn_reaction_firings, reaction_firings)
    reaction_propensities <- rbind(burn_reaction_propensities, reaction_propensities)
  }
  
  meta <- meta %>% mutate(simulation_i) %>% select(simulation_i, everything())
  
  lst(meta, counts, regulation, reaction_firings, reaction_propensities)
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
