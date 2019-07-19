#' Simulate the gold standard
#' 
#' [generate_gold_standard()] runs simulations in order to determine the gold standard
#' of the simulations.
#' [gold_standard_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the kinetics of the feature network has been generated with [generate_kinetics()].
#' @param tau The time step of the ODE algorithm used to generate the gold standard.
#' @param census_interval A granularity parameter of the gold standard time steps. Should be larger than or equal to `tau`.
#' 
#' @export
generate_gold_standard <- function(model) {
  model$gold_standard <- list()
  
  # compute changes in modules along the edges
  if (model$verbose) cat("Generating gold standard mod changes\n")
  model$gold_standard$mod_changes <- .generate_gold_standard_mod_changes(model$backbone$expression_patterns)
  
  # precompile propensity functions
  if (model$verbose) cat("Precompiling reactions for gold standard\n")
  prep_data <- .generate_gold_precompile_reactions(model)
  
  # run gold standard simulations
  if (model$verbose) cat("Running gold simulations\n")
  simulations <- .generate_gold_standard_simulations(model, prep_data)
  model$gold_standard$meta <- simulations$meta
  model$gold_standard$counts <- simulations$counts
  
  # do combined dimred
  model <- model %>% calculate_dimred()
  
  # construct simulation network from dimred
  model$gold_standard$network <- .generate_gold_standard_generate_network(model)
  
  model
}

#' @export
#' @rdname generate_gold_standard
gold_standard_default <- function(
  tau = .001,
  census_interval = .01
) {
  lst(
    tau,
    census_interval
  )
}

.generate_gold_standard_mod_changes <- function(expression_patterns) {
  expression_patterns %>% 
    mutate(
      mod_diff = module_progression %>% strsplit("\\|"),
      substate = map(mod_diff, seq_along)
    ) %>% 
    select(-module_progression) %>% 
    unnest(mod_diff, substate) %>% 
    mutate(
      mod_diff2 = strsplit(mod_diff, ","),
      mod_on = map(mod_diff2, function(x) x %>% keep(grepl("\\+", .)) %>% gsub("\\+", "", .)),
      mod_off = map(mod_diff2, function(x) x %>% keep(grepl("-", .)) %>% gsub("-", "", .))
    ) %>% 
    select(-mod_diff2) %>% 
    group_by(from, to) %>% 
    mutate(
      from_ = ifelse(row_number() == 1, from, c("", paste0(from, to, "p", (row_number() - 1)))),
      to_ = ifelse(row_number() == n(), to, from_[row_number() + 1])
    ) %>% 
    ungroup()
}

#' @importFrom gillespie compile_reactions
.generate_gold_precompile_reactions <- function(model) {
  # fetch paraneters and settings
  sim_system <- model$simulation_system
  tf_info <- model$feature_info %>% filter(is_tf)
  
  # select relevant functions
  tf_molecules <- tf_info %>% select(w, x, y) %>% gather(col, val) %>% pull(val)
  
  # filter reactions
  reactions <- sim_system$reactions %>% 
    keep(~ all(names(.$effect) %in% tf_molecules))
  
  buffer_ids <- unique(unlist(map(reactions, "buffer_ids")))
  
  comp_funs <- gillespie::compile_reactions(
    reactions = reactions,
    buffer_ids = buffer_ids,
    state_ids = tf_molecules,
    params = sim_system$parameters,
    hardcode_params = TRUE
  )
  
  lst(
    tf_molecules,
    reactions = comp_funs
  )
}

#' @importFrom pbapply timerProgressBar setTimerProgressBar
#' @importFrom gillespie ssa ode_em
.generate_gold_standard_simulations <- function(model, prep_data) {
  # fetch paraneters and settings
  mod_changes <- model$gold_standard$mod_changes
  gold_params <- model$gold_standard_params
  sim_system <- model$simulation_system
  tf_info <- model$feature_info %>% filter(is_tf)
  
  # determine ode algorithm
  algo <- gillespie::ode_em(tau = gold_params$tau, noise_strength = 0)
  
  # select relevant functions
  tf_molecules <- prep_data$tf_molecules
  reactions <- prep_data$reactions
  
  # start constructing golds
  gold_sim_outputs <- list()
  gold_sim_vectors <- list()
  gold_sim_modules <- list()
  
  start_state <- mod_changes$from[[1]]
  
  gold_sim_vectors[[start_state]] <- model$simulation_system$initial_state[tf_molecules] %>% as.matrix
  gold_sim_modules[[start_state]] <- c()
  
  timer <- pbapply::timerProgressBar(
    min = 0, 
    max = nrow(mod_changes),
    width = 50
  )
  for (i in seq_len(nrow(mod_changes))) {
    from_ <- mod_changes$from_[[i]]
    to_ <- mod_changes$to_[[i]]
    time <- mod_changes$time[[i]]
    
    # which modules are on
    mods <- 
      gold_sim_modules[[from_]] %>% 
      union(mod_changes$mod_on[[i]]) %>% 
      setdiff(mod_changes$mod_off[[i]])
    
    gold_sim_modules[[to_]] <- mods
    
    # which tfs are on
    tfs_on <- tf_info %>% filter(module_id %in% mods)
    molecules_on <- tfs_on %>% select(w, x, y) %>% gather(col, val) %>% pull(val)
    
    # fetch initial state
    new_initial_state <- rowMeans(gold_sim_vectors[[from_]])
    
    # generate nus
    new_reactions <- reactions
    rem <- setdiff(tf_molecules, molecules_on)
    if (length(rem) > 0) {
      new_reactions$state_change[match(rem, tf_molecules), ] <- 0
    }
    
    # simulation of gold standard edge
    out <- gillespie::ssa(
      initial_state = new_initial_state,
      reactions = new_reactions,
      final_time = time,
      params = sim_system$parameters,
      method = algo,
      census_interval = gold_params$census_interval,
      verbose = FALSE
    )
    
    time_out <- out$time
    state_out <- out$state
    
    meta <- tibble(from_, to_, time = time_out)
    counts <- state_out %>% Matrix::Matrix(sparse = TRUE)
      
    pbapply::setTimerProgressBar(timer, value = i)
    
    gold_sim_outputs[[i]] <- lst(meta, counts)
    
    end_state <- counts[nrow(counts), ] %>% as.matrix()
    if (!to_ %in% gold_sim_vectors) {
      gold_sim_vectors[[to_]] <- end_state
    } else {
      gold_sim_vectors[[to_]] <- cbind(gold_sim_vectors[[to_]], end_state)
    }
  }
  cat("\n")
  
  meta <- map_df(gold_sim_outputs, "meta")
  counts <- do.call(rbind, map(gold_sim_outputs, "counts")) %>% Matrix::Matrix(sparse = TRUE)
  
  meta <- meta %>% 
    left_join(mod_changes %>% select(from, to, from_, to_, substate, burn, time_per_edge = time), by = c("from_", "to_")) %>% 
    group_by(from, to) %>%
    mutate(
      simulation_i = 0,
      sim_time = time,
      time = ((substate - 1) * time_per_edge + time) / time_per_edge / max(substate)
    ) %>% 
    ungroup() %>% 
    select(-substate, -time_per_edge) %>% 
    select(simulation_i, sim_time, burn, from, to, from_, to_, time)
  
  lst(meta, counts)
}

.generate_gold_standard_generate_network <- function(model) {
  gs_dimred <- model$gold_standard$dimred
  gs_meta <- model$gold_standard$meta
  
  gs_meta %>% 
    filter(!burn) %>% 
    mutate(i = row_number()) %>% 
    group_by(from, to) %>% 
    summarise(
      length = sqrt(sum((gs_dimred[i[-1], , drop = FALSE] - gs_dimred[i[-n()], , drop = FALSE])^2)),
      directed = TRUE
    ) %>% 
    ungroup()
}
