#' @importFrom fastgssa ssa_em ssa_direct
#' @export
gold_standard_default <- function(
  ssa_algorithm = ssa_direct()
) {
  lst(
    ssa_algorithm
  )
}

#' @export
generate_gold_standard <- function(model) {
  model$gold_standard <- list()
  
  # compute changes in modules along the edges
  if (model$verbose) cat("Generating gold standard mod changes\n")
  model$gold_standard$mod_changes <- .generate_gold_standard_mod_changes(model)
  
  # precompile propensity functions
  if (model$verbose) cat("Precompiling propensity functions for gold standard\n")
  prep_data <- .generate_gold_precompile_propensity_funs(model)
  
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

.generate_gold_standard_mod_changes <- function(model) {
  mod_changes <- 
    model$backbone$expression_patterns %>% 
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

  
  mod_changes
}

#' @importFrom fastgssa compile_propensity_functions
.generate_gold_precompile_propensity_funs <- function(model) {
  # fetch paraneters and settings
  sim_system <- model$simulation_system
  tf_info <- model$feature_info %>% filter(is_tf)
  
  # select relevant functions
  tf_molecules <- tf_info %>% select(w, x, y) %>% gather(col, val) %>% pull(val)
  
  # filter propensity functions
  formulae <-
    sim_system$formulae %>% 
    filter(molecule %in% tf_molecules)
  
  # filter nus
  nus <- sim_system$nus[tf_molecules, formulae$formula_id, drop = FALSE]
  
  comp_funs <- fastgssa::compile_propensity_functions(
    propensity_funs = formulae$formula,
    buffer_ids = unlist(formulae$buffer_ids),
    state_ids = tf_molecules,
    params = sim_system$parameters,
    hardcode_params = TRUE
  )
  
  lst(
    tf_molecules,
    nus,
    propensity_funs = comp_funs
  )
}

#' @importFrom pbapply timerProgressBar setTimerProgressBar
#' @importFrom fastgssa ssa
.generate_gold_standard_simulations <- function(model, prep_data) {
  # fetch paraneters and settings
  mod_changes <- model$gold_standard$mod_changes
  sim_system <- model$simulation_system
  tf_info <- model$feature_info %>% filter(is_tf)
  
  # select relevant functions
  tf_molecules <- prep_data$tf_molecules
  nus <- prep_data$nus
  propensity_funs <- prep_data$propensity_funs
  
  # start constructing golds
  gold_sim_outputs <- list()
  gold_sim_vectors <- list()
  gold_sim_modules <- list()
  
  start_state <- mod_changes$from[[1]]
  
  gold_sim_vectors[[start_state]] <- model$simulation_system$initial_state[tf_molecules] %>% as.matrix
  gold_sim_modules[[start_state]] <- c()
  
  timer <- pbapply::timerProgressBar(min = 0, max = nrow(mod_changes))
  for (i in seq_len(nrow(mod_changes))) {
    from_ <- mod_changes$from_[[i]]
    to_ <- mod_changes$to_[[i]]
    time <- mod_changes$time[[i]]
    
    # which modules are on
    mods <- 
      gold_sim_modules[[from_]] %>% 
      union(mod_changes$mod_on[[i]]) %>% 
      setdiff(mod_changes$mod_off[[i]])
    
    # which tfs are on
    tfs_on <- tf_info %>% filter(module_id %in% mods)
    molecules_on <- tfs_on %>% select(w, x, y) %>% gather(col, val) %>% pull(val)
    
    # fetch initial state
    new_initial_state <- rowMeans(gold_sim_vectors[[from_]])
    
    # generate nus
    new_nus <- nus
    new_nus[!rownames(new_nus) %in% molecules_on] <- 0
    
    # simulation of gold standard edge
    out <- fastgssa::ssa(
      initial_state = new_initial_state,
      propensity_funs = propensity_funs,
      nu = new_nus,
      final_time = time,
      params = sim_system$parameters,
      method = model$gold_standard_params$ssa_algorithm,
      hardcode_params = TRUE,
      stop_on_neg_state = FALSE,
      verbose = FALSE,
      use_singular_optimisation = model$simulation_params$use_singular_optimisation
    )
    
    meta <- tibble(time = c(head(out$time, -1), !!time), from_, to_)
    counts <- out$state %>% Matrix::Matrix(sparse = TRUE)
    
    end_state <- out$state[nrow(out$state), ] %>% as.matrix()
    
    gold_sim_outputs[[i]] <- lst(meta, counts)
    
    gold_sim_modules[[to_]] <- mods
    
    if (!to_ %in% gold_sim_vectors) {
      gold_sim_vectors[[to_]] <- end_state
    } else {
      gold_sim_vectors[[to_]] <- cbind(gold_sim_vectors[[to_]], end_state)
    }
    
    pbapply::setTimerProgressBar(timer, value = i)
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
