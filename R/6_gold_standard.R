#' @export
gold_standard_default <- function(
  time_per_edge = 2
) {
  lst(
    time_per_edge
  )
}

#' @export
generate_gold_standard <- function(model) {
  model$gold_standard <- list()
  
  # compute changes in modules along the edges
  if (model$verbose) cat("Running gold simulations\n")
  model$gold_standard$mod_changes <- .generate_gold_standard_mod_changes(model)
  
  # run gold standard simulations
  simulations <- .generate_gold_standard_simulations(model)
  model$gold_standard$meta <- simulations %>% select(simulation_i:time)
  model$gold_standard$counts <- simulations %>% select(-simulation_i:-time) %>% as.matrix %>% Matrix::Matrix(sparse = TRUE)
  
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
      mod_off = map(mod_diff2, function(x) x %>% keep(grepl("-", .)) %>% gsub("-", "", .)),
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

.generate_gold_standard_simulations <- function(model) {
  # fetch paraneters and settings
  time_per_edge <- model$gold_standard_params$time_per_edge
  mod_changes <- model$gold_standard$mod_changes
  sim_system <- model$simulation_system
  tf_info <- model$feature_info %>% filter(is_tf)
  
  # select relevant functions
  tf_molecules <- tf_info %>% select(x, y) %>% gather(col, val) %>% pull(val)
  
  # filter propensity functions
  propensity_funs <-
    sim_system$formulae %>% 
    filter(molecule %in% tf_molecules) %>% 
    select(formula_id, formula) %>% 
    deframe
  
  # filter nus
  nus <- sim_system$nus[tf_molecules, names(propensity_funs), drop = FALSE]
  
  # start constructing golds
  gold_sim_outputs <- list()
  gold_sim_vectors <- list()
  gold_sim_modules <- list()
  
  start_state <- mod_changes$from[[1]]
  
  gold_sim_vectors[[start_state]] <- model$simulation_system$initial_state[tf_molecules] %>% as.matrix
  gold_sim_modules[[start_state]] <- c()
  
  # model$gold_standard$mod_changes %>% mutate(mod_on_str = map_chr(mod_on, paste, collapse = ","), mod_off_str = map_chr(mod_off, paste, collapse = ","))
  
  for (i in seq_len(nrow(mod_changes))) {
    from_ <- mod_changes$from_[[i]]
    to_ <- mod_changes$to_[[i]]
    
    # which modules are on
    mods <- 
      gold_sim_modules[[from_]] %>% 
      union(mod_changes$mod_on[[i]]) %>% 
      setdiff(mod_changes$mod_off[[i]])
    
    # which tfs are on
    tfs_on <- tf_info %>% filter(module_id %in% mods)
    molecules_on <- tfs_on %>% select(x, y) %>% gather(col, val) %>% pull(val)
    
    # fetch initial state
    new_initial_state <- rowMeans(gold_sim_vectors[[from_]])
    
    # generate nus
    new_nus <- nus
    new_nus[!rownames(new_nus) %in% molecules_on] <- 0
    
    # simulation of gold standard edge
    out <- fastgssa::ssa(
      initial.state = new_initial_state,
      propensity.funs = propensity_funs,
      nu = new_nus %>% as.matrix,
      final.time = time_per_edge,
      parms = sim_system$parameters,
      method = fastgssa::ssa.em(noise_strength = 0),
      recalculate.all = TRUE,
      stop.on.negstate = FALSE,
      stop.on.propensity = FALSE,
      verbose = FALSE
    )
    
    # add both burnin as normal simulation together
    simulation <-
      .generate_cells_process_ssa(out) %>%
      mutate(i, from_, to_) %>% 
      select(i, from_, to_, everything())
    
    end_state <-
      simulation %>%
      slice(n()) %>%
      select(one_of(tf_molecules)) %>%
      .[1, , drop = TRUE] %>%
      unlist() %>%
      as.matrix()
    
    gold_sim_outputs[[i]] <- simulation
    
    gold_sim_modules[[to_]] <- mods
    
    if (!to_ %in% gold_sim_vectors) {
      gold_sim_vectors[[to_]] <- end_state
    } else {
      gold_sim_vectors[[to_]] <- cbind(gold_sim_vectors[[to_]], end_state)
    }
  }
  
  
  bind_rows(gold_sim_outputs) %>% 
    left_join(mod_changes %>% select(from, to, from_, to_, substate, burn), by = c("from_", "to_")) %>% 
    group_by(from, to) %>%
    mutate(
      simulation_i = 0,
      time = ((substate - 1) * time_per_edge + t) / time_per_edge / max(substate)
    ) %>% 
    ungroup() %>% 
    select(-i, -substate) %>% 
    select(simulation_i, t, burn, from, to, from_, to_, time, everything())
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
