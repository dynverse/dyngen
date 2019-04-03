#' @export
goldstandard_default <- function(
  max_path_length = 10,
  reference_length = 100,
  smooth_window = 50
) {
  lst(
    max_path_length,
    reference_length,
    smooth_window
  )
}

generate_goldstandard <- function(model) {
  model %>% 
    .generate_goldstandard_mod_changes() %>% 
    .generate_goldstandard_simulations() %>% 
    .generate_goldstandard_predict_state()
  
}

.generate_goldstandard_mod_changes <- function(model) {
  mod_changes <- 
    model$modulenet$expression_patterns %>% 
    mutate(
      mod_diff = module_progression %>% str_split("\\|"),
      substate = map(mod_diff, seq_along)
    ) %>% 
    select(-module_progression) %>% 
    unnest(mod_diff, substate) %>% 
    mutate(
      mod_diff = str_split(mod_diff, ","),
      mod_on = map(mod_diff, function(x) x %>% keep(grepl("\\+", .)) %>% gsub("\\+", "", .)),
      mod_off = map(mod_diff, function(x) x %>% keep(grepl("-", .)) %>% gsub("-", "", .))
    ) %>% 
    group_by(from, to) %>% 
    mutate(
      from_ = ifelse(row_number() == 1, from, c("", paste0(from, to, letters[row_number() - 1]))),
      to_ = ifelse(row_number() == n(), to, from_[row_number() + 1])
    ) %>% 
    ungroup()

  # order edges
  order <- 
    mod_changes %>%
    select(from_, to_) %>% 
    igraph::graph_from_data_frame() %>% 
    igraph::make_line_graph() %>%
    igraph::topo_sort()
  
  model$goldstandard <- 
    lst(mod_changes = mod_changes[order, ])

  model
}

.generate_goldstandard_simulations <- function(model) {
  mod_changes <- model$goldstandard$mod_changes
  
  # start constructing golds
  gold_sim_outputs <- list()
  gold_sim_vectors <- list()
  gold_sim_modules <- list()
  
  start_state <- mod_changes$from[[1]]
  
  gold_sim_vectors[[start_state]] <- model$simulation_system$initial_state %>% as.matrix
  gold_sim_modules[[start_state]] <- c()
  
  for (i in seq_len(nrow(mod_changes))) {
    from_ <- mod_changes$from_[[i]]
    to_ <- mod_changes$to_[[i]]
    mod_on <- mod_changes$mod_on[[i]]
    mod_off <- mod_changes$mod_off[[i]]
    prev_mods <- gold_sim_modules[[from_]]
    mods <- union(prev_mods, mod_on) %>% setdiff(mod_off)
    
    if (model$verbose) cat("Generating gold for ", from_, "->", to_, "; modules on: ", paste(mods, collapse = ","), "\n", sep = "")
    
    # fetch initial state
    new_initial_state <- rowMeans(gold_sim_vectors[[from_]])
    
    # generate nus
    new_nus <- sim_system$nus
    filt <- colnames(new_nus) %>% purrr:::discard(grepl(paste0("_(", paste(mods, collapse = "|"), ")_"), .))
    new_nus[, filt] <- new_nus[, filt] * 0
    
    # actual simulation
    out <- fastgssa::ssa(
      initial.state = new_initial_state,
      propensity.funs = propensity_funs,
      nu = new_nus %>% as.matrix,
      final.time = sim_params$burn_time,
      parms = sim_system$parameters,
      method = fastgssa::ssa.em(noise_strength = 0),#sim_params$ssa_algorithm,
      recalculate.all = TRUE,
      stop.on.negstate = FALSE,
      stop.on.propensity = FALSE,
      verbose = FALSE
    )
    
    # add both burnin as normal simulation together
    simulation <-
      .simulate_cells_process_ssa(out) %>%
      mutate(i, from_, to_) %>% 
      select(i, from_, to_, everything())
    
    end_state <-
      simulation %>%
      slice(n()) %>%
      select(one_of(sim_system$molecule_ids)) %>%
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
  
  model$goldstandard$simulations <- 
    bind_rows(gold_sim_outputs) %>% 
    left_join(mod_changes %>% select(from, to, from_, to_, substate, burn), by = c("from_", "to_")) %>% 
    group_by(from, to) %>%
    mutate(
      simulation_i = 0,
      time = ((substate - 1) * sim_params$burn_time + t) / sim_params$burn_time / max(substate)
    ) %>% 
    ungroup() %>% 
    select(-i, -substate) %>% 
    select(simulation_i, t, burn, from, to, from_, to_, time, everything())
  
  model
}

.generate_goldstandard_predict_state <- function(model) {
  gold_sims <- model$goldstandard$simulations %>% filter(!burn)
  sims <- model$simulations %>% filter(t >= 0)
  
  meta_tr <- gold_sims %>% select(simulation_i:time) %>% mutate(edge = factor(paste0(from, "--->", to)))
  expr_tr <- gold_sims %>% select(-simulation_i:-time)
  meta_pr <- sims %>% select(-one_of(colnames(expr_tr)))
  expr_pr <- sims %>% select(one_of(colnames(expr_tr)))
  
  rownames(expr_pr) <- paste0("Step", seq_len(nrow(expr_pr)))
  
  train <- bind_cols(meta_tr %>% select(edge, time), expr_tr) %>% as.data.frame()
  
  rf <- randomForestSRC::rfsrc(Multivar(edge, time) ~ ., train)
  pred <- predict(rf, expr_pr)
  
  progressions <- 
    tibble(
      cell_id = rownames(expr_pr),
      edge = pred$classOutput$edge$class,
      time = pred$regrOutput$time$predicted,
      from = gsub("--->.*", "", edge),
      to = gsub(".*--->", "", edge)
    ) %>% 
    group_by(from, to) %>% 
    mutate(
      percentage = dynutils::scale_minmax(time)
    ) %>% 
    ungroup() %>% 
    select(cell_id, from, to, percentage)
  
  milestone_network <- 
    progressions %>% 
    transmute(from, to, length = 1, directed = TRUE) %>% 
    unique()
  
  model$goldstandard$counts <- expr_pr
  model$goldstandard$progressions <- progressions
  model$goldstandard$milestone_network <- milestone_network
  
  model
}