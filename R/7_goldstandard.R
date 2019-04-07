#' @export
goldstandard_default <- function(
  time_per_edge = 2
) {
  lst(
    time_per_edge
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
      mod_diff = module_progression %>% strsplit("\\|"),
      substate = map(mod_diff, seq_along)
    ) %>% 
    select(-module_progression) %>% 
    unnest(mod_diff, substate) %>% 
    mutate(
      mod_diff = strsplit(mod_diff, ","),
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
  time_per_edge <- model$goldstandard_params$time_per_edge
  mod_changes <- model$goldstandard$mod_changes
  sim_system <- model$simulation_system
  
  propensity_funs <-
    sim_system$formulae %>% 
    select(formula_id, formula) %>% 
    deframe
  
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
      time = ((substate - 1) * time_per_edge + t) / time_per_edge / max(substate)
    ) %>% 
    ungroup() %>% 
    select(-i, -substate) %>% 
    select(simulation_i, t, burn, from, to, from_, to_, time, everything())
  
  model
}

.generate_goldstandard_predict_state <- function(model) {
  # should the burn stages be removed earlier?
  # should the names be added to the simulations earlier?
  gold_sims <- model$goldstandard$simulations %>% filter(!burn) %>% mutate(step_id = paste0("GOLD", row_number())) %>% as.data.frame()
  sims <- model$simulations %>% filter(t >= 0) %>% mutate(step_id = paste0("Step", row_number())) %>% as.data.frame()
  
  # fetch meta data and expression data
  meta_tr <- gold_sims %>% select(simulation_i:time, step_id) %>% mutate(edge = factor(paste0(from, "--->", to)))
  expr_tr <- gold_sims %>% select(-simulation_i:-time) %>% column_to_rownames("step_id")
  meta_pr <- sims %>% select(-one_of(colnames(expr_tr)))
  expr_pr <- sims %>% select(one_of(colnames(expr_tr)), step_id) %>% column_to_rownames("step_id")
  
  # predict edges for simulation data
  # todo: use ranger random forest with probabilities and smooth it
  train <- bind_cols(meta_tr %>% select(edge, time), expr_tr) %>% as.data.frame()
  rf <- randomForestSRC::rfsrc(Multivar(edge, time) ~ ., train)
  pred <- predict(rf, expr_pr)
  
  # compute dimred
  dimred <- dyndimred::dimred_landmark_mds(rbind(expr_tr, expr_pr), distance_metric = "angular")
  
  # construct progressions
  progressions <- 
    tibble(
      simulation_i = meta_pr$simulation_i,
      step_id = rownames(expr_pr),
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
    select(simulation_i, step_id, from, to, percentage)
  
  # construct milestone network
  milestone_network <- 
    inner_join(
      meta_tr %>% group_by(from, to) %>% mutate(step_i = row_number()) %>% ungroup() %>% select(step_id_1 = step_id, from, to, step_i),
      meta_tr %>% group_by(from, to) %>% mutate(step_i = row_number() + 1) %>% ungroup() %>% select(step_id_2 = step_id, from, to, step_i),
      by = c("from", "to", "step_i")
    ) %>% 
    mutate(
      dist = sqrt(rowSums((dimred[step_id_1, , drop = FALSE] - dimred[step_id_2, , drop = FALSE])^2))
    ) %>% 
    group_by(from, to) %>% 
    summarise(length = sum(dist), directed = TRUE) %>% 
    ungroup()

  # return output  
  model$goldstandard$dimred <- dimred
  model$goldstandard$counts <- expr_pr
  model$goldstandard$progressions <- progressions
  model$goldstandard$milestone_network <- milestone_network
  
  model
}