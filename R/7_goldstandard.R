#' @export
goldstandard_default <- function(
  time_per_edge = 2
) {
  lst(
    time_per_edge
  )
}

#' @export
simulate_goldstandard <- function(model) {
  model$goldstandard <- list()
  
  # compute changes in modules along the edges
  if (model$verbose) cat("Running gold simulations\n")
  model$goldstandard$mod_changes <- .generate_goldstandard_mod_changes(model)
  
  # run gold standard simulations
  simulations <- .generate_goldstandard_simulations(model)
  gs_meta <- simulations %>% select(simulation_i:time)
  gs_counts <- simulations %>% select(-simulation_i:-time) %>% as.matrix %>% Matrix::Matrix(sparse = TRUE)
  model$goldstandard$meta <- gs_meta
  model$goldstandard$counts <- gs_counts
  
  # do combined dimred
  sim_meta <- model$simulations$meta
  sim_counts <- model$simulations$counts
  dimred <- dyndimred::dimred_landmark_mds(Matrix::rbind2(gs_counts, sim_counts), distance_metric = model$simulation_params$dimred_method)
  gs_dimred <- dimred[seq_len(nrow(gs_counts)), , drop = FALSE]
  sim_dimred <- dimred[-seq_len(nrow(gs_counts)), , drop = FALSE]
  
  model$simulations$dimred <- sim_dimred
  model$goldstandard$dimred <- gs_dimred
  
  # predict simulation edge and time with gold standard
  if (model$verbose) cat("Predicting simulation states from gold simulations\n")
  model$simulations$meta <- .generate_goldstandard_predict_state(model)
  
  # construct simulation network from dimred
  model$goldstandard$network <- .generate_goldstandard_generate_network(model)
  
  model
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
  
  mod_changes[order, ]
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
    
    # if (model$verbose) cat("Generating gold for ", from_, "->", to_, "; modules on: ", paste(mods, collapse = ","), "\n", sep = "")
    
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

.generate_goldstandard_predict_state <- function(model) {
  # fetch gold standard data
  gs_meta <- model$goldstandard$meta
  gold_ix <- !gs_meta$burn
  gs_meta <- gs_meta[gold_ix, , drop = FALSE]
  gs_counts <- model$goldstandard$counts[gold_ix, , drop = FALSE]
  gs_dimred <- model$goldstandard$dimred[gold_ix, , drop = FALSE]
  
  # fetch simulation data
  sim_meta <- model$simulations$meta
  sim_counts <- model$simulations$counts
  
  # predict edges for simulation data
  # todo: use ranger random forest with probabilities and smooth it
  train <- data.frame(gs_meta %>% transmute(edge = factor(paste0(from, "--->", to)), time), gs_counts %>% as.matrix)
  options(mc.cores = 1)
  rf <- randomForestSRC::rfsrc(Multivar(edge, time) ~ ., train)
  pred <- predict(rf, sim_counts %>% as.matrix() %>% as.data.frame())
  
  # construct progressions
  sim_meta <- 
    sim_meta %>% mutate(
      edge = pred$classOutput$edge$class,
      time = pred$regrOutput$time$predicted,
      from = gsub("--->.*", "", edge),
      to = gsub(".*--->", "", edge)
    ) %>% 
    group_by(from, to) %>% 
    mutate(
      time = dynutils::scale_minmax(time)
    ) %>% 
    ungroup()
  
  sim_meta
}

.generate_goldstandard_generate_network <- function(model) {
  gs_dimred <- model$goldstandard$dimred
  gs_meta <- model$goldstandard$meta
  
  gs_meta %>% 
    mutate(i = row_number()) %>% 
    group_by(from, to) %>% 
    summarise(
      length = sqrt(sum((gs_dimred[i[-1], , drop = FALSE] - gs_dimred[i[-n()], , drop = FALSE])^2)),
      directed = TRUE
    ) %>% 
    ungroup()
}
