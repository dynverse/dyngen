#' Simulate the gold standard
#' 
#' [generate_gold_standard()] runs simulations in order to determine the gold standard
#' of the simulations.
#' [gold_standard_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the kinetics of the feature network has been generated with [generate_kinetics()].
#' @param tau The time step of the ODE algorithm used to generate the gold standard.
#' @param census_interval A granularity parameter of the gold standard time steps. Should be larger than or equal to `tau`.
#' @param simulate_targets Also simulate the targets during the gold standard simulation
#' 
#' @export
#' 
#' @return A dyngen model.
#' 
#' @seealso [dyngen] on how to run a complete dyngen simulation
#' 
#' @examples
#' model <- 
#'   initialise_model(
#'     backbone = backbone_bifurcating(),
#'     gold_standard = gold_standard_default(tau = .01, census_interval = 1)
#'   )
#' 
#' \donttest{
#' data("example_model")
#' model <- example_model %>% generate_gold_standard()
#'   
#' plot_gold_simulations(model)
#' plot_gold_mappings(model)
#' plot_gold_expression(model)
#' }
generate_gold_standard <- function(model) {
  model$gold_standard <- list()
  
  # compute changes in modules along the edges
  if (model$verbose) cat("Generating gold standard mod changes\n")
  model <- .add_timing(model, "5_gold_standard", "generate mod changes")
  model$gold_standard$mod_changes <- .generate_gold_standard_mod_changes(model$backbone$expression_patterns)
  
  # precompile propensity functions
  if (model$verbose) cat("Precompiling reactions for gold standard\n")
  model <- .add_timing(model, "5_gold_standard", "precompiling reactions for gold standard")
  prep_data <- .generate_gold_precompile_reactions(model)
  
  # run gold standard simulations
  if (model$verbose) cat("Running gold simulations\n")
  model <- .add_timing(model, "5_gold_standard", "running gold simulations")
  simulations <- .generate_gold_standard_simulations(model, prep_data)
  model$gold_standard$meta <- simulations$meta
  model$gold_standard$counts <- simulations$counts
  
  # do combined dimred
  model <- .add_timing(model, "5_gold_standard", "compute dimred")
  model <- model %>% calculate_dimred()
  
  # construct simulation network from dimred
  model <- .add_timing(model, "5_gold_standard", "generate simulation network from dimred")
  model$gold_standard$network <- .generate_gold_standard_generate_network(model)
  
  .add_timing(model, "5_gold_standard", "end")
}

#' @export
#' @rdname generate_gold_standard
gold_standard_default <- function(
  tau = 30 / 3600,
  census_interval = 10 / 60,
  simulate_targets = FALSE
) {
  lst(
    tau,
    census_interval,
    simulate_targets
  )
}

.generate_gold_standard_mod_changes <- function(expression_patterns) {
  # satisfy r cmd check
  module_progression <- mod_diff <- substate <- mod_diff2 <- `.` <- from <- to <- from_ <- NULL
  
  expression_patterns %>% 
    mutate(
      mod_diff = module_progression %>% strsplit("\\|"),
      substate = map(mod_diff, seq_along)
    ) %>% 
    select(-module_progression) %>% 
    unnest(c(mod_diff, substate)) %>% 
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

#' @importFrom GillespieSSA2 compile_reactions
.generate_gold_precompile_reactions <- function(model) {
  # satisfy r cmd check
  is_tf <- mol_premrna <- mol_mrna <- mol_protein <- val <- NULL
  
  
  # fetch paraneters and settings
  sim_system <- model$simulation_system
  tf_info <- model$feature_info 
  if (!model$gold_standard_params$simulate_targets) {
    tf_info <- tf_info %>% filter(is_tf)
  }
  
  # select relevant functions
  tf_molecules <- tf_info %>% select(mol_premrna, mol_mrna, mol_protein) %>% gather(col, val) %>% pull(val)
  
  # filter reactions
  reactions <- sim_system$reactions %>% 
    keep(~ all(names(.$effect) %in% tf_molecules))
  
  buffer_ids <- unique(unlist(map(reactions, "buffer_ids")))
  
  comp_funs <- GillespieSSA2::compile_reactions(
    reactions = reactions,
    buffer_ids = buffer_ids,
    state_ids = tf_molecules,
    params = sim_system$parameters,
    hardcode_params = FALSE,
    fun_by = 1000L
  )
  
  lst(
    tf_molecules,
    reactions = comp_funs
  )
}

#' @importFrom pbapply timerProgressBar setTimerProgressBar
#' @importFrom GillespieSSA2 ssa ode_em
.generate_gold_standard_simulations <- function(model, prep_data) {
  # satisfy r cmd check
  is_tf <- module_id <- mol_premrna <- mol_mrna <- mol_protein <- val <- from <- to <- 
    substate <- burn <- time_per_edge <- simulation_i <- sim_time <- NULL
  
  # fetch paraneters and settings
  mod_changes <- model$gold_standard$mod_changes
  gold_params <- model$gold_standard_params
  sim_system <- model$simulation_system
  tf_info <- model$feature_info
  
  # find out whether to simulate the targets or not (default: no)
  sim_targets <- model$gold_standard_params$simulate_targets
  if (!sim_targets) {
    tf_info <- tf_info %>% filter(is_tf)
  }
  
  # determine ode algorithm
  algo <- GillespieSSA2::ode_em(tau = gold_params$tau, noise_strength = 0)
  
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
  
  if (model$verbose) {
    timer <- pbapply::timerProgressBar(
      min = 0, 
      max = nrow(mod_changes),
      width = 50
    )
  }
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
    tfs_on <- tf_info %>% filter((is.na(module_id) & sim_targets) | module_id %in% mods)
    molecules_on <- tfs_on %>% select(mol_premrna, mol_mrna, mol_protein) %>% gather(col, val) %>% pull(val)
    
    # fetch initial state
    new_initial_state <- rowMeans(gold_sim_vectors[[from_]])
    
    # generate nus
    new_reactions <- reactions
    rem <- setdiff(tf_molecules, molecules_on)
    if (length(rem) > 0) {
      new_reactions$state_change[match(rem, tf_molecules), ] <- 0
    }
    
    # simulation of gold standard edge
    out <- GillespieSSA2::ssa(
      initial_state = new_initial_state,
      reactions = new_reactions,
      final_time = time,
      params = sim_system$parameters,
      method = algo,
      census_interval = gold_params$census_interval,
      stop_on_neg_state = TRUE,
      verbose = FALSE
    )
    
    time_out <- out$time
    state_out <- out$state
    
    meta <- tibble(from_, to_, time = time_out)
    counts <- state_out %>% Matrix::Matrix(sparse = TRUE)
      
    if (model$verbose) pbapply::setTimerProgressBar(timer, value = i)
    
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
  # satisfy r cmd check
  burn <- from <- to <- time <- NULL
  
  model$backbone$expression_patterns %>%
    filter(!burn) %>% 
    transmute(
      from, 
      to,
      length = time / sum(time) * length(time),
      directed = TRUE
    )
}
