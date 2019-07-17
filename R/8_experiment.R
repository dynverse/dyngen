#' Sample cells from the simulations
#' 
#' [generate_experiment()] runs samples cells along the different simulations.
#' [experiment_snapshot()] assumes that cells are sampled from a heterogeneous pool of cells.
#' Cells will thus be sampled uniformily from the trajectory.
#' [experiment_synchronised()] assumes that all the cells are synchronised and 
#' are sampled at different timepoints.
#' 
#' @param model A dyngen intermediary model for which the simulations have been run with [generate_cells()].
#' @param weight_bw \[snapshot\] A bandwidth parameter for determining the distribution of 
#'   cells along each edge in order to perform weighted sampling.
#' @param num_timepoints \[synchronised\] The number of time points used in the experiment.
#' @param pct_between \[synchronised\] The percentage of 'unused' simulation time.
#' 
#' @rdname generate_experiment
#' @export
generate_experiment <- function(model) {
  if (model$verbose) cat("Simulating experiment\n")
  # first sample the cells from the sample, using the desired number of cells
  step_ixs <- .generate_experiment_sample_cells(model) %>% sample()
  
  sim_meta <-
    model$simulations$meta[step_ixs, , drop = FALSE] %>% 
    mutate(
      cell_id = paste0("cell", row_number()),
      step_ix = step_ixs
    ) %>% 
    select(cell_id, step_ix, simulation_i, sim_time, from, to, time)
  
  sim_wcounts <- model$simulations$counts[step_ixs, model$feature_info$w, drop = FALSE]
  rownames(sim_wcounts) <- sim_meta$cell_id
  colnames(sim_wcounts) <- model$feature_info$feature_id
  
  sim_xcounts <- model$simulations$counts[step_ixs, model$feature_info$x, drop = FALSE]
  rownames(sim_xcounts) <- sim_meta$cell_id
  colnames(sim_xcounts) <- model$feature_info$feature_id
  
  sim_ycounts <- model$simulations$counts[step_ixs, model$feature_info$y, drop = FALSE]
  rownames(sim_ycounts) <- sim_meta$cell_id
  colnames(sim_ycounts) <- model$feature_info$feature_id
  
  sim_regulation <- model$simulations$regulation[step_ixs, , drop = FALSE]
  rownames(sim_regulation) <- sim_meta$cell_id
  
  # combine into final count matrix
  model$experiment <- list(
    wcounts = sim_wcounts,
    xcounts = sim_xcounts,
    ycounts = sim_ycounts,
    feature_info =  model$feature_info,
    cell_info = sim_meta,
    regulation = sim_regulation
  )
  
  model
}

#' @export
#' @rdname generate_experiment
#' @importFrom gillespie ssa_etl
list_experiment_samplers <- function() {
  lst(
    snapshot = experiment_snapshot,
    synchronised = experiment_synchronised
  )
}

#' @rdname generate_experiment
#' @export
experiment_snapshot <- function(
  weight_bw = 0.1
) {
  lst(
    fun = .generate_experiment_snapshot,
    weight_bw
  )
}

#' @rdname generate_experiment
#' @export
experiment_synchronised <- function(
  num_timepoints = 8,
  pct_between = .75
) {
  lst(
    fun = .generate_experiment_synchronised,
    num_timepoints,
    pct_between
  )
}

.generate_experiment_sample_cells <- function(model) {
  network <- 
    model$gold_standard$network
  
  end_states <- setdiff(unique(network$to), unique(network$from))
  
  sim_meta <-
    model$simulations$meta %>%
    mutate(orig_ix = row_number()) %>% 
    filter(sim_time >= 0, !to %in% end_states | time < 1)
  
  params <- model$experiment_params
  
  params$fun(
    network = network, 
    sim_meta = sim_meta, 
    params = model$experiment_params,
    num_cells = model$numbers$num_cells
  )
}

#' @importFrom stats approxfun density
.generate_experiment_snapshot <- function(
  network,
  sim_meta,
  params,
  num_cells
) {
  network <- 
    network %>%
    mutate(
      pct = length / sum(length), 
      cum_pct = cumsum(pct),
      num_cells = diff(c(0, round(cum_pct * num_cells)))
    ) %>% 
    select(-pct, -cum_pct)
  
  map(
    seq_len(nrow(network)),
    function(i) {
      edge <- network %>% slice(i)
      meta <- 
        inner_join(sim_meta, edge %>% select(from, to), c("from", "to"))
      
      if (nrow(meta) > 1) {
        meta %>% 
          mutate(
            density = approxfun(density(time, bw = params$weight_bw))(time),
            weight = 1 / density
          ) %>% 
          {sample(.$orig_ix, size = edge$num_cells, replace = TRUE, prob = .$weight)}
      } else {
        NULL
      }
    }
  ) %>% 
    unlist()
}

.generate_experiment_synchronised <- function(
  network,
  sim_meta,
  params,
  num_cells
) {
  sim_meta2 <- 
    sim_meta %>% 
    mutate(
      t_scale = sim_time / (max(sim_time)+1e-10) * params$num_timepoints,
      timepoint_group = floor(t_scale),
      selectable = (t_scale - timepoint_group) < (1 - params$pct_between)
    ) %>% 
    filter(selectable)
  
  numbers <-
    sim_meta2 %>% 
    group_by(timepoint_group) %>% 
    summarise(n = n()) %>% 
    mutate(
      pct = n / sum(n), 
      cum_pct = cumsum(pct),
      num_cells = diff(c(0, round(cum_pct * num_cells)))
    )
  
  map2(
    numbers$timepoint_group,
    numbers$num_cells,
    function(gr, num) {
      sim_meta2 %>% filter(timepoint_group == gr) %>% pull(orig_ix) %>% sample(size = num, replace = TRUE)
    }
  ) %>% 
    unlist()
}

