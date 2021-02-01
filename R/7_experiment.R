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
#' @param realcount The name of a dataset in [realcounts]. If `NULL`, a random
#'   dataset will be sampled from [realcounts].
#' @param sample_capture_rate A function that samples values for the simulated capture rates of genes.
#' 
#' @return A dyngen model.
#' 
#' @rdname generate_experiment
#' 
#' @importFrom stats rmultinom rnorm quantile
#' @export
#' 
#' @examples
#' names(list_experiment_samplers())
#' 
#' model <- 
#'   initialise_model(
#'     backbone = backbone_bifurcating(),
#'     experiment = experiment_synchronised()
#'   )
#' \dontshow{
#' # actually use a smaller example 
#' # to reduce execution time during
#' # testing of the examples
#' model <- initialise_model(
#'   backbone = model$backbone,
#'   num_cells = 5,
#'   num_targets = 0,
#'   num_hks = 0,
#'   gold_standard_params = gold_standard_default(census_interval = 1, tau = 0.1),
#'   simulation_params = simulation_default(
#'     burn_time = 10,
#'     total_time = 10,
#'     census_interval = 1,
#'     ssa_algorithm = ssa_etl(tau = 0.1),
#'     experiment_params = simulation_type_wild_type(num_simulations = 1)
#'   )
#' )
#' }
#' \donttest{
#' model <- 
#'   model %>%
#'   generate_tf_network() %>%
#'   generate_feature_network() %>%
#'   generate_kinetics() %>%
#'   generate_gold_standard() %>%
#'   generate_cells() %>%
#'   generate_experiment() 
#'   
#' dataset <- as_dyno(model)
#' }
generate_experiment <- function(model) {
  if (model$verbose) cat("Simulating experiment\n")
  model <- .add_timing(model, "7_experiment", "sample cells")
  # first sample the cells from the sample, using the desired number of cells
  step_ixs <- .generate_experiment_sample_cells(model)
  
  cell_info <-
    model$simulations$meta[step_ixs, , drop = FALSE] %>% 
    mutate(
      step_ix = step_ixs
    )
    
  if ("group" %in% names(attributes(step_ixs))) {
    cell_info$cell_group <- attr(step_ixs, "group")
  } else {
    # only shuffle if step sampling is not used
    cell_info <- cell_info %>% sample_n(n())
  }
  
  cell_info <-
    cell_info %>% 
    mutate(
      cell_id = paste0("cell", row_number())
    ) %>% 
    select(
      .data$cell_id, 
      .data$step_ix, 
      .data$simulation_i, 
      .data$sim_time, 
      .data$from, 
      .data$to, 
      .data$time, 
      everything()
    )
  
  step_ixs <- cell_info$step_ix
  
  # collect true simulated counts of sampled cells
  tsim_counts <- model$simulations$counts[step_ixs, , drop = FALSE]
  
  # fetch real expression data
  model <- .add_timing(model, "7_experiment", "fetch realcount")
  realcount <- .generate_experiment_fetch_realcount(model)
  
  # simulate library size variation from real data
  model <- .add_timing(model, "7_experiment", "simulate library size variation")
  count_simulation <- .simulate_counts_from_realcounts(tsim_counts, realcount, cell_info, sample_capture_rate = model$experiment_params$sample_capture_rate)
  sim_counts <- count_simulation$sim_counts
  cell_info <- count_simulation$cell_info
  mol_info <- count_simulation$mol_info
  
  # split up molecules
  model <- .add_timing(model, "7_experiment", "create output")
  sim_wcounts <- sim_counts[, model$feature_info$mol_premrna, drop = FALSE]
  sim_xcounts <- sim_counts[, model$feature_info$mol_mrna, drop = FALSE]
  sim_ycounts <- sim_counts[, model$feature_info$mol_protein, drop = FALSE]
  dimnames(sim_wcounts) <- 
    dimnames(sim_xcounts) <-
    dimnames(sim_ycounts) <- 
    list(cell_info$cell_id, model$feature_info$feature_id)
  
  if (model$simulation_params$compute_cellwise_grn) {
    sim_cellwise_grn <- model$simulations$cellwise_grn[step_ixs, , drop = FALSE]
    rownames(sim_cellwise_grn) <- cell_info$cell_id
  } else {
    sim_cellwise_grn <- NULL
  }
  
  if (model$simulation_params$compute_rna_velocity) {
    sim_rna_velocity <- model$simulations$rna_velocity[step_ixs, , drop = FALSE]
    rownames(sim_rna_velocity) <- cell_info$cell_id
  } else {
    sim_rna_velocity <- NULL
  }
  
  # combine into final count matrix
  model$experiment <- list(
    counts_premrna = sim_wcounts,
    counts_mrna = sim_xcounts,
    counts_protein = sim_ycounts,
    feature_info =  model$feature_info,
    cell_info = cell_info,
    cellwise_grn = sim_cellwise_grn,
    rna_velocity = sim_rna_velocity
  )
  
  .add_timing(model, "7_experiment", "end")
}


.simulate_counts_from_realcounts <- function(
  tsim_counts, 
  realcount, 
  cell_info = tibble(cell_id = rownames(tsim_counts)), 
  sample_capture_rate = function(n) rnorm(n, 1, 0.05) %>% pmax(0)
) {
  # simulate library size variation from real data
  realsums <- Matrix::rowSums(realcount)
  dist_vals <- realsums / mean(realsums)
  lib_size <- quantile(dist_vals, runif(nrow(cell_info)))
  
  cell_info <-
    cell_info %>% 
    mutate(
      num_molecules = Matrix::rowSums(tsim_counts),
      mult = quantile(dist_vals, runif(n())),
      lib_size = sort(round(mean(.data$num_molecules) * .data$mult))[order(order(.data$num_molecules))]
    )
  
  # simulate gene capture variation
  mol_info <- 
    tibble(
      id = colnames(tsim_counts),
      capture_rate = sample_capture_rate(ncol(tsim_counts))
    )
  
  # simulate sampling of molecules
  tsim_counts_t <- Matrix::t(tsim_counts)
  new_vals <- unlist(map(seq_len(nrow(cell_info)), function(cell_i) {
    pi <- tsim_counts_t@p[[cell_i]] + 1
    pj <- tsim_counts_t@p[[cell_i + 1]]
    gene_is <- tsim_counts_t@i[pi:pj] + 1
    gene_vals <- tsim_counts_t@x[pi:pj]
    
    lib_size <- cell_info$lib_size[[cell_i]]
    cap_rates <- mol_info$capture_rate[gene_is]
    
    probs <- cap_rates * gene_vals
    probs[probs < 0] <- 0 # sometimes these can become zero due to rounding errors
    
    rmultinom(1, lib_size, probs)
  })) %>% as.numeric
  tsim_counts_t@x <- new_vals
  sim_counts <- Matrix::drop0(Matrix::t(tsim_counts_t))
  
  lst(
    sim_counts,
    cell_info,
    mol_info
  )
}

#' @export
#' @rdname generate_experiment
list_experiment_samplers <- function() {
  lst(
    snapshot = experiment_snapshot,
    synchronised = experiment_synchronised
  )
}

#' @rdname generate_experiment
#' @export
experiment_snapshot <- function(
  realcount = NULL,
  sample_capture_rate = function(n) rnorm(n, 1, .05) %>% pmax(0),
  weight_bw = 0.1
) {
  # satisfy r cmd check
  realcounts <- NULL
  
  if (is.null(realcount)) {
    data(realcounts, package = "dyngen", envir = environment())
    realcount <- sample(realcounts$name, 1)
  }
  lst(
    realcount,
    sample_capture_rate,
    fun = .generate_experiment_snapshot,
    weight_bw
  )
}

#' @rdname generate_experiment
#' @export
experiment_synchronised <- function(
  realcount = NULL,
  sample_capture_rate = function(n) rnorm(n, 1, .05) %>% pmax(0),
  num_timepoints = 8,
  pct_between = .75
) {
  # satisfy r cmd check
  realcounts <- NULL
  
  if (is.null(realcount)) {
    data(realcounts, package = "dyngen", envir = environment())
    realcount <- sample(realcounts$name, 1)
  }
  lst(
    realcount,
    sample_capture_rate,
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
    filter(.data$sim_time >= 0, !.data$to %in% end_states | .data$time < 1)
  
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
      pct = .data$length / sum(.data$length), 
      cum_pct = cumsum(.data$pct),
      num_cells = diff(c(0, round(.data$cum_pct * num_cells)))
    ) %>% 
    select(-.data$pct, -.data$cum_pct)
  
  map(
    seq_len(nrow(network)),
    function(i) {
      edge <- network %>% slice(i)
      meta <- 
        inner_join(sim_meta, edge %>% select(.data$from, .data$to), c("from", "to"))
      
      if (nrow(meta) > 1) {
        met <- meta %>% 
          mutate(
            density = approxfun(density(.data$time, bw = params$weight_bw))(.data$time),
            weight = 1 / .data$density
          )
        
        sample(met$orig_ix, size = edge$num_cells, replace = TRUE, prob = met$weight)
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
      t_scale = .data$sim_time / (max(.data$sim_time)+1e-10) * params$num_timepoints,
      timepoint_group = floor(.data$t_scale),
      selectable = (.data$t_scale - .data$timepoint_group) < (1 - params$pct_between)
    ) %>% 
    filter(.data$selectable)
  
  numbers <-
    sim_meta2 %>% 
    group_by(.data$timepoint_group) %>% 
    summarise(n = n()) %>% 
    mutate(
      pct = n / sum(n), 
      cum_pct = cumsum(.data$pct),
      num_cells = diff(c(0, round(.data$cum_pct * num_cells)))
    )
  
  map2(
    numbers$timepoint_group,
    numbers$num_cells,
    function(gr, num) {
      sim_meta2 %>% 
        filter(.data$timepoint_group == gr) %>%
        pull(.data$orig_ix) %>% 
        sample(size = num, replace = TRUE)
    }
  ) %>% 
    unlist()
}

.generate_experiment_fetch_realcount <- function(model) {
  # satisfy r cmd check
  realcounts <- NULL
  
  realcount_ <- model$experiment_params$realcount
  
  # download realcount if needed-
  realcount <- 
    if (is.character(realcount_)) {
      data(realcounts, package = "dyngen", envir = environment())
      assert_that(realcount_ %all_in% realcounts$name)
      
      url <- realcounts$url[[match(realcount_, realcounts$name)]]
      
      .download_cacheable_file(
        url = url, 
        cache_dir = model$download_cache_dir, 
        verbose = model$verbose
      )
    } else if (is.matrix(realcount_) || is_sparse(realcount_)) {
      realcount_
    } else {
      stop("realcount should be a url from dyngen::realcounts, or a sparse count matrix.")
    }
  
  realcount
}