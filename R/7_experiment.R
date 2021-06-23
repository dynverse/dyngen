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
#' @param map_reference_cpm Whether or not to try to match the CPM distribution to that of a reference dataset.
#' @param map_reference_ls Whether or not to try to match the distribution of the library sizes to that of the reference dataset.
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
#' 
#' \dontrun{
#' data("example_model")
#' model <- example_model %>% generate_experiment() 
#' 
#' plot_experiment_dimred(model)
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
  rownames(tsim_counts) <- cell_info$cell_id
  
  mol_info <- model$feature_info %>%
    gather("mol", "val", .data$mol_premrna, .data$mol_mrna, .data$mol_protein) %>% 
    select(.data$feature_id, .data$mol, .data$val)
  
  # fetch real expression data
  model <- .add_timing(model, "7_experiment", "fetch realcount")
  realcount <- .generate_experiment_fetch_realcount(model)
  
  # simulate library size variation from real data
  model <- .add_timing(model, "7_experiment", "simulate library size variation")
  
  # process mrna (both premrna and mature)
  mrna_ids <- mol_info %>% filter(.data$mol != "mol_protein") %>% pull(.data$val)
  tsim_counts_mrna <- tsim_counts[, mrna_ids, drop = FALSE]
  count_mrna_simulation <- .simulate_counts_from_realcounts(
    tsim_counts = tsim_counts_mrna,
    realcount = realcount,
    map_reference_cpm = model$experiment_params$map_reference_cpm,
    map_reference_ls = model$experiment_params$map_reference_ls
  )
  
  # process proteins
  # TODO: should use real protein dataset for this
  prot_ids <- mol_info %>% filter(.data$mol == "mol_protein") %>% pull(.data$val)
  # count_prot_simulation <- .simulate_counts_from_realcounts(tsim_counts[, prot_ids, drop = FALSE], realcount, ...)
  count_prot_simulation <- tsim_counts[, prot_ids, drop = FALSE]
  
  # split up molecules
  model <- .add_timing(model, "7_experiment", "create output")
  sim_wcounts <- count_mrna_simulation[, model$feature_info$mol_premrna, drop = FALSE]
  sim_xcounts <- count_mrna_simulation[, model$feature_info$mol_mrna, drop = FALSE]
  sim_ycounts <- count_prot_simulation[, model$feature_info$mol_protein, drop = FALSE]
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
  map_reference_cpm = TRUE,
  map_reference_ls = TRUE
) {
  # calculate lib sizes and cpms
  realcount_ls <- Matrix::rowSums(realcount)
  realcount_cpm <- realcount
  realcount_cpm@x <- realcount@x / realcount_ls[realcount@i+1] * 1e6
  
  tsim_counts_ls <- Matrix::rowSums(tsim_counts)
  tsim_counts_cpm <- tsim_counts
  tsim_counts_cpm@x <- tsim_counts@x / tsim_counts_ls[tsim_counts@i+1] * 1e6
  
  # map real density on tsim counts
  tsim_counts_cpm_new <- tsim_counts_cpm
  
  if (map_reference_cpm) {
    # tsim_counts_cpm_new@x <- stats::quantile(realcount_cpm@x, dplyr::percent_rank(tsim_counts_cpm@x))
    tsim_counts_cpm_new@x <- stats::quantile(realcount_cpm@x, sort(stats::runif(length(tsim_counts_cpm@x)))[order(order(tsim_counts_cpm@x))]) %>% unname
  }
  
  # simulate library size variation from real data
  if (map_reference_ls) {
    # lib_size_new <- round(stats::quantile(realcount_ls, dplyr::percent_rank(tsim_counts_ls)))
    lib_size_new <- round(stats::quantile(realcount_ls, sort(stats::runif(length(tsim_counts_ls)))[order(order(tsim_counts_ls))])) %>% unname
  } else {
    lib_size_new <- tsim_counts_ls
  }
  
  # simulate sampling of molecules
  tsim_counts_t <- Matrix::t(tsim_counts_cpm_new)
  new_vals <- unlist(map(seq_along(lib_size_new), function(cell_i) {
    pi <- tsim_counts_t@p[[cell_i]]
    pj <- tsim_counts_t@p[[cell_i + 1]]
    if (pi != pj) {
      pix <- seq(pi + 1, pj, 1)
      gene_is <- tsim_counts_t@i[pix] + 1
      gene_vals <- tsim_counts_t@x[pix]
      lib_size <- lib_size_new[[cell_i]]
      
      # sample 'lib_size' molecules for each of the genes, weighted by 'gene_vals'
      rmultinom(1, lib_size, gene_vals)
    } else {
      integer(0)
    }
  })) %>% as.numeric
  
  tsim_counts_t@x <- new_vals
  sim_counts <- Matrix::drop0(Matrix::t(tsim_counts_t))
  
  sim_counts
}

#' @export
#' @rdname generate_experiment
list_experiment_samplers <- function() {
  lst(
    snapshot = experiment_snapshot,
    synchronised = experiment_synchronised
  )
}

.sample_decent_realcount <- function() {
  # satisfy r cmd check
  realcounts <- NULL
  
  data(realcounts, package = "dyngen", envir = environment())
  
  realcounts %>% 
    filter(.data$qc_pass) %>%
    pull(.data$name) %>% 
    sample(1)
}

#' @rdname generate_experiment
#' @export
experiment_snapshot <- function(
  realcount = NULL,
  map_reference_cpm = TRUE,
  map_reference_ls = TRUE,
  weight_bw = 0.1
) {
  if (is.null(realcount)) {
    realcount <- .sample_decent_realcount()
  }
  lst(
    realcount,
    map_reference_cpm,
    map_reference_ls,
    fun = .generate_experiment_snapshot,
    weight_bw
  )
}

#' @rdname generate_experiment
#' @export
experiment_synchronised <- function(
  realcount = NULL,
  map_reference_cpm = TRUE,
  map_reference_ls = TRUE,
  num_timepoints = 8,
  pct_between = .75
) {
  if (is.null(realcount)) {
    realcount <- .sample_decent_realcount()
  }
  lst(
    realcount,
    map_reference_cpm,
    map_reference_ls,
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
  steps_summ <- sim_meta %>% group_by(.data$from, .data$to) %>% summarise(num_steps = n(), .groups = "drop")
  
  if (num_cells > sum(steps_summ$num_steps)) {
    stop(
      "Simulations only generated ", sum(steps_summ$num_steps), " different datapoints, whereas ", num_cells, " cells were requested.\n",
      "Increase the number of simulations (see `?simulation_default`) or decreasing the census interval (see `?simulation_default`)."
    )
  }
  if (network %>% anti_join(steps_summ, by = c("from", "to")) %>% nrow() > 0) {
    warning(
      "Certain backbone segments are not covered by any of the simulations. If this is intentional, please ignore this warning.\n",
      "  Otherwise, increase the number of simulations (see `?simulation_default`) or decreasing the census interval (see `?simulation_default`)."
    )
    steps_summ <- steps_summ %>% filter(.data$num_steps > 0)
  } 
  
  numbers_per_edge <- 
    inner_join(
      network,
      steps_summ, 
      by = c("from", "to")
    ) %>% 
    mutate(
      pct = .data$length / sum(.data$length), 
      cum_pct = cumsum(.data$pct),
      num_cells = diff(c(0, round(.data$cum_pct * !!num_cells)))
    ) %>% 
    mutate(num_steps = ifelse(is.na(.data$num_steps), 0, .data$num_steps)) %>% 
    select(-.data$pct, -.data$cum_pct)
  
  if (any(numbers_per_edge$num_cells > numbers_per_edge$num_steps)) {
    warning(
      "One of the branches did not obtain enough simulation steps.\n",
      "  Increase the number of simulations (see `?simulation_default`) or decreasing the census interval (see `?simulation_default`)."
    )
  }
  
  map(
    seq_len(nrow(numbers_per_edge)),
    function(i) {
      edge <- numbers_per_edge %>% slice(i)
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
  
  # download realcount if needed
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
