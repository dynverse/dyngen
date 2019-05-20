#' @export
list_experiment_samplers <- function() {
  lst(
    snapshot = experiment_snapshot,
    synchronised = experiment_synchronised
  )
}

#' @export
experiment_snapshot <- function(
  realcount = sample(realcounts$name, 1),
  weight_bw = 0.1
) {
  lst(
    realcount,
    sampler_type = "snapshot",
    weight_bw
  )
}

#' @export
experiment_synchronised <- function(
  realcount = sample(realcounts$name, 1),
  num_timepoints = 8,
  pct_between = .75
) {
  lst(
    realcount,
    sampler_type = "synchronised",
    num_timepoints,
    pct_between
  )
}

#' @importFrom dynutils scale_minmax
#' @importFrom stats rnorm
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
  sim_counts <- model$simulations$counts[step_ixs, model$feature_info$x, drop = FALSE]
  rownames(sim_counts) <- sim_meta$cell_id
  colnames(sim_counts) <- model$feature_info$feature_id
  
  # fetch real expression data
  realcount <- .generate_experiment_fetch_realcount(model)
  
  # map simulated counts to real counts
  trafo_counts <- map(
    sim_meta$cell_id,
    function(cid) {
      realcell <- sample.int(nrow(realcount), size = 1, replace = TRUE)
      real_expr <- sort(realcount[realcell, ])
      orig_expr <- sim_counts[cid, ]
      # add a bit of noise to the percentages to make duplicate cells a bit different
      pct_order <- 
        dynutils::scale_minmax(rank(orig_expr, ties.method = "random")) + 
        stats::rnorm(length(orig_expr), 0, .001) 
      pct_order[pct_order <= 0] <- 0
      pct_order[pct_order >= 1] <- 1
      trafo_expr <- real_expr[floor(pct_order * (length(real_expr) - 1e-10)) + 1]
      trafo_expr %>% set_names(names(orig_expr)) %>% Matrix::Matrix(sparse = TRUE) %>% Matrix::t()
    }
  )
  trafo_count <- do.call(rbind, trafo_counts)
  dimnames(trafo_count) <- dimnames(sim_counts)
  
  # generate housekeeping expression
  num_hks <- model$numbers$num_hks
  if (num_hks > 0) {
    hk_count <- 
      realcount[
        sample.int(nrow(realcount), nrow(trafo_count), replace = TRUE),
        sample.int(ncol(realcount), num_hks, replace = TRUE),
        drop = FALSE
        ]
    rownames(hk_count) <- rownames(trafo_count)
    colnames(hk_count) <- paste0("HK", seq_len(num_hks))
    hk_info <- tibble(feature_id = colnames(hk_count), is_tf = FALSE, is_hk = TRUE)
  } else {
    hk_count <- NULL
    hk_info <- NULL
  }
  
  # combine into final count matrix
  model$experiment <- list(
    counts = cbind(trafo_count, hk_count),
    feature_info = bind_rows(
      model$feature_info,
      hk_info
    ),
    cell_info = sim_meta
  )
  
  model
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
  
  if (params$sampler_type == "snapshot") {
    network <- 
      network %>%
      mutate(
        pct = length / sum(length), 
        cum_pct = cumsum(pct),
        num_cells = diff(c(0, round(cum_pct * model$numbers$num_cells)))
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
  } else if (params$sampler_type == "synchronised") {
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
        num_cells = diff(c(0, round(cum_pct * model$numbers$num_cells)))
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
}

.generate_experiment_fetch_realcount <- function(model) {
  realcount_ <- model$experiment_params$realcount
  
  # download realcount if needed-
  realcount <- 
    if (is.character(realcount_)) {
      data(realcounts, package = "dyngen", envir = environment())
      assert_that(realcount_ %all_in% realcounts$name)
      
      url <- realcounts$url[[match(realcount_, realcounts$name)]]
      
      .download_cacheable_file(url, model)
    } else if (is.matrix(realcount_) || is_sparse(realcount_)) {
      realcount_
    } else {
      stop("realcount should be a url from dyngen::realcounts, or a sparse count matrix.")
    }
  
  realcount
}
