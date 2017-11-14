#' Preprocessing of simulation prior to gold standard determination
#' 1. Smoothing of expression
preprocess_simulation_for_gs <- function(simulation, model) {
  simulation$expression_normalized <- t(dynutils::scale_quantile(t(simulation$expression), outlier_cutoff=0.05))
  simulation$expression_modules <- simulation$expression_normalized %>% t %>% as.data.frame() %>% split(model$geneinfo$module_id) %>% map(~apply(., 2, mean)) %>% do.call(rbind, .) %>% t
  
  simulation
}





divide_simulation_roughly <- function(simulation, states) {
  rough <- list()
  window <- 10
  
  max_before_breakpoint_expression <- 0.1
  min_after_breakpoint_expression <- 0.5
  
  split_min_difference <- 0.25
  convergence_max_difference <- 0.25
  
  simulation_ids <- unique(simulation$stepinfo$simulation_id)
  
  simulation_id <- simulation_ids[[1]]
  for (simulation_id in simulation_ids) {
    # extract module expression of the simulation
    expression <- simulation$expression_modules[simulation$stepinfo$simulation_id == simulation_id,]
    
    # start state
    curstate = states[[1]]
    
    await_low_expression = FALSE
    
    # when circular (at the start already), await for low expression of the given module
    if(!is.null(curstate$to) && curstate$to[[1]] == curstate$state_id) {
      await_low_expression = TRUE
      await_low_expression_module = curstate$modules
    }
    
    # initial state
    row <- list(
      simulation_id = simulation_id,
      state_id = curstate$state_id,
      prevstate_id = NA,
      start = 1
    )
    rough <- c(rough, list(row))
    
    is <- seq(1, nrow(expression)-window)
    
    i <- is[[1]]
    for (i in is) {
      next_expression <- expression[i:(i+window), curstate$modules, drop=TRUE]
      
      statechange <- FALSE # whether the state has changed
      
      if(await_low_expression) {
        low_expression_module_expression = expression[i:(i+window), curmodules, drop=TRUE]
        if(mean(low_expression_module_expression) < max_before_breakpoint_expression) {
          await_low_expression = FALSE
        }
      } else if(curstate$type == "split") {
        if (verbose) print("splitting...")
        average_next_module_expression = colMeans(next_expression)
        if(max(average_next_module_expression) - sort(average_next_module_expression, TRUE)[[2]] > split_min_difference) {
          nextstate_id <- curstate$to[which.max(average_next_module_expression)]
          statechange <- TRUE
        }
      } else if (curstate$type == "convergence") {
        # change immediatly
        nextstate_id <- curstate$to[[1]]
        statechange <- TRUE
      } else if (curstate$type == "cycle") {
        if(mean(next_expression) > min_after_breakpoint_expression) {
          nextstate_id <- curstate$to[[1]]
          statechange <- TRUE
          
          await_low_expression = TRUE
        }
      }
      
      if(statechange) {
        row <- list(
          simulation_id = simulation_id,
          state_id = nextstate_id,
          prevstate_id = curstate$state_id,
          start = i
        )
        rough <- c(rough, list(row))
        
        nextstate <- states %>% keep(~.$state_id == nextstate_id) %>% first()
        if(!is.null(nextstate)) {
          curstate <- nextstate
        } else {
          break()
        }
      }
    }
  }
  
  rough <- rough %>% bind_rows() %>% ungroup()
}





lagpad <- function(x, k=1) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  } else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}



extract_reference_from_rough_division <- function(simulation, states, rough) {
  # calculate the end of each rough piece
  # number of steps in each simulation, to calculate the end at the end of simulation
  simulation_nsteps <- simulation$stepinfo %>% group_by(simulation_id) %>% summarise(n=n()) %>% arrange(simulation_id) %>% pull(n)
  
  # actual end calculation, by lagging the start of the next piece
  rough <- rough %>% group_by(simulation_id) %>% 
    mutate(end = lagpad(start, -1)) %>% 
    mutate(end = ifelse(is.na(end), simulation_nsteps[simulation_id], end)) %>% 
    ungroup()
  
  # calculate nextstate_ids
  rough <- rough %>% group_by(simulation_id) %>% 
    mutate(nextstate_id = lagpad(state_id, -1), nextend = lagpad(end, -1)) %>% ungroup()
  
  references <- list()
  
  # now adapt start and end steps for each simulation to match those in the combined simulation expression matrix
  simulation_nsteps_cumsum <- simulation_nsteps %>% cumsum() %>% head(-1) %>% c(0, .)
  rough <- rough %>% mutate(
    start_cumulative = start + simulation_nsteps_cumsum[simulation_id],
    end_cumulative = end + simulation_nsteps_cumsum[simulation_id]
  )
  
  # now extract expression
  references <- rough %>% 
    rowwise() %>%
    mutate(
      expression = map2(start_cumulative, end_cumulative, ~simulation$expression_modules[.x:.y, ])
    ) %>% 
    ungroup()
  
  references
}


#' This function looks to make the rough estimates of splits and convergences more accurate
#' A bifurcation will be moved forward, a convergence will be moved backward
divide_simulation_precise <- function(simulation, states, rough) {
  pieces <- list()
  window <- 10
  
  split_pval_cutoff <- 0.05
  convergence_pval_cutoff <- 0.05
  
  is <- seq_len(nrow(rough))
  for (i in is) {
    roughrow <- dynutils::extract_row_to_list(rough, i)
    
    # add start if necessary
    if (is.na(roughrow$prevstate_id)) {
      piece <- list(
        start = 1,
        state_id = roughrow$state_id,
        simulation_id = roughrow$simulation_id
      )
      pieces <- c(pieces, list(piece))
    }
    
    curstate <- states %>% keep(~.$state_id == roughrow$state_id) %>% first
    if (!is.null(curstate)) {
      ## SPLITTING
      if(curstate$type == "split") {
        reference <- references %>% filter(state_id == curstate$state_id) %>% filter(simulation_id != roughrow$simulation_id)
        cur_expression <- simulation$expression_modules[roughrow$start_cumulative:roughrow$end_cumulative, curstate$modules]
        
        # move end backward, check where module expression difference was first significant
        for (j in seq(roughrow$end-roughrow$start+1, 1, -1)) {
          # get minimal distance of the expression at this step, compared with every reference
          min_distances <- map_dbl(reference$expression, function(reference_expression) {
            as.matrix(pdist(reference_expression[, curstate$modules], cur_expression[j, ])) %>% min()
          })
          
          # split according to the same nextstate_id or different nextstate_id
          same_state_distances <- min_distances[reference$nextstate_id == roughrow$nextstate_id]
          diff_state_distances <- min_distances[reference$nextstate_id != roughrow$nextstate_id]
          
          pval <- wilcox.test(same_state_distances, diff_state_distances)$p.value
          
          if(pval > split_pval_cutoff | j == 1) {
            piece <- list(
              start = roughrow$start + j,
              state_id = roughrow$nextstate_id,
              simulation_id = roughrow$simulation_id
            )
            pieces <- c(pieces, list(piece))
            break()
          }
        }
      ## CONVERGING
      } else if (curstate$type == "convergence") {
        stop("wat")
        
        reference <- references %>% filter(state_id %in% c(curstate$state_id, curstate$cross)) %>% filter(simulation_id != roughrow$simulation_id)
        
        cur_expression <- simulation$expression_modules[roughrow$start_cumulative:roughrow$end_cumulative, curstate$modules]
        
        # move start forward, check where module expression difference was first significant
        for (j in seq(1, roughrow$end - roughrow$start + 1, 1)) {
          # get minimal distance of the expression at this step, compared with every reference
          min_distances <- map_dbl(reference$expression, function(reference_expression) {
            as.matrix(pdist(reference_expression[, curstate$modules], cur_expression[j, ])) %>% min()
          })
          
          # split according to the same nextstate_id or different nextstate_id
          same_state_distances <- min_distances[reference$state_id == roughrow$state_id]
          diff_state_distances <- min_distances[reference$state_id != roughrow$state_id]
          
          pval <- wilcox.test(same_state_distances, diff_state_distances)$p.value
          
          print(pval)
          
          if(pval > convergence_pval_cutoff | j == roughrow$end - roughrow$start + 1) {
            print(glue::glue("converged {j} pval {pval}"))
            piece <- list(
              start = roughrow$start + j,
              state_id = roughrow$nextstate_id,
              simulation_id = roughrow$simulation_id
            )
            pieces <- c(pieces, list(piece))
            break()
          }
        }
      }
    }
  }
  
  pieces <- bind_rows(pieces)
}









extract_goldstandard <- function(simulation, model) {
  # preprocessing, by calculating mean expression of every module
  simulation <- preprocess_simulation_for_gs(simulation, model)
  
  
}