extract_tobecomes <- function(simulations, states) {
  quant_scale_combined = get_combined_quant_scale(simulations)
  
  ############################################
  ### This analysis is split in two parts:
  ### 1. Go over the simulation and look forward to where we are absolutely certain that a cell has chosen a certain path at the next bifurcation. Keep on looking for each bifurcation. Save the expression just before we are certain where the cell will go...
  ### 2. This expression is then used as a reference to check for every simulation step at which point the cell as almost certain it will choose a path (ie. no nearby simulations which chose an alternative path). When the cell has decide, we can put a cut.
  ###    In the second step, we also cut the expression for strictly circular trajectories, by checking whether the expression of a specified module is high again, after a period of being low
  
  tobecomes = list()
  minfoldchange = 10
  window = 10
  
  max_before_breakpoint_expression = 0.1
  min_after_breakpoint_expression = 0.5
  for (simulation_id in seq_along(simulations)) {
    print(simulation_id)
    
    expression_modules_scaled = quant_scale_combined(simulations[[simulation_id]]$expression_modules)
    
    curstate = states[[1]]
    
    nextstate_id = NA
    laststate_id = NA
    
    prevsplitstart = 1
    lastsplit = 1
    
    prevstate_id = NA
    
    prevsplits = numeric() # for each state, keep record of the previous step_id, so that we do not use the same expression twice for the same bifurcation (possible leading to a different subsequent state)
    
    await_low_expression = FALSE
    
    if(!is.null(curstate$to) && curstate$to[[1]] == curstate$state_id) {
      await_low_expression = TRUE
      await_low_expression_module = curstate$modules
    }
    
    # determine where this bifurcation will go
    for (i in seq_len(nrow(expression_modules_scaled)-window)) {
      next_expression = expression_modules_scaled[i:(i+window), curstate$modules, drop=TRUE]
      
      #print(await_low_expression)
      
      if(await_low_expression) {
        low_expression_module_expression = expression_modules_scaled[i:(i+window), await_low_expression_module, drop=TRUE]
        #print(glue::glue("{i} <<<<<<<<<<<<<<<<<< {mean(low_expression_module_expression)}"))
        if(mean(low_expression_module_expression) < max_before_breakpoint_expression) {
          await_low_expression = FALSE
          #print("switched")
        }
      } else if(curstate$type == "bifurcation") {
        ## BIFURCATION upcoming
        module_order = sort(colMeans(next_expression), decreasing = TRUE)
        
        if(module_order[[1]]-module_order[[2]] > 0.25) {
          nextstate_id = curstate$to[which(curstate$modules == names(module_order)[[1]])]
          
          if(nextstate_id == curstate$state_id) {
            await_low_expression = TRUE
            await_low_expression_module = names(module_order)[[1]]
          }
        }
      } else if (curstate$type == "convergence") {
        ## AWAITING HIGH MODULE EXPRESSION (for cycles)
        # make sure the previous expression has been low for some time, then check whether the next expression will be high
        nextstate_id = curstate$to[[1]]
        
      } else if (curstate$type == "cycle") {
        #print(glue::glue("{i} >>>>>>>>>>>>>>>> {mean(next_expression)}"))
        if(mean(next_expression) > min_after_breakpoint_expression) {
          nextstate_id = curstate$to[[1]]
          await_low_expression = TRUE
          await_low_expression_module = curstate$modules
        }
      }
      
      # if next state_id has been chosen: change current state
      if(!is.na(nextstate_id)) {
        # first add the expression to tobecomes
        
        # only get expression until the previous time this state has been included
        # avoids having expression which went left to be included in the expression which now goes right
        prevsplitstart = ifelse(as.character(curstate$state_id) %in% names(prevsplits), prevsplits[[as.character(curstate$state_id)]], prevsplitstart)
        
        tobecomes = c(tobecomes, list(list(
          simulation_id=simulation_id, 
          state_id=curstate$state_id, 
          expression=list(expression_modules_scaled[prevsplitstart:(i+window),,drop=FALSE]), 
          startstep_id = lastsplit, endstep_id = i, nextstate_id=nextstate_id, 
          piece_id=length(tobecomes) + 1,
          prevstate_id = prevstate_id
        )
        ))
        
        if(!is.na(laststate_id)) prevsplits[[as.character(laststate_id)]] = i
        
        lastsplit = i + 1
        prevstate_id = curstate$state_id
        
        nextstate = states %>% keep(~.$state_id == nextstate_id)
        #print(nextstate_id)
        laststate_id = curstate$state_id
        if(length(nextstate) == 0) {
          curstate = list(state_id=nextstate_id)
          break()
        } else {
          curstate = nextstate[[1]]
        }
        
        if(curstate$type == "convergence" && !is.null(curstate$modules)) {
          await_low_expression = TRUE
          await_low_expression_module = curstate$modules
        }
        
        nextstate_id = NA
      }
    }
    
    # add end
    
    prevsplitstart = ifelse(as.character(curstate$state_id) %in% names(prevsplits), prevsplits[[as.character(curstate$state_id)]], prevsplitstart)
    
    tobecomes = c(tobecomes, list(list(
      simulation_id=simulation_id, 
      state_id=curstate$state_id, 
      expression=list(expression_modules_scaled[prevsplitstart:nrow(expression_modules_scaled),,drop=FALSE]), 
      startstep_id = lastsplit, endstep_id = nrow(expression_modules_scaled), nextstate_id=NA, piece_id=length(tobecomes) + 1,
      prevstate_id = prevstate_id
    )
    ))
    
    simulation_idoi = simulation_id
    #pheatmap(t(expression_modules_scaled), cluster_cols=F, cluster_rows=F)
    #pheatmap(t(expression_modules_scaled), cluster_cols=F, cluster_rows=F, gaps_col = tobecomes %>% bind_rows() %>% filter(simulation_id == simulation_idoi) %>% pull(startstep_id))
  }
  #tobecomes = tobecomes %>% bind_rows() %>% mutate(piece_id=seq_len(nrow(.)))
  
  
  tobecomes
}

lagpad <- function(x, k=1) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  }
  else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}

divide_pieces <- function(simulations, tobecomes, states) {
  quant_scale_combined = get_combined_quant_scale(simulations)
  
  combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)
  
  tobecomes_combined = tobecomes %>% bind_rows() %>% mutate(piece_ids=map2(expression, piece_id, ~rep(.y, nrow(.x)))) %>% group_by(state_id) %>% summarise(expression=list(invoke(rbind, expression)), nextstate_ids=list(nextstate_id), prevstate_ids=list(prevstate_id), piece_ids=list(unlist(piece_ids))) # to becomes for each state_id
  #subsamples = map(tobecomes_combined$nextstate_ids, ~sample(seq_along(.), size=length(.) * 0.1))
  #tobecomes_combined$expression = map2(subsamples, tobecomes_combined$expression, ~.y[.x, ])
  #tobecomes_combined$nextstate_ids = map2(subsamples, tobecomes_combined$nextstate_ids, ~.y[.x])
  
  references = split(tobecomes_combined, seq_len(nrow(tobecomes_combined))) %>% map(function(x) map(x, ~.[[1]])) %>% {set_names(., map(., "state_id"))} # contains for each state associated expression values for which it is known to which next state this cell will progress
  
  pieces <- pbapply::pblapply(seq_along(simulations), function(simulation_idoi) {
    subpieces <- list()
    
    expression_modules_scaled <- quant_scale_combined(simulations[[simulation_idoi]]$expression_modules)
    subtobecomes <- tobecomes %>% keep(~.$simulation_id == simulation_idoi)
    
    prevsplit <- 1
    for (row_id in seq_along(subtobecomes)) {
      row <- subtobecomes[[row_id]]
      
      # if current state_id is not a terminal state
      if(as.character(row$state_id) %in% names(states)) {
        curstate <- states[[as.character(row$state_id)]]
        
        if(curstate$type == "bifurcation" && !is.na(row$nextstate_id)) {
          # BIFURCATING
          
          print("looking for earliest bifurcation moment")
          
          reference = references[[as.character(curstate$state_id)]]
          
          for (i in seq(row$endstep_id, row$startstep_id, -5)) {
            current_expression = expression_modules_scaled[i, ,drop=FALSE]
            distances = pdist::pdist(current_expression[, curstate$modules], reference$expression[, curstate$modules])@dist %>% tapply(reference$piece_id, min)
            
            result = wilcox.test(distances[reference$nextstate_ids == row$nextstate_id], distances[reference$nextstate_ids != row$nextstate_id], alternative="greater")
            result$p.value = ifelse(is.na(result$p.value), 1, result$p.value)
            
            if(result$p.value > 0.05) {
              print(glue::glue("split at {i}"))
              break()
            }
          }
        } else if(curstate$type == "convergence") {
          ## CHECK FOR CONVERGENCE
          
          print("looking for first convergence moment")
          
          reference = references[[as.character(curstate$to)]]
          
          for (i in seq(row$endstep_id, nrow(expression_modules_scaled), 50)) {
            current_expression = expression_modules_scaled[i, ,drop=FALSE]
            distances = pdist::pdist(current_expression, reference$expression)@dist %>% tapply(reference$piece_id, min)
            
            result = wilcox.test(distances[!(reference$prevstate_ids %in% curstate$cross)], distances[reference$prevstate_ids %in% curstate$cross])
            
            result$p.value = ifelse(is.na(result$p.value), 1, result$p.value)
            
            if(result$p.value > 0.05) {
              print(glue::glue("converged at {i}"))
              break()
            }
          }
        } else {
          i <- row$endstep_id
        }
        
        subpieces <- c(subpieces, list(tibble(
          start_step_id = prevsplit, 
          end_step_id = i, 
          state_id = curstate$state_id,
          simulation_id = simulation_idoi,
          expression = list(expression_modules_scaled[seq(prevsplit, i),,drop=F])
        )))
        prevsplit <- i + 1
        
      } else {
        if(prevsplit < row$endstep_id) {
          subpieces <- c(subpieces, list(tibble(
            start_step_id = prevsplit, 
            end_step_id = row$endstep_id, 
            state_id = row$state_id,
            simulation_id = simulation_idoi,
            expression = list(expression_modules_scaled[seq(prevsplit, row$endstep_id),,drop=F])
          )))
        }
        
        break()
      }
    }
    
    subpieces
  })
  
  pieces <- unlist(pieces, recursive = FALSE)
  
  pieces %>% bind_rows %>% 
    mutate(
      prevstate_id = ifelse(lagpad(simulation_id) == simulation_id, lagpad(state_id), NA), 
      nextstate_id = ifelse(lagpad(simulation_id, -1) == simulation_id, lagpad(state_id, -1), NA),
      piece_id = seq_along(state_id),
      nextpiece_id = ifelse(lagpad(simulation_id, -1) == simulation_id, lagpad(state_id, -1), NA),
      cell_id = map(expression, ~rownames(.)),
      step_id = map2(start_step_id, end_step_id, ~seq(.x, .y)),
      piecestep_id = map2(start_step_id, end_step_id, ~seq_len(.y-.x+1))
    )
}


extract_piece_times <- function(pieces, states) {
  map(unique(pieces$state_id), function(state_idoi) {
    state <- states %>% keep(~.$state_id == state_idoi) 
    if(length(state)) {iscycle = state[[1]] %>% {.$state_id %in% .$to}} else {iscycle=FALSE} # check whether this state is a cycle, for circular averaging of times
    
    expressions <- pieces %>% filter(state_id==state_idoi) %>% .$expression 
    expressions_filtered <- expressions %>% map(~.[as.integer(seq(1, nrow(.), length.out=min(20, nrow(.)))), ,drop=FALSE])
    
    combined <- expressions_filtered %>% invoke(rbind, .)
    
    init.traj = expressions_filtered[[expressions %>% map_int(nrow) %>% which.max]]
    fit <- princurve::principal.curve(combined, start = init.traj, plot.true = F, trace = F, stretch = 0)
    times <- fit$lambda / max(fit$lambda)
    
    splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
    
    times <- splitAt(times, (map_int(expressions_filtered, nrow) %>% cumsum) + 1) # split according to piece
    
    map(seq_along(expressions), function(piece_id) {
      expressionoi <- expressions[[piece_id]]
      expressionoi_filtered <- expressions_filtered[[piece_id]]
      timesoi <- times[[piece_id]]
      
      #linear interpolation, only if we estimated time of a subsample of all steps
      if(length(timesoi) < nrow(expressionoi)) {
        approx(which(rownames(expressionoi) %in% rownames(expressionoi_filtered)), timesoi, seq_len(nrow(expressionoi)))$y %>% tibble(time=., cell_id=rownames(expressionoi), scaletime = time/max(time))
      } else {
        tibble(time=timesoi, cell_id=rownames(expressionoi), scaletime = timesoi/max(timesoi))
      }
    }) %>% bind_rows()
  }) %>% bind_rows()
}

# due to a bug in smoothing, R will halt somewhere in C code if the smoothing is run in parallel. If you want to extract the gold standard in parallel, smooth beforehand
extract_goldstandard <- function(simulation, verbose=FALSE, smooth=FALSE) {
  gs <- list()
  
  # smoothing + calculating module expression
  if(is.null(simulation$simulations[[1]]$expression_modules)) simulation$simulations <- smoothe_simulations(simulation$simulations, simulation$model)
  
  if(verbose) {print("raw division")}
  tobecomes <- extract_tobecomes(simulation$simulations, simulation$model$states)
  
  if(verbose) {print("exact division")}
  pieces <- divide_pieces(simulation$simulations, tobecomes, simulation$model$states)
  
  cellinfo <- pieces %>% unnest(cell_id, step_id, piecestep_id) %>% select(-start_step_id, -end_step_id)
  
  if(verbose) {print("extracting times")}
  cellinfo <- cellinfo %>% left_join(extract_piece_times(pieces, simulation$model$states), by="cell_id")
  gs$cellinfo <- cellinfo
  
  if(verbose) {print("converting to milestones")}
  gs <- c(gs, get_milestones(simulation$simulations, simulation$model, gs))
  
  gs
}

#plot_goldstandard(experiment, gs)


get_combined_quant_scale <- function(simulations) {
  combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)
  quant_scale_combined
}

plot_goldstandard <- function(gs,experiment=NULL, simulations=NULL, plot_what="progression") {
  ##
  source("scripts/evaluation/methods/dimred.R")
  if(!is.null(simulations)) {
    combined = map(simulations, ~.$expression) %>% do.call(rbind, .)
    expression = combined[gs$cellinfo$cell_id, ]
    
    dimred_methods <- c(ica)
  } else {
    expression = experiment$expression %>% set_rownames(experiment$cellinfo$step_id)
    
    dimred_methods <- c(ica, tsne, mds, mds_smacof)
  }
  
  expression = log2(expression + 1)
  
  
  map(dimred_methods, function(space_func) {
    space = space_func(expression, ndim = 2)
    #space = mds_smacof(expression, ndim = 2)
    #space = dambiutils::mds_withlandmarks(expression, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S
    
    
    if(plot_what=="progression") {
      plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cell_id=rownames(.)) %>% left_join(gs$cellinfo, by="cell_id")
      
      list(
        ggplot(plotdata %>% mutate(state_id=factor(state_id))) + geom_point(aes(Comp1, Comp2, color=state_id, group=simulation_id))
        ,
        ggplot(plotdata %>% mutate(state_id=factor(state_id))) + geom_point(aes(Comp1, Comp2, color=time, group=simulation_id)) + viridis::scale_color_viridis()
      )
    } else if (plot_what == "percentages") {
      plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cell_id=rownames(.)) %>% left_join(gs$cellinfo, by="cell_id") %>% right_join(gs$milestone_percentages, by="cell_id")
      
      list(
        ggplot(plotdata %>% mutate(state_id=factor(state_id))) + geom_point(aes(Comp1, Comp2, color=percentage)) + facet_wrap(~milestone) + viridis::scale_color_viridis(option="B", direction = -1)
      )
    }
  })
}


plot_goldstandard_percentages <- function(experiment, gs) {
  ##
  combined = map(experiment$simulations, ~.$expression) %>% do.call(rbind, .)
  expression = combined[gs$cellinfo$cell_id, ]
  
  expression = experiment$expression %>% set_rownames(experiment$cellinfo$simulationstep_id)
  ##
  
  expression = log2(expression + 1)
  
  source("scripts/evaluation/methods/dimred.R")
  
  map(c(ica, tsne, mds, mds_smacof), function(space_func) {
    space = space_func(expression, ndim = 2)
    #space = mds_smacof(expression, ndim = 2)
    #space = dambiutils::mds_withlandmarks(expression, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S
    
    plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cell_id=rownames(.)) %>% left_join(gs$cellinfo, by="cell_id") %>% right_join(gs$percentages, by="cell_id")
    
    list(
      ggplot(plotdata %>% mutate(state_id=factor(state_id))) + geom_point(aes(Comp1, Comp2, color=percentage)) + facet_wrap(~milestone) + viridis::scale_color_viridis(option="B", direction = -1)
    )
  })
}



# needs: model$states, gs$cellinfo
get_milestones <- function(simulations, model, gs) {
  # construct the milestone_network
  # similar to statenet, but with deduplication of circular edges
  # and added start-end milestones
  
  statenet <- map(model$states, function(x) {if(is.null(x$to)){x$to = NA};data.frame(from=x$state_id, to=x$to)}) %>% bind_rows()
  milestone_network <- statenet
  cellinfo <- gs$cellinfo %>% mutate(state_id = as.numeric(state_id))
  
  maxtimes <- gs$cellinfo %>% group_by(state_id) %>% summarise(maxtime=max(time)) %>% {set_names(.$maxtime, .$state_id)}
  
  milestone2state_id <- unique(c(statenet$from, statenet$to)) %>% keep(~!is.na(.)) %>% {set_names(., .)}
  
  
  ## deduplication of cycles: from 1->1 to 1->2->3->1
  if(nrow(milestone_network) > 0) {
    maxpiece_id <- (milestone_network %>% select(from, to) %>% max) + 1 # next name for a piece_id
    
    # divide each circular path in three subpaths (three milestones)
    for (state_id in milestone_network %>% filter(from == to) %>% .$from) {
      newpieces <- c(state_id, seq(maxpiece_id, length.out=2))
      
      timepartition <- maxtimes[[state_id]] / 3
      
      cellinfo[cellinfo$state_id == state_id, ] <- cellinfo[cellinfo$state_id == state_id, ] %>% mutate(
        state_id = newpieces[pmax(1, ceiling(time/timepartition))],
        time = time %% timepartition
      )
      
      milestone2state_id[as.character(newpieces)] <- state_id
      
      maxtimes <- c(maxtimes, set_names(rep(timepartition, 2), as.character(newpieces[c(2, 3)])))
      maxtimes[[as.character(state_id)]] <- timepartition
      
      milestone_network <- bind_rows(milestone_network, tibble(from=c(newpieces), to=c(newpieces[c(2, 3,1)])))
      
      maxpiece_id <- maxpiece_id + 2
    }
    milestone_network <- milestone_network %>% filter(is.na(to) | (from != to))
  }
  
  # add terminal nodes
  terminalnodes <- c(milestone_network$from[is.na(milestone_network$to)], milestone_network$to[!(milestone_network$to %in% milestone_network$from)]) %>% keep(~!is.na(.))
  milestone_network <- milestone_network %>% bind_rows(tibble(from = terminalnodes, to=seq(max(c(milestone_network$from, milestone_network$to) %>% keep(~!is.na(.))) + 1, length.out=length(terminalnodes)))) %>% filter(!is.na(to))
  
  milestone_network$length <- maxtimes[milestone_network$from]
  
  # join with milestone_network, find from percentage, use this to get to percentage(s) for every to
  milestone_percentages <- cellinfo %>% left_join(milestone_network, by=c("state_id"="from")) %>% 
    mutate(from=state_id, from_percentage=1-time/maxtimes[as.character(from)]) %>% 
    group_by(cell_id) %>% 
    mutate(to_percentage=(1-from_percentage)/length(to)) %>% 
    summarise(milestone=list(c(from[[1]], to)), percentage=list(c(from_percentage[[1]], to_percentage))) %>% 
    unnest(milestone, percentage)
  
  progression <- cellinfo %>% left_join(milestone_network, by=c("state_id"="from")) %>% 
    mutate(from=state_id, from_percentage=1-time/maxtimes[as.character(from)]) %>% 
    group_by(cell_id) %>% 
    mutate(to_percentage=(1-from_percentage)/length(to)) %>% select(from, to, to_percentage) %>% rename(percentage=to_percentage)
  
  tibble::lst(milestone_percentages, milestone_network, progression)
}