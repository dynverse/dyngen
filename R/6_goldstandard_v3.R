extract_tobecomes <- function(simulations, piecestates) {
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
  for (simulationid in seq_along(simulations)) {
    print(simulationid)
    
    expression_modules_scaled = quant_scale_combined(simulations[[simulationid]]$expression_modules)
    
    curpiecestate = piecestates[[1]]
    
    nextpiecestateid = NA
    lastpiecestateid = NA
    
    prevsplitstart = 1
    lastsplit = 1
    
    prevpiecestateid = NA
    
    prevsplits = numeric() # for each piecestate, keep record of the previous stepid, so that we do not use the same expression twice for the same bifurcation (possible leading to a different subsequent state)
    
    await_low_expression = FALSE
    
    if(!is.null(curpiecestate$to) && curpiecestate$to[[1]] == curpiecestate$piecestateid) {
      await_low_expression = TRUE
      await_low_expression_module = curpiecestate$modules
    }
    
    # determine where this bifurcation will go
    for (i in seq_len(nrow(expression_modules_scaled)-window)) {
      next_expression = expression_modules_scaled[i:(i+window), curpiecestate$modules, drop=TRUE]
      
      #print(await_low_expression)
      
      if(await_low_expression) {
        low_expression_module_expression = expression_modules_scaled[i:(i+window), await_low_expression_module, drop=TRUE]
        #print(glue::glue("{i} <<<<<<<<<<<<<<<<<< {mean(low_expression_module_expression)}"))
        if(mean(low_expression_module_expression) < max_before_breakpoint_expression) {
          await_low_expression = FALSE
          #print("switched")
        }
      } else if(curpiecestate$type == "bifurcation") {
        ## BIFURCATION upcoming
        module_order = sort(colMeans(next_expression), decreasing = TRUE)
        
        if(module_order[[1]]-module_order[[2]] > 0.25) {
          nextpiecestateid = curpiecestate$to[which(curpiecestate$modules == names(module_order)[[1]])]
          
          if(nextpiecestateid == curpiecestate$piecestateid) {
            await_low_expression = TRUE
            await_low_expression_module = names(module_order)[[1]]
          }
        }
      } else if (curpiecestate$type == "convergence") {
        ## AWAITING HIGH MODULE EXPRESSION (for cycles)
        # make sure the previous expression has been low for some time, then check whether the next expression will be high
        nextpiecestateid = curpiecestate$to[[1]]
        
      } else if (curpiecestate$type == "cycle") {
        #print(glue::glue("{i} >>>>>>>>>>>>>>>> {mean(next_expression)}"))
        if(mean(next_expression) > min_after_breakpoint_expression) {
          nextpiecestateid = curpiecestate$to[[1]]
          await_low_expression = TRUE
          await_low_expression_module = curpiecestate$modules
        }
      }
      
      # if next piecestateid has been chosen: change current piecestate
      if(!is.na(nextpiecestateid)) {
        # first add the expression to tobecomes
        
        # only get expression until the previous time this piecestate has been included
        # avoids having expression which went left to be included in the expression which now goes right
        prevsplitstart = ifelse(as.character(curpiecestate$piecestateid) %in% names(prevsplits), prevsplits[[as.character(curpiecestate$piecestateid)]], prevsplitstart)
        
        tobecomes = c(tobecomes, list(list(
          simulationid=simulationid, 
          piecestateid=curpiecestate$piecestateid, 
          expression=list(expression_modules_scaled[prevsplitstart:(i+window),,drop=FALSE]), 
          startstepid = lastsplit, endstepid = i, nextpiecestateid=nextpiecestateid, 
          pieceid=length(tobecomes) + 1,
          prevpiecestateid = prevpiecestateid
        )
        ))
        
        if(!is.na(lastpiecestateid)) prevsplits[[as.character(lastpiecestateid)]] = i
        
        lastsplit = i + 1
        prevpiecestateid = curpiecestate$piecestateid
        
        nextpiecestate = piecestates %>% keep(~.$piecestateid == nextpiecestateid)
        #print(nextpiecestateid)
        lastpiecestateid = curpiecestate$piecestateid
        if(length(nextpiecestate) == 0) {
          curpiecestate = list(piecestateid=nextpiecestateid)
          break()
        } else {
          curpiecestate = nextpiecestate[[1]]
        }
        
        if(curpiecestate$type == "convergence" && !is.null(curpiecestate$modules)) {
          await_low_expression = TRUE
          await_low_expression_module = curpiecestate$modules
        }
        
        nextpiecestateid = NA
      }
    }
    
    # add end
    
    prevsplitstart = ifelse(as.character(curpiecestate$piecestateid) %in% names(prevsplits), prevsplits[[as.character(curpiecestate$piecestateid)]], prevsplitstart)
    
    tobecomes = c(tobecomes, list(list(
      simulationid=simulationid, 
      piecestateid=curpiecestate$piecestateid, 
      expression=list(expression_modules_scaled[prevsplitstart:nrow(expression_modules_scaled),,drop=FALSE]), 
      startstepid = lastsplit, endstepid = nrow(expression_modules_scaled), nextpiecestateid=NA, pieceid=length(tobecomes) + 1,
      prevpiecestateid = prevpiecestateid
    )
    ))
    
    simulationidoi = simulationid
    #pheatmap(t(expression_modules_scaled), cluster_cols=F, cluster_rows=F)
    #pheatmap(t(expression_modules_scaled), cluster_cols=F, cluster_rows=F, gaps_col = tobecomes %>% bind_rows() %>% filter(simulationid == simulationidoi) %>% pull(startstepid))
  }
  #tobecomes = tobecomes %>% bind_rows() %>% mutate(pieceid=seq_len(nrow(.)))
  
  
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

divide_pieces <- function(simulations, tobecomes, piecestates) {
  quant_scale_combined = get_combined_quant_scale(simulations)
  
  combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)
  
  tobecomes_combined = tobecomes %>% bind_rows() %>% mutate(pieceids=map2(expression, pieceid, ~rep(.y, nrow(.x)))) %>% group_by(piecestateid) %>% summarise(expression=list(invoke(rbind, expression)), nextpiecestateids=list(nextpiecestateid), prevpiecestateids=list(prevpiecestateid), pieceids=list(unlist(pieceids))) # to becomes for each piecestateid
  #subsamples = map(tobecomes_combined$nextpiecestateids, ~sample(seq_along(.), size=length(.) * 0.1))
  #tobecomes_combined$expression = map2(subsamples, tobecomes_combined$expression, ~.y[.x, ])
  #tobecomes_combined$nextpiecestateids = map2(subsamples, tobecomes_combined$nextpiecestateids, ~.y[.x])
  
  references = split(tobecomes_combined, seq_len(nrow(tobecomes_combined))) %>% map(function(x) map(x, ~.[[1]])) %>% {set_names(., map(., "piecestateid"))} # contains for each piecestate associated expression values for which it is known to which next piecestate this cell will progress
  
  pieces <- pbapply::pblapply(seq_along(simulations), function(simulationidoi) {
    subpieces <- list()
    
    expression_modules_scaled <- quant_scale_combined(simulations[[simulationidoi]]$expression_modules)
    subtobecomes <- tobecomes %>% keep(~.$simulationid == simulationidoi)
    
    prevsplit <- 1
    for (rowid in seq_along(subtobecomes)) {
      row <- subtobecomes[[rowid]]
      
      # if current piecestateid is not a terminal state
      if(as.character(row$piecestateid) %in% names(piecestates)) {
        curpiecestate <- piecestates[[as.character(row$piecestateid)]]
        
        if(curpiecestate$type == "bifurcation" && !is.na(row$nextpiecestateid)) {
          # BIFURCATING
          
          print("looking for earliest bifurcation moment")
          
          reference = references[[as.character(curpiecestate$piecestateid)]]
          
          for (i in seq(row$endstepid, row$startstepid, -5)) {
            current_expression = expression_modules_scaled[i, ,drop=FALSE]
            distances = pdist::pdist(current_expression[, curpiecestate$modules], reference$expression[, curpiecestate$modules])@dist %>% tapply(reference$pieceid, min)
            
            result = wilcox.test(distances[reference$nextpiecestateids == row$nextpiecestateid], distances[reference$nextpiecestateids != row$nextpiecestateid], alternative="greater")
            result$p.value = ifelse(is.na(result$p.value), 1, result$p.value)
            
            if(result$p.value > 0.05) {
              print(glue::glue("split at {i}"))
              break()
            }
          }
        } else if(curpiecestate$type == "convergence") {
          ## CHECK FOR CONVERGENCE
          
          print("looking for first convergence moment")
          
          reference = references[[as.character(curpiecestate$to)]]
          
          for (i in seq(row$endstepid, nrow(expression_modules_scaled), 50)) {
            current_expression = expression_modules_scaled[i, ,drop=FALSE]
            distances = pdist::pdist(current_expression, reference$expression)@dist %>% tapply(reference$pieceid, min)
            
            result = wilcox.test(distances[!(reference$prevpiecestateids %in% curpiecestate$cross)], distances[reference$prevpiecestateids %in% curpiecestate$cross])
            
            result$p.value = ifelse(is.na(result$p.value), 1, result$p.value)
            
            if(result$p.value > 0.05) {
              print(glue::glue("converged at {i}"))
              break()
            }
          }
        } else {
          i <- row$endstepid
        }
        
        subpieces <- c(subpieces, list(tibble(
          start_stepid = prevsplit, 
          end_stepid = i, 
          piecestateid = curpiecestate$piecestateid,
          simulationid = simulationidoi,
          expression = list(expression_modules_scaled[seq(prevsplit, i),,drop=F])
        )))
        prevsplit <- i + 1
        
      } else {
        if(prevsplit < row$endstepid) {
          subpieces <- c(subpieces, list(tibble(
            start_stepid = prevsplit, 
            end_stepid = row$endstepid, 
            piecestateid = row$piecestateid,
            simulationid = simulationidoi,
            expression = list(expression_modules_scaled[seq(prevsplit, row$endstepid),,drop=F])
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
      prevpiecestateid = ifelse(lagpad(simulationid) == simulationid, lagpad(piecestateid), NA), 
      nextpiecestateid = ifelse(lagpad(simulationid, -1) == simulationid, lagpad(piecestateid, -1), NA),
      pieceid = seq_along(piecestateid),
      nextpieceid = ifelse(lagpad(simulationid, -1) == simulationid, lagpad(piecestateid, -1), NA),
      cellid = map(expression, ~rownames(.)),
      stepid = map2(start_stepid, end_stepid, ~seq(.x, .y)),
      piecestepid = map2(start_stepid, end_stepid, ~seq_len(.y-.x+1))
    )
}


extract_piece_times <- function(pieces, piecestates) {
  map(unique(pieces$piecestateid), function(piecestateidoi) {
    piecestate <- piecestates %>% keep(~.$piecestateid == piecestateidoi) 
    if(length(piecestate)) {iscycle = piecestate[[1]] %>% {.$piecestateid %in% .$to}} else {iscycle=FALSE} # check whether this piecestate is a cycle, for circular averaging of times
    
    expressions <- pieces %>% filter(piecestateid==piecestateidoi) %>% .$expression 
    expressions_filtered <- expressions %>% map(~.[as.integer(seq(1, nrow(.), length.out=min(20, nrow(.)))), ,drop=FALSE])
    
    combined <- expressions_filtered %>% invoke(rbind, .)
    
    init.traj = expressions_filtered[[expressions %>% map_int(nrow) %>% which.max]]
    fit <- princurve::principal.curve(combined, start = init.traj, plot.true = F, trace = F, stretch = 0)
    times <- fit$lambda / max(fit$lambda)
    
    splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
    
    times <- splitAt(times, (map_int(expressions_filtered, nrow) %>% cumsum) + 1) # split according to piece
    
    map(seq_along(expressions), function(pieceid) {
      expressionoi <- expressions[[pieceid]]
      expressionoi_filtered <- expressions_filtered[[pieceid]]
      timesoi <- times[[pieceid]]
      
      #linear interpolation, only if we estimated time of a subsample of all steps
      if(length(timesoi) < nrow(expressionoi)) {
        approx(which(rownames(expressionoi) %in% rownames(expressionoi_filtered)), timesoi, seq_len(nrow(expressionoi)))$y %>% tibble(time=., cellid=rownames(expressionoi), scaletime = time/max(time))
      } else {
        tibble(time=timesoi, cellid=rownames(expressionoi), scaletime = timesoi/max(timesoi))
      }
    }) %>% bind_rows()
  }) %>% bind_rows()
}

# due to a bug in smoothing, R will halt somewhere in C code if the smoothing is run in parallel. If you want to extract the gold standard in parallel, smooth beforehand
extract_goldstandard <- function(experiment, verbose=FALSE, smooth=FALSE) {
  gs <- list()
  
  # smoothing + calculating module expression
  if(smooth) experiment$simulations <- smoothe_simulations(experiment$simulations, experiment$model)
  
  if(verbose) {print("raw division")}
  tobecomes <- extract_tobecomes(experiment$simulations, experiment$model$piecestates)
  
  if(verbose) {print("exact division")}
  pieces <- divide_pieces(experiment$simulations, tobecomes, experiment$model$piecestates)
  
  cellinfo <- pieces %>% unnest(cellid, stepid, piecestepid) %>% select(-start_stepid, -end_stepid)
  
  if(verbose) {print("extracting times")}
  cellinfo <- cellinfo %>% left_join(extract_piece_times(pieces, experiment$model$piecestates), by="cellid")
  gs$cellinfo <- cellinfo
  
  if(verbose) {print("converting to milestones")}
  gs <- c(gs, get_milestones(experiment, gs))
  
  gs
}

#plot_goldstandard(experiment, gs)


get_combined_quant_scale <- function(simulations) {
  combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)
  quant_scale_combined
}

plot_goldstandard <- function(experiment, gs, simulations=NULL) {
  ##
  if(!is.null(simulations)) {
    combined = map(experiment$simulations, ~.$expression) %>% do.call(rbind, .)
    expression = combined[gs$cellinfo$cellid, ]
  } else {
    expression = experiment$expression %>% set_rownames(experiment$cellinfo$stepid)
  }
  
  expression = log2(expression + 1)
  
  source("scripts/evaluation/methods/dimred.R")
  
  
  map(c(ica, tsne, mds, mds_smacof), function(space_func) {
    space = space_func(expression, ndim = 2)
    #space = mds_smacof(expression, ndim = 2)
    #space = dambiutils::mds_withlandmarks(expression, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S
    
    plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cellid=rownames(.)) %>% left_join(gs$cellinfo, by="cellid")
    
    list(
      ggplot(plotdata %>% mutate(piecestateid=factor(piecestateid))) + geom_point(aes(Comp1, Comp2, color=piecestateid, group=simulationid))
      ,
      ggplot(plotdata %>% mutate(piecestateid=factor(piecestateid))) + geom_point(aes(Comp1, Comp2, color=time, group=simulationid)) + viridis::scale_color_viridis()
    )
  })
}


plot_goldstandard_percentages <- function(experiment, gs) {
  ##
  combined = map(experiment$simulations, ~.$expression) %>% do.call(rbind, .)
  expression = combined[gs$cellinfo$cellid, ]
  
  expression = experiment$expression %>% set_rownames(experiment$cellinfo$simulationstepid)
  ##
  
  expression = log2(expression + 1)
  
  source("scripts/evaluation/methods/dimred.R")
  
  map(c(ica, tsne, mds, mds_smacof), function(space_func) {
    space = space_func(expression, ndim = 2)
    #space = mds_smacof(expression, ndim = 2)
    #space = dambiutils::mds_withlandmarks(expression, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S
    
    plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cellid=rownames(.)) %>% left_join(gs$cellinfo, by="cellid") %>% right_join(gs$percentages, by="cellid")
    
    list(
      ggplot(plotdata %>% mutate(piecestateid=factor(piecestateid))) + geom_point(aes(Comp1, Comp2, color=percentage)) + facet_wrap(~milestone) + viridis::scale_color_viridis(option="B", direction = -1)
    )
  })
}




get_milestones <- function(experiment, gs) {
  # construct the milestonenet
  # similar to piecestatenet, but with deduplication of circular edges
  # and added start-end milestones
  
  piecestatenet <- map(experiment$model$piecestates, function(x) {if(is.null(x$to)){x$to = NA};data.frame(from=x$piecestateid, to=x$to)}) %>% bind_rows()
  milestonenet <- piecestatenet
  cellinfo <- gs$cellinfo %>% mutate(piecestateid = as.numeric(piecestateid))
  
  maxtimes <- gs$cellinfo %>% group_by(piecestateid) %>% summarise(maxtime=max(time)) %>% {set_names(.$maxtime, .$piecestateid)}
  
  milestone2piecestateid <- unique(c(piecestatenet$from, piecestatenet$to)) %>% keep(~!is.na(.)) %>% {set_names(., .)}
  
  
  ## deduplication of cycles: from 1->1 to 1->2->3->1
  if(nrow(milestonenet) > 0) {
    maxpieceid <- (milestonenet %>% select(from, to) %>% max) + 1 # next name for a pieceid
    
    # divide each circular path in three subpaths (three milestones)
    for (piecestateid in milestonenet %>% filter(from == to) %>% .$from) {
      newpieces <- c(piecestateid, seq(maxpieceid, length.out=2))
      
      timepartition <- maxtimes[[piecestateid]] / 3
      
      cellinfo[cellinfo$piecestateid == piecestateid, ] <- cellinfo[cellinfo$piecestateid == piecestateid, ] %>% mutate(
        piecestateid = newpieces[pmax(1, ceiling(time/timepartition))],
        time = time %% timepartition
      )
      
      milestone2piecestateid[as.character(newpieces)] <- piecestateid
      
      maxtimes <- c(maxtimes, set_names(rep(timepartition, 2), as.character(newpieces[c(2, 3)])))
      maxtimes[[as.character(piecestateid)]] <- timepartition
      
      milestonenet <- bind_rows(milestonenet, tibble(from=c(newpieces), to=c(newpieces[c(2, 3,1)])))
      
      maxpieceid <- maxpieceid + 2
    }
    milestonenet <- milestonenet %>% filter(is.na(to) | (from != to))
  }
  
  # add terminal nodes
  terminalnodes <- c(milestonenet$from[is.na(milestonenet$to)], milestonenet$to[!(milestonenet$to %in% milestonenet$from)]) %>% keep(~!is.na(.))
  milestonenet <- milestonenet %>% bind_rows(tibble(from = terminalnodes, to=seq(max(c(milestonenet$from, milestonenet$to) %>% keep(~!is.na(.))) + 1, length.out=length(terminalnodes)))) %>% filter(!is.na(to))
  
  milestonenet$length <- maxtimes[milestonenet$from]
  
  # join with milestonenet, find from percentage, use this to get to percentage(s) for every to
  percentages <- cellinfo %>% left_join(milestonenet, by=c("piecestateid"="from")) %>% 
    mutate(from=piecestateid, from_percentage=1-time/maxtimes[as.character(from)]) %>% 
    group_by(cellid) %>% 
    mutate(to_percentage=(1-from_percentage)/length(to)) %>% 
    summarise(milestone=list(c(from, to)), percentage=list(c(from_percentage, to_percentage))) %>% 
    unnest(milestone, percentage)
  
  tibble::lst(percentages, milestonenet)
}

