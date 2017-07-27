simulations <- smoothe_simulations(experiment$simulations, model)

combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)

piecestates <- model$piecestates

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
  
  if(curpiecestate$to[[1]] == curpiecestate$piecestateid) {
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
      
      if(module_order[[1]]-module_order[[2]] > 0.5) {
        nextpiecestateid = curpiecestate$to[which(curpiecestate$modules == names(module_order)[[1]])]
        
        if(nextpiecestateid %in% c(3, 4)) {
          print(paste0(nextpiecestateid, "   ", module_order[[1]],"    ", module_order[[2]]))
        }
        
        if(nextpiecestateid == curpiecestate$piecestateid) {
          await_low_expression = TRUE
          await_low_expression_module = names(module_order)[[1]]
        }
      }
    } else if (curpiecestate$type == "convergence") {
      ## AWAITING HIGH MODULE EXPRESSION (for cycles)
      # make sure the previous expression has been low for some time, then check whether the next expression will be high
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

tobecomes_combined = tobecomes %>% bind_rows() %>% mutate(pieceids=map2(expression, pieceid, ~rep(.y, nrow(.x)))) %>% group_by(piecestateid) %>% summarise(expression=list(invoke(rbind, expression)), nextpiecestateids=list(nextpiecestateid), prevpiecestateids=list(prevpiecestateid), pieceids=list(unlist(pieceids))) # to becomes for each piecestateid
#subsamples = map(tobecomes_combined$nextpiecestateids, ~sample(seq_along(.), size=length(.) * 0.1))
#tobecomes_combined$expression = map2(subsamples, tobecomes_combined$expression, ~.y[.x, ])
#tobecomes_combined$nextpiecestateids = map2(subsamples, tobecomes_combined$nextpiecestateids, ~.y[.x])

references = split(tobecomes_combined, seq_len(nrow(tobecomes_combined))) %>% map(function(x) map(x, ~.[[1]])) %>% {set_names(., map(., "piecestateid"))} # contains for each piecestate associated expression values for which it is known to which next piecestate this cell will progress

# now check for every simulation at what particular points its bifurcation happens (ie. at what point the cell is sure it will become a particular state)

confidence <- 0.95
bifurcation_mindist <- 0.4
convergence_mindist <- 1

step_size <- 10

p_value_bifurcation <- 1e-3
p_value_convergence <- 0.05
distance_until_possible_convergence <- 1

simulationid = 40
pieces = mclapply(seq_along(simulations), function(simulationid) {
  pieces = list()
  
  print(glue::glue("############### {simulationid}"))
  
  expression_modules_scaled = quant_scale_combined(simulations[[simulationid]]$expression_modules)
  
  #pheatmap(expression_modules_scaled, cluster_cols=F, cluster_rows=F, gaps_row=600)
  
  curpiecestate = piecestates[[1]]
  nextpiecestateid = curpiecestate$piecestateid
  
  prevsplit = 1
  
  tmplist = tibble()
  
  for (i in seq(1, nrow(expression_modules_scaled), step_size)) {
    current_expression = expression_modules_scaled[i, ,drop=FALSE]
    #current_expression = expression_modules_scaled[i:min(i+window, nrow(expression_modules_scaled)), ,drop=FALSE]
    
    if(length(curpiecestate$modules) > 1) {
      ## BIFURCATION upcoming
      reference = references[[as.character(curpiecestate$piecestateid)]]
      
      distances = pdist::pdist(current_expression[, curpiecestate$modules], reference$expression[, curpiecestate$modules])@dist %>% tapply(reference$pieceid, min)
      
      mean_distances = tapply(distances, reference$nextpiecestateids, mean)
      result = t.test(distances[reference$nextpiecestateids == reference$nextpiecestateids[[1]]], distances[reference$nextpiecestateids != reference$nextpiecestateids[[1]]])
      result$p.value = ifelse(is.na(result$p.value), 1, result$p.value)
      
      print(mean_distances)
      
      tmplist = bind_rows(tmplist, data.frame(t = result$statistic, p=result$p.value, split="split", i=i))
      
      #tmplist = bind_rows(tmplist, tibble(distance=distances, nextpiecestateid=reference$nextpiecestateids, i=i, pieceid=names(distances)))
      
      print(simulationid)
      
      if(result$p.value < p_value_bifurcation) {
        nextpiecestateid = as.numeric(names(mean_distances))[which.min(mean_distances)]
      }
      
      # to check for convergence, expression won't be split initially, wait until complete split to start checking for convergence
      split = FALSE
      split_expression = current_expression
      
      prevpiecestate = curpiecestate
    } else if(!is.null(curpiecestate$cross)) {
      ## CONVERGENCE CHECK upcoming, first wait for split, then wait for convergence again
      if(length(curpiecestate$cross) > 1) stop("multiple crosses not supported yet")
      reference_cross = references[[as.character(curpiecestate$cross)]]
      close_cross = sum(pdist::pdist(current_expression, reference_cross$expression)@dist < convergence_mindist)
      distances_cross = pdist::pdist(current_expression, reference_cross$expression)@dist %>% tapply(reference_cross$pieceid, min)
      
      reference_same = references[[as.character(curpiecestate$piecestateid)]]
      close_same = sum(pdist::pdist(current_expression, reference_same$expression)@dist < convergence_mindist)
      distances_same = pdist::pdist(current_expression, reference_same$expression)@dist %>% tapply(reference_same$pieceid, min)
      
      result = t.test(distances_cross, distances_same, exact=FALSE)
      result$p.value = ifelse(is.na(result$p.value), 1, result$p.value) # if variances are 0
      
      tmplist = bind_rows(tmplist, data.frame(t = result$statistic, p=result$p.value, split="convergence", i=i, dist=pdist(current_expression, split_expression)@dist))
      
      
      # check previous bifurcation
      reference = references[[as.character(prevpiecestate$piecestateid)]]
      
      distances = pdist::pdist(current_expression[, prevpiecestate$modules], reference$expression[, prevpiecestate$modules])@dist %>% tapply(reference$pieceid, min)
      
      mean_distances = tapply(distances, reference$nextpiecestateids, mean)
      result2 = t.test(distances[reference$nextpiecestateids == reference$nextpiecestateids[[1]]], distances[reference$nextpiecestateids != reference$nextpiecestateids[[1]]])
      result2$p.value = ifelse(is.na(result2$p.value), 1, result2$p.value)
      
      
      if(split) {
        # already split, check for convergence now
        if(result$p.value > p_value_convergence && result2$p.value > p_value_bifurcation) {
          nextpiecestateid = curpiecestate$to[[1]]
        }
      } else {
        # not yet split, check for split
        dist_from_split = pdist(current_expression, split_expression)@dist
        if(dist_from_split > distance_until_possible_convergence) {
          split = TRUE
        }
      }
      
      #print(glue::glue("{close_cross} ___ {close_same}"))
      
    } else if(length(curpiecestate$modules) == 1) {
      ## MODULE CHECK upcoming
      prev_expression = expression_modules_scaled[max(i-window, 0):i, curpiecestate$modules, drop=TRUE]
      if(mean(prev_expression) < max_before_breakpoint_expression) {
        if(mean(current_expression[, curpiecestate$modules]) > min_after_breakpoint_expression) {
          nextpiecestateid = curpiecestate$to[[1]]
        }
      }
    }
    
    if(nextpiecestateid != curpiecestate$piecestateid) {
      pieces = c(pieces, list(tibble(
        simulationid=simulationid, 
        piecestateid=curpiecestate$piecestateid, 
        expression=list(expression_modules_scaled[prevsplit:i,,drop=FALSE]), 
        start_stepid=prevsplit,
        end_stepid=i,
        nextpiecestateid=nextpiecestateid)
      ))
      
      print(paste0(curpiecestate$piecestateid, " -> ", nextpiecestateid, "         ", i))
      
      prevsplit = i
      
      nextpiecestate = piecestates %>% keep(~.$piecestateid == nextpiecestateid)
      
      if(length(nextpiecestate) == 0) {
        pieces = c(pieces, list(tibble(
          simulationid=simulationid, 
          piecestateid=nextpiecestateid, 
          expression=list(expression_modules_scaled[i:nrow(expression_modules_scaled),,drop=FALSE]), 
          start_stepid=i,
          end_stepid=nrow(expression_modules_scaled),
          nextpiecestateid=NA)
        ))
        break
      } else {
        curpiecestate = nextpiecestate[[1]]
      }
    }
  }
  
  # add latest piece if not added yet
  last_stepid = ifelse(length(pieces) > 0, pieces[[length(pieces)]]$end_stepid, 1)
  if(last_stepid < nrow(expression_modules_scaled)) {
    pieces = c(pieces, list(tibble(
      simulationid=simulationid, 
      piecestateid=nextpiecestateid, 
      expression=list(expression_modules_scaled[last_stepid:nrow(expression_modules_scaled),,drop=FALSE]), 
      start_stepid=last_stepid,
      end_stepid=nrow(expression_modules_scaled),
      nextpiecestateid=NA)
    ))
  }
  
  pieces = pieces %>% bind_rows()
  
  #cellinfo = pieces %>% rowwise() %>% mutate(piecestateids = list(rep(piecestateid, end_stepid - start_stepid))) %>% ungroup() %>% select(piecestateids) %>% unnest() %>% mutate(piecestateids=factor(piecestateids)) %>% as.data.frame() %>% set_rownames(rownames(expression_modules_scaled))
  #pheatmap(t(expression_modules_scaled), cluster_cols=F, cluster_rows=F, annotation_col = cellinfo)
  
  pieces
}, mc.cores = 8)
pieces = pieces %>% bind_rows
pieces = pieces %>% mutate(ncells=map_int(expression, nrow))


pieces %>% filter(piecestateid == 2) %>% {pheatmap(t(invoke(rbind, .$expression)), cluster_cols=F, cluster_rows=F, gaps_col=cumsum(map_int(.$expression, ~nrow(.))))}


simulationstepinfo = pieces %>% mutate(simulationstepid = map(expression, ~rownames(.)), stepid = map2(start_stepid, end_stepid, ~seq(.x, .y))) %>% unnest(simulationstepid, stepid)
