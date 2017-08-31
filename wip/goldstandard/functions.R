extract_tobecomes <- function(simulations, piecestates) {
  combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)
  
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
  tobecomes_combined = tobecomes %>% bind_rows() %>% mutate(pieceids=map2(expression, pieceid, ~rep(.y, nrow(.x)))) %>% group_by(piecestateid) %>% summarise(expression=list(invoke(rbind, expression)), nextpiecestateids=list(nextpiecestateid), prevpiecestateids=list(prevpiecestateid), pieceids=list(unlist(pieceids))) # to becomes for each piecestateid
  #subsamples = map(tobecomes_combined$nextpiecestateids, ~sample(seq_along(.), size=length(.) * 0.1))
  #tobecomes_combined$expression = map2(subsamples, tobecomes_combined$expression, ~.y[.x, ])
  #tobecomes_combined$nextpiecestateids = map2(subsamples, tobecomes_combined$nextpiecestateids, ~.y[.x])
  
  references = split(tobecomes_combined, seq_len(nrow(tobecomes_combined))) %>% map(function(x) map(x, ~.[[1]])) %>% {set_names(., map(., "piecestateid"))} # contains for each piecestate associated expression values for which it is known to which next piecestate this cell will progress
  
  
  pieces <- list()
  
  for (simulationidoi in seq_along(simulations)) {
    expression_modules_scaled <- quant_scale_combined(simulations[[simulationidoi]]$expression_modules)
    subtobecomes <- tobecomes %>% keep(~.$simulationid == simulationidoi)
    
    prevsplit <- 1
    for (rowid in seq_along(subtobecomes)) {
      row <- subtobecomes[[rowid]]
      
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
          
          for (i in seq(row$endstepid, nrow(expression_modules_scaled), 5)) {
            current_expression = expression_modules_scaled[i, ,drop=FALSE]
            distances = pdist::pdist(current_expression, reference$expression)@dist %>% tapply(reference$pieceid, min)
            
            result = wilcox.test(distances[!(reference$prevpiecestateids %in% curpiecestate$cross)], distances[reference$prevpiecestateids %in% curpiecestate$cross])
            
            result$p.value = ifelse(is.na(result$p.value), 1, result$p.value)
            
            pvals = c(pvals, result$p.value)
            
            print(result$p.value)
            
            if(result$p.value > 0.05) {
              print(glue::glue("converged at {i}"))
              break()
            }
          }
        } else {
          i <- row$endstepid
        }
        
        pieces <- c(pieces, list(tibble(
          start_stepid = prevsplit, 
          end_stepid = i, 
          piecestateid = curpiecestate$piecestateid,
          simulationid = simulationidoi,
          expression = list(expression_modules_scaled[seq(prevsplit, i),,drop=F])
        )))
        prevsplit <- i + 1
        
      } else {
        pieces <- c(pieces, list(tibble(
          start_stepid = prevsplit, 
          end_stepid = row$endstepid, 
          piecestateid = row$piecestateid,
          simulationid = simulationidoi,
          expression = list(expression_modules_scaled[seq(prevsplit, row$endstepid),,drop=F])
        )))
        
        break()
      }
    }
  }
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
pieces %>% bind_rows %>% {hist(.$start_stepid)}

pieces %>% filter(piecestateid == 4) %>% sample_n(5) %>% {pheatmap(t(invoke(rbind, .$expression)), cluster_cols=F, cluster_rows=F, gaps_col=cumsum(map_int(.$expression, ~nrow(.))))}










pheatmap(t(expression_modules_scaled), cluster_cols=F, cluster_rows=F, gaps_col = pieces %>% filter(simulationid == simulationidoi) %>% pull(end_stepid))




##
source("~/thesis/projects/dyngen/scripts/evaluation/methods/dimred.R")
combined = map(experiment$simulations, ~.$expression) %>% do.call(rbind, .)
combined = combined[simulationstepinfo$simulationstepid, ]
space = ica(combined, ndim = 2)
space = tsne(combined, ndim = 2)
space = SCORPIUS::reduce.dimensionality.landmarked(combined, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S

plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% bind_cols(simulationstepinfo)
ggplot(plotdata %>% mutate(piecestateid=factor(piecestateid))) + geom_path(aes(Comp1, Comp2, color=piecestateid, group=simulationid))




###
expression = experiment$expression %>% {log2(.+1)}
rownames(expression) = experiment$cellinfo$simulationstepid
source("~/thesis/projects/dyngen/scripts/evaluation/methods/dimred.R")
space = ica(expression, ndim = 2)
space = tsne(expression, ndim = 2)
space = mds(expression, ndim = 2)
space = mds_smacof(expression, ndim = 2)

plotdata = bind_cols(space %>% as.data.frame, simulationstepinfo %>% slice(match(rownames(space), simulationstepid)))
ggplot(plotdata) + geom_point(aes(Comp2, Comp1, color=factor(piecestateid)))# + viridis::scale_color_viridis(option="A")
