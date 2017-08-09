get_bias_knn = function(bifurcation, experiment, cellinfo, gs, k=10) {
  model = experiment$model
  endmodules = bifurcation$endmodules
  
  piecesoi = gs$pieces %>% filter(stateid %in% bifurcation$stateid) %>% filter(!is.na(nextstateid)) # take those pieces which are currently bifurcating, then filter those with no next state
  piecesoi$nextstateidsoi = piecesoi$nextstateids %>% map(~.[. %in% bifurcation$nextstateids]) # get the next stateids of interest
  piecesoi = piecesoi %>% rowwise() %>% filter(length(nextstateidsoi) > 0) %>% ungroup() # filter those pieces which have any stateidoi in their next states
  piecesoi$nextstateidoi = map_int(piecesoi$nextstateidsoi, first)
  piecesoi = piecesoi %>% select(-nextstateids, -nextstateidsoi)
  piecesoi$original_expression = map2(piecesoi$simulationid, piecesoi$cells, ~experiment$simulations[[.x]]$expression[.y, ])
  piecesoi$module_expression = map(piecesoi$original_expression, function(expression) {
    modulecounts <- get_module_counts(expression, model$modulemembership)[, endmodules]
    modulecounts %>% as.data.frame()
  })
  observations = piecesoi %>% select(-cells, -expression, -original_expression) %>% unnest() %>% mutate(simulationid=factor(simulationid), nextstateidoi=factor(nextstateidoi), nextstateid=factor(nextstateid))  # contains different simulation steps, and at each point the expression of the modules of interest, together with their nextstateid, which can be used to find for each "real" cell what the probability will be to get to a next stateid
  
  plot(
    ggplot(observations) +
      geom_path(aes_string(endmodules[[1]], endmodules[[2]], group="simulationid", color="nextstateidoi")) + coord_equal()
  )
  
  observations_expression = observations[, endmodules]
  
  cellsoi = cellinfo %>% filter(stateid %in% bifurcation$stateids)
  cellsoi_expression = experiment$expression_modules[cellsoi$cell, endmodules]
  
  nearest_neighbors = lapply(seq_along(rownames(cellsoi_expression)), function(rowid) {
    row = cellsoi_expression[rowid, ]
    differences = sweep(observations_expression, 2, row) %>% abs() %>% rowMeans()
    ordered = order(differences, decreasing=F)
    data.frame(cell=rownames(cellsoi_expression)[rowid], nextstateidoi=observations$nextstateidoi[ordered[seq_len(k)]])
  }) %>% bind_rows()
  cellbiases = nearest_neighbors %>% group_by(cell, nextstateidoi) %>% summarise(n=n()) %>% mutate(p=n/sum(n)) %>% ungroup() %>% select(-n) %>% spread(nextstateidoi, p, drop=FALSE, fill=0)
  
  if(!all(bifurcation$nextstateids %in% colnames(cellbiases))) stop("not all biases recovered!!")
  
  #cellsoi_expression %>% as.data.frame() %>% rownames_to_column("cell") %>% left_join(cellbiases, by="cell") %>% 
  #ggplot() + geom_point(aes(M3, M4, color=`2`)) + scale_color_gradientn(colors=c("red", "black", "blue")) + coord_equal()
  
  cellbiases
}

# get bias using local density
get_bias_density = function(bifurcation, experiment, cellinfo, gs, k=10) {
  model = experiment$model
  endmodules = bifurcation$endmodules
  
  piecesoi = gs$pieces %>% filter(stateid %in% bifurcation$stateid) %>% filter(!is.na(nextstateid)) # take those pieces which are currently bifurcating, then filter those with no next state
  piecesoi$nextstateidsoi = piecesoi$nextstateids %>% map(~.[. %in% bifurcation$nextstateids]) # get the next stateids of interest
  piecesoi = piecesoi %>% rowwise() %>% filter(length(nextstateidsoi) > 0) %>% ungroup() # filter those pieces which have any stateidoi in their next states
  piecesoi$nextstateidoi = map_int(piecesoi$nextstateidsoi, first)
  piecesoi = piecesoi %>% select(-nextstateids, -nextstateidsoi)
  piecesoi$original_expression = map2(piecesoi$simulationid, piecesoi$cells, ~experiment$simulations[[.x]]$expression[.y, ])
  piecesoi$module_expression = map(piecesoi$original_expression, function(expression) {
    modulecounts <- get_module_counts(expression, model$modulemembership)[, endmodules]
    modulecounts %>% as.data.frame()
  })
  observations = piecesoi %>% select(-cells, -expression, -original_expression) %>% unnest() %>% mutate(simulationid=factor(simulationid), nextstateidoi=factor(nextstateidoi), nextstateid=factor(nextstateid))  # contains different simulation steps, and at each point the expression of the modules of interest, together with their nextstateid, which can be used to find for each "real" cell what the probability will be to get to a next stateid
  
  plot(
    ggplot(observations) +
      geom_path(aes_string(endmodules[[1]], endmodules[[2]], group="simulationid", color="nextstateidoi")) + coord_equal()
  )
  
  observations_expression = observations[, endmodules]
  
  cellsoi = cellinfo %>% filter(stateid %in% bifurcation$stateids)
  cellsoi_expression = experiment$expression_modules[cellsoi$cell, endmodules]
  
  H = matrix(c(0.1, 0, 0, 0.1)*0.05, ncol=2)
  observations_expression = observations %>% filter(nextstateidoi == bifurcation$nextstateids[[1]]) %>% .[, endmodules]
  density1 = ks::kde(observations_expression, H=H)
  observations_expression = observations %>% filter(nextstateidoi == bifurcation$nextstateids[[2]]) %>% .[, endmodules]
  density2 = ks::kde(observations_expression, H=H)
  
  pheatmap(log2((density1$estimate+0.001)/(density2$estimate+0.001)), cluster_cols=F, cluster_rows=F)
  
  x1s = ((predict(density1, x=cellsoi_expression) + 0.00001)/(predict(density2, x=cellsoi_expression) + 0.00001))
  x2s = 1
  perc1 = x1s / (1+x1s)
  perc2 = 1-perc1
  cellbiases = tibble(cell=cellsoi$cell)
  cellbiases[, as.character(bifurcation$nextstateids[[1]])] = perc1
  cellbiases[, as.character(bifurcation$nextstateids[[2]])] = perc2
  
  cellbiases
}

get_bias = get_bias_knn

# get_milestones = function(experiment, gs) {
#   # construct the milestonenet
#   # similar to statenet, but with deduplication of circular edges
#   # and added start-end milestones
#   
#   ## deduplication of cycles: from 1->1 to 1->2->3->1
#   maxprogressions = gs$reference$cellinfo %>% group_by(stateid) %>% summarise(maxprogression=max(progression)) %>% {set_names(.$maxprogression, .$stateid)}
#   milestonenet = gs$statenet %>% select(from, to)
#   cellinfo = gs$cellinfo %>% mutate(stateid = as.numeric(stateid))
#   
#   milestone2stateid = unique(c(gs$statenet$from, gs$statenet$to)) %>% {set_names(., .)}
#   
#   if(nrow(milestonenet) > 0) {
#     maxpieceid = (milestonenet %>% select(from, to) %>% max) + 1 # next name for a pieceid
#     for (stateid in milestonenet %>% filter(from == to) %>% .$from) {
#       newpieces = c(stateid, seq(maxpieceid, length.out=2))
#       
#       progressionpartition = maxprogressions[[stateid]] / 3
#       
#       cellinfo[cellinfo$stateid == stateid, ] <- cellinfo[cellinfo$stateid == stateid, ] %>% mutate(
#         stateid = newpieces[ceiling(progression/progressionpartition)],
#         progression = progression %% progressionpartition
#       )
#       
#       milestone2stateid[as.character(newpieces)] = stateid
#       
#       maxprogressions <- c(maxprogressions, set_names(rep(progressionpartition, 2), as.character(newpieces[c(2, 3)])))
#       maxprogressions[[as.character(stateid)]] <- progressionpartition
#       
#       milestonenet <- bind_rows(milestonenet, tibble(from=c(newpieces), to=c(newpieces[c(2, 3,1)])))
#       
#       maxpieceid = maxpieceid + 2
#     }
#     milestonenet <- milestonenet %>% filter(from != to)
#   }
#   
#   # add terminal nodes
#   terminalnodes = gs$piecenet$piece %>% keep(!(. %in% milestonenet$from))
#   milestonenet = milestonenet %>% bind_rows(tibble(from = terminalnodes, to=seq(max(c(milestonenet$from, milestonenet$to)) + 1, length.out=length(terminalnodes))))
#   milestonenet$length = maxprogressions[milestonenet$from]
#   
#   
#   # stateids: where to get the reference cells, to which states can these belong? nextstateids: possible follow up stateids.; tentstateid: the starting stateid for the tent
#   bifurcations = lapply(gs$statenet %>% count(from) %>% filter(n>1) %>% .$from, function(bifurcation_stateid) {
#     nextstateids = gs$statenet %>% filter(from==bifurcation_stateid) %>% .$to
#     nextstate_firststates = map_int(nextstateids, ~gs$states[[.]][[1]])
#     endmodules = experiment$model$modulenodes$module[match(nextstate_firststates, experiment$model$modulenodes$state)]
#     
#     list(
#       stateids = bifurcation_stateid,
#       nextstateids = nextstateids,
#       tentstateid = bifurcation_stateid,
#       endmodules = endmodules
#     )
#   })
#   
#   print(bifurcations)
#   
#   cellbiases = bifurcations %>% map(get_bias, experiment=experiment, gs=gs, cellinfo=cellinfo) %>% c(list(tibble(cell=cellinfo$cell)), .) %>% plyr::join_all(type="full", by="cell")
#   ## add end states
#   
#   percentages = lapply(seq_len(nrow(cellinfo)), function(cellid) {
#     cellinfo = cellinfo[cellid, ]
#     nextstateids = milestonenet %>% filter(from == cellinfo$stateid) %>% .$to
#     
#     frompercentage = cellinfo$progression/maxprogressions[[as.character(cellinfo$stateid)]]
#     
#     if(length(nextstateids) > 1) {
#       # if bifurcating
#       topercentages = cellbiases %>% filter(cell == cellinfo$cell) %>% .[, as.character(milestone2stateid[as.character(nextstateids)])] %>% as.list() %>% unlist()
#     } else {
#       # if just progressing forward
#       topercentages = 1 %>% set_names(nextstateids)
#     }
#     
#     percentages = (1-frompercentage) * topercentages
#     percentages[[as.character(cellinfo$stateid)]] = frompercentage
#     
#     tibble(cell=cellinfo$cell, percentage=percentages, milestone=as.numeric(names(percentages)))
#   }) %>% bind_rows()
#   percentages = percentages %>% spread(milestone, percentage, fill=0)
#   
#   milestonenet = milestonenet %>% mutate(from=as.character(from), to=as.character(to))
#   colnames(percentages) = as.character(colnames(percentages))
#   
#   named_list(
#     milestone_names = percentages %>% select(-cell) %>% colnames, 
#     milestone_net = milestonenet, 
#     milestone_percentages = percentages %>% rename(id=cell)#,
#     #milestone_percentages_notent = untent_percentages(milestonenet, percentages) %>% rename(id=cell)
#   )
# }



untent_percentages = function(milestone_net, milestone_percentages) {
  bifurcations = milestone_net %>% group_by(from) %>% summarise(to=list(to), nto=n()) %>% filter(nto>1)
  if(nrow(bifurcations) > 0) {
    bifurcations = bifurcations %>% {split(., 1:nrow(.))}
    for(bifurcation in bifurcations) {
      rowids = milestone_percentages[, c(bifurcation$from, bifurcation$to[[1]])] %>% {(. > 0) & (.<1)} %>% apply(1, all)
      milestone_percentages[rowids,bifurcation$to[[1]]] = ((1-milestone_percentages[rowids, bifurcation$from])/2) %>% as.matrix() %>% as.numeric()
    }
  }
  milestone_percentages
}