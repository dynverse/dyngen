get_bias = function(bifurcation, experiment, gs, k=50) {
  endmodules = bifurcation$endmodules
  
  piecesoi = gs$pieces %>% filter(piecestateid %in% bifurcation$piecestateid) %>% filter(!is.na(nextpiecestateid)) # take those pieces which are currently bifurcating, then filter those with no next state
  piecesoi$nextpiecestateidsoi = piecesoi$nextpiecestateids %>% map(~.[. %in% bifurcation$nextpiecestateids]) # get the next piecestateids of interest
  piecesoi = piecesoi %>% rowwise() %>% filter(length(nextpiecestateidsoi) > 0) %>% ungroup() # filter those pieces which have any piecestateidoi in their next states
  piecesoi$nextpiecestateidoi = map_int(piecesoi$nextpiecestateidsoi, first)
  piecesoi = piecesoi %>% select(-nextpiecestateids, -nextpiecestateidsoi)
  piecesoi$original_expression = map2(piecesoi$simulationid, piecesoi$cells, ~experiment$simulations[[.x]]$expression[.y, ])
  piecesoi$module_expression = map(piecesoi$original_expression, function(expression) {
    modulecounts <- get_module_counts(expression, model$modulemembership)[, endmodules]
    modulecounts %>% as.data.frame()
  })
  observations = piecesoi %>% select(-cells, -expression, -original_expression) %>% unnest() %>% mutate(simulationid=factor(simulationid), nextpiecestateidoi=factor(nextpiecestateidoi), nextpiecestateid=factor(nextpiecestateid))  # contains different simulation steps, and at each point the expression of the modules of interest, together with their nextpiecestateid, which can be used to find for each "real" cell what the probability will be to get to a next piecestateid
  
  plot(
    ggplot(observations) +
      geom_path(aes_string(endmodules[[1]], endmodules[[2]], group="simulationid", color="nextpiecestateidoi")) + coord_equal()
  )
  
  observations_expression = observations[, endmodules]
  
  cellsoi = gs$cellinfo %>% filter(piecestateid %in% bifurcation$piecestateids)
  cellsoi_expression = experiment$expression_modules[cellsoi$cell, endmodules]
  
  nearest_neighbors = lapply(seq_along(rownames(cellsoi_expression)), function(rowid) {
    row = cellsoi_expression[rowid, ]
    differences = sweep(observations_expression, 2, row) %>% abs() %>% rowMeans()
    ordered = order(differences, decreasing=F)
    data.frame(cell=rownames(cellsoi_expression)[rowid], nextpiecestateidoi=observations$nextpiecestateidoi[ordered[seq_len(k)]])
  }) %>% bind_rows()
  cellbiases = nearest_neighbors %>% group_by(cell, nextpiecestateidoi) %>% summarise(n=n()) %>% mutate(p=n/sum(n)) %>% ungroup() %>% select(-n) %>% spread(nextpiecestateidoi, p, fill=0)

  #cellsoi_expression %>% as.data.frame() %>% rownames_to_column("cell") %>% left_join(cellbiases, by="cell") %>% 
  #ggplot() + geom_point(aes(M3, M4, color=`2`)) + scale_color_gradientn(colors=c("red", "black", "blue")) + coord_equal()
  
  cellbiases
}

get_milestones = function(experiment) {
  gs = experiment$gs
  
  # piecestateids: where to get the reference cells, to which piecestates can these belong? nextpiecestateids: possible follow up piecestateids.; tentpiecestateid: the starting piecestateid for the tent
  bifurcations = lapply(gs$piecestatenet %>% count(from) %>% filter(n>1) %>% .$from, function(bifurcation_piecestateid) {
    nextpiecestateids = gs$piecestatenet %>% filter(from==bifurcation_piecestateid) %>% .$to
    nextpiecestate_firststates = map_int(nextpiecestateids, ~gs$piecestates[[.]][[1]])
    endmodules = experiment$model$modulenodes$module[match(nextpiecestate_firststates, experiment$model$modulenodes$state)]
    
    list(
      piecestateids = bifurcation_piecestateid,
      nextpiecestateids = nextpiecestateids,
      tentpiecestateid = bifurcation_piecestateid,
      endmodules = endmodules
    )
  })
  
  
  cellbiases = bifurcations %>% map(get_bias,) %>% join_all(type="full", by="cell")
  ## add end states
  maxprogressions = gs$reference$cellinfo %>% group_by(piecestateid) %>% summarise(maxprogression=max(progression)) %>% {set_names(.$maxprogression, .$piecestateid)}
  
  milestonenet = gs$piecestatenet %>% select(from, to)
  terminalnodes = unique(c(milestonenet$from, milestonenet$to)) %>% keep(!(. %in% milestonenet$from))
  milestonenet = milestonenet %>% bind_rows(tibble(from = terminalnodes, to=seq(max(c(milestonenet$from, milestonenet$to)) + 1, length.out=length(terminalnodes))))
  milestonenet$length = maxprogressions[milestonenet$from]
  
  
  percentages = lapply(seq_len(nrow(gs$cellinfo)), function(cellid) {
    cellinfo = gs$cellinfo[cellid, ]
    nextpiecestateids = milestonenet %>% filter(from == cellinfo$piecestateid) %>% .$to
    
    frompercentage = cellinfo$progression/maxprogressions[[as.character(cellinfo$piecestateid)]]
    
    if(length(nextpiecestateids) > 1) {
      # if bifurcating
      topercentages = cellbiases %>% filter(cell == cellinfo$cell) %>% .[, as.character(nextpiecestateids)] %>% as.list() %>% unlist()
    } else {
      # if just progressing forward
      topercentages = 1 %>% set_names(nextpiecestateids)
    }
    
    percentages = (1-frompercentage) * topercentages
    percentages[[as.character(cellinfo$piecestateid)]] = frompercentage
    
    tibble(cell=cellinfo$cell, percentage=percentages, milestone=as.numeric(names(percentages)))
  }) %>% bind_rows()
  percentages = percentages %>% spread(milestone, percentage, fill=0)
  
  milestonenet = milestonenet %>% mutate(from=as.character(from), to=as.character(to))
  colnames(percentages) = as.character(colnames(percentages))
  
  named_list(
    milestone_names = percentages %>% select(-cell) %>% colnames, 
    milestone_net = milestonenet, 
    milestone_percentages = percentages %>% rename(id=cell))
}