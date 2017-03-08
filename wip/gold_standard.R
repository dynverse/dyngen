simulations = qsub_lapply(1:24, function(i) {
  simulate_one_cell(model, burntime, totaltime, ncells=NULL)
})
# for (i in seq_len(length(simulations))) {
#   x = simulations[[i]]$cellinfo$simulationtime
#   simulations[[i]]$expression_smooth = apply(simulations[[i]]$expression, 2, function(y) {
#     lo = loess(y~x, span=0.5)
#     lo$fitted
#   }) %>% set_rownames(rownames(simulations[[i]]$expression))
# }

get_piecenet = function(statenet) {
  piecenet = statenet %>% select(from, to) %>% mutate(contains = map(seq_len(nrow(statenet)), ~integer()))
  for (node in unique(c(statenet$from, statenet$to))) {
    froms = piecenet %>% filter(from == node)
    tos = piecenet %>% filter(to == node)
    
    if ((nrow(froms) == 1) && (nrow(tos) == 1) && froms$to != node) {
      newfrom = tos$from
      newto = froms$to
      contains = c(froms$contains[[1]], tos$contains[[1]], node)
      last = node
      piecenet = piecenet %>% filter(from != node) %>% filter(to != node) %>% bind_rows(tibble(from=newfrom, to=newto, contains=list(contains), previous=last))
    }
  }
  piecenet %<>% mutate(piece=seq_len(nrow(.)))
  piecenet %>% graph_from_data_frame() %>% plot
  
  piecenet
}

get_piecestates = function(piecenet) {
  piecestates <- lapply(seq_len(nrow(piecenet)), function(rowid){
    piece = c()
    if(piecenet$from[[rowid]] == 1) {
      piece = c(1)
    }
    #piece = c(piecenet$from[[rowid]]) # remove if not include first interaction (except state 1)
    c(piece, piecenet$contains[[rowid]], piecenet$to[[rowid]])
  })
  
  piecestates
}


# first get the smoothed module expression of all simulations
smoothe_simulations = function(simulations) {
  for (simulationid in seq_len(length(simulations))) {
    expression_smooth = simulations[[simulationid]]$expression %>% zoo::rollmean(50, c("extend", "extend", "extend")) %>% set_rownames(rownames(simulations[[simulationid]]$expression))
    expression_modules = get_module_counts(expression_smooth, model$modulemembership)
    simulations[[simulationid]]$expression_smooth = expression_smooth
    simulations[[simulationid]]$expression_modules = expression_modules
  }
  simulations
}

# the linear state ordering of every piece
# based on the state progression, divide the expression data into linear pieces
divide_simulation = function(progressioninfo, piecestates, expression_smooth) {
  statesunique = integer()
  statesunique_windows = list()
  curstate = 1
  states = progressioninfo$state %>% c(NA)
  curstart = 1
  for(i in seq_len(length(states))) {
    if(is.na(states[[i]]) || states[[i]] != curstate) {
      statesunique = c(statesunique, curstate)
      statesunique_windows = c(statesunique_windows, list(rownames(expression_smooth)[seq(curstart, i-1)]))
      curstart = i
      curstate = states[[i]]
    }
  }
  
  print(statesunique)
  
  pieces = list()
  # check for every possible piecestate
  for (piecestateid in seq_len(length(piecestates))) {
    piecestate = piecestates[[piecestateid]]
    # check for every state sequence
    for (statesuniqueid in seq_len(length(statesunique))) {
      start = statesuniqueid
      end = statesuniqueid+length(piecestate)-1
      if(all(statesunique[seq(start, end)] %>% keep(~!is.na(.)) == piecestate)) {
        piece_cells = statesunique_windows[seq(start, end)] %>% unlist
        piece = 
          {tibble(
            cells=list(piece_cells), 
            piecestateid=piecestateid, 
            expression = list(expression_smooth[piece_cells, ])
          )}
        pieces = c(pieces, list(piece))
      }
    }
  }
  pieces <- pieces %>% bind_rows
  
  pieces
}



# now use this to scale all simulations
divide_simulations = function(simulations, piecestates, model) {
  combined = map(simulations, "expression_modules") %>% do.call(rbind, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.01)
  quant_scale_combined = function(x) SCORPIUS::apply.quant.scale(x, attributes(combined_scaled)$center, attributes(combined_scaled)$scale)
  
  pieces = list()
  for (simulationid in seq_len(length(simulations))) {
    expression_modules_scaled = quant_scale_combined(simulations[[simulationid]]$expression_modules)
    expression_modules_scaled = simulations[[simulationid]]$expression_modules
    
    progressioninfo = get_states(expression_modules_scaled, model, minexpression = 0.3, minratio = 2) %>% select(-cell) %>% bind_cols(simulations[[simulationid]]$cellinfo)
    ggplot(progressioninfo) + geom_area(aes(simulationtime, group=state, fill=factor(state)), position="fill", stat="bin", bins=50)
    ggplot(progressioninfo) + geom_point(aes(simulationtime, state))
    
    pheatmap(expression_modules_scaled %>% t, cluster_cols=F, cluster_rows=F, annotation_col = progressioninfo %>% mutate(state=factor(state)) %>% as.data.frame() %>% {set_rownames(., .$cell)} %>% select(state))
    
    expression_pieces = divide_simulation(progressioninfo, piecestates, simulations[[simulationid]]$expression_smooth)
    #expression_pieces[[2]] %>% t %>% pheatmap(cluster_cols=F, scale="row", cluster_rows=T)
    
    expression_pieces$simulationid = simulationid
    
    print(expression_pieces$piecestateid)
    print("---")
    
    pieces = c(pieces, list(expression_pieces))
  }
  pieces = pieces %>% bind_rows()
  
  pieces
}

# determine average expression, mapped to the first piece
average_pieces = function(piecesoi, model) {
  total = tibble()
  piece1id = piecesoi$expression %>% map_int(nrow) %>% order() %>% {.[round(length(.)/2)]}
  piece1 = piecesoi$expression[[piece1id]]
  piece1_times = seq_len(nrow(piece1))
  for (i in seq_len(nrow(piecesoi)-1)+1) {
    piece2 = piecesoi[i, ]$expression[[1]]
    
    celldistances = SCORPIUS::euclidean.distance(piece1, piece2)
    
    total = data.frame(pieceid = i, time=piece1_times[apply(celldistances, 2, which.min)], piece2 %>% reshape2::melt(varnames=c("cell", "gene"), value.name="expression")) %>% as_tibble() %>% bind_rows(total)
  }
  meanexpression = total %>% reshape2::acast(time~gene, mean, value.var = "expression")
  meanexpression_times = as.numeric(rownames(meanexpression))
  
  window = 30
  meanexpression = meanexpression %>% zoo::rollmean(window, c("extend", "extend", "extend")) %>% set_rownames(rownames(meanexpression))
  
  
  #meanexpression %>% t %>% pheatmap(cluster_rows=T, cluster_cols=F, scale="row")
  get_module_counts(meanexpression, model$modulemembership) %>% t %>% pheatmap(cluster_rows=F, cluster_cols=F)
  
  rownames(meanexpression) = seq_len(nrow(meanexpression))/nrow(meanexpression) * map_int(piecesoi$expression, nrow) %>% mean # get realtime estimate based on average number of cells in each piece
  
  meanexpression
}

get_reference_expression = function(pieces, model) {
  reference_list = lapply(unique(pieces$piecestateid) %>% sort, function(piecestateidoi) {
    piecesoi = pieces %>% filter(piecestateid == piecestateidoi)
    
    meanexpression = average_pieces(piecesoi, model)
    
    print(piecesoi %>% nrow)
    
    list(cellinfo = tibble(piecestateid=piecestateidoi, progression=as.numeric(rownames(meanexpression))), meanexpression=meanexpression)
  })
  reference_expression = map(reference_list, "meanexpression") %>% do.call(rbind, .)
  reference_cellinfo = map(reference_list, "cellinfo") %>% bind_rows() %>% as.data.frame()
  
  rownames(reference_expression) = seq_len(nrow(reference_expression))
  rownames(reference_cellinfo) = rownames(reference_expression)
  reference_cellinfo$piecestateid = factor(reference_cellinfo$piecestateid)
  
  #reference_expression %>% t %>% pheatmap(cluster_cols=F, cluster_rows=T, annotation_col=reference_cellinfo, gaps_col = which(diff(as.numeric(reference_cellinfo$piecestateid)) != 0))
  reference_expression %>% get_module_counts(., model$modulemembership) %>% t %>% pheatmap(cluster_cols=F, cluster_rows=T, annotation_col=reference_cellinfo, gaps_col = which(diff(as.numeric(reference_cellinfo$piecestateid)) != 0))
  
  list(expression=reference_expression, cellinfo=reference_cellinfo)
}

piecenet = get_piecenet(model$statenet)
piecestates = get_piecestates(piecenet)

simulations = smoothe_simulations(simulations)

pieces = divide_simulations(simulations, piecestates, model)

reference = get_reference_expression(pieces, model)



for (i in seq_len(length(simulations))) {
  cellcors = SCORPIUS::euclidean.distance(reference_expression, simulations[[i]]$expression)
  bestreferencecell = cellcors %>% apply(2, which.min)
  simulations[[i]]$cellinfo$progression = reference_cellinfo$progression[bestreferencecell]
  simulations[[i]]$cellinfo$piecestateid = reference_cellinfo$piecestateid[bestreferencecell]
}


# for plotting: subsample
for (i in seq_len(length(simulations))) {
  sample = sample(c(T, F), size=nrow(simulations[[i]]$expression), T, c(9999999, 1))
  
  sample[simulations[[i]]$expression[sample,] %>% apply(1, mean) %>% {. == 0}] = FALSE
  
  print(sum(sample))
  
  simulations[[i]]$subcellinfo = simulations[[i]]$cellinfo[sample,]
  simulations[[i]]$subexpression = simulations[[i]]$expression[sample,]
  
  simulations[[i]]$expression[sample,] %>% apply(1, sd) %>% {. == 0} %>% sum %>% print
  simulations[[i]]$subexpression %>% {apply(., 1, sd) == 0} %>% sum %>% print
  
  print("---")
}

overallexpression = map(simulations, "subexpression") %>% do.call(rbind, .)
overallexpression = overallexpression %>% set_rownames(1:nrow(overallexpression))
#overallcellinfo = map(simulations, ~.$subcellinfo) %>% bind_rows()
overallcellinfo = assign_progression(overallexpression, reference)
overallcellinfo$simulationid = map(seq_len(length(simulations)), ~rep(., nrow(simulations[[.]]$expression))) %>% unlist %>% factor()
overallcellinfo$observed = T

#overallexpression = rbind(overallexpression, reference_expression)
#overallcellinfo = bind_rows(overallcellinfo, tibble(simulationid=rep(NA, nrow(reference_cellinfo)), observed=FALSE))

space =  lmds(overallexpression, 3) %>% as.data.frame() %>% bind_cols(overallcellinfo)
space %>% filter(observed) %>% ggplot() + geom_path(aes(Comp1, Comp2, group=simulationid, color=progression)) + geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3)
space %>% filter(observed) %>% ggplot() + geom_path(aes(Comp1, Comp2, group=simulationid, color=progression)) + geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3) + viridis::scale_color_viridis(option="A")
space %>% filter(observed) %>% ggplot() + geom_path(aes(V1, V2, group=simulationid, color=piecestateid)) + geom_path(aes(V1, V2), data=space %>% filter(!observed), size=3)

for (i in unique(as.numeric(space$simulationid))) rgl::lines3d(space %>% filter(simulationid == i), col = rainbow(length(unique(space$simulationid)))[[i]])

aveRDS(space, file="space2.rds")


space = readRDS(file="space.rds")



#####

assign_progression = function(expression, reference) {
  cellcors = SCORPIUS::correlation.distance(reference$expression, expression)
  bestreferencecell = cellcors %>% apply(2, which.min)
  tibble(
    cell = colnames(cellcors),
    progression = reference$cellinfo$progression[bestreferencecell],
    piecestateid = reference$cellinfo$piecestateid[bestreferencecell]
  )
}


mat1 = matrix(runif(9), ncol=3)
mat2 = matrix(runif(12), ncol=3)

SCORPIUS::euclidean.distance(mat1, mat2) %>% mean
SCORPIUS::euclidean.distance(mat1[,3:1], mat2[, 3:1]) %>% mean
