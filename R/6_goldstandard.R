#' @import dplyr
get_module_counts = function(counts, modulemembership) {
  found = lapply(modulemembership, function(module) {all(module %in% colnames(counts))}) %>% unlist
  if(any(!found)) {
    stop(list(modulemembership, found, colnames(counts)) %>% str)
  }
  
  lapply(modulemembership, function(module) apply(counts[,as.character(module),drop=F], 1, mean)) %>% do.call(cbind, .)
}

#' @import dplyr
get_states_and_progression = function(expression_modules, modulenodes, expression_cutoff=2, expression_max=4) {
  modulesoi = modulenodes %>% filter(!is.na(state))
  
  expression_maxs = apply(expression_modules[, modulesoi$module], 2, quantile, probs=0.95)
  
  expression_modules[, modulesoi$module] %>% apply(1, function(x) {
    maxmodule = (1:length(x))[x>expression_cutoff] %>% last
    if (is.na(maxmodule)) {
      warning("no good state, picked max expression")
      maxmodule = which.max(x)
    }
    state = modulesoi$state[maxmodule]
    expression_max = expression_maxs[[maxmodule]]
    progression = pmax(pmin((x[state] - expression_cutoff)/(expression_max-expression_cutoff), 1), 0)
    
    tibble(state=state, progression=progression)
  }) %>% bind_rows %>% mutate(cell=rownames(expression_modules))
}

## Get cell states and progression

# shouldn't we start from the end and retrace our steps? What if due to stochsticity a primed cell still choses for another path?
#' @import dplyr
get_states = function(expression_modules, model, minexpression = 0.8, minratio = 1.5, state=1) {
  statechanged = T
  states = state
  for (cell in rownames(expression_modules)[2:nrow(expression_modules)]) {
    if (statechanged) {
      nextstates = model$statenet %>% filter(from == state) %>% .$to
      state2module = model$modulenodes %>% filter(!is.na(state)) %>% {set_names(.$module, .$state)}
      nextmodules = state2module[nextstates %>% as.character]
    }
    
    statechanged = F
    if(length(nextmodules) == 1) {
      if(expression_modules[cell,nextmodules] > minexpression) {
        state = nextstates[[1]]
        statechanged = T
      }
    } else if (length(nextmodules) > 1) {
      expression_modules_cell = expression_modules[cell,nextmodules]
      if(any(expression_modules_cell > minexpression)) {
        ratios = expression_modules_cell[which.max(expression_modules_cell)]/expression_modules_cell[-which.max(expression_modules_cell)]
        if(all(ratios > minratio)) {
          state = nextstates[[which.max(expression_modules_cell)]]
          statechanged = T
        }
      }
    }
    states = c(states, state)
  }
  
  tibble(state=states, cell=rownames(expression_modules))
}
















## Get cell distances in gold standard
# Goal: get a pairwise distance matrix between all cells, based on their branch and progression in the gold standard
#' @import dplyr
get_line_graph = function(net) {
  g = net %>% igraph::graph_from_data_frame() %>% igraph::make_line_graph()
  if(length(igraph::V(g)) == 1) {
    return(tibble(from=1, to=-1)) # if there is only one state, add a "pseudo-state", so that igraph stops whining
  }
  g %>% igraph::as_data_frame()
}


#' @import dplyr
get_branchconnected_statenet = function(statenet, statenodes) {
  # network between states
  # at bifurcation points, the two states after a bifurcation are connected through a head2head edge
  statenet = statenet %>% mutate(type = "head2tail")
  statenet = statenet %>% group_by(from) %>% summarise(n=n()) %>% filter(n>1) %>% .$from %>% map(function(fromoi) {
    statesoi = statenet %>% filter(from==fromoi) %>% .$to
    gtools::combinations(length(statesoi), r=2, v=statesoi) %>% as.data.frame() %>% rename(from=V1, to=V2)
  }) %>% bind_rows() %>% mutate(type="head2head") %>% bind_rows(statenet)
  print(statenet$type)
  if(nrow(statenet) > 1) igraph::graph_from_data_frame(statenet[,c("from", "to")]) %>% plot(edge.color=factor(statenet$type))
  snet = statenet[,c("from", "to")] %>% igraph::graph_from_data_frame(vertices=unique(c(statenet$from, statenet$to)) %>% sort)
  list(snet=snet, statenet=statenet, statenodes=statenodes)
}


# head2head = weight 2 (double the path)
# head2tail = weight 1 ()
# start and end state: check 

# the distance if the states are the same is just the difference in progression time within this stage
traj_distance_samestate = function(time1, time2, verbose=F) abs(time1-time2)

#' @import dplyr
get_state_edge = function(state1, state2, statenet) {
  statenet %>% filter((from == state1 & to==state2) | (from == state2 & to == state1)) %>% as.list()
}

# processessing the shortest path between 2 states through the state network
# returns a function which can be used to get the distance of two cells in these two respective states
# we calculate these functions here so that they do not need to be recalculated for every pair of cells in the same two states
#' @import dplyr
get_maxprogressions = function(snet, states) {
  snet$statenodes %>% filter(state %in% states) %>% .$maxprogression %>% sum
}

#' @import dplyr
process_sp = function(sp, snet) {
  statenet = snet$statenet
  extracost = 0
  if(length(sp) == 1) {
    traj_distance_samestate
  } else {
    edge = get_state_edge(sp[[1]], sp[[2]], snet$statenet)
    if(edge$type == "head2head" || edge$from == sp[[2]]) {
      start = "from0"
    } else {
      start = "from1"
    }
    if(edge$type == "head2head" && length(sp) > 2) {
      #extracost = 1
    }
    
    edge = get_state_edge(sp[[length(sp)-1]], sp[[length(sp)]], snet$statenet)
    if(edge$type == "head2head" || edge$from == sp[[length(sp)-1]]) {
      end = "from0"
    } else {
      end = "from1"
    }
    
    extracost = extracost + get_maxprogressions(snet, sp[c(-1, -length(sp))])
    start_maxprogression = get_maxprogressions(snet, sp[[1]])
    end_maxprogression = get_maxprogressions(snet, sp[[length(sp)]])
    
    function(time1, time2, verbose=F) {
      #if(verbose) print(time1);print(time2);print(extracost)
      totalcost = extracost
      if(start == "from0") {
        totalcost = totalcost + time1
      } else {
        totalcost = totalcost + start_maxprogression - time1
      }
      if(end == "from0") {
        totalcost = totalcost + time2
      } else {
        totalcost = totalcost + end_maxprogression - time2
      }
      totalcost
    }
  }
}
get_cell_distance_func = function(snet, state1, state2) {
  state1 = as.character(state1) # because igraph works with character names and order of vertices is not numerical
  state2 = as.character(state2)
  igraph::shortest_paths(snet$snet, state1, state2, mode = "all")$vpath[[1]] %>% names() %>% as.integer() %>% process_sp(snet)
}
#func = get_cell_distance_func(snet, 1, 2)
#func(0.9, 0.1)
#func(0.2, 0.3)

#' @import dplyr
#' @import purrr
get_cell_distances = function(cellinfo, snet) {
  cellinfo$state = factor(cellinfo$state) # make sure this is a factor, because if some states are missing this will mess up the indexing
  allstates = unique(cellinfo$state) %>% sort
  cell_distance_funcs = map(allstates, function(state1) map(allstates, ~get_cell_distance_func(snet, state1, .)))
  
  get_cell_distance = function(cell_distance_funcs, cell1, cell2) {
    time1 = cellinfo$progression[[cell1]]
    state1 = cellinfo$state[[cell1]]
    time2 = cellinfo$progression[[cell2]]
    state2 = cellinfo$state[[cell2]]
    
    #print(list(time1, state1, time2, state2))
    
    cell_distance_funcs[[state1]][[state2]](time1, time2)
  }
  
  allcells = 1:nrow(cellinfo)
  cell_distances = lapply(allcells, function(cell1) {
    time1 = cellinfo$progression[[cell1]]
    state1 = cellinfo$state[[cell1]]
    lapply(allcells, function(cell2) {
      time2 = cellinfo$progression[[cell2]]
      state2 = cellinfo$state[[cell2]]
      cell_distance_funcs[[state1]][[state2]](time1, time2)
    }) %>% as.double
  }) %>% {matrix(unlist(.), ncol=length(allcells))}
  
  cellsoi = cellinfo$state %in% c(1,2)
  cell_distances[] %>% pheatmap::pheatmap(cluster_cols=F, cluster_rows=F)
  
  dimnames(cell_distances) = list(cellinfo$cell, cellinfo$cell)
  
  cell_distances
}
