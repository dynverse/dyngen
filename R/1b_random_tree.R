generate_random_tree = function(seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  dig_deeper = function(net=tibble(from=c(1), to=c(2)), level=0, parentstate=2, decay=1.4) {
    if(decay <= 1) stop("decay has to be > 1")
    for (i in c(1,2)) {
      maxstate = ifelse(nrow(net), max(net$to), 0)
      if(runif(1) < decay^(-level)) {
        net = net %>% bind_rows(tibble(from=parentstate, to=maxstate+1))
        net = dig_deeper(net, level+1, maxstate+1, decay)
      }
    }
    return(net)
  }
  
  states = dig_deeper()
  states = states %>% group_by(from) %>% summarize(singular=n() == 1) %>% right_join(states, by="from")
  igraph::graph_from_data_frame(states[, c("from", "to")]) %>% plot
  
  states
}

#' @importFrom utils combn
from_states_to_modulenet = function(states) {
  dig_deeper = function(states, curstate=1, visitedstates=c(), modules=list(modulenet=tibble(), modulenodes=tibble(module=1, state=1)), connecting_module=1) {
    modulenet = modules$modulenet
    modulenodes = modules$modulenodes
    to = states %>% filter(from == curstate) %>% .$to
    curmodule = ifelse(nrow(modulenodes), max(modulenodes$module), 0)
    if(length(to) > 1) {
      if(is.null(connecting_module)) stop("need a connecting module before a split...")
      receiving_modules = seq(curmodule+1, length.out=length(to)) %>% set_names(to)
      submodulenet = tibble(from=connecting_module, to=receiving_modules, strength=1, effect=1) %>% 
        bind_rows(
          combn(receiving_modules, 2) %>% cbind2(utils::combn(receiving_modules %>% rev, 2)) %>% t %>% as.data.frame %>% rename(from=V1, to=V2) %>% mutate(strength=100, effect=-1)
        )
      submodulenodes = tibble(module=receiving_modules, state=to)
    } else {
      receiving_modules = (curmodule+1)
      if (length(to) > 0) {receiving_modules = set_names(receiving_modules, to)}
      if(!is.null(connecting_module)) {
        submodulenet = tibble(from=connecting_module, to=receiving_modules, strength=1, effect=1)
      } else {
        submodulenet = tibble()
      }
      submodulenodes = tibble(module=receiving_modules, state=to)
    }
    
    modulenet %<>% bind_rows(submodulenet)
    modulenodes %<>% bind_rows(submodulenodes)
    
    modules = list(modulenet=modulenet, modulenodes=modulenodes)
    
    visitedstates = c(visitedstates, curstate)
    for(nextstate in to) {
      if(nextstate %in% visitedstates) {
        print("reached this state already!!")
        stop("cycles don't work yet... you can only go forward!")
      } else if (nrow(states %>% filter(from == nextstate))){
        modules = dig_deeper(states, nextstate, visitedstates, modules, receiving_modules[as.character(nextstate)])
        
      }
    }
    return(modules)
  }
  
  modules = dig_deeper(states, 1)
  modules$modulenodes %<>% mutate(a0=ifelse(state==1, 1, 0))
  
  modules$modulenet = modules$modulenet %>% mutate(cooperativity = 2, from=paste0("M", from), to=paste0("M", to))
  modules$modulenodes$module = paste0("M", modules$modulenodes$module)
  modules$modulenodes$celltype = 1
  modules$modulenodes$burn = c(1, seq_len(nrow(modules$modulenodes)-1))
  c(list(celltypes=tibble(celltype=1, dies=FALSE), statenet=states %>% select(-singular)), modules)
}
