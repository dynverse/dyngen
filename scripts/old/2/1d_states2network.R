states = tibble(from=c(1, 2), to=c(2, 3))
states = states %>% group_by(from) %>% summarize(singular=n() == 1) %>% right_join(states, by="from")
  



dig_deeper = function(states, curstate=1, visitedstates=c(), modules=list(modulenet=tibble(), modulenodes=tibble()), connecting_module=NULL) {
  modulenet = modules$modulenet
  modulenodes = modules$modulenodes
  to = states %>% filter(from == curstate) %>% .$to
  curmodule = ifelse(nrow(modulenodes), max(modulenodes$module), 0)
  if(length(to) > 1) {
    if(is.null(connecting_module)) stop("need a connecting module before a split...")
    receiving_modules = seq(curmodule+1, length.out=length(to)) %>% set_names(to)
    submodulenet = tibble(from=connecting_module, to=receiving_modules, strength=1, effect=1) %>% 
      bind_rows(
        combn(receiving_modules, 2) %>% cbind2(combn(receiving_modules %>% rev, 2)) %>% t %>% as.data.frame %>% rename(from=V1, to=V2) %>% mutate(strength=100, effect=-1)
      )
    submodulenodes = tibble(module=receiving_modules, state=curstate)
  } else {
    receiving_modules = (curmodule+1)
    if (length(to) > 0) {receiving_modules = set_names(receiving_modules, to)}
    if(!is.null(connecting_module)) {
      submodulenet = tibble(from=connecting_module, to=receiving_modules, strength=1, effect=1)
    } else {
      submodulenet = tibble()
    }
    submodulenodes = tibble(module=receiving_modules, state=curstate)
  }
  
  modulenet %<>% bind_rows(submodulenet)
  modulenodes %<>% bind_rows(submodulenodes)
  
  modules = list(modulenet=modulenet, modulenodes=modulenodes)
  
  visitedstates = c(visitedstates, curstate)
  for(nextstate in to) {
    if(nextstate %in% visitedstates) {
      print("reached this state already!!")
      stop("cycles don't work yet... you can only go forward!")
    } else {
      modules = dig_deeper(states, nextstate, visitedstates, modules, receiving_modules[as.character(nextstate)])
      
    }
  }
  return(modules)
}
  
modules = dig_deeper(states, 1)
modules$modulenodes %<>% mutate(a0=ifelse(state==1, 1, 0))
list2env(modules, .GlobalEnv)
