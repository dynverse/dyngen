#' Generate a random tree
#' 
#' @param treeseed A seed to set
#' @param decay The decay rate
#' 
#' @importFrom stats runif
generate_random_tree <- function(treeseed=NULL, decay=1.4) {
  if(!is.null(treeseed)) set.seed(treeseed)
  
  recurse <- function(net, level=0, parentstage=2, decay=1.4, parentstate = "S1") {
    if (decay <= 1) stop("decay has to be > 1")
    for (i in c(1,2)) {
      maxstage <- ifelse(nrow(net), max(net$to), 0)
      if (stats::runif(1) < decay^(-level)) {
        net <- net %>% 
          add_row(from = parentstage, to = maxstage+1) %>% 
          recurse(level + 1, maxstage+1, decay)
      }
    }
    return(net)
  }
  
  initial_net <- data_frame(from = 1, to = 2)
  
  recurse(initial_net) %>% 
    mutate_all(as.character) %>% 
    group_by(from) %>% 
    mutate(singular=n() == 1)
}

#' Convert tree to states list, used later by the gold standard to determine where there are bifurcations
#' The modules of each bifurcation will be added later by `from_stages_to_modulenet`
#' 
#' @param stagenet The state network to be converted.
generate_tree_states <- function(stagenet) {
  bifurcating_froms <- table(stagenet$from) %>% keep(~.>1) %>% names
  stage2state <- paste0("S", seq_along(bifurcating_froms)) %>% set_names(bifurcating_froms)
  states <- list()
  for (fromoi in bifurcating_froms) {
    tos <- stagenet %>% filter(from == fromoi) %>% pull(to)
    stage2state <- c(stage2state, paste0("S", seq(length(stage2state)+1, length.out = length(tos))) %>% set_names(tos))
    states <- c(states, list(list("state_id"=stage2state[[fromoi]], "type"="bifurcation", to=stage2state[tos])))
  }
  states <- states %>% set_names(map_chr(states, "state_id"))
  lst(stage2state, states)
}

#' Convert a stage network to a module network
#' 
#' @param stagenet A network of stages
#' 
#' @importFrom utils combn
from_stages_to_modulenet <- function(stagenet) {
  statesinfo <- generate_tree_states(stagenet)
  stage2state <- statesinfo$stage2state
  states <- statesinfo$states
  
  modulename_generator <- function(modulenodes, n) {paste0("M", seq(nrow(modulenodes)+1, length.out=n))}
  
  # recursive function to generate the module network
  dig_deeper <- function(
    stagenet, 
    curstage="1", 
    visitedstages=c(), 
    modules=list(
      modulenet=tibble(), 
      modulenodes=tibble(module_id="M1", stage="1"),
      states = states
    ), 
    connecting_module="M1"
  ) {
    modulenet <- modules$modulenet
    modulenodes <- modules$modulenodes
    states <- modules$states
    
    to <- stagenet %>% filter(from == curstage) %>% .$to
    
    receiving_modules <- modulename_generator(modulenodes, length(to)) %>% set_names(to)
    # if splitting
    if(length(to) > 1) {
      if(is.null(connecting_module)) stop("need a connecting module before a split...")
      # repression between modules of the different stages
      submodulenet <- tibble(from=connecting_module, to=receiving_modules, strength=1, effect=1) %>% 
        bind_rows(
          combn(receiving_modules, 2) %>% cbind2(utils::combn(receiving_modules %>% rev, 2)) %>% t %>% as.data.frame(stringsAsFactors=FALSE) %>% rename(from=V1, to=V2) %>% mutate(strength=100, effect=-1)
        )
      submodulenodes <- tibble(module_id=receiving_modules, stage=to)
      
      # add module info to the states
      states[[stage2state[[curstage]]]]$modules <- receiving_modules
      
    # if not splitting
    } else {
      if (length(to) > 0) {receiving_modules <- set_names(receiving_modules, to)}
      if(!is.null(connecting_module)) {
        submodulenet <- tibble(from=connecting_module, to=receiving_modules, strength=1, effect=1)
      } else {
        submodulenet <- tibble()
      }
      submodulenodes <- tibble(module_id=receiving_modules, stage=to)
    }
    
    modulenet <- modulenet %>% bind_rows(submodulenet)
    modulenodes <- modulenodes %>% bind_rows(submodulenodes)
    
    modules <- lst(modulenet, modulenodes, states)
    
    # now dig deeper into the next stages
    visitedstages = c(visitedstages, curstage)
    for(nextstage in to) {
      if(nextstage %in% visitedstages) {
        print("reached this stage already!!")
        stop("cycles don't work yet... you can only go forward!")
      } else if (nrow(stagenet %>% filter(from == nextstage))){
        modules = dig_deeper(stagenet, nextstage, visitedstages, modules, receiving_modules[as.character(nextstage)])
        
      }
    }
    return(modules)
  }
  
  modules <- dig_deeper(stagenet, "1")
  modulenodes <- modules$modulenodes
  modulenet <- modules$modulenet
  states <- modules$states
  
  # process module information
  modulenodes <- modulenodes %>% 
    mutate(
      a0=ifelse(stage == 1, 1, 0),
      cell_id = 1,
      burn=ifelse(stage == 1, TRUE, FALSE)
    )
  modulenet <- modulenet %>% 
    mutate(
      cooperativity = 2,
      randomize=TRUE
    )
  
  cells <- tibble(cell_id=1, dies=FALSE)
  lst(cells, modulenodes, modulenet, states)
}
