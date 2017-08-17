#' @import dplyr
#' @import stringr
#' @import dambiutils
generate_kinetics = function(vargroups, variables, nus.changes) {
  ## generating the (initial) parameters of the system
  
  R = sapply(vargroups$r, function(r) 20)
  D = sapply(vargroups$d, function(d) 5)
  
  P = sapply(vargroups$p, function(p) 2)
  Q = sapply(vargroups$q, function(q) 2*ifelse(!is.null(variables[[q]]$strength) && !is.na(variables[[q]]$strength), variables[[q]]$strength, 1))
  
  K = sapply(vargroups$k, function(k) ifelse(is.na(variables[[k]]$strength),2, 2/variables[[k]]$strength))
  
  g = c(1)
  #P[paste0("p_", g)] = P[paste0("p_", g)] /5
  #Q[paste0("q_", g)] = Q[paste0("q_", g)] /5
  
  A0 = map_dbl(variables[vargroups$a0], ~ifelse(is.null(.$value), 0, .$value))
  
  #A0[[1]] = 1
  #A0[[4]] = 0
  #A0[paste0("a0_", ldtfs)] = 1 # IF CYCLE
  #A0[modules[[1]]] = 1
  
  A = sapply(vargroups$a, function(a) ifelse(variables[[a]]$effect==1, 1, 0))
  
  initial.state = numeric()
  initial.state[vargroups$x] = 0
  initial.state[vargroups$y] = 0
  initial.state[vargroups$u] = 1
  initial.state[vargroups$b] = 0
  initial.state[vargroups$rg] = 1
  
  formulae.nus = nus_to_matrix(nus.changes, names(initial.state))
  
  params = c(R, D, K, P, Q, A0,A, rg=1,dg=1,kg=1,pg=1,qg=1, a1=1)
  
  if(!is.null(vargroups$assoc)) {params = c(params, sapply(vargroups$assoc, function(assoc) 1 * ifelse(!is.null(variables[[assoc]]$strength) && !is.na(variables[[assoc]]$strength), variables[[assoc]]$strength, 1)))}
  if(!is.null(vargroups$deassoc)) {params = c(params, sapply(vargroups$deassoc, function(x) 1))}
  
  params = params[unique(names(params))] # remove duplicates, retain first
  
  if(any(is.na(params))) stop("Some parameters are NA!")
  
  tibble::lst(formulae.nus, params, initial.state)
}

## determine start state genes (active during burn-in)
#' @import dplyr
#' @importFrom purrr %>% map map_int
determine_burngenes = function(model) {
  variables_genes = map(model$variables, "gene") %>% keep(~!is.null(.)) %>% unlist()
  (variables_genes %in% (model$geneinfo %>% filter(as.logical(burn)) %>% .$gene)) %>% variables_genes[.]
}



nus_to_matrix <- function(nus.changes, allmolecules) {
  formulae.nus = sapply(nus.changes, function(nu.changes) {
    nu = sapply(allmolecules, function(r) 0)
    nu[names(nu.changes)] = nu.changes
    nu
  })
  if(class(formulae.nus) != "matrix") stop("formulae.nus shoudl become matrix")
  
  formulae.nus
}





#' Merge different models
#' @export
merge_models <- function(model1, model2) {
  model <- list()
  for(element in c("formulae.nus", "formulae", "formulae.strings", "params", "variables", "initial.state", "nus.changes")) {
    model[[element]] = model[[element]] <- c(model1[[element]], model2[[element]][setdiff(names(model2[[element]]), names(model1[[element]]))])
  }
  
  model$net <- bind_rows(model1$net, model2$net)
  model$geneinfo <- bind_rows(model2$geneinfo, model1$geneinfo) %>% group_by(gene) %>% filter(duplicated(gene) | n()==1) %>% ungroup()
  model$celltypes <- model1$celltypes
  model$modulemembership <- map(unique(names(model1$modulemembership), names(model2$modulemembership)), ~unique(c(model1$modulemembership[[.]], model2$modulemembership[[.]]))) %>% set_names(unique(names(model1$modulemembership), names(model2$modulemembership)))
  model$vargroups <- map(unique(names(model1$vargroups), names(model2$vargroups)), ~unique(c(model1$vargroups[[.]], model2$vargroups[[.]]))) %>% set_names(unique(names(model1$vargroups), names(model2$vargroups)))
  model$statenet <- bind_rows(model1$statenet, model2$statenet)
  model$net <- bind_rows(model1$net, model2$net)
  model$formulae.nus <- nus_to_matrix(model$nus.changes, names(model$initial.state))
  
  model$burngenes = determine_burngenes(model)
  
  model
}