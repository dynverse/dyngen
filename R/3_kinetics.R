#' @import stringr
generate_kinetics = function(vargroups, variables, nus_df, formulae) {
  ## generating the (initial) parameters of the system
  
  # parameters
  R = map_dbl(vargroups$r, ~20) %>% set_names(vargroups$r)
  D = map_dbl(vargroups$d, ~5) %>% set_names(vargroups$d)
  
  P = map_dbl(vargroups$p, ~2) %>% set_names(vargroups$p)
  Q = map_dbl(vargroups$q, function(q) 2*ifelse(!is.null(variables[[q]]$strength) && !is.na(variables[[q]]$strength), variables[[q]]$strength, 1)) %>% set_names(vargroups$q)
  
  K = map_dbl(vargroups$k, function(k) ifelse(is.na(variables[[k]]$strength),2, 2/variables[[k]]$strength)) %>% set_names(vargroups$k)
  
  RG = map_dbl(vargroups$rg, ~1) %>% set_names(vargroups$rg)
  DG = map_dbl(vargroups$dg, ~1) %>% set_names(vargroups$dg)
  KG = map_dbl(vargroups$kg, ~1) %>% set_names(vargroups$kg)
  PG = map_dbl(vargroups$pg, ~1) %>% set_names(vargroups$pg)
  QG = map_dbl(vargroups$qg, ~1) %>% set_names(vargroups$qg)
  
  A0 = map_dbl(variables[vargroups$a0], ~ifelse(is.null(.$value), 0, .$value)) %>% set_names(vargroups$a0)
  
  A = map_dbl(vargroups$a, function(a) ifelse(variables[[a]]$effect==1, 1, 0)) %>% set_names(vargroups$a)
  
  params = c(RG, DG, KG, PG, QG, R, D, K, P, Q, A0, A)
  
  # molecules
  initial_state = numeric()
  initial_state[vargroups$x] = 0
  initial_state[vargroups$y] = 0
  # initial_state[vargroups$u] = 1
  # initial_state[vargroups$b] = 0
  
  # fix order, to match initial state and formulae
  nus_df$molecule <- factor(nus_df$molecule, levels=names(initial_state))
  nus_df$formula_id <- factor(nus_df$formula_id, levels=names(formulae))
  
  nus <- nus_df %>% reshape2::acast(molecule~formula_id, fill=0, fun.aggregate = mean, value.var = "effect", drop=FALSE)
  
  if (any((nus %>% abs %>% apply(2, sum))==0)) {warning("Some molecules never change")}
  
  if(!is.null(vargroups$assoc)) {params = c(params, sapply(vargroups$assoc, function(assoc) 1 * ifelse(!is.null(variables[[assoc]]$strength) && !is.na(variables[[assoc]]$strength), variables[[assoc]]$strength, 1)))}
  if(!is.null(vargroups$deassoc)) {params = c(params, sapply(vargroups$deassoc, function(x) 1))}
  
  params = params[unique(names(params))] # remove duplicates, retain first
  if(any(is.na(params))) stop("Some parameters are NA!")
  
  lst(nus, params, initial_state)
}

## determine start state genes (active during burn-in)
determine_burn_variables = function(model) {
  variables_burn <- map(model$variables, "gene_id") %>% keep(~!is.null(.)) %>% unlist()
  variables_burn_genes <- names(variables_burn)[variables_burn %in% (model$geneinfo %>% filter(as.logical(burn)) %>% pull(gene_id))]
  # add other variables?
  intersect(variables_burn_genes, names(model$initial_state)) # only retain variables that change
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




# -------------------------- DELETE
#' Merge different models
#' 
#' @param model1 The first model
#' @param model2 The second model
#' 
#' @export
merge_models <- function(model1, model2) {
  model <- list()
  for(element in c("formulae.nus", "formulae", "formulae.strings", "params", "variables", "initial_state", "nus.changes")) {
    model[[element]] = model[[element]] <- c(model1[[element]], model2[[element]][setdiff(names(model2[[element]]), names(model1[[element]]))])
  }
  
  model$net <- bind_rows(model1$net, model2$net)
  model$geneinfo <- bind_rows(model2$geneinfo, model1$geneinfo) %>% group_by(gene) %>% filter(duplicated(gene) | n()==1) %>% ungroup()
  model$celltypes <- model1$celltypes
  model$modulemembership <- map(unique(names(model1$modulemembership), names(model2$modulemembership)), ~unique(c(model1$modulemembership[[.]], model2$modulemembership[[.]]))) %>% set_names(unique(names(model1$modulemembership), names(model2$modulemembership)))
  model$vargroups <- map(unique(names(model1$vargroups), names(model2$vargroups)), ~unique(c(model1$vargroups[[.]], model2$vargroups[[.]]))) %>% set_names(unique(names(model1$vargroups), names(model2$vargroups)))
  model$statenet <- bind_rows(model1$statenet, model2$statenet)
  model$net <- bind_rows(model1$net, model2$net)
  model$formulae.nus <- nus_to_matrix(model$nus.changes, names(model$initial_state))
  
  model$burngenes = determine_burngenes(model)
  
  model
}