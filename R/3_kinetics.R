# helper functions to get configuration ids (of TF binding)
# get the configuration of a binding based on the configuration_id
get_binding_configuration <- function(configuration_id, n_regulators) {
  number2binary <- function(number, noBits) {
    binary_vector = rev(as.numeric(intToBits(number)))
    if(missing(noBits)) {
      return(binary_vector)
    } else {
      binary_vector[-(1:(length(binary_vector) - noBits))]
    }
  }
  rev(number2binary(configuration_id, n_regulators))
}

# get all configuration ids, eg. with 2 regulators -> 1,2,3, with one regulator -> 1, with no regulators -> numeric()
get_configuration_ids <- function(n) {
  if(n == 0) {
    numeric()
  } else {
    seq(1, 2**(n)-1)
  }
}

#' Randomize kinetics of genes
#' @param geneinfo
#' @param net
randomize_gene_kinetics <- function(geneinfo, net) {
  # sample r, d, p, q and k ----------------------
  
  sample_r <- function(n) runif(n, 10, 20)
  sample_d <- function(n) runif(n, 2, 8)
  sample_p <- function(n) runif(n, 2, 8)
  sample_q <- function(n) runif(n, 1, 5)
  sample_c <- function(n) runif(n, 1, 2)
  
  sample_k <- function(n, r, d, p, q) r/d * p/q / (runif(n, 2, 3))
  
  geneinfo <- geneinfo %>% mutate(
    r = sample_r(n()),
    d = sample_d(n()),
    p = sample_p(n()),
    q = sample_q(n()),
    c = sample_c(n()),
    k = sample_k(n(), r, d, p, q)
  )
  
  # sample a0 and a ---------------------------------------
  # include a list of effects in the geneinfo
  geneinfo <- geneinfo %>% 
    select(-matches("effects"), -matches("regulator_ids")) %>% # remove previously added effects
    left_join(
    net %>% group_by(to) %>% 
      summarise(effects = list(effect), regulator_ids = list(from))
    , by=c("gene_id"="to")
  )
  
  # for a0
  calculate_a0 <- function(effects) {
    if(length(effects) == 0) {
      1
    } else if(sum(effects == -1) == 0) {
      0.0001
    } else if(sum(effects == 1) == 0) {
      1
    } else {
      0.5
    }
  }
  
  # for calculating a
  # calculate a of a configuration based on the incoming active effects
  calculate_a <- function(configuration_id, effects) {
    bound <- get_binding_configuration(configuration_id, length(effects))
    
    if(any(effects[bound] == -1)) {
      0
    } else {
      1
    }
  }
  
  # calculate all a based on effects
  calculate_all_a <- function(effects) {
    # for every binding configuration
    configuration_ids <- get_configuration_ids(length(effects))
    map_dbl(configuration_ids, calculate_a, effects=effects) %>% set_names(configuration_ids)
  }
  
  geneinfo <- geneinfo %>% 
    mutate(
      a0 = ifelse(!is.na(a0), a0, map_dbl(effects, calculate_a0)),
      configuration_ids = map(effects, ~get_configuration_ids(length(.))),
      as = map(effects, calculate_all_a)
    )
  
  lst(geneinfo, net)
}

#' Randomize kinetics of cells
#' @param cells
randomize_cell_kinetics <- function(cells) {
  sample_kg <- function(n) 1
  sample_rg <- function(n) 1
  sample_pg <- function(n) 1
  sample_qg <- function(n) 1
  sample_dg <- function(n) 1
  
  cells %>% 
    mutate(
      kg = sample_kg(n()),
      rg = sample_kg(n()),
      pg = sample_kg(n()),
      qg = sample_kg(n()),
      dg = sample_kg(n())
    )
}

extract_params <- function(geneinfo, net, cells) {
  params <- c()
  # extract r, d, p, q, a0, k
  params <- geneinfo %>% select(r, d, p, q, a0, k, c, gene_id) %>% 
    gather("param_type", "param_value", -gene_id) %>% 
    mutate(param_id = glue::glue("{param_type}_{gene_id}")) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  params <- geneinfo %>% 
    unnest(as, configuration_ids) %>% 
    mutate(param_id=glue::glue("a_{.$gene_id}_{.$configuration_ids}"), param_value=as) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  params <- cells %>% select(kg, rg, pg, qg, dg, cell_id) %>% 
    gather("param_type", "param_value", -cell_id) %>% 
    mutate(param_id = glue::glue("{param_type}_{cell_id}")) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  params
}

generate_system <- function(net, geneinfo, cells) {
  # Randomize
  randomized_gene_kinetics <- randomize_gene_kinetics(model$geneinfo, model$net)
  geneinfo <- randomized_gene_kinetics$geneinfo
  net <- randomized_gene_kinetics$net
  
  cells <- randomize_cell_kinetics(cells)
  
  # Create variables
  geneinfo$x <- glue("x_{geneinfo$gene_id}")
  geneinfo$y <- glue("y_{geneinfo$gene_id}")
  molecule_ids <- c(geneinfo$x, geneinfo$y)
  
  # Extract formulae & kinetics
  formulae_changes <- generate_formulae(net, geneinfo, cells)
  initial_state <- rep(0, length(molecule_ids)) %>% setNames(molecule_ids)
  params <- extract_params(geneinfo, net, cells)
  
  # Extract nus
  formulae <- formulae_changes$formula %>% set_names(formulae_changes$formula_id)
  formulae_changes$molecule <- factor(formulae_changes$molecule, levels=molecule_ids) # fix order
  formulae_changes$formula_id <- factor(formulae_changes$formula_id, levels=names(formulae)) # fix order
  nus <- reshape2::acast(formulae_changes, molecule~formula_id, value.var="effect", fun.aggregate = first, fill = 0, drop = F)
  
  # Burn in
  burn_variables <- geneinfo %>% filter(as.logical(burn)) %>% select(x, y) %>% unlist() %>% unname()
  
  lst(formulae, initial_state, params, nus, burn_variables, molecule_ids)
}

## determine start state genes (active during burn-in)
determine_burn_variables = function(model) {
  variables_burn <- map(model$variables, "gene_id") %>% keep(~!is.null(.)) %>% unlist()
  variables_burn_genes <- names(variables_burn)[variables_burn %in% (model$geneinfo %>% filter(as.logical(burn)) %>% pull(gene_id))]
  # add other variables?
  intersect(variables_burn_genes, names(model$initial_state)) # only retain variables that change
}