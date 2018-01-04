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
  as.logical(rev(number2binary(configuration_id, n_regulators)))
}

# get all configuration ids, eg. with 2 regulators -> 1,2,3, with one regulator -> 1, with no regulators -> numeric()
get_configuration_ids <- function(n) {
  if(n == 0) {
    numeric()
  } else {
    seq(1, 2**(n)-1)
  }
}

get_default_kinetics_samplers <- function() {
  list(
    sample_r = function(n) runif(n, 10, 200), 
    sample_d = function(n) runif(n, 2, 8), 
    sample_p = function(n) runif(n, 2, 8), 
    sample_q = function(n) runif(n, 1, 5),
    calculate_a0 = function(effects) {
      if(length(effects) == 0) {
        1
      } else if(sum(effects == -1) == 0) {
        0.0001
      } else if(sum(effects == 1) == 0) {
        1
      } else {
        0.5
      }
    },
    calculate_a = function(configuration_id, effects) {
      bound <- get_binding_configuration(configuration_id, length(effects))
      if(any(effects[bound] == -1)) {
        0
      } else {
        1
      }
    },
    sample_strength = function(n) runif(n, 1, 20),
    calculate_k = function(max_protein, strength) {
      max_protein/2/strength
    },
    sample_cooperativity = function(n) runif(n, 1, 4)
  )
}

#' Randomize kinetics of genes
#' @rdname generate_system
#' @export
randomize_gene_kinetics <- function(
  geneinfo, 
  net, 
  samplers = get_default_kinetics_samplers()) {
  # sample r, d, p, q and k ----------------------
  
  geneinfo <- geneinfo %>% mutate(
    r = samplers$sample_r(n()),
    d = samplers$sample_d(n()),
    p = samplers$sample_p(n()),
    q = samplers$sample_q(n())
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
  
  # for a0 and a
  
  # calculate all a based on effects
  calculate_all_a <- function(effects) {
    # for every binding configuration
    configuration_ids <- get_configuration_ids(length(effects))
    map_dbl(configuration_ids, samplers$calculate_a, effects=effects) %>% set_names(configuration_ids)
  }
  
  geneinfo <- geneinfo %>% 
    mutate(
      a0 = ifelse(!is.na(a0), a0, map_dbl(effects, samplers$calculate_a0)),
      configuration_ids = map(effects, ~get_configuration_ids(length(.))),
      as = map(effects, calculate_all_a)
    )
  
  # calculate interaction cooperativity and binding strength
  # calculate k
  geneinfo <- geneinfo %>% mutate(max_protein = r/d * p/q) # calculate maximal protein
  net <- net %>% 
    left_join(geneinfo %>% select(max_protein, gene_id), by=c("from" = "gene_id")) %>%  # add maximal protein of regulator
    mutate(
      strength = ifelse(is.na(strength), samplers$sample_strength(n()), strength),
      k = samplers$calculate_k(max_protein, strength)
    )
  
  # calculate c
  net <- net %>% mutate(
    c = ifelse(is.na(cooperativity), samplers$sample_cooperativity(n()), cooperativity)
  )
  
  lst(geneinfo, net)
}

#' Randomize kinetics of cells
#' @rdname generate_system
#' @export
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
  # extract r, d, p, q, a0
  params <- geneinfo %>% select(r, d, p, q, a0, gene_id) %>% 
    gather("param_type", "param_value", -gene_id) %>% 
    mutate(param_id = glue::glue("{param_type}_{gene_id}")) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  # extract a
  params <- geneinfo %>% 
    unnest(as, configuration_ids) %>% 
    mutate(param_id=glue::glue("a_{.$gene_id}_{.$configuration_ids}"), param_value=as) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  # extract k and c
  params <- net %>% select(c, k, from, to) %>% 
    gather("param_type", "param_value", -from, -to) %>% 
    mutate(param_id = glue::glue("{param_type}_{to}_{from}")) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  # extract cell parameters
  params <- cells %>% select(kg, rg, pg, qg, dg, cell_id) %>% 
    gather("param_type", "param_value", -cell_id) %>% 
    mutate(param_id = glue::glue("{param_type}_{cell_id}")) %>% 
    {set_names(.$param_value, .$param_id)} %>% 
    c(params)
  
  params
}

#' Generate a system from a network using 
#' 
#' @param net Regulatory network dataframe
#' @param geneinfo Geneinfo dataframe
#' @param cells Cells dataframe
#' @param samplers The samplers for the kinetics parameters
#' @export
generate_system <- function(net, geneinfo, cells, samplers) {
  # Randomize
  randomized_gene_kinetics <- randomize_gene_kinetics(geneinfo, net, samplers)
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
  
  lst(formulae, initial_state, params, nus, burn_variables, molecule_ids, geneinfo, net, cells)
}