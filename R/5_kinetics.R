#' @export
simulation_setup_default <- function(
  sample_r = function(n) runif(n, 10, 200),
  sample_d = function(n) runif(n, 2, 8),
  sample_p = function(n) runif(n, 2, 8),
  sample_q = function(n) runif(n, 1, 5),
  
  sample_effect = function(n) sample(c(-1, 1), n, replace = TRUE, prob = c(.25, .75)),
  sample_strength = function(n) runif(n, 1, 20),
  sample_cooperativity = function(n) runif(n, 0.5, 2)
) {
  lst(
    sample_r,
    sample_d,
    sample_p,
    sample_q,
    
    sample_effect,
    sample_strength,
    sample_cooperativity
  )
}

#' @export
simulation_setup_custom <- function(
  sample_r = function(n) rnorm(n, 100, 20) %>% pmax(10),
  sample_d = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_p = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_q = function(n) rnorm(n, 3, .5) %>% pmax(1),
  
  sample_effect = function(n) sample(c(-1, 1), n, replace = TRUE, prob = c(.25, .75)),
  sample_strength = function(n) rexp(n, rate = .1),
  sample_cooperativity = function(n) runif(n, 0.5, 2)
) {
  lst(
    sample_r,
    sample_d,
    sample_p,
    sample_q,
    
    sample_effect,
    sample_strength,
    sample_cooperativity
  )
}

# why rlang %|%
`%|||%` <- function(x, y) {
  ifelse(is.na(x), y, x)
}

#' @export
generate_simulation_setup <- function(model) {
  # generate kinetics params
  model <- .kinetics_generate_gene_kinetics(model)
  
  # generate formulae
  formulae <- .kinetics_generate_formulae(model)
  
  # create variables
  model$feature_info$x <- paste0("x_", model$feature_info$feature_id)
  model$feature_info$y <- paste0("y_", model$feature_info$feature_id)
  molecule_ids <- c(model$feature_info$x, model$feature_info$y)

  initial_state <- set_names(
    rep(0, length(molecule_ids)),
    molecule_ids
  )
  
  # extract params
  parameters <- .kinetics_extract_parameters(model)
  
  # extract nus
  nus <- 
    formulae %>% 
    {Matrix::sparseMatrix(
      i = match(.$molecule, molecule_ids),
      j = seq_len(nrow(.)),
      x = .$effect, 
      dims = c(length(molecule_ids), nrow(.)),
      dimnames = list(molecule_ids, .$formula_id)
    )}
  
  # determine variables to be used during burn in
  burn_variables <- 
    model$feature_info %>% 
    filter(burn) %>% 
    select(x, y) %>% 
    gather(col, val) %>% 
    pull(val)
    
  # return system
  model$simulation_system <- lst(
    formulae, 
    molecule_ids,
    initial_state,
    parameters,
    nus,
    burn_variables
  )
  
  model
}

.kinetics_generate_gene_kinetics <- function(model) {
  if (model$verbose) cat("Generating kinetics for ", nrow(model$feature_info), " features\n", sep = "")
  params <- model$simulation_setup_params
  
  feature_info <-
    model$feature_info %>% 
    mutate(
      r = params$sample_r(n()),
      d = params$sample_d(n()),
      p = params$sample_p(n()),
      q = params$sample_q(n()),
      max_protein = r / d * p / q
    )
  
  # remove previous 'as' if present
  if (feature_info %has_name% "as") {
    feature_info <- feature_info %>% select(-as)
  }
  
  # sample effect, cooperativity and strength
  # calculate k
  feature_network <- 
    model$feature_network %>% 
    left_join(feature_info %>% select(from = feature_id, max_protein), by = "from") %>% 
    mutate(
      effect = effect %|||% params$sample_effect(n()),
      cooperativity = cooperativity %|||% params$sample_cooperativity(n()),
      strength = strength %|||% params$sample_strength(n()),
      k = .kinetics_calculate_k(max_protein, strength)
    )
  
  # calculate a0 and a
  feature_info <- 
    left_join(
      feature_info,
      feature_network %>% 
        rename(feature_id = to) %>% 
        group_by(feature_id) %>% 
        summarise(
          a0_2 = .kinetics_calculate_a0(effect),
          as = .kinetics_calculate_a(effect)
        ),
      by = "feature_id"
    ) %>% 
    mutate(
      a0 = a0 %|||% a0_2
    ) %>% 
    select(-a0_2)
  
  model$feature_info <- feature_info
  model$feature_network <- feature_network
  
  model
}

.kinetics_generate_formulae <- function(model) {
  if (model$verbose) cat("Generating formulae\n")
  
  # add helper information to feature info
  feature_info <- 
    model$feature_info %>% 
    left_join(
      model$feature_network %>% 
        group_by(feature_id = to) %>% 
        do({tibble(regulators = list(.))}) %>% 
        ungroup(),
      by = "feature_id"
    ) %>% 
    left_join(
      model$feature_network %>% 
        group_by(feature_id = from) %>% 
        summarise(num_targets = n()),
      by = "feature_id"
    )
  
  # generate formula per feature
  bind_rows(pbapply::pblapply(
    seq_len(nrow(feature_info)),
    cl = model$num_cores,
    function(i) {
      formulae <- list()
      
      info <- feature_info %>% extract_row_to_list(i)
      
      fid <- info$feature_id
      
      x <- paste0("x_", fid)
      y <- paste0("y_", fid)
      p <- paste0("p_", fid)
      q <- paste0("q_", fid)
      r <- paste0("r_", fid)
      d <- paste0("d_", fid)
      
      a0 <- paste0("a0_", fid)
      
      activation_function <- 
        if (!is.null(info$regulators)) {
          rid <- info$regulators$from
          reg_ys <- paste0("y_", rid)
          reg_ks <- paste0("k_", fid, "_", rid)
          reg_cs <- paste0("c_", fid, "_", rid)
          
          reg_affinities <- paste0("((", reg_ys, "/", reg_ks, ")**", reg_cs, ")")
          
          ix <- seq(1, 2 ^ length(reg_affinities) - 1)
          num_states <- length(ix)
          config_as <- paste0("a_", fid, "_", ix)
          
          config_affinities <- 
            map(
              seq_along(reg_affinities),
              function(i) {
                mask <- 2 ^ (i - 1)
                ifelse(bitwAnd(ix, mask) != 0, reg_affinities[[i]], NA)
              }
            ) %>% 
            do.call(cbind, .) %>% 
            apply(1, function(x) {
              na.omit(x) %>% paste(collapse = "*")
            })
          
          activation <- paste0(config_as, " * ", config_affinities, collapse = " + ")
          activation_denominator <- paste0(config_affinities, collapse = "+")
          
          # now make the activation function, by summing the activations
          activation_function <- paste0("(", a0, " + ", activation, ")/(1 + ", activation_denominator, ")")
          
        } else {
          a0
        }
      
      
      # mRNA production
      formulae[[length(formulae)+1]] <- 
        tibble(
          formula_id = paste0("mrna_production_", fid),
          formula = paste0(r, " * ", activation_function),
          effect = 1,
          molecule = x
        )
      
      # mRNA degradation
      formulae[[length(formulae)+1]] <- 
        tibble(
          formula_id = paste0("mrna_degradation_", fid),
          formula = paste0(d, " * ", x),
          effect = -1,
          molecule = x
        )
      
      # add protein synthesis if target gene also regulates other genes
      if (!is.na(info$num_targets) && info$num_targets > 0) {
        # protein production
        formulae[[length(formulae)+1]] <- 
          tibble(
            formula_id = paste0("protein_production_", fid),
            formula = paste0(p, " * ", x),
            effect = 1,
            molecule = y
          )
        
        # protein degradation
        formulae[[length(formulae)+1]] <- 
          tibble(
            formula_id = paste0("protein_degradation_", fid),
            formula = paste0(q, " * ", y),
            effect = -1,
            molecule = y
          )
      }
      
      bind_rows(formulae)
    }
  ))
}

.kinetics_extract_parameters <- function(model) {
  # extract r, d, p, q, and a0
  feature_params <- 
    model$feature_info %>% 
    select(feature_id, r, d, p, q, a0) %>% 
    gather(param, value, -feature_id) %>% 
    mutate(id = paste0(param, "_", feature_id)) %>% 
    select(id, value) %>% 
    deframe()
  
  # extract as
  as_params <- 
    model$feature_info %>% 
    select(feature_id, as) %>% 
    mutate(as = map(as, function(x) if (is.null(x)) numeric(0) else x), configs = map(as, seq_along)) %>% 
    unnest(as, configs) %>% 
    mutate(id = paste0("a_", feature_id, "_", configs)) %>% 
    select(id, as) %>% 
    deframe()
  
  # extract k and c
  edge_params <- 
    model$feature_network %>% 
    select(from, to, k, c = cooperativity) %>% 
    gather(param, value, -from, -to) %>% 
    mutate(id = paste0(param, "_", to, "_", from)) %>% 
    select(id, value) %>% 
    deframe()
  
  c(feature_params, as_params, edge_params)
}

.kinetics_calculate_k <- function(max_protein, strength) {
  max_protein / 2 / strength
}

.kinetics_calculate_a0 <- function(effects) {
  case_when(
    all(effects != 1) ~ 1,
    all(effects != -1) ~ 0.0001,
    TRUE ~ 0
  )
}

.kinetics_calculate_a <- function(effects) {
  ix <- seq(1, 2 ^ length(effects) - 1)
  num_states <- length(ix)
  a <- rep(1, num_states)
  
  for (i in seq_along(effects)) {
    if (effects[[i]] == -1) {
      mask <- 2 ^ (i-1)
      a[bitwAnd(ix, mask) != 0] <- 0
    }
  }
  list(a)
}