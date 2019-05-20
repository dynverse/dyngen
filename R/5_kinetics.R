#' @export
kinetics_default <- function(
  sample_wpr = function(n) runif(n, 10, 200),
  sample_wdr = function(n) runif(n, 2, 8),
  sample_xpr = function(n) runif(n, 2, 8),
  sample_xdr = function(n) runif(n, 2, 8),
  sample_ydr = function(n) runif(n, 2, 8),
  sample_ypr = function(n) runif(n, 1, 5),
  
  sample_effect = function(n) sample(c(-1, 1), n, replace = TRUE, prob = c(.25, .75)),
  sample_strength = function(n) runif(n, 1, 20),
  sample_cooperativity = function(n) runif(n, 0.5, 2)
) {
  lst(
    sample_wpr,
    sample_wdr,
    sample_xpr,
    sample_xdr,
    sample_ypr,
    sample_ydr,
    
    sample_effect,
    sample_strength,
    sample_cooperativity
  )
}

#' @export
kinetics_custom <- function(
  sample_wpr = function(n) rnorm(n, 100, 20) %>% pmax(10),
  sample_wdr = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_xpr = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_xdr = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_ypr = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_ydr = function(n) rnorm(n, 3, .5) %>% pmax(1),
  
  sample_effect = function(n) sample(c(-1, 1), n, replace = TRUE, prob = c(.25, .75)),
  sample_strength = function(n) runif(n, 10, 40),
  sample_cooperativity = function(n) runif(n, 0.5, 2)
) {
  lst(
    sample_wpr,
    sample_wdr,
    sample_xpr,
    sample_xdr,
    sample_ypr,
    sample_ydr,
    
    sample_effect,
    sample_strength,
    sample_cooperativity
  )
}

#' @export
kinetics_fixed <- function(
  sample_wpr = function(n) rep(20, n),
  sample_wdr = function(n) rep(10, n),
  sample_xpr = function(n) rep(10, n),
  sample_xdr = function(n) rep(10, n),
  sample_ypr = function(n) rep(10, n),
  sample_ydr = function(n) rep(1, n),
  
  sample_effect = function(n) rep(1, n),
  sample_strength = function(n) rep(1, n),
  sample_cooperativity = function(n) rep(2, n)
) {
  lst(
    sample_wpr,
    sample_wdr,
    sample_xpr,
    sample_xdr,
    sample_ypr,
    sample_ydr,
    
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
generate_kinetics <- function(model) {
  # generate kinetics params
  model <- .kinetics_generate_gene_kinetics(model)
  
  # generate formulae
  formulae <- .kinetics_generate_formulae(model)
  
  # create variables
  fid <- model$feature_info$feature_id
  model$feature_info$w <- paste0("w_", fid)
  model$feature_info$x <- paste0("x_", fid)
  model$feature_info$y <- paste0("y_", fid)
  
  molecule_ids <- c(model$feature_info$w, model$feature_info$x, model$feature_info$y)

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
    select(w, x, y) %>% 
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
  params <- model$kinetics_params
  
  feature_info <-
    model$feature_info %>% 
    mutate(
      wpr = params$sample_wpr(n()),
      wdr = params$sample_wdr(n()),
      xpr = params$sample_xpr(n()),
      xdr = params$sample_xdr(n()),
      ypr = params$sample_ypr(n()),
      ydr = params$sample_ydr(n()),
      max_protein = wpr / wdr * xpr / xdr * ypr / ydr
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
      k = max_protein / 2 / strength
    )
  
  # calculate a0 and a
  feature_info <- 
    left_join(
      feature_info,
      feature_network %>% 
        rename(feature_id = to) %>% 
        group_by(feature_id) %>% 
        summarise(
          a0_2 = .kinetics_calculate_a0(effect)
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
    function(i) {
      formulae <- list()
      
      info <- feature_info %>% extract_row_to_list(i)
      
      fid <- info$feature_id
      
      w <- paste0("w_", fid)
      x <- paste0("x_", fid)
      y <- paste0("y_", fid)
      
      wpr <- paste0("wpr_", fid)
      wdr <- paste0("wdr_", fid)
      xpr <- paste0("xpr_", fid)
      xdr <- paste0("xdr_", fid)
      ypr <- paste0("ypr_", fid)
      ydr <- paste0("ydr_", fid)
      
      a0 <- paste0("a0_", fid)
      
      wpr_function <- 
        if (!is.null(info$regulators)) {
          rid <- info$regulators$from
          eff <- info$regulators$effect
          reg_ys <- paste0("y_", rid)
          reg_ks <- paste0("k_", fid, "_", rid)
          reg_cs <- paste0("c_", fid, "_", rid)
          buffer_var <- paste0("buffer_", fid, "_", rid)
          
          reg_affinity_calc <- paste(paste0(buffer_var, " = 1 + pow(", reg_ys, "/", reg_ks, ", ", reg_cs, "); "), collapse = "")
          
          numerator <-
            if (sum(eff > 0) > 0) {
              paste0(a0, " - 1 + ", paste(buffer_var[eff > 0], collapse = " * "))
            } else {
              a0
            }
          denominator <- paste(buffer_var, collapse = " * ")
          
          paste0(reg_affinity_calc, wpr, " * (", numerator, ")/(", denominator, ")")
        } else {
          paste0(wpr, " * ", a0)
        }
      
      
      # pre-mRNA production
      formulae[[length(formulae)+1]] <- 
        tibble(
          formula_id = paste0("premrna_production_", fid),
          formula = wpr_function,
          effect = 1,
          molecule = w
        )
      
      # pre-mRNA degradation
      formulae[[length(formulae)+1]] <- 
        tibble(
          formula_id = paste0("premrna_degradation_", fid),
          formula = paste0(wdr, " * ", w),
          effect = -1,
          molecule = w
        )
      
      # mRNA production
      formulae[[length(formulae)+1]] <- 
        tibble(
          formula_id = paste0("mrna_production_", fid),
          formula = paste0(xpr, " * ", w),
          effect = 1,
          molecule = x
        )
      
      # mRNA degradation
      formulae[[length(formulae)+1]] <- 
        tibble(
          formula_id = paste0("mrna_degradation_", fid),
          formula = paste0(xdr, " * ", x),
          effect = -1,
          molecule = x
        )
      
      # add protein synthesis if target gene also regulates other genes
      if (!is.na(info$num_targets) && info$num_targets > 0) {
        # protein production
        formulae[[length(formulae)+1]] <- 
          tibble(
            formula_id = paste0("protein_production_", fid),
            formula = paste0(ypr, " * ", x),
            effect = 1,
            molecule = y
          )
        
        # protein degradation
        formulae[[length(formulae)+1]] <- 
          tibble(
            formula_id = paste0("protein_degradation_", fid),
            formula = paste0(ydr, " * ", y),
            effect = -1,
            molecule = y
          )
      }
      
      bind_rows(formulae)
    }
  ))
}

.kinetics_extract_parameters <- function(model) {
  # extract m to qr, d, p, q, and a0
  feature_params <- 
    model$feature_info %>% 
    select(feature_id, wpr, wdr, xpr, xdr, ypr, ydr, a0) %>% 
    gather(param, value, -feature_id) %>% 
    mutate(id = paste0(param, "_", feature_id)) %>% 
    select(id, value) %>% 
    deframe()
  
  # extract k and c
  edge_params <- 
    model$feature_network %>% 
    select(from, to, k, c = cooperativity) %>% 
    gather(param, value, -from, -to) %>% 
    mutate(id = paste0(param, "_", to, "_", from)) %>% 
    select(id, value) %>% 
    deframe()
  
  c(feature_params, edge_params)
}

.kinetics_calculate_a0 <- function(effects) {
  case_when(
    all(effects == -1) ~ 1,
    all(effects == 1) ~ 0.0001,
    TRUE ~ 0.5
  )
}