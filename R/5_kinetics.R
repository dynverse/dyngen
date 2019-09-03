#' Determine the kinetics of the feature network 
#' 
#' [generate_kinetics()] samples the kinetics of genes in the feature network for which 
#'   the kinetics have not yet been defined.
#' [kinetics_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the feature network has been generated with [generate_feature_network()].
#' @param sample_wpr A function specifying the distribution from which to sample the pre-mRNA production rate.
#' @param sample_wsr A function specifying the distribution from which to sample the splicing rate.
#' @param sample_xdr A function specifying the distribution from which to sample the mRNA decay rate.
#' @param sample_ypr A function specifying the distribution from which to sample the protein production rate.
#' @param sample_ydr A function specifying the distribution from which to sample the protein decay rate.
#' @param sample_ind A function specifying the distribution from which to sample the regulator independence factor.
#' @param sample_effect A function specifying the distribution from which to sample the effect of an interaction.
#' @param sample_strength A function specifing the distribution from which to sample the strength of an interaction.
#' @param sample_cooperativity A function specifing the distribution from which to sample the cooperativity of an interaction from.
#' 
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
  
  # determine variables to be used during burn in
  burn_variables <- 
    model$feature_info %>% 
    filter(burn) %>% 
    select(w, x, y) %>% 
    gather(col, val) %>% 
    pull(val)
    
  # return system
  model$simulation_system <- lst(
    reactions = formulae, 
    molecule_ids,
    initial_state,
    parameters,
    burn_variables
  )
  
  model
}

#' @export
#' @rdname generate_kinetics
#' @importFrom stats rnorm runif
kinetics_default <- function(
  sample_wpr = function(n) rnorm(n, 100, 20) %>% pmax(10),
  sample_wsr = function(n) rnorm(n, 10, 2) %>% pmax(2),
  sample_xdr = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_ypr = function(n) rnorm(n, 5, 1) %>% pmax(2),
  sample_ydr = function(n) rnorm(n, 3, .5) %>% pmax(1),
  sample_ind = function(n) runif(n, 0, 1),
  
  sample_effect = function(n) sample(c(-1, 1), n, replace = TRUE, prob = c(.25, .75)),
  sample_strength = function(n) 10 ^ runif(n, log10(1), log10(100)),
  sample_cooperativity = function(n) runif(n, 0.5, 2)
) {
  lst(
    sample_wpr,
    sample_wsr,
    sample_xdr,
    sample_ypr,
    sample_ydr,
    sample_ind,
    
    sample_effect,
    sample_strength,
    sample_cooperativity
  )
}

.kinetics_generate_gene_kinetics <- function(model) {
  if (model$verbose) cat("Generating kinetics for ", nrow(model$feature_info), " features\n", sep = "")
  params <- model$kinetics_params
  
  feature_info <-
    model$feature_info %>% 
    mutate(
      wpr = params$sample_wpr(n()),
      wsr = params$sample_wsr(n()),
      xdr = params$sample_xdr(n()),
      ypr = params$sample_ypr(n()),
      ydr = params$sample_ydr(n()),
      max_protein = wpr / xdr * ypr / ydr,
      ind = ind %|% params$sample_ind(n())
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
      effect = effect %|% params$sample_effect(n()),
      cooperativity = cooperativity %|% params$sample_cooperativity(n()),
      strength = strength %|% params$sample_strength(n()),
      k = max_protein / 2 / strength
    )
  
  # calculate ba and a
  feature_info <- 
    left_join(
      feature_info,
      feature_network %>% 
        rename(feature_id = to) %>% 
        group_by(feature_id) %>% 
        summarise(
          ba_2 = .kinetics_calculate_ba(effect)
        ),
      by = "feature_id"
    ) %>% 
    mutate(
      # 1 for genes that are not being regulated by any other genes,
      # yet did not already have an ba defined
      ba = ba %|% ba_2 %|% 1 
    ) %>% 
    select(-ba_2)
  
  model$feature_info <- feature_info
  model$feature_network <- feature_network
  
  model
}

#' @importFrom GillespieSSA2 reaction
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
  out <- pbapply::pblapply(
    seq_len(nrow(feature_info)),
    cl = model$num_cores,
    function(i) {
      info <- feature_info %>% extract_row_to_list(i)
      
      fid <- info$feature_id
      
      w <- paste0("w_", fid)
      x <- paste0("x_", fid)
      y <- paste0("y_", fid)
      
      wpr <- paste0("wpr_", fid)
      wsr <- paste0("wsr_", fid)
      xdr <- paste0("xdr_", fid)
      ypr <- paste0("ypr_", fid)
      ydr <- paste0("ydr_", fid)
      
      ba <- paste0("ba_", fid)
      ind <- paste0("ind_", fid)
      
      if (!is.null(info$regulators)) {
        rid <- info$regulators$from
        eff <- info$regulators$effect
        reg_ys <- paste0("y_", rid)
        reg_ks <- paste0("k_", rid, "_", fid)
        reg_cs <- paste0("c_", rid, "_", fid)
        regulation_var <- paste0("regulation_", rid, "_", fid)
        
        reg_affinity_calc <- paste(paste0(regulation_var, " = pow(", reg_ys, "/", reg_ks, ", ", reg_cs, "); "), collapse = "")
        
        # Several optimisations have been applied.
        #
        # original:
        #   [ba + x0 + x0x1 + x1] / [x0 + x0x1 + x1 + 1],
        #   with xi = (yi / ki) ^ ci
        #
        # factorise:
        #   [ba + (x0 + 1) * (x1 + 1) - 1] / (x0 + 1) / (x1 + 1) / (x2 + 1)
        #
        # use buffer to remember calculations:
        #   [ba - 1 + buf0 * buf1] / buf0 / buf1 / buf2,
        # with buf0 = x0 + 1, buf1 = x1 + 1, buf2 = x2 + 1
        
        numerator <-
          if (sum(eff > 0) > 0) {
            paste0(ba, " - pow(", ind, ",", sum(eff > 0), ") + ", paste("(", regulation_var[eff > 0], " + ", ind, ")", collapse = " * ", sep = ""))
          } else {
            ba
          }
        denominator <- paste0("1 + ", paste("(", regulation_var, " + 1)", collapse = " * ", sep = ""))
        
        wpr_function <- paste0(reg_affinity_calc, wpr, " * (", numerator, ")/(", denominator, ")")
      } else {
        wpr_function <- paste0(wpr, " * ", ba)
        regulation_var <- character()
      }
      
      formulae <- list(
        # pre-mRNA production
        reaction(
          name = paste0("transcription_", fid),
          effect = set_names(1, w),
          propensity = wpr_function
        ),
        # splicing
        reaction(
          name = paste0("splicing_", fid),
          effect = set_names(c(1, -1), c(x, w)),
          propensity = paste0(wsr, " * ", w)
        ),
        # mRNA degradation
        reaction(
          name = paste0("mrna_degradation_", fid),
          effect = set_names(-1, x),
          propensity = paste0(xdr, " * ", x)
        ),
        # protein production
        reaction(
          name = paste0("translation_", fid), 
          effect = set_names(1, y),
          propensity = paste0(ypr, " * ", x)
        ),
        # protein degradation
        reaction(
          name = paste0("protein_degradation_", fid),
          effect = set_names(-1, y),
          propensity = paste0(ydr, " * ", y)
        )
      )
      
      formulae[[1]]$buffer_ids <- regulation_var
      
      formulae
    }
  )
  
  unlist(out, recursive = FALSE)
}

.kinetics_extract_parameters <- function(model) {
  # extract m to qr, d, p, q, and ba
  feature_params <- 
    model$feature_info %>% 
    select(feature_id, wpr, wsr, xdr, ypr, ydr, ba, ind) %>% 
    gather(param, value, -feature_id) %>% 
    mutate(id = paste0(param, "_", feature_id)) %>% 
    select(id, value) %>% 
    deframe()
  
  # extract k and c
  edge_params <- 
    model$feature_network %>% 
    select(from, to, k, c = cooperativity) %>% 
    gather(param, value, -from, -to) %>% 
    mutate(id = paste0(param, "_", from, "_", to)) %>% 
    select(id, value) %>% 
    deframe()
  
  c(feature_params, edge_params)
}

.kinetics_calculate_ba <- function(effects) {
  case_when(
    all(effects == -1) ~ 1,
    all(effects == 1) ~ 0.0001,
    TRUE ~ 0.5
  )
}