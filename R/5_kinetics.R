#' Determine the kinetics of the feature network 
#' 
#' [generate_kinetics()] samples the kinetics of genes in the feature network for which 
#'   the kinetics have not yet been defined.
#' [kinetics_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the feature network has been generated with [generate_feature_network()].
#' @param sample_wpr A function specifying the distribution from which to sample the pre-mRNA production rate.
#' @param sample_whl A function specifying the distribution from which to sample the pre-mRNA half-life.
#' @param sample_wsr A function specifying the distribution from which to sample the splicing rate.
#' @param sample_xhl A function specifying the distribution from which to sample the pre-mRNA half-life.
#' @param sample_ypr A function specifying the distribution from which to sample the protein production rate.
#' @param sample_yhl A function specifying the distribution from which to sample the pre-mRNA half-life.
#' @param sample_independence A function specifying the distribution from which to sample the regulator independence factor.
#' @param sample_effect A function specifying the distribution from which to sample the effect of an interaction.
#' @param sample_strength A function specifying the distribution from which to sample the strength of an interaction.
#' @param sample_hill A function specifying the distribution from which to sample the hill coefficient of an interaction from.
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
  parameters <- .kinetics_extract_parameters(
    model$feature_info, 
    model$feature_network
  )
  
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

# df <- bind_rows(
#   tibble(lab = "before", x = log(2) / rnorm_bounded(100000, 3, .5, min = 1)),
#   tibble(lab = "after", x = rnorm_bounded(100000, .24, .043, min = .13)),
# )
# 
# df <- bind_rows(
#   tibble(lab = "before", x = log(2) / rnorm_bounded(100000, 5, 1, min = 2)),
#   tibble(lab = "after", x = rnorm_bounded(100000, .145, .032, min = .075)),
# )
# df %>% group_by(lab) %>% summarise(mean = mean(x), sd = sd(x), min = min(x))
# ggplot(df) + geom_density(aes(x, colour = lab))

#' @export
#' @rdname generate_kinetics
#' @importFrom stats rnorm runif
kinetics_default <- function(
  sample_wpr = function(n) rnorm_bounded(n, 50, 10, min = 10),
  sample_whl = function(n) rnorm_bounded(n, .15, .03, min = .05),
  sample_wsr = function(n) rnorm_bounded(n, 5, 1, min = 1),
  sample_xhl = function(n) rnorm_bounded(n, .15, .03, min = .05),
  sample_ypr = function(n) rnorm_bounded(n, 5, 1, min = 1),
  sample_yhl = function(n) rnorm_bounded(n, .25, .05, min = .1),
  sample_independence = function(n) runif(n, 0, 1),
  
  sample_effect = function(n) sample(c(-1, 1), n, replace = TRUE, prob = c(.25, .75)),
  sample_strength = function(n) 10 ^ runif(n, log10(1), log10(100)),
  sample_hill = function(n) rnorm_bounded(n, 2, 2, min = 1, max = 10)
) {
  lst(
    sample_wpr,
    sample_whl,
    sample_wsr,
    sample_xhl,
    sample_ypr,
    sample_yhl,
    sample_independence,
    
    sample_effect,
    sample_strength,
    sample_hill
  )
}


.kinetics_generate_gene_kinetics <- function(model) {
  if (model$verbose) cat("Generating kinetics for ", nrow(model$feature_info), " features\n", sep = "")
  params <- model$kinetics_params
  
  feature_info <-
    model$feature_info %>% 
    mutate(
      wpr = params$sample_wpr(n()),
      whl = params$sample_whl(n()),
      wdr = log(2) / whl,
      wsr = params$sample_wsr(n()),
      xhl = params$sample_xhl(n()),
      xdr = log(2) / xhl,
      ypr = params$sample_ypr(n()),
      yhl = params$sample_yhl(n()),
      ydr = log(2) / yhl,
      independence = independence %|% params$sample_independence(n())
    )
  
  # sample effect, hill and strength
  # calculate k
  feature_network <- 
    model$feature_network %>% 
    mutate(
      effect = effect %|% params$sample_effect(n()),
      hill = hill %|% params$sample_hill(n()),
      strength = strength %|% params$sample_strength(n())
    )
  out <- .kinetics_calculate_dissociation(feature_info, feature_network)
  feature_info <- out$feature_info
  feature_network <- out$feature_network
  
  # calculate ba and a
  feature_info <- 
    left_join(
      feature_info,
      feature_network %>% 
        rename(feature_id = to) %>% 
        group_by(feature_id) %>% 
        summarise(
          basal_2 = .kinetics_calculate_basal(effect)
        ),
      by = "feature_id"
    ) %>% 
    mutate(
      # 1 for genes that are not being regulated by any other genes,
      # yet did not already have an ba defined
      basal = basal %|% basal_2 %|% 1 
    ) %>% 
    select(-basal_2)
  
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
      wdr <- paste0("wdr_", fid)
      wsr <- paste0("wsr_", fid)
      xdr <- paste0("xdr_", fid)
      ypr <- paste0("ypr_", fid)
      ydr <- paste0("ydr_", fid)
      
      basal <- paste0("bas_", fid)
      independence <- paste0("ind_", fid)
      
      if (!is.null(info$regulators)) {
        rid <- info$regulators$from
        eff <- info$regulators$effect
        str <- info$regulators$strength
        reg_ys <- paste0("y_", rid)
        reg_diss <- paste0("dis_", rid, "_", fid)
        reg_hills <- paste0("hill_", rid, "_", fid)
        regulation_var <- paste0("chi_", rid, "_", fid)
        reg_strs <- paste0("str_", rid, "_", fid)
        
        reg_affinity_calc <- paste(paste0(regulation_var, " = ", reg_strs, " * pow(", reg_ys, "/", reg_diss, ", ", reg_hills, "); "), collapse = "")
        
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
            paste0(basal, " - pow(", independence, ",", sum(eff > 0), ") + ", paste("(", regulation_var[eff > 0], " + ", independence, ")", collapse = " * ", sep = ""))
          } else {
            basal
          }
        denominator <- paste("(", regulation_var, " + 1)", collapse = " * ", sep = "")
        
        wpr_function <- paste0(reg_affinity_calc, wpr, " * (", numerator, ")/(", denominator, ")")
      } else {
        wpr_function <- paste0(wpr, " * ", basal)
        regulation_var <- character()
      }
      
      formulae <- list(
        # pre-mRNA production
        reaction(
          name = paste0("transcription_", fid),
          effect = set_names(1, w),
          propensity = wpr_function
        ),
        # pre-mRNA degradation
        reaction(
          name = paste0("premrna_degradation_", fid),
          effect = set_names(-1, w),
          propensity = paste0(wdr, " * ", w)
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

.kinetics_extract_parameters <- function(feature_info, feature_network) {
  # extract production / degradation rates, ind and bas
  feature_params <- 
    feature_info %>% 
    select(feature_id, wpr, wdr, wsr, xdr, ypr, ydr, bas = basal, ind = independence) %>% 
    gather(param, value, -feature_id) %>% 
    mutate(id = paste0(param, "_", feature_id)) %>% 
    select(id, value) %>% 
    deframe()
  
  # extract dis, hill, str
  edge_params <- 
    feature_network %>% 
    select(from, to, dis = dissociation, hill, str = strength) %>% 
    gather(param, value, -from, -to) %>% 
    mutate(id = paste0(param, "_", from, "_", to)) %>% 
    select(id, value) %>% 
    deframe()
  
  c(feature_params, edge_params)
}

.kinetics_calculate_basal <- function(effects) {
  case_when(
    all(effects == -1) ~ 1,
    all(effects == 1) ~ 0.0001,
    TRUE ~ 0.5
  )
}

.kinetics_calculate_dissociation <- function(feature_info, feature_network) {
  remove <- c("max_w", "max_x", "max_y", "dissociation", "k", "max_protein")
  
  feature_info <- feature_info[, !colnames(feature_info) %in% remove]
  feature_network <- feature_network[, !colnames(feature_network) %in% remove]
  
  feature_info <- 
    feature_info %>%
    mutate(
      max_w = wpr / (wdr + wsr),
      max_x = wsr / xdr * max_w,
      max_y = ypr / ydr * max_x
    )
  
  feature_network <- 
    feature_network %>% 
    left_join(feature_info %>% select(from = feature_id, max_y), by = "from") %>% 
    mutate(
      dissociation = max_y / 2
    )
  
  lst(feature_info, feature_network)
}