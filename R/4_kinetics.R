#' Determine the kinetics of the feature network 
#' 
#' [generate_kinetics()] samples the kinetics of genes in the feature network for which 
#'   the kinetics have not yet been defined.
#' [kinetics_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen intermediary model for which the feature network has been generated with [generate_feature_network()].
#' 
#' @details To write different kinetics settings, you need to write three functions
#' with interface `function(feature_info, feature_network, cache_dir, verbose)`. 
#' Described below are the default kinetics samplers.
#' 
#' `sampler_tfs()` mutates the `feature_info` data frame by adding the following columns:
#'  * `transcription_rate`: the rate at which pre-mRNAs are transcribed, 
#'     in pre-mRNA / hour. Default distribution: U(1, 2).
#'  * `translation_rate`:  the rate at which mRNAs are translated into proteins,
#'     in protein / mRNA / hour. Default distribution: U(100, 150).
#'  * `mrna_halflife`: the half-life of (pre-)mRNA molecules, in hours. 
#'     Default distribution: U(2.5, 5).
#'  * `protein_halflife`: the half-life of proteins, in hours. 
#'     Default distribution: U(5, 10).
#'  * `splicing_rate`: the rate at which pre-mRNAs are spliced into mRNAs, 
#'     in reactions / hour. Default value: log(2) / (10/60), which corresponds to a half-life of 10 minutes.
#'  * `independence`: the degree to which all regulators need to be bound for transcription to occur (0), or 
#'     whether transcription can occur if only one of the regulators is bound (1).
#'     
#' `sampler_nontfs()` samples the `transcription_rate`, `translation_rate`, 
#'   `mrna_halflife` and `protein_halflife` from a supplementary file of Schwannhäusser et al., 
#'   2011, doi.org/10.1038/nature10098. `splicing_rate` is by default the same as in `sampler_tfs()`. 
#'   `independence` is sampled from U(0, 1).
#'   
#' `sampler_interactions()` mutates the `feature_network` data frame by adding the following columns.
#'  * `effect`: the effect of the interaction; upregulating = +1, downregulating = -1.
#'    By default, sampled from {-1, 1} with probabilities {.25, .75}.
#'  * `strength`: the strength of the interaction. Default distribution: 10^U(0, 2).
#'  * `hill`: the hill coefficient. Default distribution: N(2, 2) with a minimum of 1 and a maximum of 10.
#' 
#' @export
#' 
#' @return A dyngen model.
#' 
#' @seealso [dyngen] on how to run a complete dyngen simulation
#' 
#' @examples
#' model <- 
#'   initialise_model(
#'     backbone = backbone_bifurcating(),
#'     kinetics_params = kinetics_default()
#'   )
#' 
#' \donttest{
#' data("example_model")
#' model <- example_model %>%
#'   generate_kinetics()
#' }
generate_kinetics <- function(model) {
  model <- .add_timing(model, "4_kinetics", "checks")
  
  # satisfy r cmd check
  burn <- mol_premrna <- mol_mrna <- mol_protein <- val <- NULL
  
  assert_that(
    !is.null(model$feature_info),
    !is.null(model$feature_network)
  )
  
  # generate kinetics params
  model <- .add_timing(model, "4_kinetics", "generate kinetics")
  model <- .kinetics_generate_gene_kinetics(model)
  
  # generate formulae
  model <- .add_timing(model, "4_kinetics", "generate formulae")
  formulae <- .kinetics_generate_formulae(model)
  
  # create variables
  model <- .add_timing(model, "4_kinetics", "create variables")
  fid <- model$feature_info$feature_id
  model$feature_info$mol_premrna <- paste0("mol_premrna_", fid)
  model$feature_info$mol_mrna <- paste0("mol_mrna_", fid)
  model$feature_info$mol_protein <- paste0("mol_protein_", fid)
  
  molecule_ids <- c(
    model$feature_info$mol_premrna, 
    model$feature_info$mol_mrna,
    model$feature_info$mol_protein
  )

  initial_state <- set_names(
    rep(0, length(molecule_ids)),
    molecule_ids
  )
  
  # extract params
  model <- .add_timing(model, "4_kinetics", "extract parameters")
  parameters <- .kinetics_extract_parameters(
    model$feature_info, 
    model$feature_network
  )
  
  # determine variables to be used during burn in
  burn_variables <- 
    model$feature_info %>% 
    filter(burn) %>% 
    select(mol_premrna, mol_mrna, mol_protein) %>% 
    gather(col, val) %>% 
    pull(val)
    
  # return system
  model <- .add_timing(model, "4_kinetics", "create output")
  model$simulation_system <- lst(
    reactions = formulae, 
    molecule_ids,
    initial_state,
    parameters,
    burn_variables
  )
  
  .add_timing(model, "4_kinetics", "end")
}

#' @export
#' @rdname generate_kinetics
#' @importFrom stats runif
kinetics_default <- function() {
  # satisfy r cmd check
  transcription_rate <- translation_rate <- mrna_halflife <- protein_halflife <- 
    independence <- splicing_rate <- effect <- strength <- hill <- NULL
  
  sampler_tfs <-
    function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
      feature_info %>% mutate(
        transcription_rate = transcription_rate %|% runif(n(), 10, 20),
        translation_rate = translation_rate %|% runif(n(), 100, 150),
        mrna_halflife = mrna_halflife %|% runif(n(), 2.5, 5),
        protein_halflife = protein_halflife %|% runif(n(), 5, 10),
        independence = independence %|% 1,
        splicing_rate = splicing_rate %|% (log(2) / 2)
      )
    }
  
  sampler_interactions <-
    function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
      feature_network %>% 
        mutate(
          effect = effect %|% sample(c(-1L, 1L), n(), replace = TRUE, prob = c(.25, .75)),
          strength = strength %|% 10 ^ runif(n(), log10(1), log10(100)),
          hill = hill %|% rnorm_bounded(n(), 2, 2, min = 1, max = 10)
        )
    }
  
  lst(sampler_tfs, sampler_nontfs = sampler_tfs, sampler_interactions)
}


.kinetics_generate_gene_kinetics <- function(model) {
  # satisfy r cmd check
  is_tf <- mrna_halflife <- protein_halflife <- to <- feature_id <- effect <- 
    basal <- basal_2 <- num_molecules <- mult <- id <- NULL
  
  if (model$verbose) cat("Generating kinetics for ", nrow(model$feature_info), " features\n", sep = "")
  params <- model$kinetics_params
  
  # fetch feature info and network
  feature_info <- model$feature_info %>%
    .kinetics_add_columns(c("transcription_rate", "splicing_rate", "translation_rate", "mrna_halflife", "protein_halflife", "independence"), NA_real_)
  feature_network <- model$feature_network %>%
    .kinetics_add_columns(c("effect", "strength", "hill"), NA_real_)
  
  # generate relatively stable kinetics for TFs
  feature_info_tf <- params$sampler_tfs(
    feature_info %>% filter(is_tf),
    feature_network
  )
  
  # sample kinetics from dataset for non-tfs
  feature_info_nontf <- params$sampler_tfs(
    feature_info %>% filter(!is_tf),
    feature_network, 
    cache_dir = model$download_cache_dir,
    verbose = model$verbose
  )

  # combine feature info
  feature_info <- 
    bind_rows(feature_info_tf, feature_info_nontf) %>% 
    mutate(
      mrna_decay_rate = log(2) / mrna_halflife,
      protein_decay_rate = log(2) / protein_halflife
    )
  
  # sample network kinetics from distributions
  feature_network <- params$sampler_interactions(feature_info, feature_network)
  
  # calculate k
  dis_out <- .kinetics_calculate_dissociation(feature_info, feature_network)
  feature_info <- dis_out$feature_info
  feature_network <- dis_out$feature_network
  
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
      # yet did not already have a value for 'basal' defined
      basal = basal %|% basal_2 %|% 1 
    ) %>% 
    select(-basal_2)
  
  model$feature_info <- feature_info
  model$feature_network <- feature_network
  
  model
}

#' @importFrom GillespieSSA2 reaction
.kinetics_generate_formulae <- function(model) {
  # satisfy r cmd check
  from <- to <- `.` <- NULL
  
  
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
      
      w <- paste0("mol_premrna_", fid)
      x <- paste0("mol_mrna_", fid)
      y <- paste0("mol_protein_", fid)
      
      transcription_rate <- paste0("transcription_rate_", fid)
      splicing_rate <- paste0("splicing_rate_", fid)
      translation_rate <- paste0("translation_rate_", fid)
      mrna_decay_rate <- paste0("mrna_decay_rate_", fid)
      protein_decay_rate <- paste0("protein_decay_rate_", fid)
      
      basal <- paste0("bas_", fid)
      independence <- paste0("ind_", fid)
      
      if (!is.null(info$regulators)) {
        rid <- info$regulators$from
        eff <- info$regulators$effect
        str <- info$regulators$strength
        reg_ys <- paste0("mol_protein_", rid)
        reg_diss <- paste0("dis_", rid, "_", fid)
        reg_hills <- paste0("hill_", rid, "_", fid)
        reg_strs <- paste0("str_", rid, "_", fid)
        regulation_var <- paste0("chi_", rid, "_", fid)
        
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
        
        act_function <- paste0(reg_affinity_calc, transcription_rate, " * (", numerator, ")/(", denominator, ")")
      } else {
        act_function <- paste0(transcription_rate, " * ", basal)
        regulation_var <- character()
      }
      
      formulae <- list(
        # pre-mRNA production
        reaction(
          name = paste0("transcription_", fid),
          effect = set_names(1, w),
          propensity = paste0(act_function)
        ),
        # splicing
        reaction(
          name = paste0("splicing_", fid),
          effect = set_names(c(1, -1), c(x, w)),
          propensity = paste0(splicing_rate, " * ", w)
        ),
        # protein production
        reaction(
          name = paste0("translation_", fid), 
          effect = set_names(1, y),
          propensity = paste0(translation_rate, " * ", x)
        ),
        # pre-mRNA degradation
        reaction(
          name = paste0("premrna_degradation_", fid),
          effect = set_names(-1, w),
          propensity = paste0(mrna_decay_rate, " * ", w)
        ),
        # mRNA degradation
        reaction(
          name = paste0("mrna_degradation_", fid),
          effect = set_names(-1, x),
          propensity = paste0(mrna_decay_rate, " * ", x)
        ),
        # protein degradation
        reaction(
          name = paste0("protein_degradation_", fid),
          effect = set_names(-1, y),
          propensity = paste0(protein_decay_rate, " * ", y)
        )
      )
      
      formulae[[1]]$buffer_ids <- regulation_var
      
      formulae
    }
  )
  
  unlist(out, recursive = FALSE)
}

.kinetics_extract_parameters_as_df <- function(feature_info, feature_network) {
  # satisfy r cmd check
  feature_id <- transcription_rate <- splicing_rate <- translation_rate <- mrna_decay_rate <- protein_decay_rate <-
    basal <- independence <- param <- value <- id <- from <- to <- dissociation <- hill <- strength <- `.` <- from <- NULL
  
  # extract production / degradation rates, ind and bas
  feature_params <- 
    feature_info %>% 
    select(feature_id, transcription_rate, splicing_rate, translation_rate, mrna_decay_rate, protein_decay_rate, bas = basal, ind = independence) %>% 
    gather(param, value, -feature_id) %>% 
    mutate(id = paste0(param, "_", feature_id), type = "feature_info")
  
  # extract dis, hill, str
  edge_params <- 
    feature_network %>% 
    select(from, to, dis = dissociation, hill, str = strength) %>% 
    gather(param, value, -from, -to) %>% 
    mutate(id = paste0(param, "_", from, "_", to), type = "feature_network")
  
  bind_rows(feature_params, edge_params)
}


.kinetics_extract_parameters <- function(feature_info, feature_network) {
  .kinetics_extract_parameters_as_df(feature_info, feature_network) %>% 
    select(.data$id, .data$value) %>% 
    deframe()
}

.kinetics_calculate_basal <- function(effects) {
  case_when(
    all(effects == -1) ~ 1,
    all(effects == 1) ~ 0.0001,
    TRUE ~ 0.5
  )
}

.kinetics_calculate_dissociation <- function(feature_info, feature_network) {
  # satisfy r cmd check
  transcription_rate <- mrna_decay_rate <- splicing_rate <- max_premrna <- translation_rate <- 
    protein_decay_rate <- max_mrna <- feature_id <- max_protein <- NULL
  
  remove <- c("max_premrna", "max_mrna", "max_protein", "dissociation", "k", "max_protein")
  
  feature_info <- feature_info[, !colnames(feature_info) %in% remove]
  feature_network <- feature_network[, !colnames(feature_network) %in% remove]
  
  feature_info <- 
    feature_info %>%
    mutate(
      max_premrna = transcription_rate / (mrna_decay_rate + splicing_rate),
      max_mrna = splicing_rate / mrna_decay_rate * max_premrna,
      max_protein = translation_rate / protein_decay_rate * max_mrna
    )
  
  feature_network <- 
    feature_network %>% 
    left_join(feature_info %>% select(from = feature_id, max_protein), by = "from") %>% 
    mutate(
      dissociation = max_protein / 2
    )
  
  lst(feature_info, feature_network)
}

.kinetics_add_columns <- function(df, colnames, fill = NA) {
  for (colname in colnames) {
    if (!colname %in% colnames(df)) {
      df[[colname]] <- fill
    }
  }
  df
}