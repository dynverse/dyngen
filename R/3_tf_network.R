#' Generate a transcription factor network from the backbone
#' 
#' [generate_tf_network()] generates the transcription factors (TFs) that 
#' drive the dynamic process a cell undergoes. 
#' [tf_network_default()] is used to configure parameters pertaining this process.
#' 
#' @param model A dyngen initial model created with [initialise_model()].
#' @param min_tfs_per_module The number of TFs to generate per module in the backbone.
#' @param sample_num_regulators A function to generate the number of TFs per module each TF will be regulated by.
#' @param weighted_sampling When determining what TFs another TF is regulated by, whether to perform weighted
#'  sampling (by rank) or not.
#' 
#' @export
#' @rdname generate_tf_network
#' 
#' @return A dyngen model.
#' 
#' @examples
#' model <- 
#'   initialise_model(
#'     backbone = backbone_bifurcating()
#'   )
#' \dontshow{
#' # actually use a smaller example 
#' # to reduce execution time during
#' # testing of the examples
#' model <- initialise_model(
#'   backbone = model$backbone,
#'   num_cells = 5,
#'   num_targets = 0,
#'   num_hks = 0,
#'   gold_standard_params = gold_standard_default(census_interval = 1, tau = 0.1),
#'   simulation_params = simulation_default(
#'     burn_time = 10,
#'     total_time = 10,
#'     census_interval = 1,
#'     ssa_algorithm = ssa_etl(tau = 0.1),
#'     experiment_params = simulation_type_wild_type(num_simulations = 1)
#'   )
#' )
#' }
#' \donttest{
#' model <- model %>%
#'   generate_tf_network() %>%
#'   generate_feature_network() %>%
#'   generate_kinetics() %>%
#'   generate_gold_standard() %>%
#'   generate_cells() %>%
#'   generate_experiment()
#'   
#' plot_feature_network(model)
#'   
#' dataset <- as_dyno(model)
#' }
generate_tf_network <- function(
  model
) {
  if (model$verbose) cat("Generating TF network\n")
  
  model %>% 
    .generate_tf_info() %>% 
    .generate_tf_network()
}

#' @export
#' @rdname generate_tf_network
tf_network_default <- function(
  min_tfs_per_module = 1L,
  sample_num_regulators = function() 2,
  weighted_sampling = FALSE
) {
  lst(
    min_tfs_per_module,
    sample_num_regulators,
    weighted_sampling
  )
}

.generate_tf_info <- function(model) {
  # satisfy r cmd check
  module_id <- num_tfs <- feature_id <- NULL
  
  module_info <- model$backbone$module_info
  module_network <- model$backbone$module_network
  numbers <- model$numbers
  
  module_info <- 
    module_info %>% mutate(
      num_tfs = .generate_partitions(
        num_elements = max(numbers$num_modules, numbers$num_tfs),
        num_groups = numbers$num_modules, 
        min_elements_per_group = model$tf_network_params$min_tfs_per_module
      ),
      feature_id = map2(module_id, num_tfs, function(module_id, num_tfs) paste0(module_id, "_TF", seq_len(num_tfs))),
      is_tf = TRUE,
      is_hk = FALSE
    )
  
  model$feature_info <- 
    module_info %>% 
    unnest(feature_id) %>% 
    select(feature_id, everything(), -num_tfs)
  
  model
}

.generate_tf_network <- function(model) {
  # satisfy r cmd check
  to <- from <- module_id <- feature_id <- NULL
  
  module_network <- model$backbone$module_network
  tf_info <- model$feature_info
  
  # initialise model structures
  tf_network <- map(tf_info$feature_id, ~list()) %>% set_names(tf_info$feature_id)
  num_targets <- rep(0, length(tf_network)) %>% set_names(tf_info$feature_id)
  
  # go over each tf to find their regulators
  for (i in sample.int(nrow(tf_info))) {
    fi <- tf_info$feature_id[[i]]
    mi <- tf_info$module_id[[i]]
    
    # what are the regulating modules, make sure that each regulating
    # module has at least one regulator
    regulating_modules <- module_network %>% filter(to == mi) %>% pull(from)
    
    # generate regulating tfs
    edges <- 
      map_df(
        regulating_modules,
        function(mreg) {
          # what tfs are in this module
          candidate_regulating_tfs <- tf_info %>% filter(module_id == mreg) %>% pull(feature_id)
          
          # how many regulators will we sample?
          num_regulating_tfs <-
            model$tf_network_params$sample_num_regulators() %>% 
            min(length(candidate_regulating_tfs)) %>%
            max(1L)
          
          weights <- 
            if (model$tf_network_params$weighted_sampling) {
              # do weighted sampling based on the number of
              # targes the candidate regulator already has
              num_targets[candidate_regulating_tfs] + 1
            } else {
              NULL
            }
          regulating_tfs <- sample(candidate_regulating_tfs, num_regulating_tfs, replace = FALSE, prob = weights)
          
          tibble(
            from = regulating_tfs,
            to = fi,
            from_module = mreg,
            to_module = mi
          )
        }
      )
    
    if (nrow(edges) > 0) {
      num_targets[edges$from] <- num_targets[edges$from] + 1
      tf_network[[fi]] <- edges
    }
  }
  
  model$feature_network <- 
    bind_rows(unname(tf_network)) %>% 
    left_join(module_network %>% rename(from_module = from, to_module = to), by = c("from_module", "to_module"))
  
  model
}
