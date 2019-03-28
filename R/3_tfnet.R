#' @export
tfgen_random <- function(
  percentage_tfs = 0.05,
  min_tfs_per_module = 1L
) {
  lst(
    percentage_tfs,
    min_tfs_per_module
  )
}

#' @export
generate_tfnet <- function(
  model
) {
  if (model$verbose) cat("Generating TF network\n")
  
  model %>% 
    .generate_tf_info() %>% 
    .generate_tf_network()
}

.generate_tf_info <- function(model) {
  module_info <- model$modulenet$module_info
  module_network <- model$modulenet$module_network
  
  # number of genes in main based on number of genes in platform
  num_traj_features <- round(model$platform$num_features * model$platform$pct_trajectory_features)
  num_tfs <- round(num_traj_features * model$tfgen_params$percentage_tfs)
  num_targets <- num_traj_features - num_tfs
  num_modules <- nrow(module_info)
  
  model$feature_numbers <- lst(
    num_features = model$platform$num_features,
    num_traj_features,
    num_tfs,
    num_targets,
    num_modules
  )
  
  module_info <- 
    module_info %>% mutate(
      num_tfs = .generate_partitions(
        num_elements = max(num_modules, num_tfs),
        num_groups = num_modules, 
        min_elements_per_group = model$tfgen_params$min_tfs_per_module
      ),
      feature_id = map2(module_id, num_tfs, function(module_id, num_tfs) paste0(module_id, "_TF", seq_len(num_tfs))),
      is_tf = TRUE,
      is_trajectory_driven = TRUE
    )
  
  model$tf_info <- 
    module_info %>% 
    unnest(feature_id) %>% 
    select(feature_id, everything(), -num_tfs)
  
  model
}

.generate_tf_network <- function(model) {
  module_network <- model$modulenet$module_network
  tf_info <- model$tf_info
  
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
    regulating_tfs <- 
      map(
        regulating_modules,
        function(mreg) {
          # what tfs are in this module
          candidate_regulating_tfs <- tf_info %>% filter(module_id == mreg) %>% pull(feature_id)
          
          # how many regulators will we sample?
          num_regulating_tfs <- sample.int(min(3L, length(candidate_regulating_tfs)))
          
          # do weighted sampling based on the number of 
          # targes the candidate regulator already has
          weights <- num_targets[candidate_regulating_tfs] + 1
          
          sample(candidate_regulating_tfs, num_regulating_tfs, replace = FALSE, prob = weights)
        }
      ) %>% 
      unlist()
    
    if (length(regulating_tfs) > 0) {
      num_targets[regulating_tfs] <- num_targets[regulating_tfs] + 1
      tf_network[[fi]] <- tibble(regulator = regulating_tfs, target = fi)
    }
  }
  
  model$tf_network <- bind_rows(tf_network)
  
  model
}