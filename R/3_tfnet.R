#' @export
tf_random <- function(
  percentage_tfs = 0.05,
  min_tfs_per_module = 1L,
  edge_retainment = function(n) round(n / 2) %>% max(1)
) {
  lst(
    percentage_tfs,
    min_tfs_per_module,
    edge_retainment
  )
}

#' @export
generate_tfnet <- function(
  data
) {
  module_info <- data$modulenet$module_info
  module_network <- data$modulenet$module_network
  
  ###########################
  ## GENERATING TF NETWORK ##
  ###########################
  if (data$verbose) cat("Generating TF network\n")
  
  # number of genes in main based on number of genes in platform
  num_traj_features <- round(data$platform$num_features * data$platform$pct_trajectory_features)
  num_tfs <- round(num_traj_features * data$tfgen_params$percentage_tfs)
  num_targets <- num_traj_features - num_tfs
  
  num_modules <- nrow(module_info)
  
  config$tf_info <- 
    module_info %>% mutate(
      num_tfs = .num_tf_per_module_sampler(
        num_tfs = max(num_modules, num_tfs),
        num_modules = num_modules, 
        min_tfs_per_module = config$tfgen_params$min_tfs_per_module
      ),
      feature_id = map2(module_id, num_tfs, function(module_id, num_tfs) paste0(module_id, "_TF", seq_len(num_tfs))),
      is_tf = TRUE,
      is_tractory_driven = TRUE
    ) %>% 
    unnest(feature_id) %>% 
    select(feature_id, everything(), -num_tfs)
  
  data$tf_regnet <- .generate_tf_regnet(module_network, data$tf_info)
  
  data  
}

.num_tf_per_module_sampler <- function(num_tfs, num_modules, min_tfs_per_module) {
  assert_that(
    min_tfs_per_module >= 1, 
    num_modules * min_tfs_per_module <= num_tfs
  )
  
  sample(
    seq(0, num_tfs - num_modules * min_tfs_per_module),
    num_modules - 1,
    replace = TRUE
  ) %>% 
    sort() %>% 
    c(0, ., num_tfs - num_modules) %>% 
    as.integer() %>% 
    diff() %>% 
    { . + min_tfs_per_module}
}

.generate_tf_regnet <- function(module_network, tf_info) {
  # initialise data structures
  regnet <- map(tf_info$feature_id, ~c()) %>% set_names(tf_info$feature_id)
  num_targets <- rep(0, length(regnet)) %>% set_names(tf_info$feature_id)
  
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
    
    num_targets[regulating_tfs] <- num_targets[regulating_tfs] + 1
    regnet[[fi]] <- tibble(regulator = regulating_tfs, target = fi)
  }
  
  bind_rows(regnet)
}