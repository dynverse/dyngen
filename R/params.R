#' Base params
#' 
#' @export
base_params = list(
  model = list(
    # modulenet
    modulenet_name = "linear",
    # treeseed = 1,
    
    # reference
    platform = readRDS(paste0(find.package("dyngen"), "/ext_data/platforms/small.rds")),
    
    # network between tfs
    ngenes_per_module_sampler = function(n_genes, n_modules) {
      diff(c(0, sort(sample(1:(n_genes-1), n_modules-1, replace = FALSE)), n_genes)) # randomly distribute the tfs over all modules
    },
    edge_retainment = function(n) sample(seq(1, min(1, n)), 1),
    main_targets_ratio = 0.05,
    # edge_retainment = function(n) 1,
    
    # extra targets
    target_adder_name = "realnet",
    realnet_name = "regulatorycircuits",
    damping = 0.05,
    ntargets_sampler = function(n_genes, n_regulators) {
      diff(c(0, sort(sample(1:(n_genes-1), n_regulators-1, replace = FALSE)), n_genes)) # randomly distribute the targets over all modules
    },
    
    #system
    samplers = list(
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
  ),
  simulation = list(
    totaltime = 20,
    burntime = 2,
    nsimulations = 32,
    ssa_algorithm = fastgssa::ssa.em(noise_strength=8),
    local = TRUE
  ),
  experiment = list(
    sampler = sample_snapshot,
    platform = readRDS(paste0(find.package("dyngen"), "/ext_data/platforms/small.rds"))
  ),
  gs = list(
    max_path_length = 20,
    reference_length = 50,
    smooth_window = 50
  ),
  normalisation = list(
    nmads = 3
  )
)# %>% list2env(.GlobalEnv)

#' @rdname base_params
#' @export
simple_params = list(
  model = list(
    # modulenet
    modulenet_name = "linear",
    
    #system
    samplers = list(
      sample_r = function(n) 20, 
      sample_d = function(n) 10, 
      sample_p = function(n) 10,
      sample_q = function(n) 1,
      calculate_a0 = function(effects) {
        if(length(effects) == 0) {
          1
        } else if(sum(effects == -1) == 0) {
          0
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
      sample_strength = function(n) 1,
      calculate_k = function(max_protein, strength) {
        max_protein/2/strength
      },
      sample_cooperativity = function(n) 2
    )
  ),
  simulation = list(
    nsimulations = 6
  )
) %>% modifyList(base_params, .)