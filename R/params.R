base_params = list(
  model = list(
    # modulenet
    modulenet_name = "linear",
    # treeseed = 1,
    
    # reference
    platform = readRDS(paste0(find.package("dyngen"), "/ext_data/platforms/cell-cycle_leng.rds")),
    
    # network between tfs
    ngenes_per_module_generator = function(ngenes_per_module_mean) {
      function(n) {
        sample(1:(ngenes_per_module_mean * 2), n, replace=TRUE)
      }
    },
    edge_retainment = function(n) sample(seq(1, min(1, n)), 1),
    main_targets_ratio = 0.05,
    # edge_retainment = function(n) 1,
    
    # extra targets
    target_adder_name = "realnet",
    realnet_name = "regulatorycircuits",
    damping = 0.05,
    ntargets_sampler_generator = function(ntargets_mean) {
      function() {sample(1:(ntargets_mean*2), 1)}
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
    local = 8,
    nsimulations = 32,
    ssa_algorithm = fastgssa::ssa.em(noise_strength=8)
  ),
  experiment = list(
    # experiment setting
    sampler = sample_snapshot
    # add_housekeeping = FALSE,
    # n_housekeeping_genes = 500,
    # housekeeping_reference_means = readRDS(paste0(find.package("dyngen"), "/ext_data/housekeeping_reference_means.rds"))[[1]]
  ),
  gs = list(
    max_path_length = 20,
    reference_length = 50,
    smooth_window = 50
  ),
  normalization = list(
    nmads = 3
  )
)# %>% list2env(.GlobalEnv)



simple_params = list(
  model = list(
    # modulenet
    modulenet_name = "linear",
    
    # network between tfs
    ngenes_per_module= function(n) 1, 
    edge_retainment = function(n) 1,
    
    # extra targets
    ntargets_sampler = function() 0,
    
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
    nsimulations = 12
  )
) %>% modifyList(base_params, .)