library(tidyverse)
library(dyngen)
library(fastgssa)

base_model <- function(num_cells, num_tfs, num_targets, num_hks) {
  initialise_model(
    num_cells = num_cells,
    num_tfs = num_tfs,
    num_targets = num_targets,
    num_hks = num_hks,
    backbone = backbone_linear(),
    tf_network_params = tf_network_random(min_tfs_per_module = max(floor(num_tfs / 20), 1)),
    feature_network_params = feature_network_default(target_resampling = 5000),
    simulation_params = simulation_default(num_simulations = 100, census_interval = .025),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8
  )
}

backbone_modifiers <- list(
  bifurcating = function(model) {
    model$backbone <- backbone_bifurcating()
    model
  },
  bifurcating_converging = function(model) {
    model$backbone <- backbone_bifurcating_converging()
    model
  },
  bifurcating_cycle = function(model) {
    model$backbone <- backbone_bifurcating_cycle()
    model$simulation_params$total_time <- 20
    model
  },
  bifurcating_loop = function(model) {
    model$backbone <- backbone_bifurcating_loop()
    model
  },
  binary_tree = function(model) {
    model$backbone <- backbone_binary_tree()
    model
  },
  branching = function(model) {
    model$backbone <- backbone_branching()
    model
  },
  consecutive_bifurcating = function(model) {
    model$backbone <- backbone_consecutive_bifurcating()
    model
  },
  converging = function(model) {
    model$backbone <- backbone_converging()
    model$simulation_params$burn_time <- 5
    model
  },
  cycle = function(model) {
    model$backbone <- backbone_cycle()
    model$simulation_params$total_time <- 20
    model
  },
  disconnected = function(model) {
    model$backbone <- backbone_disconnected()
    model$simulation_params$burn_time <- 6
    model
  },
  linear = function(model) {
    model$backbone <- backbone_linear()
    model
  },
  trifurcating = function(model) {
    model$backbone <- backbone_trifurcating()
    model
  }
)

size_cats <- tribble(
  ~size_cat, ~num_cells, ~num_tfs, ~num_targets, ~num_hks,
  "pico", 100, 20, 40, 40,
  "nano", 200, 30, 85, 85,
  "micro", 500, 40, 230, 230,
  "milli", 1000, 50, 475, 475,
  "one", 5000, 100, 1450, 1450,
  "kilo", 10000, 200, 4900, 4900,
  "mega", 20000, 200, 4900, 4900,
  "giga", 50000, 200, 4900, 4900
)

step_funs <- list(
  base = function(...) {
    model <- 
      base_model(num_cells = num_cells, num_tfs = num_tfs, num_targets = num_targets, num_hks = num_hks) %>% 
      backbone_modifiers[[backbone_name]]()
    
    # fixing num_modules value for new backbone
    model$numbers$num_modules <- nrow(model$backbone$module_info)
    
    cat("Params:\n")
    cat(paste0("  ", names(params), " = ", params, collapse = "\n"), "\n", sep = "")
    
    model
  },
  features = function(...) {
    model <- model %>% 
      generate_tf_network() %>% 
      generate_feature_network() %>% 
      generate_kinetics()
    
    g1 <- plot_backbone(model)
    g2 <- plot_feature_network(model, show_targets = FALSE)
    g3 <- plot_feature_network(model)
    g <- patchwork::wrap_plots(
      patchwork::wrap_plots(g1, g2, ncol = 1),
      g3, nrow = 1
    )  
    ggsave(paste0(output_dir, "/step2_plot_networks.pdf"), g, width = 20, height = 12)
    
    model
  },
  gs = function(...) {
    model <- read_rds(model2_file) %>% 
      generate_gold_standard()
    
    g4 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
    g5 <- plot_gold_expression(model)
    ggsave(paste0(output_dir, "/step3_plot_gold.pdf"), g4, width = 8, height = 8)
    ggsave(paste0(output_dir, "/step3_plot_gold_expression.pdf"), g5, width = 12, height = 12)
    
    model
  },
  simulations = function(...) {
    model <- read_rds(model3_file) %>% 
      generate_cells()
    
    g6 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
    g7 <- plot_gold_mappings(model) + scale_colour_brewer(palette = "Set3")
    g8 <- plot_simulations(model)
    g9 <- plot_simulation_expression(model)
    g <- patchwork::wrap_plots(g6, g8, nrow = 1)
    ggsave(paste0(output_dir, "/step4_plot_simulations.pdf"), g, width = 12, height = 6)
    ggsave(paste0(output_dir, "/step4_plot_gold_mappings.pdf"), g7, width = 12, height = 10)
    
    pdf(paste0(output_dir, "/step4_simulation_expression.pdf"), width = 12, height = 12)
    tryCatch({
      for (sim_i in seq_len(model$simulation_params$num_simulations)) {
        g <- 
          plot_simulation_expression(model, simulation_i = sim_i) +
          labs(title = paste0("Simulation ", sim_i))
        print(g)
      }
    }, finally = {
      dev.off()
    })
    
    model
  },
  experiment = function(...) {
    model <- read_rds(model4_file) %>% 
      generate_experiment()
    model
  },
  traj = function(...) {
    model <- read_rds(model5_file)
    traj <- model %>% 
      wrap_dyngen_dataset()
    
    # :scream:
    g1 <- plot_backbone(model)
    g2 <- plot_feature_network(model, show_targets = FALSE)
    g3 <- plot_feature_network(model)
    g4 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
    g5 <- plot_gold_expression(model)
    g6 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
    g8 <- plot_simulations(model)
    g9 <- plot_simulation_expression(model)
    g <- patchwork::wrap_plots(
      g1, g2, g3, g4, g5, g6, g8, g9,
      byrow = FALSE,
      ncol = 4
    ) +
      patchwork::plot_annotation(tag_levels = "A")
    ggsave(paste0(output_dir, "/step6_plot_all.pdf"), g, width = 40, height = 30)
    
    traj
  }
)

cross_df <- 
  crossing(
    size_cats,
    step = seq_along(step_funs),
    rep = 1:3,
    backbone_name = names(backbone_modifiers)
  ) %>% 
  mutate(
    name = paste0(backbone_name, "_rep", rep, "_", size_cat),
    seed = row_number()
  ) %>% 
  filter(size_cat == "pico", step <= 4)

catt <- function(..., sep = "") {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., sep = sep)
}

walk(seq_len(nrow(cross_df)), function(i) {
  params <- cross_df %>% dynutils::extract_row_to_list(i)
  
  params$output_dir <- paste0("scripts_tmp/generate_all/", params$name)
  if (!dir.exists(params$output_dir)) dir.create(params$output_dir, recursive = TRUE)
  
  params$step_files <- paste0(params$output_dir, "/step", seq_along(step_funs), "_model_", names(step_funs), ".rds")
  
  if (file.exists(params$step_files[[step]])) {
    return()
  }
  
  if (params$step > 1) {
    params$model <- read_rds(params$step_files[[params$step - 1]])
  }
  
  catt("Processing ", params$name, ", step ", params$step, "\n", sep = "")
  sink(paste0(params$output_dir, "/log.txt"), append = params$step != 1)
  on.exit(sink())
  
  set.seed(params$seed)
  out <- do.call(step_funs[[params$step]], params)
  write_rds(out, params$step_files[[params$step]])
})
