library(tidyverse)
library(dyngen)
library(fastgssa)

backbone_models <- list_backbones()

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
    list2env(list(...), environment())
    
    cat("Params:\n")
    orig_params <- cross_df %>% dynutils::extract_row_to_list(i)
    cat(paste0("  ", names(orig_params), " = ", orig_params, collapse = "\n"), "\n", sep = "")
    
    # generate backbone
    backbone <- backbone_models[[backbone_name]]()
    
    # determine burn time
    burn_time <- 
      case_when(
        backbone_name == "converging" ~ 5,
        backbone_name == "disconnected" ~ 6,
        TRUE ~ 2
      )
    cat("  burn_time=", burn_time, "\n", sep = "")
    
    # determine final time
    exp_pat <- backbone$expression_patterns
    sim_time_sum <- exp_pat %>% filter(start) %>% pull(from) %>% set_names(rep(0, length(.)), .)
    for (i in seq_len(nrow(exp_pat))) {
      sim_time_sum[[exp_pat$to[[i]]]] <- sim_time_sum[[exp_pat$from[[i]]]] + exp_pat$time[[i]]
    }
    total_time <- max(sim_time_sum * 1.2)
    cat("  total_time=", total_time, "\n", sep = "")
    
    # generate model parameters
    initialise_model(
      num_cells = num_cells,
      num_tfs = num_tfs,
      num_targets = num_targets,
      num_hks = num_hks,
      backbone = backbone,
      tf_network_params = tf_network_random(
        min_tfs_per_module = max(floor(num_tfs / 20), 1)
      ),
      feature_network_params = feature_network_default(target_resampling = 5000),
      simulation_params = simulation_default(
        num_simulations = 100, 
        census_interval = .025, 
        burn_time = burn_time, 
        total_time = total_time
      ),
      verbose = TRUE,
      download_cache_dir = "~/.cache/dyngen",
      num_cores = 8
    )
  },
  features = function(...) {
    list2env(list(...), environment())
    
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
    list2env(list(...), environment())
    
    model <- model %>% 
      generate_gold_standard()
    
    g4 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
    g5 <- plot_gold_expression(model)
    ggsave(paste0(output_dir, "/step3_plot_gold.pdf"), g4, width = 8, height = 8)
    ggsave(paste0(output_dir, "/step3_plot_gold_expression.pdf"), g5, width = 12, height = 12)
    
    model
  },
  simulations = function(...) {
    list2env(list(...), environment())
    
    model <- model %>% 
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
    list2env(list(...), environment())
    
    model %>% 
      generate_experiment()
  },
  traj = function(...) {
    list2env(list(...), environment())
    
    traj <- model %>% 
      wrap_dyngen_dataset()
    
    # :scream:
    g1 <- plot_backbone(model)
    g2 <- plot_feature_network(model, show_targets = FALSE)
    g3 <- plot_feature_network(model)
    g4 <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Set3")
    g5 <- plot_gold_expression(model)
    g6 <- plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Set3")
    g8 <- plot_simulations(model)
    g9 <- plot_simulation_expression(model)
    g <- patchwork::wrap_plots(
      g1, g2, g3, g4, g5, g6, g8, g9,
      byrow = FALSE,
      ncol = 4
    ) +
      patchwork::plot_annotation(tag_levels = "A")
    ggsave(paste0(output_dir, "/step6_plot_all.pdf"), g, width = 40, height = 20)
    
    traj
  }
)

cross_df <- 
  crossing(
    size_cats,
    step = seq_along(step_funs),
    rep = 1:3,
    backbone_name = names(backbone_models)
  ) %>% 
  mutate(
    name = paste0(backbone_name, "_rep", rep, "_", size_cat),
    seed = row_number()
  ) %>% 
  filter(size_cat %in% c("pico", "nano", "micro", "milli"))

catt <- function(..., sep = "") {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., sep = sep)
}

walk(seq_len(nrow(cross_df)), function(i) {
  params <- cross_df %>% dynutils::extract_row_to_list(i)
  
  params$i <- i
  
  params$output_dir <- paste0("scripts_tmp/generate_all/", params$name)
  if (!dir.exists(params$output_dir)) dir.create(params$output_dir, recursive = TRUE)
  
  params$step_files <- paste0(params$output_dir, "/step", seq_along(step_funs), "_model_", names(step_funs), ".rds")
  
  if (file.exists(params$step_files[[params$step]])) {
    return()
  }
  
  if (params$step > 1) {
    params$model <- read_rds(params$step_files[[params$step - 1]])
  }
  
  catt("Processing ", params$name, ", step ", params$step, "\n", sep = "")
  sink(paste0(params$output_dir, "/log.txt"), append = params$step != 1)
  on.exit(sink())
  
  catt("Processing ", params$name, ", step ", params$step, "\n", sep = "")
  set.seed(params$seed)
  out <- do.call(step_funs[[params$step]], params)
  write_rds(out, params$step_files[[params$step]])
  catt("Done!\n", sep = "")
})
