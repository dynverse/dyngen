library(tidyverse)
library(dynutils)
#library(dyngen)

updates <- tribble(
  ~modulenet_name, ~totaltime,
  "linear", 5,
  "bifurcating", 5,
  "linear_long", 30,
  "cycle", 30,
  "consecutive_bifurcating", 10,
  "bifurcating_converging", 15,
  "trifurcating", 10,
  "converging", 10,
  "bifurcating_loop", 30
)
nreplicates <- 2
updates <- updates[rep(1:nrow(updates), each=nreplicates), ] %>% mutate(replicate_id = rep(1:nreplicates, nrow(updates)))

settings <- pmap(
  updates,
  function(modulenet_name, totaltime, replicate_id) {
    modifyList(dyngen:::simple_params, list(model = list(modulenet_name = modulenet_name), simulation=list(totaltime = totaltime)))
  }
)

settings <- settings[15]

ncores <- 8

params <- settings[[5]]

folder <- "~/thesis/projects/dynverse/dynalysis/analysis/data/datasets/synthetic/v5/"
# unlink(folder)
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

models <- map(seq_along(settings), function(settings_i) {
  print(glue::glue("{settings_i} / {length(settings)} ======================================"))
  params <- settings[[settings_i]]
  options(ncores = ncores)
  
  model <- invoke(dyngen:::generate_model_from_modulenet, params$model)
  model
})
saveRDS(models, paste0(folder, "models.rds"))

simulations <- map(seq_along(settings), function(settings_i) {
  print(glue::glue("{settings_i} / {length(settings)} ======================================"))
  params <- settings[[settings_i]]
  model <- models[[settings_i]]
  options(ncores = ncores)

  simulation <- invoke(dyngen:::simulate_multiple, params$simulation, model$system)
  simulation
})
saveRDS(simulations, paste0(folder, "simulations.rds"))
simulations <- readRDS(paste0(folder, "simulations.rds"))

goldstandards <- map(seq_along(settings), function(settings_i) {
  print(glue::glue("{settings_i} / {length(settings)} ======================================"))
  params <- settings[[settings_i]]
  model <- models[[settings_i]]
  simulation <- simulations[[settings_i]]
  options(ncores = ncores)
  
  print("Preprocessing")
  simulation <- dyngen:::preprocess_simulation_for_gs(simulation, model, params$gs$smooth_window) # do preprocessing separate, otherwise zoo will stay in an infinite loop in case of later error
  gs <- invoke(dyngen:::extract_goldstandard, params$gs, simulation, model, preprocess=FALSE)
  gs$checks <- dyngen:::check_goldstandard(gs)
  gs
})
saveRDS(goldstandards, paste0(folder, "goldstandards.rds"))

settings_i <- 1
walk(seq_along(settings), function(settings_i) {
  print(glue::glue("{settings_i} / {length(settings)} ======================================"))
  params <- settings[[settings_i]]
  model <- models[[settings_i]]
  simulation <- simulations[[settings_i]]
  gs <- goldstandards[[settings_i]]
  
  pdf(paste0(folder, "/", settings_i, ".pdf"), width=12, height=12)
  dyngen:::plot_net(model, label=FALSE, main_only = FALSE)
  dyngen:::plot_modulenet(model)
  dyngen:::plot_goldstandard(simulation, model, gs)
  dev.off()
})



