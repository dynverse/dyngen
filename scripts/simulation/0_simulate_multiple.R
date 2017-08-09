library(tidyverse)

.version <- 4
.datasets_location <- paste0("/home/wouters/thesis/projects/dyngen/results/", .version, "/")

## Model settings----
nreplicates <- 4
modelgenerators <- list(
  list(modulenetname = "linear", totaltime=4),
  list(modulenetname = "linear_long", totaltime=10),
  list(modulenetname = "consecutive_bifurcating", totaltime=6),
  list(modulenetname = "bifurcating_convergence", totaltime=8),
  list(modulenetname = "trifurcating", totaltime=8),
  list(modulenetname = "bifurcating_cycle", totaltime=20),
  list(modulenetname = "bifurcating_loop", totaltime=20)
)

takesettings <- list(
  list(type="snapshot", ncells=500),
  list(type="synchronized", ntimepoints=10)
)

expand_lists <- function(...) {
  dots <- list(...)
  map(dots, ~seq_len(length(.))) %>% 
    expand.grid() %>% 
    set_names(names(dots)) %>% 
    {split(., seq_len(nrow(.)))} %>% 
    map(function(x){
      map2(seq_len(length(x)), x, ~dots[[.x]][[.y]]) %>% set_names(names(x))
    }) %>% 
    setNames(NULL)
}

## Generate models-------------
modelsettings <- expand_lists(modelgenerator=modelgenerators, replicate=seq_len(nreplicates))
models <- map(modelsettings, function(modelsetting) {
  generate_model(modulenetname = modelsetting$modelgenerator$modulenetname)
})
models <- map2(models, modelsettings, function(model, modelsetting) {
  model$info <- list(
    id = paste0(modelsetting$modelgenerator$modulenetname, "_", modelsetting$replicate),
    version = .version
  )
  class(model) <- "dyngen::model"
  model$modelsetting <- modelsetting
  model
})
saver(models, "models")

models <- loader(overviewer("models")$id, "models")

## Simulate models----------------
simulations <- parallel::mclapply(models, function(model) {
  simulations <-run_simulations(model, model$modelsetting$modelgenerator$totaltime, 2, 32, local=FALSE)
}, mc.cores=8)
simulations <- map(simulations, function(x) {x$simulations <- smoothe_simulations(x$simulations, x$model); x})
simulations <- map2(simulations, models, function(simulation, model) {
  simulation$info <- list(
    id = paste0(model$info$id),
    model_id = model$info$id,
    version = .version
  )
  simulation
})
saver(simulations, "simulations")

## Extract gold standards-------------------
# only retain smoothed module data, as this is the only data needed for gold standards
# this makes it easier to transfer data to the cluster
sub = simulations %>% map(function(simulations) {
  simulations$simulations = map(simulations$simulations, ~list(expression_modules = .$expression_modules))
  simulations
})

goldstandards <- pbapply::pblapply(sub, function(x){
  print("----")
  library(magrittr)
  library(tidyverse)
  library(dambiutils)
  dyngen:::extract_goldstandard(x, smooth=FALSE, verbose = TRUE)
})#, mc.cores=8)#, qsub_environment = list2env(list()))
goldstandards <- map2(goldstandards, simulations, function(gs, simulation) {
  gs$info <- list(
    id=simulation$info$id, 
    simulation_id = simulation$info$id, 
    model_id = simulation$info$model_id,
    version = .version
  )
  gs
})

saver(goldstandards, "goldstandards")

## Extract experiments-------------------
experimentsettings <- expand.grid(simulation=seq_along(simulations), takesetting=seq_along(takesettings))
experiments <- map2(experimentsettings$simulation, experimentsettings$takesetting, function(simulation_id, takesetting_id) {
  experiment <- run_experiment(simulations[[simulation_id]], takesettings[[takesetting_id]])
  experiment$info$model_id <- models[[simulation_id]]$info$id
  experiment$info$goldstandard_id <- goldstandards[[simulation_id]]$info$id
  experiment$info$simulation_id <- simulations[[simulation_id]]$info$id
  experiment$info$id <- paste0(models[[simulation_id]]$info$id, "_", takesettings[[takesetting_id]]$type)
  experiment$info$version <- .version
  experiment
})
saver(experiments, "experiments")

# plot gold standards-------------------
walk(experiments, function(experiment) {
  print(experiment$info$id)
  gs <- goldstandards %>% keep(~(.$info$id == experiment$info$goldstandard_id)) %>% first
  
  plot = plot_goldstandard(experiment=experiment, gs) %>% unlist(recursive=FALSE) %>% 
    cowplot::plot_grid(plotlist = ., ncol=2) %>% cowplot::add_sub(experiment$info$id)
  
  dir.create(dirname(paste0("images/", experiment$info$id, ".png")), recursive = TRUE)
  ggsave(paste0("images/", experiment$info$id, ".png"), plot, width = 10, height=10)
})


## Run scRNAseq----------------------
platforms = readr::read_tsv("data/platforms.tsv")
# platforms2 = expand.grid(sequencerate = c(seq_len(10)/10), cellcapturerate = seq_len(10)/10) %>% as.data.frame()
# platforms = platforms[1, ] %>% select(-sequencerate, -cellcapturerate) %>% data.frame(platforms2)
# platforms$platformname = seq_len(nrow(platforms)) %>% as.character
datasets = lapply(experiments, function(experiment) {
  lapply(seq_len(nrow(platforms)), function(platform_id) {
    platform = platforms[platform_id, ] %>% as.list
    dataset = dyngen:::run_scrnaseq(experiment, platform)
    dataset$info = experiment$info
    dataset$info$experiment_id <- experiment$info$id
    dataset$info$id <- paste0(experiment$info$id, "_", platform$platformname)
    dataset$info$version <- .version
    dataset
  })
}) %>% unlist(recursive=F)
saver(datasets, "datasets")

## Sync to prism
PRISM:::rsync_remote("", "~/thesis/projects/dyngen/results", "prism", "/group/irc/shared/dyngen_results")

## Download from prism
PRISM:::rsync_remote("prism", "/group/irc/shared/dyngen_results/results", "", "~/Workspace/papers/ti_eval/dyngen")
