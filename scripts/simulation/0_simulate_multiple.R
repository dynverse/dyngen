library(dplyr)

.datasets_location = "/home/wouters/thesis/projects/dyngen/results/"
.version = "4/"

## Model settings----
nreplicates <- 4
modelgenerators = list(
  list(modulenetname = "linear", totaltime=4),
  list(modulenetname = "linear_long", totaltime=10),
  list(modulenetname = "consecutive_bifurcating", totaltime=6),
  list(modulenetname = "bifurcating_convergence", totaltime=8),
  list(modulenetname = "trifurcating", totaltime=8),
  list(modulenetname = "bifurcating_cycle", totaltime=20),
  list(modulenetname = "bifurcating_loop", totaltime=20)
)
takesettings = list(
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
    id = paste0(.version, modelsetting$modelgenerator$modulenetname, "_", modelsetting$replicate)
  )
  model$modelsetting <- modelsetting
  model
})
saver(models, "models")

## Simulate models----------------
simulations <- mclapply(models, function(model) {
  simulations <-run_simulations(model, model$modelsetting$modelgenerator$totaltime, 2, 32, local=FALSE)
}, mc.cores=8)
simulations <- map(simulations, function(x) {x$simulations <- smoothe_simulations(x$simulations, x$model); x})
simulations <- map(simulations, function(simulation) {
  simulation$info <- list(
    id = paste0(simulation$model$info$id),
    modelid = simulation$model$info$id
  )
  simulation
})
saver(simulations, "simulations")

sub = simulations %>% map(function(simulations) {
  simulations$simulations = map(simulations$simulations, ~list(expression_modules = .$expression_modules))
  simulations
})

## Extract gold standards-------------------
goldstandards <- pbapply::pblapply(sub, function(x){
  print("----")
  library(magrittr)
  library(tidyverse)
  library(dambiutils)
  dyngen:::extract_goldstandard(x, smooth=FALSE)
})#, mc.cores=8)#, qsub_environment = list2env(list()))
goldstandards <- map2(goldstandards, simulations, function(gs, simulation) {
  gs$info <- list(id=simulation$info$id, simulationid = simulation$info$id, modelid = simulation$info$modelid)
  gs
})

saver(goldstandards, "goldstandards")

## Extract experiments-------------------
experimentsettings <- expand.grid(simulation=seq_along(simulations), takesetting=seq_along(takesettings))
experiments <- map2(experimentsettings$simulation, experimentsettings$takesetting, function(simulationid, takesettingid) {
  experiment <- run_experiment(simulations[[simulationid]], takesettings[[takesettingid]])
  experiment$info$modelid <- models[[simulationid]]$info$id
  experiment$info$goldstandardid <- goldstandards[[simulationid]]$info$id
  experiment$info$simulationid <- simulations[[simulationid]]$info$id
  experiment$info$id <- paste0(models[[simulationid]]$info$id, "_", takesettings[[takesettingid]]$type)
  experiment
})
saver(experiments, "experiments")

# plot gold standards-------------------
walk(experiments, function(experiment) {
  print(experiment$info$id)
  gs <- goldstandards %>% keep(~(.$info$id == experiment$info$goldstandardid)) %>% first
  
  plot = plot_goldstandard(experiment, gs) %>% unlist(recursive=FALSE) %>% 
    cowplot::plot_grid(plotlist = ., ncol=2) %>% add_sub(experiment$info$id)
  
  dir.create(dirname(paste0("images/", experiment$info$id, ".png")), recursive = TRUE)
  ggsave(paste0("images/", experiment$info$id, ".png"), plot, width = 10, height=10)
})


## Run scRNAseq----------------------
platforms = readr::read_tsv("data/platforms.tsv")
# platforms2 = expand.grid(sequencerate = c(seq_len(10)/10), cellcapturerate = seq_len(10)/10) %>% as.data.frame()
# platforms = platforms[1, ] %>% select(-sequencerate, -cellcapturerate) %>% data.frame(platforms2)
# platforms$platformname = seq_len(nrow(platforms)) %>% as.character
datasets = lapply(experiments, function(experiment) {
  lapply(seq_len(nrow(platforms)), function(platformid) {
    platform = platforms[platformid, ] %>% as.list
    dataset = dyngen:::run_scrnaseq(experiment, platform)
    dataset$info = experiment$info
    dataset$info$experimentid <- experiment$info$id
    dataset$info$id <- paste0(experiment$info$id, "_", platform$platformname)
    dataset
  })
}) %>% unlist(recursive=F)
saver(datasets, "datasets")

## Sync to prism
PRISM:::rsync_remote("", "~/thesis/projects/dyngen/results", "prism", "/group/irc/shared/dyngen_results")

## Download from prism
PRISM:::rsync_remote("prism", "/group/irc/shared/dyngen_results/results", "", "~/Workspace/papers/ti_eval/dyngen")
