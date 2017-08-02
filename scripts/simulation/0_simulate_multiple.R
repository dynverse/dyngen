library(dplyr)

.datasets_location = "/home/wouters/thesis/projects/dyngen/results/"
.version = 2


refresh_datasets();refresh_experiments();refresh_goldstandards();refresh_models()


list_datasets()

## Run experiments
# Experiment settings
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

.version = "3/"

modelsettings <- expand_lists(modelgenerator=modelgenerators, replicate=seq_len(nreplicates))
models <- map(modelsettings, function(modelsetting) {
  model <- generate_model(modulenetname = modelsetting$modelgenerator$modulenetname)
  model$info$id <- paste0(.version, modelsetting$modelgenerator$modulenetname, "_", modelsetting$replicate)
  model$modelsetting <- modelsetting
  model
})
models %>% walk(save_model)

simulations <- mclapply(models, function(model) {
  run_simulations(model, model$modelsetting$modelgenerator$totaltime, 2, 32, local=FALSE)
}, mc.cores=16)


goldstandards = map(simulations[1], function(x){
  print("----")
  library(magrittr)
  library(tidyverse)
  library(dambiutils)
  dyngen:::extract_goldstandard(x, smooth=TRUE)
})#, qsub_environment = list2env(list()))




experimentsettings <- expand.grid(simulation=seq_along(simulations), takesetting=seq_along(takesettings))
experiments <- map2(experimentsettings$simulation, experimentsettings$takesetting, function(simulationid, takesettingid) {
  experiment <- run_experiment(simulations[[simulationid]], takesettings[[takesettingid]])
  experiment$info$modelid <- models[[simulationid]]$info$id
  experiment$info$id <- paste0(models[[simulationid]]$info$id, "_", takesettings[[takesettingid]]$type)
  experiment
})

###

experimentsettings <- list(
  tibble(treeseed = c(10, 100, 1000, 10000, 100000, 1000000, 20, 200, 2000, 20000, 200000, 2000000), totaltime=20)
) %>% bind_rows() %>% mutate(ncells=500, experimentname=paste0("tree_", seq_along(treeseed)))
models = map(experimentsettings$treeseed, ~generate_model(treeseed = .))
models %>% walk(save_model)

###
experiments %>% walk(save_experiment)

experimentids = readRDS("results/experiments.rds")$id[grepl("2017_08_01", readRDS("results/experiments.rds")$id)]
experiments = map(experimentids, load_experiment, contents_experiment(T, T, T, T, T, T, T))

#dyngen:::extract_goldstandard(experiments[[3]])

## Extract gold standard

experiments = map(experiments, function(x) {x$simulations <- smoothe_simulations(x$simulations, x$model); x})

goldstandards = mclapply(experiments, function(x){
  library(magrittr)
  library(tidyverse)
  library(dambiutils)
  dyngen:::extract_goldstandard(x, smooth = FALSE)
}, mc.cores = 8)#, qsub_environment = list2env(list()))
walk(goldstandards, save_goldstandard)
remove_duplicate_goldstandards()


plots = map2(experiments, goldstandards, function(experiment, gs) {
  cowplot::plot_grid(plotlist=unlist(plot_goldstandard(experiment, gs), recursive=FALSE), ncol=2)
})


experiment = load_experiment(datasetsinfo$experimentid[[10]], contents_experiment(T, T, T, T, T, T, T))
gs = extract_goldstandard(experiment)
gs %>% save_goldstandard()
remove_duplicate_goldstandards()


## Run scRNAseq
platforms = readr::read_tsv("data/platforms.tsv")
# platforms2 = expand.grid(sequencerate = c(seq_len(10)/10), cellcapturerate = seq_len(10)/10) %>% as.data.frame()
# platforms = platforms[1, ] %>% select(-sequencerate, -cellcapturerate) %>% data.frame(platforms2)
# platforms$platformname = seq_len(nrow(platforms)) %>% as.character
datasets = lapply(experiments, function(experiment) {
  lapply(seq_len(nrow(platforms)), function(platformid) {
    platform = platforms[platformid, ] %>% as.list
    dataset = dyngen:::run_scrnaseq(experiment, platform)
  })
}) %>% unlist(recursive=F)
datasets %>% walk(save_dataset)

## Sync to prism
PRISM:::rsync_remote("", "~/thesis/projects/dyngen/results", "prism", "/group/irc/shared/dyngen_results")


## Download from prism
PRISM:::rsync_remote("prism", "/group/irc/shared/dyngen_results/results", "", "~/Workspace/papers/ti_eval/dyngen")
