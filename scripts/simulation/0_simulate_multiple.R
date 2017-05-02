library(dplyr)

.datasets_location = "/home/wouters/thesis/projects/dyngen/results/"

## Run experiments
# Experiment settings
nreplicates <- 2
experimentsettings <- list(
  tibble(modulenetname = "linear", totaltime=4, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "cycle", totaltime=20, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "consecutive_bifurcating", totaltime=6, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "bifurcating_convergence", totaltime=8, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "trifurcating", totaltime=8, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "bifurcating_cycle", totaltime=20, replicate=seq_len(nreplicates))
) %>% bind_rows() %>% mutate(ncells=500, experimentname=paste0(modulenetname, "_", replicate))

models = map(experimentsettings$modulenetname, generate_model)
models %>% walk(save_model)

run_settingid <- function(experimentid) {
  setting <- experimentsettings[experimentid,] %>% as.list()
  experiment <- dyngen:::run_experiment(models[[experimentid]], setting$totaltime, ncells=setting$ncells)
}



#experimentid = first(which(experimentsettings$modulenetname == "bifurcating_cycle"))
#experiment = run_settingid(experimentid)
#experiments = list(experiment)

experiments = parallel::mclapply(seq_len(nrow(experimentsettings)), run_settingid, mc.cores = 8)
experiments %>% walk(save_experiment)

experiments = map(readRDS("results/experiments.rds")$id, load_experiment, contents_experiment(T, T, T, T, T, T, T))

## Extract gold standard
goldstandards = qsub_lapply(experiments, function(x){
  library(magrittr)
  library(tidyverse)
  library(dambiutils)
  dyngen:::extract_goldstandard(x)
}, qsub_environment = list2env(list()))
walk(goldstandards, save_goldstandard)
remove_duplicate_goldstandards()


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