library(dplyr)

source("scripts/simulation/0_simulation.R")

## Run experiments
# Experiment settings
nreplicates = 2
experimentsettings = list(
  tibble(modulenetname = "linear", totaltime=4, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "cycle", totaltime=20, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "consecutive_bifurcating", totaltime=6, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "bifurcating_convergence", totaltime=8, replicate=seq_len(nreplicates))
) %>% bind_rows() %>% mutate(ncells=500, experimentname=paste0(modulenetname, "_", replicate))

run_settingid = function(experimentid) {
  setting = experimentsettings[experimentid,] %>% as.list()
  print(setting)
  do.call(wrapper, setting)
}

wrapper = function(modulenetname, totaltime, replicate, ncells, experimentname) {
  run_experiment(modulenetname, totaltime, ncells=ncells)
}

experiments = mclapply(seq_len(nrow(experimentsettings)), run_settingid, mc.cores = 8)

experiments = mclapply(experiments, function(experiment) {
  snet = get_branchconnected_statenet(experiment$model$statenet)
  experiment$cell_distances = get_cell_distances(experiment$cellinfo, snet$snet, snet$statenet)
  experiment
}, mc.cores = 8)

## Run scRNAseq
platforms = read_tsv("data/platforms.tsv")
# platforms2 = expand.grid(sequencerate = c(seq_len(10)/10), cellcapturerate = seq_len(10)/10) %>% as.data.frame()
# platforms = platforms[1, ] %>% select(-sequencerate, -cellcapturerate) %>% data.frame(platforms2)
# platforms$platformname = seq_len(nrow(platforms)) %>% as.character
datasets = lapply(experiments, function(experiment) {
  lapply(seq_len(nrow(platforms)), function(platformid) {
    platform = platforms[platformid, ] %>% as.list
    dataset = run_scrnaseq(experiment, platform)
  })
}) %>% unlist(recursive=F)

## Extract gold standard
datasets = lapply(datasets, function(dataset) {
  dataset$gs = dataset$gs = extract_goldstandard(dataset, verbose=T)
  dataset
})

saveRDS(datasets, file="results/datasets.rds")

##

list2env(datasets[[12]], .GlobalEnv)
list2env(dataset, .GlobalEnv)
