library(dplyr)

source("scripts/simulation/0_simulation.R")

## Run experiments

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
  qsub.conf = qsub.configuration(exec.before = c("module unload python", "module load python", "module load R/x86_64/3.2.2"), memory = "1G", name=experimentname)
  
  run_experiment(modulenetname, totaltime, ncells=ncells, qsub.conf=qsub.conf)
}

experiments = mclapply(seq_len(nrow(experimentsettings)), run_settingid, mc.cores = 8)

experiments = mclapply(experiments, function(experiment) {
  snet = get_branchconnected_statenet(experiment$model$statenet)
  experiment$cell_distances = get_cell_distances(experiment$cellinfo, snet$snet, snet$statenet)
  experiment
}, mc.cores = 8)

## Run scRNAseq
datasets = lapply(experiments, function(experiment) {
  lapply(seq_len(nrow(platforms)), function(platformid) {
    platform = platforms[platformid, ] %>% as.list
    dataset = run_scrnaseq(experiment, platform)
  })
}) %>% unlist(recursive=F)

saveRDS(datasets, file="results/datasets.rds")

##

list2env(datasets[[12]], .GlobalEnv)
list2env(dataset, .GlobalEnv)
