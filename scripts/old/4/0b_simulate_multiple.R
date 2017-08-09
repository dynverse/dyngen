library(dplyr)

datasettings = list(
  tibble(modulenetname = "linear", totaltime=4, replicate=seq_len(5)),
  tibble(modulenetname = "cycle", totaltime=20, replicate=seq_len(5)),
  tibble(modulenetname = "consecutive_bifurcating", totaltime=6, replicate=seq_len(5)),
  tibble(modulenetname = "bifurcating_convergence", totaltime=8, replicate=seq_len(5))
) %>% bind_rows() %>% mutate(ncells=, datasetname=paste0(modulenetname, "_", replicate))

run_settingid = function(datasettingid) {
  setting = datasettings[datasettingid,] %>% as.list()
  print(setting)
  do.call(wrapper, setting)
}

wrapper = function(modulenetname, totaltime, replicate, ncells, datasetname) {
  source("scripts/4/0a_simulate_all.R", local=TRUE)
  
  list(E=E, Emrna=Emrna, counts=counts, modulenet=modulenet, modulenodes=modulenodes, genes=genes, net=net)
}

datasets = mclapply(seq_len(nrow(datasettings)), run_settingid, mc.cores = 2)
saveRDS(datasets, file="results/datasets.rds")

##

datasets = readRDS("results/datasets.rds")

list2env(datasets[[12]], .GlobalEnv)
list2env(dataset, .GlobalEnv)
