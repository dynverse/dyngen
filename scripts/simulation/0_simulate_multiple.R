library(dplyr)

## Run experiments
# Experiment settings
nreplicates <- 2
experimentsettings <- list(
  tibble(modulenetname = "linear", totaltime=4, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "cycle", totaltime=20, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "consecutive_bifurcating", totaltime=6, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "bifurcating_convergence", totaltime=8, replicate=seq_len(nreplicates)),
  tibble(modulenetname = "trifurcating", totaltime=8, replicate=seq_len(nreplicates))
) %>% bind_rows() %>% mutate(ncells=500, experimentname=paste0(modulenetname, "_", replicate))

run_settingid <- function(experimentid) {
  setting <- experimentsettings[experimentid,] %>% as.list()
  print(setting)
  do.call(wrapper, setting)
}

wrapper <- function(modulenetname, totaltime, replicate, ncells, experimentname) {
  dyngen::run_experiment(modulenetname, totaltime, ncells=ncells)
}

experiments <- parallel::mclapply(seq_len(nrow(experimentsettings)), run_settingid, mc.cores = 8)

## Extract gold standard
gss <- mclapply(experiments, function(experiment) {
  gs <- extract_goldstandard(experiment, verbose=T)
  gs
}, mc.cores=4)

experiments <- map2(experiments, gss, function(experiment, gs) {
  experiment$gs <- gs
  experiment
})

experiment <- experiments[[5]]

## Run scRNAseq
platforms <- read_tsv("data/platforms.tsv")
# platforms2 <- expand.grid(sequencerate = c(seq_len(10)/10), cellcapturerate = seq_len(10)/10) %>% as.data.frame()
# platforms <- platforms[1, ] %>% select(-sequencerate, -cellcapturerate) %>% data.frame(platforms2)
# platforms$platformname <- seq_len(nrow(platforms)) %>% as.character
datasets <- lapply(experiments, function(experiment) {
  lapply(seq_len(nrow(platforms)), function(platformid) {
    platform <- platforms[platformid, ] %>% as.list
    dataset <- run_scrnaseq(experiment, platform)
  })
}) %>% unlist(recursive=F)



saveRDS(datasets, file="results/datasets.rds")

datasets <- readRDS("results/datasets.rds")

##

list2env(datasets[[12]], .GlobalEnv)
list2env(dataset, .GlobalEnv)
