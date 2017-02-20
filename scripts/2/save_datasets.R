for (i in 1:10) {
  source("./scripts/2/0a_simulate_all.R")
  
  saveRDS(list(expression_sequenced=Emrna2, expression_cell=Emrna, net=net), file=paste0("results/expression/linear/", i, ".rds"))
}

datasets = lapply(1:10, function(i) {readRDS(file=paste0("results/expression/linear/", i, ".rds"))})
datasetinfo = tibble(datasetid=1:10, modulenetname=modulenetname)


saveRDS(list(datasets=datasets, datasetinfo=datasetinfo), file="results/expression/linear.rds")