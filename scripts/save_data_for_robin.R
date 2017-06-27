datasetids = readRDS("results/datasets.rds") %>% filter(stringr::str_detect(experimentid, "05_04")) %>% .$id
for (i in seq_along(datasetids)) {
  datasetid = datasetids[[i]]
  dataset = dyngen:::load_dataset(datasetid, dyngen:::contents_dataset(goldstandard=NULL))
  
  folder = file.path("results/data_robin/170504/", i)
  dir.create(folder, recursive=T)
  
  piecestatenet = get_piecestatenet(get_piecenet(dataset$model$statenet))
  write.csv(piecestatenet, file.path(folder, "statenet.csv"))
  
  distances = dataset$counts %>% SCORPIUS::correlation.distance()
  write.csv(distances, file.path(folder, "distances.csv"))
  
  write.csv(dataset$counts, file.path(folder, "counts.csv"))
  
  SCORPIUS::quant.scale(log2(dataset$counts + 1)) %>% write.csv(file.path(folder, "counts2.csv"))
  
  if(nrow(piecestatenet)) {
    pdf(file.path(folder, "statenet.pdf"))
    piecestatenet %>% graph_from_data_frame() %>% plot()
    dev.off()
  }
}