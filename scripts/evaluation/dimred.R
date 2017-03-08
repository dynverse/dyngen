settings = expand.grid(datasetid=seq_len(length(datasets)), dimredname=names(dimreds)) %>% split(1:nrow(.))

spaces = settings %>% keep(~T) %>% 
  mclapply2(function(x) dimreds[[x$dimredname]](datasets[[x$datasetid]]$counts))


get_knn = function(distances, k=100) {
  diag(distances) = max(distances)+0.1
  knn = apply(distances, 1, function(row) order(row)[1:k])
  split(knn, rep(1:ncol(knn), each=nrow(knn)))
}


## predict simulation time
simulationtime_comparison = map(settings, ~tibble(datasets[[.$datasetid]]$cellinfo$simulationtime, )
scores = simulationtime_comparison %>% group_by(datasetid, dimredname) %>% summarise(simcor = cor(observedtime, simulationtime))

datasetinfo = tibble(datasetid=seq_len(length(datasets)), modulenetname=map_chr(datasets, ~.$model$modulenetname))
scores %<>% left_join(datasetinfo, by="datasetid")

ggplot(scores) + geom_point(aes(dimredname, simcor, color=factor(modulenetname)))


## overlap of nearby cells
jaccard = function(x1, x2) length(intersect(x1, x2))/length(union(x1, x2))
cellcomparison = map2(settings, spaces, function(setting, space) {
  dataset = datasets[[setting$datasetid]]
  
  ## first get the distances between every pair of cells
  original = dataset$expression[rownames(dataset$counts), ]
  simulationtimes = dataset$cellinfo$simulationtime[rownames(dataset$counts)]
  
  distances = SCORPIUS::euclidean.distance(space)
  
  k = 200
  knn_observed = get_knn(distances, k)
  knn = get_knn(dataset$gs$celldistances[rownames(dataset$counts), rownames(dataset$counts)])
  
  #> simulation time prediction
  observedtimes = map_dbl(knn_observed, function(neighborhood) {mean(simulationtimes[neighborhood])})
  
  #> overlap of nearby cells
  jaccards = map_dbl(1:length(knn), function(i) jaccard(knn[[i]], knn_observed[[i]]))
  
  data.frame(setting, observedtime=observedtimes, simulationtime=simulationtimes, jaccard=jaccards)
}) %>% bind_rows()

scores = cellcomparison %>% group_by(datasetid, dimredname) %>% summarise(simcor = cor(observedtime, simulationtime), avjaccard=mean(jaccard))
scores %<>% left_join(datasetinfo, by="datasetid")
ggplot(scores) + geom_point(aes(dimredname, avjaccard, color=factor(modulenetname))) + facet_wrap(~modulenetname)
