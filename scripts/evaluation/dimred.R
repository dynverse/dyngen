settings = expand.grid(datasetid=seq_len(length(datasets)), dimredname=names(dimreds))%>% mutate(settingid=seq_len(length(datasetid))) %>% split(1:nrow(.)) %>% map(as.list)

spaces = settings %>% keep(~T) %>% map(~list(counts=datasets[[.$datasetid]]$counts, dimred=dimreds[[.$dimredname]])) %>% 
  qsub_lapply(
    function(x) x$dimred(x$counts), 
    override_qsub_config(wait=F, stop_on_error = F),
    qsub_environment=c("process_dimred")
  )
spaces = qsub_retrieve(spaces)
settings = map2(settings, spaces, function(setting, space) {setting$space = space;setting})
map_lgl(settings, ~is.na(.$space)[[1]])
settings = map_if(settings, ~is.na(.$space)[[1]], function(x) {
  counts = datasets[[x$datasetid]]$counts
  x$space=replicate(nrow(counts), rnorm(3)) %>% t
  rownames(x$space) = rownames(counts)
  x
})
settings = map(settings, function(x) {
  counts = datasets[[x$datasetid]]$counts
  rownames(x$space) = rownames(counts)
  x
})

get_knn = function(distances, k=100) {
  diag(distances) = max(distances)+0.1
  knn = apply(distances, 1, function(row) names(row)[order(row)[1:k]])
  knn = split(knn, rep(1:ncol(knn), each=nrow(knn)))
  names(knn) = rownames(distances)
  knn
}


## overlap of nearby cells
jaccard = function(x1, x2) length(intersect(x1, x2))/length(union(x1, x2))
cellcomparison = map(settings, function(setting) {
  dataset = datasets[[setting$datasetid]]
  space = setting$space
  
  ## first get the distances between every pair of cells
  original = dataset$expression[rownames(dataset$counts), ]
  simulationtimes = dambiutils::pick(dataset$cellinfo, cell=rownames(dataset$counts))$simulationtime %>% set_names(rownames(dataset$counts))
  
  distances = SCORPIUS::euclidean.distance(space)
  
  k = 10
  knn_observed = get_knn(distances, k)
  knn = get_knn(dataset$gs$celldistances[rownames(dataset$counts), rownames(dataset$counts)], k)
  
  #> simulation time prediction
  observedtimes = map_dbl(knn_observed, function(neighborhood) {mean(simulationtimes[neighborhood])})
  
  #> overlap of nearby cells
  jaccards = map_dbl(1:length(knn), function(i) jaccard(knn[[i]], knn_observed[[i]]))
  
  data.frame(setting, observedtime=observedtimes, simulationtime=simulationtimes, jaccard=jaccards, cell=names(observedtimes))
}) %>% bind_rows()

scores = cellcomparison %>% group_by(datasetid, dimredname) %>% summarise(simcor = cor(observedtime, simulationtime, method="kendall"), avjaccard=mean(jaccard))
datasetinfo = tibble(datasetid=seq_len(length(datasets)), modulenetname=map_chr(datasets, ~.$model$modulenetname), platform=map_chr(datasets, ~.$platform$platformname))
scores %<>% left_join(datasetinfo, by="datasetid")
ggplot(scores) + geom_boxplot(aes(modulenetname, avjaccard, color=factor(modulenetname))) + facet_grid(platform~dimredname)
ggplot(scores) + geom_boxplot(aes(dimredname, simcor, color=factor(dimredname))) + facet_grid(platform~modulenetname)
