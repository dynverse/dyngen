## Get cell distances in gold standard
# Goal: get a pairwise distance matrix between all cells, based on their branch and progression in the gold standard

# network between states
# at bifurcation points, the two states after a bifurcation are connected through a head2head edge
statenet = modulenet %>% filter(effect == 1) %>% filter(from != to) %>% select(from, to) %>% mutate(type="head2tail")
statenet = statenet %>% group_by(from) %>% summarise(n=n()) %>% filter(n>1) %>% .$from %>% map(function(fromoi) {
  statesoi = statenet %>% filter(from==fromoi) %>% .$to
  gtools::combinations(length(statesoi), r=2, v=statesoi) %>% as.data.frame() %>% rename(from=V1, to=V2)
}) %>% bind_rows() %>% mutate(type="head2head") %>% bind_rows(statenet)
igraph::graph_from_data_frame(statenet) %>% plot(edge.color=factor(statenet$type))
snet = statenet %>% graph_from_data_frame(vertices=unique(statenet$from, statenet$to) %>% order)

# head2head = weight 2 (double the path)
# head2tail = weight 1 ()
# start and end state: check 

# the distance if the states are the same is just the difference in progression time within this stage
traj_distance_samestate = function(time1, time2, verbose=F) abs(time1-time2)


get_state_edge = function(state1, state2) {
  statenet %>% filter((from == state1 & to==state2) | (from == state2 & to == state1)) %>% as.list()
}

# processessing the shortest path between 2 states through the state network
# returns a function which can be used to get the distance of two cells in these two respective states
# we calculate these functions here so that they do not need to be recalculated for every pair of cells in the same two states
process_sp = function(sp) {
  extracost = 0
  if(length(sp) == 1) {
    traj_distance_samestate
  } else {
    edge = get_state_edge(sp[[1]], sp[[2]])
    if(edge$type == "head2head" || edge$from == sp[[2]]) {
      start = "from0"
    } else {
      start = "from1"
    }
    if(edge$type == "head2head" && length(sp) > 2) {
      #extracost = 1
    }
    
    edge = get_state_edge(sp[[length(sp)-1]], sp[[length(sp)]])
    if(edge$type == "head2head" || edge$from == sp[[length(sp)-1]]) {
      end = "from0"
    } else {
      end = "from1"
    }
    
    extracost = extracost + length(sp) - 2
    
    function(time1, time2, verbose=F) {
      #if(verbose) print(time1);print(time2);print(extracost)
      totalcost = extracost
      if(start == "from0") {
        totalcost = totalcost + time1
      } else {
        totalcost = totalcost + 1 - time1
      }
      if(end == "from0") {
        totalcost = totalcost + time2
      } else {
        totalcost = totalcost + 1 - time2
      }
      totalcost
    }
  }
}
get_cell_distance_func = function(snet, state1, state2) {
  state1 = as.character(state1) # because igraph works with character names and order of vertices is not numerical
  state2 = as.character(state2)
  shortest_paths(snet, state1, state2, mode = "all")$vpath[[1]] %>% names() %>% as.integer() %>% process_sp()
}
func = get_cell_distance_func(snet, 1, 2)
func(0.9, 0.1)
func(0.2, 0.3)

phenoData(E)$state = factor(phenoData(E)$state) # make sure this is a factor, because if some states are missing this will mess up the indexing
allstates = unique(phenoData(E)$state) %>% sort
cell_distance_funcs = map(allstates, function(state1) map(allstates, ~get_cell_distance_func(snet, state1, .)))

get_cell_distance = function(cell_distance_funcs, cell1, cell2) {
  time1 = phenoData(E)$progression[[cell1]]
  state1 = phenoData(E)$state[[cell1]]
  time2 = phenoData(E)$progression[[cell2]]
  state2 = phenoData(E)$state[[cell2]]
  
  #print(list(time1, state1, time2, state2))
  
  cell_distance_funcs[[state1]][[state2]](time1, time2)
}

allcells = 1:ncol(E)
cell_distances = lapply(allcells, function(cell1) {
  time1 = phenoData(E)$progression[[cell1]]
  state1 = phenoData(E)$state[[cell1]]
  lapply(allcells, function(cell2) {
    time2 = phenoData(E)$progression[[cell2]]
    state2 = phenoData(E)$state[[cell2]]
    cell_distance_funcs[[state1]][[state2]](time1, time2)
  }) %>% as.double
}) %>% {matrix(unlist(.), ncol=length(allcells))}

cellsoi = phenoData(E)$state %in% c(1,2)
cell_distances[] %>% pheatmap(cluster_cols=F, cluster_rows=F)
cell_distances

#### TI evaluation ####
# Evaluate TI methods based on a given cell_distances_observed, all pairwise distances between cells generated from trajectory model
## rank correlation: for every cell, how well do the rankings of cell distances correspond?
cellcors = map_dbl(1:ncol(Efiltered), function(cellid) cor(cell_distances[,cellid], cell_distances_observed[,cellid], method="spearman"))
cellcors %>% hist
cellcors %>% mean

## nearest neighbours: for every cell, how well do its neigbours correspond?
k = 100
knn = apply(cell_distances, 1, function(row) order(row)[1:k])
knn_observed = apply(cell_distances_observed, 1, function(row) order(row)[1:k])
knn %>% pheatmap(cluster_cols=F, cluster_rows=F)
knn = split(knn, rep(1:ncol(knn), each=nrow(knn)))
knn_observed = split(knn_observed, rep(1:ncol(knn_observed), each=nrow(knn_observed)))
knn_random = lapply(1:length(knn), function(i) sample(1:ncol(cell_distances), k))

jaccard = function(x1, x2) length(intersect(x1, x2))/length(union(x1, x2))
jaccards = map_dbl(1:length(knn), function(i) jaccard(knn[[i]], knn_observed[[i]]))
jaccards_random = map_dbl(1:length(knn), function(i) jaccard(knn[[i]], knn_random[[i]]))

plotdata = bind_rows(
  tibble(jaccard=jaccards, type="observed"),
  tibble(jaccard=jaccards_random, type="random")
)
ggplot(plotdata) + geom_histogram(aes(jaccard, group=type, fill=type), alpha=0.7, position="identity")


#### Dimensionality reduction evaluation ####
get_knn = function(distances, k=100) {
  diag(distances) = max(distances)+0.1
  knn = apply(distances, 1, function(row) order(row)[1:k])
  split(knn, rep(1:ncol(knn), each=nrow(knn)))
}

## first get the distances between every pair of cells
original = exprs(Efiltered)
celltimes = phenoData(Efiltered)$time

space_distances = lapply(spaces, euclidean.distance)
space_distances$original = correlation.distance(t(original))
space_distances$random = space_distances$original %>% sample %>% matrix(nrow=nrow(space_distances$original), dimnames = dimnames(space_distances$original))

scores = tibble(space=character())

## predict simulation/progression time
k = 20
# get for every cell its knn and calculate the average time
celltimes_observeds = lapply(space_distances, function(distances) {
  knn = get_knn(distances, k)
  map_dbl(knn, function(knn) {mean(celltimes[knn])})
})
celltimes_observeds = names(celltimes_observeds) %>% map(~tibble(space=., cell=names(celltimes_observeds[[.]]), time_observed=celltimes_observeds[[.]])) %>% bind_rows %>% mutate(cell=as.integer(cell))

celltimes_observeds = tibble(time=celltimes, cell=1:length(celltimes)) %>% right_join(celltimes_observeds, on="cell")

ggplot(celltimes_observeds) + geom_point(aes(time, time_observed)) + facet_wrap(~space)

scores = celltimes_observeds %>% group_by(space) %>% summarise(score=cor(time, time_observed, method="spearman")) %>% mutate(scorename = "simtimecor") %>% left_join(scores, by="space")

# nearby cells
k = 25
knn_observeds = lapply(space_distances, function(distances) {
  get_knn(distances, k)
})
knn = get_knn(cell_distances, k)
knn_observeds$random = lapply(1:length(knn), function(i) sample(1:ncol(cell_distances), k))
knn_observeds$state = get_knn(cell_distances, k)

jaccards = lapply(names(knn_observeds), function(spacename) {
  knn_observed = knn_observeds[[spacename]]
  jaccards = map_dbl(1:length(knn), function(i) jaccard(knn[[i]], knn_observed[[i]]))
  tibble(jaccard=jaccards, cells=1:length(knn_observed), space=spacename)
}) %>% bind_rows()
scores = jaccards %>% group_by(space) %>% summarise(score=mean(jaccard)) %>% mutate(scorename="meanjaccard") %>% bind_rows(scores)



scores %>% filter(space!="state") %>% ggplot() + geom_bar(aes(space, score), stat="identity") + facet_wrap(~scorename) + coord_flip()
