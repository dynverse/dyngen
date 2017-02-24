# intermediate results
# timing of steps


mds = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(Efiltered))),ndim = 3)



## SLICER
k = SLICER::select_k(exprs(Efiltered) %>% t, kmin=5)
traj_lle = lle::lle(exprs(Efiltered) %>% t, m=2, k)$Y
traj_graph = SLICER::conn_knn_graph(traj_lle,5)
ends = SLICER::find_extreme_cells(traj_graph, traj_lle)
start = 1
cells_ordered = SLICER::cell_order(traj_graph, start)
branches = SLICER::assign_branches(traj_graph,start) %>% factor
draw.trajectory.plot(mds, branches)

pseudotime = cells_ordered %>% order %>% {./max(.)} %>% {tibble(time=., cellid=1:length(.), branch=branches)}
rgl::plot3d(space, col=branches)


cells_ordered


## DPT

# setRepositories(ind = 1:2)
# devtools::install_url(paste0('https://www.helmholtz-muenchen.de/fileadmin/ICB/software/',
#                              c('destiny/destiny_1.3.4.tar.gz', 'dpt_0.6.0.tar.gz')))
library(dpt)
ts = dpt::Transitions(t(exprs(Efiltered)))
pt = dpt::dpt(ts, branching = TRUE)
ev = eigen(as.matrix(ts@transitions), TRUE)$vectors
dm = as.data.frame(ev[, -1])
colnames(dm) = paste0('DC', seq_len(ncol(dm)))
qplot(DC1, DC2, data = dm, colour = pt$Branch)
draw.trajectory.plot(mds, pt$Branch)




## MPATH
write.table(exprs(Efiltered), "tmp.tsv", sep="\t")
Mpath::landmark_designation("tmp.tsv")





## Slingshot
#BiocInstaller::biocLite("kstreet13/slingshot")
space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(Efiltered))))
clus = factor(kmeans(space, 10)$cluster)
startcluster = clus %>% {tibble(cell=names(.), cluster=., simulationtime=phenoData(E)$time)} %>% group_by(cluster) %>% summarise(meantime=mean(simulationtime)) %>% {.$cluster[which.min(.$meantime)]}
lineages = slingshot::get_lineages(space, clus, start.clus=startcluster)
slingshot::plot_tree(space, clus, lineages)

SCORPIUS::draw.trajectory.plot(space, clus)
crv = slingshot::get_curves(space, clus, lineages)
slingshot::plot_curves(space, clus, crv)

bind_cols(map(crv, "w")) %>% View



## Monocle

space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(Efiltered))),ndim = 3)
spaces$mds = space## Monocle
num.genes = 1000
num.paths = 3
library(monocle)

counts = Efiltered %>% exprs %>% t
group.index = phenoData(Efiltered)$time %>% cut(5, labels=F)

requireNamespace("dplyr")
counts.mat <- as.matrix(counts) %>% t
pd <- new("AnnotatedDataFrame", data = data.frame(row.names=rownames(counts), num.genes=rowSums(counts!=min(counts))))
fd <- new("AnnotatedDataFrame", data = data.frame(row.names=colnames(counts), times.expressed=colSums(counts!=min(counts))))
cds <- newCellDataSet(counts.mat, phenoData = pd, featureData = fd)

cds <- detectGenes(cds, min_expr = min(counts.mat)+1e-20)

pct.cells.expressed <- fData(cds)$num_cells_expressed / nrow(counts)
var <- apply(counts, 2, var)
range <- apply(counts, 2, function(x) max(x) - min(x))
score <- pct.cells.expressed * var * range

ordering_genes <- names(score)[order(score, decreasing = T)[seq_len(min(num.genes, length(score)))]]

cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, 3)
cds <- orderCells(cds, num_paths=num.paths, reverse=F)

space2 <- t(cds@reducedDimS)
colnames(space2) <- paste0("Comp", seq_len(ncol(space2)))

mst <- minSpanningTree(cds)
diam <- as.data.frame(as.vector(V(mst)[get.diameter(mst, weights = NA)]$name))
colnames(diam) <- c("sample_name")
path <- space2[diam$sample_name,]

value <- dplyr::percent_rank(pData(cds)[,"Pseudotime"])

path <- path[order(value[rownames(path)]),]
rownames(path) <- NULL

statemeantimes = tapply(phenoData(Efiltered)$time, phenoData(cds)$State, mean)
root_state = which.min(statemeantimes)
cds = monocle::orderCells(cds, root_state = root_state)

monocle::plot_cell_trajectory(cds)
SCORPIUS::draw.trajectory.plot(mds, phenoData(cds)$State)
SCORPIUS::draw.trajectory.plot(space2, phenoData(cds)$State)
SCORPIUS::draw.trajectory.plot(space2, phenoData(E)$state)


cds@cellPairwiseDistances %>% pheatmap(cluster_cols=F, cluster_rows=F)
cell_distances_observed = cds@cellPairwiseDistances
