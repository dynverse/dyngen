Emarques = read_tsv("~/GSE75330_Marques_et_al_mol_counts2.tab")
Emarques %>% rownames

genes = Emarques$cellid
Emarques = Emarques[,2:ncol(Emarques)] %>% as.matrix
rownames(Emarques) = genes

apply(Emarques, 1, sd) %>% log %>% qplot

genesoi = (apply(Emarques, 1, max) > 1) & (apply(Emarques, 1, sd) > quantile(apply(Emarques, 1, sd), 0.8))
Emarques2 = Emarques[genesoi,]


library(SCORPIUS)

cellcors = correlation.distance(t(Emarques2))
clustering = hclust(as.dist(cellcors), method="ward.D2")
cellabels = cutree(clustering, k=14)
Eclusters = limma::avearrays(Emarques, cellabels)
wrongclusters = colnames(Eclusters)[Eclusters["Lum",] > 0.5]
correctcells = !(cellabels %in% wrongclusters)
Emarques = Emarques[,correctcells]
cellcors = cellcors[correctcells,correctcells]
cellabels = cellabels[correctcells]

## MDS
space = reduce.dimensionality(cellcors, ndim = 3)

## TSNE
space = tsne::tsne(as.dist(cellcors), k = 2)
colnames(space) = paste0("Comp", 1:ncol(space))

save(Emarques, cellcors, cellabels, space, file="results/marques.RData")

load(file="results/marques.RData")

seed = set.seed(7)
traj = infer.trajectory(space)
SCORPIUS::draw.trajectory.plot(space, factor(cellabels), path = traj$path)
SCORPIUS::draw.trajectory.plot(space, apply(Emarques, 2, sum), path = traj$path) + viridis::scale_color_viridis(trans="log", option="A", direction=-1)
SCORPIUS::draw.trajectory.plot(space, Emarques[gimp$gene[1] %>% as.character, ] + 1, path = traj$path) + viridis::scale_color_viridis(trans="log", option="A", direction=-1)
SCORPIUS::draw.trajectory.plot(space, Emarques["Sox10", ] + 1, path = traj$path) + viridis::scale_color_viridis(trans="log", option="A", direction=-1)

#gimp <- gene.importances(Emarques[(apply(Emarques, 1, sd) > quantile(apply(Emarques, 1, sd), 0.95)),] %>% t, traj$time, num.permutations = 0, ntree=100)
Erf = Emarques[(apply(Emarques, 1, sd) > quantile(apply(Emarques, 1, sd), 0.9)),] %>% t
sample = sample(rownames(Erf), 200)
result = randomForest::randomForest(Erf[sample,], traj$time[sample], ntree=500, importance=T)
gimp = result$importance %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% rename(importance=IncNodePurity) %>% arrange(-importance)

Goi = c("Id2", "Sox10", "Ascl1", "Tcf7l2", "Id4", "Mal", "Mog", "Plp1", "Myrf", "Olig2")
SCORPIUS::draw.trajectory.heatmap(Emarques[Goi,] %>% t, traj$time, scale.features = T, show.labels.row = T)
Goi = c("Ngfr", "Sox10")
Emarques[Goi,] %>% t %>% scale %>% t %>%reshape2::melt(varnames=c("gene", "cell"), value.name="expression") %>% 
  mutate(time=traj$time[cell]) %>% 
  ggplot(aes(time, expression, color=gene)) + geom_smooth()



expr <- Emarques[gimp$gene[1:100] %>% as.character,] %>% t %>% quant.scale
modules <- extract.modules(expr)
SCORPIUS::draw.trajectory.heatmap(expr, traj$time, modules = modules, scale.features = F)
