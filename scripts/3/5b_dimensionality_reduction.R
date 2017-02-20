library(SCORPIUS)

Efiltered =  Emrna2[apply(exprs(Emrna2), 1, sd) > 0,apply(exprs(Emrna2), 2, sd) > 0]

# MDS
space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(Efiltered))),ndim = 3)

# TSNE
space = tsne::tsne(as.dist(correlation.distance(t(exprs(Efiltered)))), k = 2)
colnames(space) = paste0("Comp", 1:ncol(space))

# DIFFUSION MAP
space = diffusionMap::diffuse(as.dist(correlation.distance(t(exprs((Efiltered))))), neigen=3)
space = space$X[,c(1:2)]
colnames(space) = paste0("Comp", 1:ncol(space))


rgl::plot3d(space, col=viridis::magma(100)[cut(phenoData(Efiltered)$time, 100)], box=F, size=4)


## Trajectory inference
plotdata = as.data.frame(space)
plotdata$progression = phenoData(Efiltered)$time
ggplot(plotdata) + geom_point(aes(Comp1, Comp2, color=progression)) +  scale_colour_distiller(palette = "RdYlBu") + theme_classic()

trajectory = infer.trajectory(space)
draw.trajectory.plot(space, phenoData(Efiltered)$time, trajectory$final.path) + scale_colour_distiller(palette = "RdYlBu") + theme_classic()
rownames(E) = c(1:nrow(E))
draw.trajectory.heatmap(t(exprs(Efiltered)), trajectory$time, as.factor(cut(phenoData(Efiltered)$time, breaks=length(phenoData(Efiltered)$time), labels=F)), show.labels.row = T)
draw.trajectory.heatmap(t(exprs(Efiltered)), phenoData(Efiltered)$time, as.factor(cut(phenoData(Efiltered)$time, breaks=length(phenoData(Efiltered)$time), labels=F)), show.labels.row = T)

plotdata = data.frame(observed_time = trajectory$time, known_time=phenoData(Efiltered)$time)
ggplot(plotdata) + geom_point(aes(trajectory$time, known_time))
