library(SCORPIUS)

Efiltered = counts
Efiltered =  Efiltered[apply(exprs(Efiltered), 1, sd) > 0,apply(exprs(Efiltered), 2, sd) > 0]


spaces = list()
# PCA
space = prcomp(exprs(Efiltered))$rotation[,1:3]
colnames(space) = paste0("Comp", 1:ncol(space))
spaces$pca = space

# MDS
space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(Efiltered))),ndim = 3)
spaces$mds = space

# TSNE
space = tsne::tsne(as.dist(correlation.distance(t(exprs(Efiltered)))), k = 2)
colnames(space) = paste0("Comp", 1:ncol(space))
spaces$tsne = space

# DIFFUSION MAP
space = diffusionMap::diffuse(as.dist(correlation.distance(t(exprs((Efiltered))))), neigen=2)
space = space$X[,c(1:2)]
colnames(space) = paste0("Comp", 1:ncol(space))
spaces$diffusion_map = space

# ICA
space = fastICA::fastICA(t(scale(exprs(Efiltered))), 3)$S
colnames(space) = paste0("Comp", 1:ncol(space))
spaces$ica = space

# LLE
ndim = 2
k = lle::calc_k(t(scale(exprs(Efiltered))), ndim)
k = k$k[which.min(k$rho)]
space = lle::lle(t(exprs(Efiltered)), ndim, k)$Y
colnames(space) = paste0("Comp", 1:ncol(space))
spaces$lle = space

# ForceAtlas2
distances = correlation.distance(t(exprs((Efiltered)))) %>% reshape2::melt(varnames=c("from", "to"), value.name="distance")
distances %>% .$distance %>% hist
distancecutoff = 0.5
#graph = distances %>% filter(distance > distancecutoff) %>% graph_from_data_frame()
graph = distances %>% group_by(from) %>% top_n(3, distance) %>% graph_from_data_frame()
layout = ForceAtlas2::layout.forceatlas2(graph, directed=F, plotstep=9999999, iterations=3000)
plot(graph, layout=layout)
space = layout
colnames(space) = paste0("Comp", 1:ncol(space))
spaces$forceatlasv2 = space


lapply(names(spaces), function(spacename) {
  space = spaces[[spacename]]
  plotdata = as.data.frame(space)
  plotdata$progression = phenoData(Efiltered)$time
  plotdata$state = phenoData(E)$state
  plots = list(
    ggplot(plotdata) + geom_jitter(aes(Comp1, Comp2, color=progression), width=.01, height=.01) + viridis::scale_color_viridis(option = "A", direction = -1, end=0.8) + theme_classic() + guides(color=FALSE) + ggtitle(spacename)
  )
  
  if ("state" %in% colnames(plotdata)) {
    plots = c(plots, list(ggplot(plotdata) + geom_jitter(aes(Comp1, Comp2, color=factor(state)), width=.01, height=.01) + theme_classic() + guides(color=FALSE)))
  }

  plots
}) %>% unlist(recursive=F) %>% cowplot::plot_grid(plotlist=., ncol=2)


rgl::plot3d(space, col=viridis::magma(20)[phenoData(Efiltered)$time %>% cut(20, labels=F)])


## Trajectory inference
plotdata = as.data.frame(space)
plotdata$progression = phenoData(Efiltered)$time
ggplot(plotdata) + geom_jitter(aes(Comp1, Comp2, color=progression), width=.01, height=.01) + viridis::scale_color_viridis(option = "A", direction = -1, end=0.8) + theme_classic()

trajectory = infer.trajectory(space)
draw.trajectory.plot(space, phenoData(Efiltered)$time, trajectory$path) + scale_colour_distiller(palette = "RdYlBu") + theme_classic()
rownames(E) = c(1:nrow(E))
draw.trajectory.heatmap(t(exprs(Efiltered)), trajectory$time, as.factor(cut(phenoData(Efiltered)$time, breaks=length(phenoData(Efiltered)$time), labels=F)), show.labels.row = T)
draw.trajectory.heatmap(t(exprs(Efiltered)), phenoData(Efiltered)$time, as.factor(cut(phenoData(Efiltered)$time, breaks=length(phenoData(Efiltered)$time), labels=F)), show.labels.row = T)

plotdata = data.frame(observed_time = trajectory$time, known_time=phenoData(Efiltered)$time)
ggplot(plotdata) + geom_point(aes(trajectory$time, known_time))


# GNG
mds = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(Efiltered))),ndim = 20)
result = GNG::gng(mds)
spaces$gng = result$space.proj

edges = result$edges %>% bind_cols(result$node.space[result$edges$i,] %>% as.data.frame %>% set_names(paste0("i_", colnames(result$node.space)))) %>% bind_cols(result$node.space[result$edges$j,] %>% as.data.frame %>% set_names(paste0("j_", colnames(result$node.space))))
draw.trajectory.plot(space, phenoData(Efiltered)$time)  + viridis::scale_color_viridis(option = "A", direction = -1, end=0.8) + theme_classic() + guides(color=FALSE) + geom_point(aes(Comp1, Comp2), data=result$node.space %>% as.data.frame, fill="#3793d6", color="black",pch=21, size=5) + geom_segment(aes(x=i_Comp1, y=i_Comp2, xend=j_Comp1, yend=j_Comp2), data=edges)


edges = result$edges %>% bind_cols(result$node.proj[result$edges$i,] %>% as.data.frame %>% set_names(paste0("i_", colnames(result$node.proj)))) %>% bind_cols(result$node.proj[result$edges$j,] %>% as.data.frame %>% set_names(paste0("j_", colnames(result$node.proj))))
space.proj = result$space.proj
ggplot() + geom_point(aes(GNG_X, GNG_Y, colour = time), data.frame(space.proj, time = phenoData(Efiltered)$time)) + scale_colour_distiller(palette = "RdYlBu") + theme_classic() + geom_point(aes(GNG_X, GNG_Y), data=result$node.proj %>% as.data.frame) + geom_segment(aes(x=i_GNG_X, y=i_GNG_Y, xend=j_GNG_X, yend=j_GNG_Y), data=edges) + coord_equal()

