library(tidyverse)
library(cowplot)

expression1 = read_tsv("~/Downloads/GSE70240_Gmp.txt") %>% as.data.frame
rownames(expression1) = expression1$uid
expression1 = as.matrix(expression1[, 2:ncol(expression1)])

expression2 = read_tsv("~/Downloads/GSE70236_Cmp.txt") %>% as.data.frame
rownames(expression2) = expression2$uid
expression2 = as.matrix(expression2[, 2:ncol(expression2)])

expression3 = read_tsv("~/Downloads/GSE70243_LK.CD34+.txt") %>% as.data.frame
rownames(expression3) = expression3$uid
expression3 = as.matrix(expression3[, 2:ncol(expression3)])

expression4 = read_tsv("~/Downloads/GSE70244_Lsk.txt") %>% as.data.frame
rownames(expression4) = expression4$uid
expression4 = as.matrix(expression4[, 2:ncol(expression4)])

expression = gdata::cbindX(expression1, expression2, expression3, expression4)

sampleinfo = data.frame(gate=c(rep("Gmp", ncol(expression1)), rep("Cmp", ncol(expression2)), rep("LK.CD34+", ncol(expression3)), rep("Lsk", ncol(expression4))), row.names=colnames(expression))

fig1b = read_tsv("~/Downloads/nature19348-f1_fig1b.csv")
clusters = fig1b[1, 3:ncol(fig1b)] %>% as.numeric %>% setNames(colnames(fig1b)[3:ncol(fig1b)])
names(clusters) = gsub(".*:(.*)", "\\1", names(clusters))
cellcluster_names = c("HSCP1", "HSCP2", "Meg", "Eryth", "Multi-Lin", "MDP", "Mono", "Gran", "Myelocyte")
clusters = clusters[rownames(sampleinfo)]
clusters = cellcluster_names[clusters] %>% factor(levels=cellcluster_names)
sampleinfo$cluster = clusters

expression = ExpressionSet(expression, phenoData = AnnotatedDataFrame(sampleinfo))

library(SCORPIUS)

expression_filtered = expression[,phenoData(expression)$cluster %in% c("Multi-Lin", "MDP", "Mono")]
plot(exprs(expression_filtered)["Irf8",], exprs(expression_filtered)["Zeb2",])
Goi = fig1b$UID[fig1b$`row_clusters-flat`[2:nrow(fig1b)] == 5]
expression_filtered = expression_filtered[,apply(expression_filtered[Goi,] %>% exprs, 2, mean) < 2]
#expression_filtered = expression_filtered[apply(expression_filtered, 1, sd) > 0,]

outliers = SCORPIUS::outlier.filter(correlation.distance(t(exprs(expression_filtered))))
expression_filtered = expression_filtered[,outliers]

space = reduce.dimensionality(correlation.distance(t(exprs(expression_filtered))), ndim = 2)

#space = tsne::tsne(as.dist(euclidean.distance(t(exprs(expression_filtered)))), k = 2)
#colnames(space) = paste0("Comp", 1:ncol(space))

traj = infer.trajectory(space)
draw.trajectory.plot(space, progression.group = phenoData(expression_filtered)$cluster, path = traj$path)
draw.trajectory.plot(space, progression.group = exprs(expression_filtered)["Klf4",], path = traj$path)+ viridis::scale_color_viridis(option="B")
draw.trajectory.plot(space, progression.group = exprs(expression_filtered)["Irf8",], path = traj$path)+ viridis::scale_color_viridis(option="B")
draw.trajectory.plot(space, progression.group = exprs(expression_filtered)["Cdk2",], path = traj$path)+ viridis::scale_color_viridis(option="B")

time = traj$time
time = -space[, "Comp1"]

phenoData(expression_filtered)$time = time

importances = SCORPIUS::gene.importances(t(exprs(expression_filtered)), time = time, ntree = 100)
draw.trajectory.heatmap(t(exprs(expression_filtered))[,importances %>% top_n(100, importance) %>% .$gene %>% as.character], time, show.labels.row = T)

Goi = c("Irf8", "Klf4")
Goi = c("Zeb2", "Klf4")
Goi = c("Irf8", "Zeb2")

plotdata = data.frame(time=time, bind_rows(lapply(Goi, function(g) {data.frame(gene=g, expression=exprs(expression_filtered)[g,])})), phenoData(expression_filtered) %>% as("data.frame"))
plot1 = plotdata %>% ggplot() + geom_smooth(aes(time, expression, group=gene, color=gene), span=0.7) + geom_point(aes(time, expression, group=gene, color=gene), alpha=0.3)
plotdata = data.frame(time=time, from=exprs(expression_filtered)[Goi[[1]],], to=exprs(expression_filtered)[Goi[[2]],], phenoData(expression_filtered) %>% as("data.frame"))
plotdata = plotdata %>% mutate(from_smooth=loess(from~time, data=plotdata, control=loess.control(surface="direct"), span=5) %>% predict(plotdata$time), to_smooth=loess(to~time, data=plotdata, control=loess.control(surface="direct"), span=5) %>% predict(plotdata$time))
plotdata = plotdata %>% mutate(from_smooth=plotdata %>% arrange(time) %$% smooth.spline(time, from, nknots=8) %>% predict(plotdata$time) %>% .$y, to_smooth=plotdata %>% arrange(time) %$% smooth.spline(time, to, nknots=8) %>% predict(plotdata$time) %>% .$y)
plot2 = plotdata %>% ggplot() + geom_point(aes(from, to, color=time)) + geom_path(aes(from_smooth, to_smooth), data=plotdata %>% arrange(time), arrow=arrow(), size=2) + viridis::scale_color_viridis(option="A") + xlab(Goi[[1]]) + ylab(Goi[[2]])
plot_grid(plotlist=list(plot1, plot2))

expression_smooth2 = lapply(featureNames(expression_filtered), function(g) {
  #smooth.spline(phenoData(expression_filtered)$time, exprs(expression_filtered)[g, ]) %>% predict(phenoData(expression_filtered)$time) %>% .$y
  loess(expression~time, data.frame(time=phenoData(expression_filtered)$time, expression=exprs(expression_filtered)[g, ]) ) %>% predict(phenoData(expression_filtered)$time)
})
expression_smooth = matrix(unlist(expression_smooth2), ncol=length(expression_smooth2[[1]]), byrow=T)
dimnames(expression_smooth) = dimnames(expression_filtered)
expression_smooth = ExpressionSet(expression_smooth, phenoData(expression_filtered))

plot(phenoData(expression_smooth)$time,exprs(expression_smooth)["Irf8",])
plot(phenoData(expression_filtered)$time,exprs(expression_filtered)["Spi1",])


g = "Irf5"

loess(expression~time, data.frame(time=phenoData(expression_filtered)$time, expression=exprs(expression_filtered)[g, ])) %>% predict(phenoData(expression_filtered)$time) %>% plot(phenoData(expression_filtered)$time, .)


# check delays
Goi = c("Irf8", "Klf4")
Goi = c("Irf8", "Csf1r")
Goi = c("Irf8", "Zeb2")
Goi = c("Irf8", "Csf1r")
Goi = c("Irf8", "Tgfbi")
Goi = c("Irf8", "Klf6")
Goi = c("Irf8", "Atf7")
Goi = c("Irf8", "Atf3")
Goi = c("Irf8", "Spi1")
Goi = c("Spi1", "Irf8")
Goi = c("Irf8", "Fce1g")
from = Goi[[1]]
to = Goi[[2]]

delays = seq(-0.9, 0.9, 0.01)
correlation_shifts = map_dbl(delays, function(delay) {
  cor(exprs(shift_expression(expression_smooth, delay))[from,],
      exprs(expression_smooth)[to,]
      )
})
plot(delays, correlation_shifts)
delay = delays[which.max(correlation_shifts)]
delay
plot(exprs(shift_expression(expression_smooth, delay))[from,],exprs(expression_smooth)[to,])
plot(exprs(shift_expression(expression_smooth, 0))[from,],exprs(expression_smooth)[to,])


##

expression_smooth_filtered = expression_smooth[apply(expression_smooth, 1, sd) > 0,]
gene_cors = cor(t(exprs(expression_smooth_filtered)))
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38825

library(GEOquery)
geo = getGEO("GSE38810")
gse = geo[[1]]
geo_expression = exprs(gse)
gene_symbols = fData(gse)$GENE_SYMBOL
gene_symbols[is.na(gene_symbols)] = ""
geo_expression = t(limma::avearrays(t(exprs(gse)), gene_symbols))
geo_expression["Ccr5",]

lfcs = (-log2(apply(geo_expression[,c("GSM949931", "GSM949932")], 1, mean)) + log2(apply(geo_expression[,c("GSM949933", "GSM949934")], 1, mean)))

lfcs %>% sort(T) %>% data.frame(lfc=.) %>% View
