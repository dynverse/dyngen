# this paper seems to do something similar: http://bioinformatics.oxfordjournals.org/content/26/18/2281.full.pdf
# but not really, as they shift the expression profiles of a gene vertically between time series (we want to align in the other dimension) => to cluster genes (modules)

library(dtw)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(pheatmap)

Emrna1 = readRDS("../1.rds")
Emrna2 = readRDS("../2.rds")

distances = dist(t(exprs(Emrna1)), t(exprs(Emrna2)))
result = dtw(distances)
plot(result)
pheatmap(distances, cluster_cols=F, cluster_rows=F, show_rownames = F, show_colnames = F)

dtwPlotThreeWay(result, exprs(Emrna1)["x_3",],  exprs(Emrna2)["x_3",])

# postprocessing to a more useful format
get_links = function(result) {
  idx1 = cumsum(result$stepsTaken %in% c(1, 3))+1
  idx2 = cumsum(result$stepsTaken %in% c(1, 2))+1
  
  matrix(c(idx1, idx2), nrow=2, byrow = T)
}

links = get_links(result)
pheatmap(links, cluster_cols=F, cluster_rows=F)

# visualization 
plotdata = bind_rows(
  data.frame(
    t(exprs(Emrna1))[links[1,],],
    time=phenoData(Emrna1)$time[links[1,]],
    sample=1,
    link=1:ncol(links)
  ),
  data.frame(
    t(exprs(Emrna2))[links[2,],],
    time=phenoData(Emrna2)$time[links[2,]],
    sample=2,
    link=1:ncol(links)
  )
)
plotdata$sample = factor(plotdata$sample)

linkplotdata = data.frame(
  linkleft=exprs(Emrna1)[,links[1,]] %>% t, 
  linkright=exprs(Emrna2)[, links[2,]] %>% t,
  timeleft=phenoData(Emrna1)$time[links[1,]],
  timeright=phenoData(Emrna2)$time[links[2,]],
  row.names=NULL # ignore rownames, otherwise he's gonna whine that there are duplicates
)

linkplotdata = data.frame(
  exprs(Emrna1)[,links[1,]] %>% t %>% melt(value.name="expression", varnames=c("cellid", "gene"))  %>% rename(expression1=expression, cellid1=cellid),
  exprs(Emrna2)[,links[2,]] %>% t %>% melt(value.name="expression", varnames=c("cellid", "gene")) %>% select(expression, cellid) %>% rename(expression2=expression, cellid2=cellid),
  time1=phenoData(Emrna1)$time[links[1,]],
  time2=phenoData(Emrna2)$time[links[2,]]
) %>% filter(gene %in% c("x_1", "x_2", "x_3", "x_4", "x_5"))

ggplot(linkplotdata) + 
  geom_line(aes(time1, expression1), color="red") +
  geom_line(aes(time2, expression2), color="blue") + 
  geom_segment(aes(time1, expression1, xend=time2, yend=expression2, color=time1>time2), alpha=0.5) +
  facet_wrap(~gene)

# combined expression
Emrna1_linked = Emrna1[,links[1,]]
featureNames(Emrna1_linked) = paste0(featureNames(Emrna1_linked), "_1")
Emrna2_linked = Emrna2[,links[2,]]
featureNames(Emrna2_linked) = paste0(featureNames(Emrna2_linked), "_2")

Elinked = rbind2(exprs(Emrna1_linked), exprs(Emrna2_linked))
pheatmap(Elinked, cluster_rows = T, cluster_cols = F, annotation_row = data.frame(sample=c(rep(1, nrow(Elinked)/2), rep(2, nrow(Elinked)/2)), row.names = rownames(Elinked)))

