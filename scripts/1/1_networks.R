library(fastgssa)
library(plyr)
library(tidyverse)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)
qsub.conf = qsub.configuration(exec.before = c("module unload python", "module load python"), memory = "10G")
library(readr)
library(Biobase)
library(igraph)
library(dyngen)
library(cowplot)

source("scripts/testje4_functions.R")

net = read_tsv("data/networks/linear.tsv")

ldtfs = sort(unique(c(net$from, net$to)))

### ONLY RUN IF YOU WANT EXTRA TARGET GENES
subtfs = c(last(ldtfs)+1:10)
net = bind_rows(net, data.frame(from=sample(ldtfs, length(subtfs), replace = T), to=subtfs, effect=1, strength=1, cooperativity=1))

tfs = c(ldtfs, subtfs)

targets = c(last(subtfs)+1:100)
net = bind_rows(net, data.frame(from=sample(tfs, length(targets), replace = T), to=targets, effect=1, strength=1, cooperativity=1))
###

G = sort(unique(union(net$from, net$to)))
tfs = sort(unique(net$from))


graph = graph_from_data_frame(net)
layout <- layout_with_fr(graph)
plot.igraph(graph, edge.color = c("red", "blue")[as.numeric(factor(net$effect, levels = c(-1,1)))], layout=layout, vertex.size = c(5,20)[as.numeric(factor(G %in% ldtfs, levels = c(F, T)))], edge.arrow.size=0.5, vertex.label=lapply(G, function(g) ifelse(g%in%ldtfs, g, "")), edge.loop.angle=0.1)

jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
pheatmap(sapply(G, function(i) sapply(G, function(j) jaccard(net$from[net$to==i], net$from[net$to==j]))))
