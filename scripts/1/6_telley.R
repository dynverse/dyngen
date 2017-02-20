Etelley = read_tsv("~/Downloads/telley_govindan_science2016/counts.tsv")
genes = Etelley$X1
Etelley = Etelley[,2:ncol(Etelley)] %>% as.matrix
rownames(Etelley) = genes

cell.info = tibble(cell=colnames(Etelley), hour=Etelley %>% colnames %>% gsub("X(.*)H_.*", "\\1", .) %>% as.numeric)

Etelley = Etelley[((apply(Etelley, 1, sd) %>% log) > 5) & ((apply(Etelley, 1, sd) %>% log) > 5), apply(Etelley, 2, sum) > 5e06]

cellcors = correlation.distance(t(Etelley))

#
space = reduce.dimensionality(cellcors, ndim = 1)
space = cbind(space, space)
colnames(space) = c("Comp1", "Comp2")

#
space = tsne::tsne(as.dist(cellcors), k = 2)
colnames(space) = paste0("Comp", 1:ncol(space))

##

seed = set.seed(1)
traj = infer.trajectory(space)
SCORPIUS::draw.trajectory.plot(space,cell.info[match(colnames(Etelley), cell.info$cell),]$hour %>% factor, path = traj$path)

library(org.Mm.eg.db)
symbol2entrez = map_chr(as.list(org.Mm.egSYMBOL2EG), ~.[[1]])
entrez2ensembl = map_chr(as.list(org.Mm.egENSEMBL), ~.[[1]])
symbol2ensembl = entrez2ensembl[symbol2entrez] %>% setNames(names(symbol2entrez)) %>% keep(~!is.na(.))
ensembl2symbol = setNames(names(symbol2ensembl), symbol2ensembl)

SCORPIUS::draw.trajectory.plot(space, Etelley[symbol2ensembl["Sox9"], ] + 1, path = traj$path) + viridis::scale_color_viridis(trans="log", option="A", direction=-1)
SCORPIUS::draw.trajectory.plot(space, apply(Etelley, 2, sum), path = traj$path) + viridis::scale_color_viridis(trans="log", option="A", direction=-1)

Goi = c("Foxp1")
SCORPIUS::draw.trajectory.heatmap(Etelley[symbol2ensembl[Goi],,drop=F] %>% t, traj$time, scale.features = T, show.labels.row = T)

gimp <- gene.importances(Etelley[(apply(Etelley, 1, sd) > quantile(apply(Etelley, 1, sd), 0.5)),] %>% t, traj$time, num.permutations = 0, ntree=100)

expr <- (t(Etelley[gimp$gene[1:500] %>% as.character,]) / apply(Etelley, 2, sum)) %>% quant.scale
colnames(expr) = ensembl2symbol[colnames(expr)]
modules <- extract.modules(expr)
SCORPIUS::draw.trajectory.heatmap(expr, traj$time, modules = modules, scale.features = F, show.labels.row = T)


Goi = c("Sox2", "Neurog2", "Tbr1", "Neurod2", "Tcf20")
Etelley[symbol2ensembl[Goi],,drop=F] %>% t %>% scale %>% t %>% reshape2::melt(varnames=c("gene", "cell"), value.name="expression") %>% mutate(gene=ensembl2symbol[gene %>% as.character]) %>% 
  mutate(time=traj$time[cell]) %>% 
  ggplot(aes(time, expression, color=gene)) + geom_smooth()
