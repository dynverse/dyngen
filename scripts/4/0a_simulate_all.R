scriptfolder = "./scripts/4/"
source(file.path(scriptfolder, "1b_regulatorymodules.R"))
source(file.path(scriptfolder, "2_formulae.R"))
source(file.path(scriptfolder, "3_kinetics.R"))
source(file.path(scriptfolder, "4_simulate_functions.R"))

modulenetname = "linear";celltypes = tibble(celltype=c(1), dies=c(F))
modulenetname = "linear_intercell";celltypes = tibble(celltype=c(1, 2), dies=c(T, F))
modulenetname = "cycle";celltypes = tibble(celltype=c(1), dies=c(F))
modulenetname = "bifurcating";celltypes = tibble(celltype=c(1), dies=c(F))
#modulenetname = "bifurcating_cycle"
#modulenetname = "trifurcating"
modulenetname = "consecutive_bifurcating";celltypes = tibble(celltype=c(1), dies=c(F))
#modulenetname = "bifurcating_convergence"

# load module net
modulenodes = read_tsv(paste0("data/networks/", modulenetname, "_nodes.tsv"))
modulenodes$celltype = 1
modulenet = read_tsv(paste0("data/networks/", modulenetname, ".tsv"))

# generate gene network
modulenet_to_modules(modulenet, modulenodes) %>% list2env(.GlobalEnv)
allgenes = ldtfs = sort(unique(union(net$from, net$to)))
#add_targets_individually(net) %>% list2env(.GlobalEnv)
add_targets_shared(net) %>% list2env(.GlobalEnv)
net$effect[net$effect == 0] = sample(c(1, -1), sum(net$effect == 0), replace=T)

net = bind_rows(net, tibble(from=last(ldtfs), to=max(net$to)+1, effect=1, strength=1, cooperativity=2)) # add death marker


G = allgenes = sort(unique(union(net$from, net$to)))
tfs = sort(unique(net$from))
length(tfs)/length(G)

## list genes
gene2module = unlist(lapply(1:length(modulemembership), function(moduleid) setNames(rep(moduleid, length(modulemembership[[moduleid]])), modulemembership[[moduleid]])))
gene2module = c(gene2module, net %>% group_by(to) %>% summarize(firstf = first(from)) %>% mutate(module=gene2module[as.character(firstf)]) %>% select(-firstf) %>% dplyr::rename(gene=to) %>% filter(!(gene %in% names(gene2module))) %>% {setNames(.$module, .$gene)}) # assign genes to its first parent module
genes = gene2module %>% 
{tibble(gene=as.numeric(names(.)), module=.)}
genes = genes %>% bind_rows(tibble(gene=allgenes[!(allgenes %in% genes$gene)])) %>% # add genes not in one of the modules
  mutate(tf=gene %in% net$from) %>% 
  left_join(modulenodes, by="module")
genes$a0 = ifelse(genes$gene %in% ldtfs, genes$a0, NA) # target genes will adapt the a0, not use the a0 of the module

## visualizing the networks
graph = graph_from_data_frame(modulenet, vertices = modulenodes)
layout <- layout_with_fr(graph)
#png(file.path(imagefolder, "net_consecutive_bifurcating.png"), pointsize = 30, width=1000, height=1000)
plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "green")[as.numeric(factor(modulenet$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = 20, edge.arrow.size=0.5, edge.loop.angle=0.1, vertex.color=c("#222222", "#662222")[modulenodes$a0+1], vertex.label.color="white", vertex.label.size=20)
#dev.off()

graph = graph_from_data_frame(net)
layout <- layout_with_fr(graph)
ldtf_filter = as.numeric(factor(G %in% ldtfs, levels = c(F, T)))
plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "#7cd637")[as.numeric(factor(net$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = c(1,5)[ldtf_filter], edge.arrow.size=0.5, vertex.label=lapply(G, function(g) ifelse(g%in%ldtfs, "", "")), edge.loop.angle=0.1, vertex.color=c("black", "white")[ldtf_filter])

jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
pheatmap(sapply(G, function(i) sapply(G, function(j) jaccard(net$from[net$to==i], net$from[net$to==j]))))


generate_formulae(net, genes, celltypes) %>% list2env(.GlobalEnv)
generate_kinetics(vargroups, variables, nus.changes) %>% list2env(.GlobalEnv)

## determine start state genes
variables_genes = map(variables, "gene") %>% keep(~!is.null(.)) %>% unlist()
modulesoi = modulenodes %>% filter(celltype == 2) %>% .$module
modulesoi = 1
burngenes = (variables_genes %in% (genes %>% filter(module %in% modulesoi | is.na(module)) %>% .$gene)) %>% variables_genes[.]

newtime = 5

##
totaltime = 8
burntime = 2
newtime = estimate_convergence(8, verbose=T, totaltime=15) %>% max() %>% {.+1}

E = expression_multiple_cells(burntime, totaltime, 1000)
E = expression_one_cell(burntime, totaltime)
E = expression_multiple_cells_split_local(burntime, totaltime)

E = E[,!(exprs(E) %>% apply(2, function(x) any(is.na(x))))]
E = E[,phenoData(E)$time %>% order]
Emrna = E[str_detect(featureNames(E), "x_")]

pheatmap(SCORPIUS::quant.scale(exprs(Emrna[,phenoData(Emrna)$time %>% order]) %>% t) %>% t, cluster_cols = F, cluster_rows=T,scale="none")
exprs(E)["rg_1",] %>% plot
exprs(E)[paste0("x_", (net$to %>% last)),] %>% plot

counts = Emrna %>% exprs %>% {.*100} %>% round() %>% abs()
Emrna2 = libprep(counts, amplifyrate = c(0.01, 0.03), verbose=T, verbose_plot_cell= 10, verbose_follow_gene=c("x_3","x_10", "x_50"))
colnames(Emrna2) = sampleNames(Emrna)
Emrna2 = Emrna2 %>% ExpressionSet(phenoData(Emrna))
pheatmap(exprs(Emrna2)[,phenoData(Emrna2)$time %>% order], cluster_cols = F, scale="none", cluster_row=F, show_rownames = F, show_colnames = F)

pheatmap(SCORPIUS::quant.scale(exprs(Emrna2[,phenoData(Emrna2)$time %>% order]) %>% t, 0.05) %>% t, cluster_cols = F, scale="none", cluster_row=T, show_rownames = F, show_colnames = F)


Emodules = lapply(modulemembership, function(module) apply(exprs(E)[intersect(featureNames(E), paste0("x_", module)),,drop=F], 2, mean)) %>% do.call(rbind, .)
pheatmap(SCORPIUS::quant.scale(t(Emodules[,phenoData(Emrna)$time %>% order]), 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, labels_row = 1:nrow(Emodules))
base = 2
basemax = 4
cellstateinfo = Emodules %>% apply(2, function(x) {
  state = (1:length(x))[x>base] %>% last
  if (is.na(state)) {
    warning("no good state")
    state = 1
  }
  progression = pmin((x[state] - base)/(basemax-base), 1)
  
  tibble(state=state, progression=progression)
}) %>% bind_rows
phenoData(E)$state = cellstateinfo$state
phenoData(E)$progression = cellstateinfo$progression
ggplot(phenoData(E) %>% as("data.frame")) + geom_area(aes(time, group=state, fill=factor(state)), position="fill", stat="bin", bins=30)




exprs(E)[vargroups$rg,] %>% reshape2::melt(varnames=c("variable", "cell"), value.name="expression") %>% ggplot() + geom_line(aes(cell, expression, group=variable, color=variable))
