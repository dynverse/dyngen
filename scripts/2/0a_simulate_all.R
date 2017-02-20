scriptfolder = "./scripts/3/"
source(file.path(scriptfolder, "1b_regulatorymodules.R"))
source(file.path(scriptfolder, "2_formulae.R"))
source(file.path(scriptfolder, "3_kinetics.R"))
source(file.path(scriptfolder, "4_simulate_functions.R"))

modulenetname = "linear"
#modulenetname = "cycle"
#modulenetname = "bifurcating"
#modulenetname = "bifurcating_cycle"
#modulenetname = "trifurcating"
modulenetname = "consecutive_bifurcating"
#modulenetname = "bifurcating_convergence"

# load module net
modulenodes = read_tsv(paste0("data/networks/", modulenetname, "_nodes.tsv"))
modulenet = read_tsv(paste0("data/networks/", modulenetname, ".tsv"))

# generate gene network
modulenet_to_modules(modulenet, modulenodes) %>% list2env(.GlobalEnv)
allgenes = ldtfs = sort(unique(union(net$from, net$to)))
#add_targets_individually(net) %>% list2env(.GlobalEnv)
add_targets_shared(net) %>% list2env(.GlobalEnv)

#net = bind_rows(net, tibble(from=last(ldtfs), to=max(net$to)+1, effect=1, strength=1, cooperativity=2)) # add death marker


G = allgenes = sort(unique(union(net$from, net$to)))
tfs = sort(unique(net$from))
length(tfs)/length(G)

## list genes
genes = unlist(lapply(1:length(modulemembership), function(moduleid) setNames(rep(moduleid, length(modulemembership[[moduleid]])), modulemembership[[moduleid]]))) %>% 
{tibble(gene=as.numeric(names(.)), module=.)}
genes = genes %>% bind_rows(tibble(gene=allgenes[!(allgenes %in% genes$gene)])) %>% # add genes not in one of the modules
  mutate(tf=gene %in% net$from) %>% 
  left_join(modulenodes, by="module")

## visualizing the networks
graph = graph_from_data_frame(modulenet, vertices = modulenodes)
layout <- layout_with_fr(graph)
plot.igraph(graph, edge.color = c("red", "blue", "green")[as.numeric(factor(modulenet$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = 20, edge.arrow.size=0.5, edge.loop.angle=0.1, vertex.color=modulenodes$a0)

graph = graph_from_data_frame(net)
layout <- layout_with_fr(graph)
ldtf_filter = as.numeric(factor(G %in% ldtfs, levels = c(F, T)))
plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "#7cd637")[as.numeric(factor(net$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = c(1,5)[ldtf_filter], edge.arrow.size=0.5, vertex.label=lapply(G, function(g) ifelse(g%in%ldtfs, "", "")), edge.loop.angle=0.1, vertex.color=c("black", "white")[ldtf_filter])

jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
pheatmap(sapply(G, function(i) sapply(G, function(j) jaccard(net$from[net$to==i], net$from[net$to==j]))))


generate_formulae(net) %>% list2env(.GlobalEnv)
generate_kinetics(vargroups, variables, nus.changes) %>% list2env(.GlobalEnv)

## determine start state genes
variables_genes = map(variables, "gene") %>% keep(~!is.null(.)) %>% unlist()
burngenes = (variables_genes %in% (genes %>% filter(!(module > 1) | is.na(module)) %>% .$gene)) %>% variables_genes[.]

newtime = 5

##
totaltime = 15
burntime = 2
newtime = estimate_convergence(8, verbose=T, totaltime=15) %>% max() %>% {.+1}

E = expression_multiple_cells(burntime, 8)
E = expression_one_cell(burntime, 15)
Emrna = E[str_detect(featureNames(E), "x_")]

pheatmap(SCORPIUS::quant.scale(exprs(Emrna[,phenoData(Emrna)$time %>% order]) %>% t) %>% t, cluster_cols = F, cluster_rows=T,scale="none")
exprs(E)["rg",] %>% plot
exprs(E)["x_120",] %>% plot

counts = Emrna %>% exprs %>% {.*100} %>% round() %>% abs()
Emrna2 = libprep(counts) %>% ExpressionSet(phenoData(Emrna))
pheatmap(exprs(Emrna2)[,phenoData(Emrna2)$time %>% order], cluster_cols = F, scale="none", cluster_row=F, show_rownames = F, show_colnames = F)

pheatmap(SCORPIUS::quant.scale(exprs(Emrna2[,phenoData(Emrna2)$time %>% order]), 0.05), cluster_cols = F, scale="none", cluster_row=T, show_rownames = F, show_colnames = F)
