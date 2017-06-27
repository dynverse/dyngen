wnet = read_tsv("data/realnetworks/preproc/regulatory_circuits/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/01_neurons_fetal_brain.txt.gz", col_names = c("from", "to", "weight"), col_types=cols(col_character(), col_character(), col_double()))

net = wnet %>% top_n(250000, weight) %>% select(-weight)
adj = net %>% mutate(weight=1) %>% reshape2::acast(from~to, fill=0) %>% {.>0}

#net2 = tibble(from=c(1, 1, 1, 1, 2, 2, 2, 2, 2), to=c(3, 4, 2,5,  3, 2, 1, 4, 5))
#adj = net2 %>% mutate(weight=1) %>% reshape2::acast(from~to, fill=0) %>% {.>0}

#distances = dist(t(adj))
tfs = intersect(net$from, net$to)
distances = vegan::vegdist(t(adj[tfs, tfs]), "jaccard") # jaccard distance
distances[is.na(distances)] = 1
clusters = hclust(distances)
plot(clusters)

labels = clusters %>% cutree(h = 0.6)

labels2modules = function(labels, labelnames=names(labels)) {
  map(unique(labels), function(label) labelnames[labels == label])
}
allmodules = labels2modules(labels)
adj[tfs, allmodules[[10]],drop=F] %>% apply(1, sum) %>% hist

subnet = net %>% filter(to %in% tfs)
modules_sent = map(allmodules, ~filter(subnet, from %in% .) %>% .$to %>% unique)
modules_received = map(allmodules, ~filter(subnet, to %in% .) %>% .$from %>% unique)

jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
strength = function(x, z) {length(intersect(x, z))/length(z)}
sentstrengths = map(allmodules, function(module) (map(modules_sent, strength, z=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
receivestrengths = map(allmodules, function(module) (map(modules_received, strength, z=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
sentjaccards = map(allmodules, function(module) (map(modules_sent, jaccard, y=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))


mnet = (sentstrengths >= 1) %>% reshape2::melt(varnames=c("from", "to"), value.name="weight") %>% filter(weight) %>% filter(from!=to)
mgraph = mnet %>% graph_from_data_frame()# %>% as.undirected()
plot(mgraph)

mnet_wanted = tibble(from=1:4, to=2:5)
mnet_wanted = read_tsv(paste0("data/networks/consecutive_bifurcating//", "modulenet.tsv")) %>% filter(from!=to) %>% filter(!(to %in% c(10, 8, 9)))
mnet_wanted = model$modulenet

if(is.numeric(mnet_wanted$from)) {
  mnet_wanted$from = paste0("M", mnet_wanted$from)
  mnet_wanted$to = paste0("M", mnet_wanted$to)
}

for (perc in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) {
  for (i in 1:10) {
    print(perc)
    mgraph = mnet %>% sample_frac(perc) %>% graph_from_data_frame()# %>% as.undirected()
    mgraph_wanted = mnet_wanted %>% graph_from_data_frame() %>% simplify()
    
    results = subgraph_isomorphisms(mgraph_wanted, mgraph)
    
    if(length(results) > 0) {
      break
    }
  }
  if(length(results) > 0) {
    break
  }
}

edgedifferences = (map_int(results, function(modulevertices) {
  modulenamesoi = modulevertices %>% names %>% as.integer()
  subgraph = mnet %>% filter((from %in% modulenamesoi) & (to %in% modulenamesoi))
  nrow(subgraph)
}) - mnet_wanted %>% nrow)

modulenamesoi = results[[edgedifferences %>% which.min()]] %>% names %>% as.integer()
mnet_got = mnet %>% filter((from %in% modulenamesoi) & (to %in% modulenamesoi))
mgraph_got = mnet_got %>% graph_from_data_frame()
mgraph_got %>% plot

wanted_modulenames_mapper = set_names(V(mgraph_wanted) %>% names, modulenamesoi)
mnet_got$from_mapped = wanted_modulenames_mapper[mnet_got$from %>% as.character]
mnet_got$to_mapped = wanted_modulenames_mapper[mnet_got$to %>% as.character]

gene2module = map(seq_along(allmodules), ~(rep(., length(allmodules[[.]])))) %>% unlist() %>% set_names(unlist(allmodules))

genesoi = allmodules[modulenamesoi] %>% unlist
netoi = filter(net, (from %in% genesoi) & (to %in% genesoi)) %>% mutate(from_module=gene2module[from], to_module=gene2module[to])
netoi = mnet_got %>% left_join(mnet_wanted %>% mutate(present=T), by=c("from_mapped"="from", "to_mapped"="to")) %>% right_join(netoi, by=c("from"="from_module", "to"="to_module"), suffix=c("_module", "_gene")) %>% rename(from_module=from, to_module=to, from=from_gene, to=to_gene)

netoi = netoi %>% filter(!is.na(present))
graphoi = netoi[, c("from", "to")] %>% graph_from_data_frame()

genesoi_modulenames = gene2module[names(V(graphoi))] %>% as.character() %>% wanted_modulenames_mapper[.] %>% as.character %>% factor(levels=names(V(mgraph_wanted))) %>% set_names(names(V(graphoi)))

colors = rainbow(length(levels(genesoi_modulenames)))
V(graphoi)$color = colors[genesoi_modulenames %>% as.numeric()]
E(graphoi)$color = netoi$effect %>% factor(levels=c(1, -1)) %>% as.numeric() %>% c("blue", "red")[.]
plot(graphoi)

weight.community=function(row,membership,weigth.within=10,weight.between=1){
  if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
    weight=weigth.within
  }else{
    weight=weight.between
  }
  return(weight)
}
E(graphoi)$weight = apply(get.edgelist(graphoi), 1, weight.community, membership=genesoi_modulenames %>% set_names(names(V(graphoi))))
graphoi$layout = layout.fruchterman.reingold(graphoi, weights=E(graphoi)$weight)
plot(graphoi, vertex.size=6, arrow.size=0.1)


colors = rainbow(length(levels(genesoi_modulenames)))
V(mgraph_wanted)$color = colors
plot(mgraph_wanted)



notfs = setdiff(net$to, net$from)
tfsoi = unique(c(netoi$from, netoi$to))

netoi_targets = net %>% filter(to %in% notfs) %>% filter(from %in% tfsoi) %>% 
  mutate(effect=sample(c(1, -1), length(from), T), strength=1, cooperativity = 2)
netoi_targets = netoi_targets %>% sample_n(1000)
netoi_all = bind_rows(netoi, netoi_targets)
netoi_all %>% mutate(weight=1) %>% reshape2::acast(from~to, fill=0, value.var="weight") %>% pheatmap()
unique(c(netoi_all$to, netoi_all$from)) %>% length
