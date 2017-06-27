wnet = read_tsv("data/realnetworks/regulatory_circuits/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/01_neurons_fetal_brain.txt.gz", col_names = c("from", "to", "weight"), col_types=cols(col_character(), col_character(), col_double()))

net = wnet %>% top_n(250000, weight) %>% select(-weight)
graph = net %>% graph_from_data_frame()

communities = walktrap.community(graph)
communities = edge.betweenness.community(graph)
communities = fastgreedy.community(graph %>% as.undirected())
communities = label.propagation.community(graph %>% as.undirected())
communities = infomap.community(graph)

adj = graph %>% as_adjacency_matrix()
adj = adj[communities$membership %>% order, communities$membership %>% order]

labels = communities$membership
labels = cut_at(communities, no=length(communities$names) / 10) %>% set_names(communities$names)

allmodules = list()
labels2modules = function(labels, labelnames=names(labels)) {
  map(unique(labels), function(label) labelnames[labels == label])
}

modules = list(names(V(graph)))
allmodules = list()

allmodules = c(modules %>% keep(~between(length(.), 2, 50)), allmodules)
print(length(allmodules))
nextmodules = modules %>% keep(~length(.) > 50)
print(length(nextmodules))
modules = map(nextmodules, ~split_net(filter(net, (from %in% .) & (to %in% .)))) %>% unlist(recursive=F)
print(length(modules))

split_net = function(net) {
  graph = net %>% graph_from_data_frame()
  communities = walktrap.community(graph)
  labels = cut_at(communities, no=length(communities$names) / 10) %>% set_names(communities$names)
  labels2modules(labels)
}

adj = graph %>% as_adjacency_matrix()
allmodules %>% keep(~length(.)>1) %>% map_dbl(~sum(adj[., .])/(length(.) * length(.) - length(.))) %>% plot


modules_sent = map(allmodules, ~filter(net, from %in% .) %>% .$to %>% unique)
modules_received = map(allmodules, ~filter(net, to %in% .) %>% .$from %>% unique)

jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
strength = function(x, z) {length(intersect(x, z))/length(z)}
sentstrengths = map(allmodules, function(module) (map(modules_sent, strength, z=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
receivestrengths = map(allmodules, function(module) (map(modules_received, strength, z=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
sentjaccards = map(allmodules, function(module) (map(modules_sent, jaccard, y=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))



mnet = sentstrengths %>% apply(1, function(row) {
  #order(row, decreasing=T)[1:5]
  which(row > 0.3)
}) %>% {map2(., seq_along(.), ~tibble(from=.x, to=.y))} %>% bind_rows() %>% filter(from!=to)

mnet_wanted = tibble(from=1:4, to=2:5)
mnet_wanted = read_tsv(paste0("data/networks/consecutive_bifurcating//", "modulenet.tsv")) %>% filter(from!=to) %>% filter(!(to %in% c(10, 8, 9)))

mgraph = mnet %>% graph_from_data_frame()# %>% as.undirected()
mgraph_wanted = mnet_wanted %>% graph_from_data_frame()# %>% as.undirected()
plot(mgraph_wanted)

# Subgraph isomorphism problem
# https://en.wikipedia.org/wiki/Subgraph_isomorphism_problem

results = subgraph_isomorphisms(mgraph_wanted, mgraph)

n = length(E(mgraph_wanted))
possible_deletions = lapply(seq_len(n-1), function(x) combn(seq_len(n), x, simplify=F)) %>% unlist(F)
possible_deletions = c(list(numeric()), possible_deletions)
possible_graphs = map(possible_deletions, ~delete.edges(mgraph_wanted, .))
for(graph in possible_graphs) {
  results = subgraph_isomorphisms(graph, mgraph)
  if(length(results) > 0) {
    break
  }
}

##

edgedifferences = (map_int(results, function(modulevertices) {
  modulenamesoi = modulevertices %>% names %>% as.integer()
  subgraph = mnet %>% filter((from %in% modulenamesoi) & (to %in% modulenamesoi))
  nrow(subgraph)
}) - mnet_wanted %>% nrow)

modulenamesoi = results[[edgedifferences %>% which.max()]] %>% names %>% as.integer()
mnet_got = mnet %>% filter((from %in% modulenamesoi) & (to %in% modulenamesoi))
subgraph = mnet_got %>% graph_from_data_frame()
subgraph %>% plot

wanted_modulenames_mapper = set_names(V(mgraph_wanted) %>% names %>% as.integer(), modulenamesoi)
mnet_got$from_mapped = wanted_modulenames_mapper[mnet_got$from %>% as.character]
mnet_got$to_mapped = wanted_modulenames_mapper[mnet_got$to %>% as.character]

genesoi = allmodules[modulenamesoi] %>% unlist
netoi = filter(net, (from %in% genesoi) & (to %in% genesoi))
graphoi = netoi %>% graph_from_data_frame()

gene2module = map(seq_along(allmodules), ~(rep(., length(allmodules[[.]])))) %>% unlist() %>% set_names(unlist(allmodules))

genesoi_modulenames = gene2module[names(V(graphoi))] %>% as.character() %>% factor(levels=names(V(subgraph)))

colors = rainbow(length(levels(genesoi_modulenames)))
V(graphoi)$color = colors[genesoi_modulenames]
plot(graphoi)

V(subgraph)$color = colors
plot(subgraph)


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
