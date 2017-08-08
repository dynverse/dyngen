realnet = read_csv("data/realnetworks/regulatorycircuits.csv")


get_modulenetwork_from_network = function(realnet) {
  adj = realnet %>% mutate(weight=1) %>% reshape2::acast(from~to, fill=0) %>% {.>0}
  
  #net2 = tibble(from=c(1, 1, 1, 1, 2, 2, 2, 2, 2), to=c(3, 4, 2,5,  3, 2, 1, 4, 5))
  #adj = net2 %>% mutate(weight=1) %>% reshape2::acast(from~to, fill=0) %>% {.>0}
  
  #distances = dist(t(adj))
  tfs = intersect(realnet$from, realnet$to)
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
  
  subnet = realnet %>% filter(to %in% tfs)
  modules_sent = map(allmodules, ~filter(subnet, from %in% .) %>% .$to %>% unique)
  modules_received = map(allmodules, ~filter(subnet, to %in% .) %>% .$from %>% unique)
  
  jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
  strength = function(x, z) {length(intersect(x, z))/length(z)}
  sentstrengths = map(allmodules, function(module) (map(modules_sent, strength, z=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
  receivestrengths = map(allmodules, function(module) (map(modules_received, strength, z=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
  sentjaccards = map(allmodules, function(module) (map(modules_sent, jaccard, y=module) %>% unlist())) %>% unlist() %>% matrix(ncol=length(modules_sent))
  
  mnet = ((sentstrengths >= 1)) %>% reshape2::melt(varnames=c("from", "to"), value.name="weight") %>% filter(weight) %>% filter(from!=to)
  mgraph = mnet %>% graph_from_data_frame()# %>% as.undirected()
  
  tibble::lst(mnet, allmodules)
}

extract_network_from_modulenet = function(real_modulenet, mnet_wanted) {
  mnet = real_modulenet$mnet
  allmodules = real_modulenet$allmodules
  
  if(is.numeric(mnet_wanted$from)) {
    mnet_wanted$from = paste0("M", mnet_wanted$from)
    mnet_wanted$to = paste0("M", mnet_wanted$to)
  }
  
  for (perc in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) {
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
  genesoi_modulenames = gene2module[genesoi] %>% as.character() %>% wanted_modulenames_mapper[.] %>% as.character %>% factor(levels=names(V(mgraph_wanted))) %>% set_names(genesoi)
  netoi = filter(realnet, (from %in% genesoi) & (to %in% genesoi)) %>% mutate(from_module=gene2module[from], to_module=gene2module[to])
  netoi = mnet_wanted %>% filter(from == to) %>% left_join(tibble(gene=names(genesoi_modulenames), module=genesoi_modulenames), by=c("from"="module")) %>% mutate(from=gene, to=gene, from_module=gene2module[from], to_module=gene2module[to]) %>% select(from, to, from_module, to_module) %>% bind_rows(netoi) # add missing self links
  netoi %<>% mutate(from_module = wanted_modulenames_mapper[as.character(from_module)], to_module = wanted_modulenames_mapper[as.character(to_module)])
  netoi = mnet_wanted %>% mutate(present=T) %>% right_join(netoi, by=c("from"="from_module", "to"="to_module"), suffix=c("_module", "_gene")) %>% rename(from_module=from, to_module=to, from=from_gene, to=to_gene)
  
  netoi = netoi %>% filter(!is.na(present)) # remove edges not present in module network
  graphoi = netoi[, c("from", "to")] %>% graph_from_data_frame()
  
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
  
  
  notfs = setdiff(realnet$to, realnet$from)
  tfsoi = unique(c(netoi$from, netoi$to))
  
  netoi_targets = realnet %>% filter(to %in% notfs) %>% filter(from %in% tfsoi) %>% 
    mutate(effect=sample(c(1, -1), length(from), T), strength=1, cooperativity = 2)
  netoi_targets = netoi_targets %>% sample_n(100)
  netoi_all = bind_rows(netoi, netoi_targets)
  netoi_all %>% mutate(weight=1) %>% reshape2::acast(from~to, fill=0, value.var="weight") %>% pheatmap()
  unique(c(netoi_all$to, netoi_all$from)) %>% length
  
  netoi_all = netoi_all %>% group_by(from, to) %>% filter(row_number()==1) %>% ungroup() # remove duplicates
  
  process_extracted_net(netoi_all, tfsoi, genesoi_modulenames)
}


process_extracted_net = function(netoi_all, tfsoi, genesoi_modulenames) {
  realnet = netoi_all
  
  allgenes = sort(unique(union(realnet$from, realnet$to)))
  genemapper = paste0("G", seq_along(allgenes)) %>% set_names(allgenes)
  realnet$from = genemapper[realnet$from]
  realnet$to = genemapper[realnet$to]
  allgenes = sort(unique(union(realnet$from, realnet$to)))
  
  ldtfs = tfsoi
  geneinfo = tibble(gene=allgenes) %>% mutate(
    tf = gene %in% sort(unique(realnet$from)),
    ldtf = gene %in% paste0("G", ldtfs)
  )
  #geneinfo$celltype = 1
  
  modulemembership = map(levels(genesoi_modulenames), ~genemapper[names(genesoi_modulenames)[genesoi_modulenames == .]])
  modulemembership = map(modulemembership, ~intersect(., allgenes)) # remove genes which are not regulated by any other gene
  names(modulemembership) = levels(genesoi_modulenames)
  
  tibble::lst(net=realnet %>% select(from, to, effect, strength, cooperativity), geneinfo, modulemembership, allgenes)
}

real_modulenet = get_modulenetwork_from_network(realnet)
