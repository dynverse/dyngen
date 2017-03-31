# load module net
#' @import dplyr
#' @import readr
load_modulenet = function(modulenetname) {
  modulenodes = read_tsv(paste0("data/networks/", modulenetname, "/modulenodes.tsv"), col_types=cols()) %>% 
    mutate(module=paste0("M", module))
  if(!("celltype" %in% colnames(modulenodes))) {
    modulenodes$celltype = 1
  }
  modulenet = read_tsv(paste0("data/networks/", modulenetname, "/modulenet.tsv"), col_types=cols()) %>% 
    mutate(from=paste0("M", from), to=paste0("M", to))
  celltypes = read_tsv(paste0("data/networks/", modulenetname, "/celltypes.tsv"), col_types=cols())
  
  statenet = read_tsv(paste0("data/networks/", modulenetname, "/statenet.tsv"), col_types=cols())
  
  dambiutils::named_list(modulenodes, modulenet, celltypes, modulenetname, statenet)
}

## add extra target genes to every ldtf
#' @import dplyr
add_targets_individually = function(net) {
  allgenes = ldtfs = sort(unique(union(net$from, net$to)))
  for(ldtf in ldtfs) {
    nnewtargets = sample(1:2, 1)
    if (nnewtargets > 10) {
      subnet = dyngen::generate.ba.with.modules(nnewtargets, nnewtargets*2, 0.05, 0.05)$data.frame %>% rename(from=i, to=j) %>% mutate(effect=0, strength=1, cooperativity=1) %>% mutate(from=from+max(allgenes), to=to+max(allgenes)) %>% mutate(from=replace(from, from==max(allgenes)+1, ldtf))
    } else if (nnewtargets > 0) {
      subnet = tibble(from=ldtf, to=(max(allgenes)+1):(max(allgenes) + nnewtargets + 1), effect=0, strength=1, cooperativity=1)
    } else {
      subnet= tibble()
    }
    
    net = bind_rows(net, subnet)
    allgenes = sort(unique(union(net$from, net$to)))
  }
  
  net
}

## add extra target genes to every ldtf in a modular nature, where some target modules are regulated by multiple ldtf modules
#' @import dplyr
add_targets_shared = function(net) {
  allgenes = ldtfs = sort(unique(union(net$from, net$to)))
  for(i in 1:length(ldtfs)) {
    nnewtargets = sample(1:6, 1)
    newtargets = (max(allgenes)+1):(max(allgenes) + nnewtargets + 1)
    
    subnet = expand.grid(from=sample(ldtfs, size = sample(1:3, size=1, prob=c(0.4, 0.4, 0.2))), to=newtargets) %>% as.data.frame()
    
    subnet = subnet %>% mutate(effect=0, strength=1, cooperativity=1)
    
    net = bind_rows(net, subnet)
    allgenes = sort(unique(union(net$from, net$to)))
  }
  net
}

# convert modulenet to modules including lineage determining transcription factors (ldtfs)
#' @import dplyr
#' @import tibble
modulenet_to_modules = function(modulenet, modulenodes, ngenespermodule=4) {
  modulenames = c(modulenet$from, modulenet$to) %>% unique
  nmodules = ngenespermodule * (length(modulenames))
  modulemembership = split(seq_len(nmodules), ceiling(seq_len(nmodules)/ngenespermodule)) %>% set_names(modulenames)
  
  net = lapply(seq_len(nrow(modulenet)), function(i) {
    from_module = modulenet[i,]$from
    to_module = modulenet[i,]$to
    from = modulemembership[[from_module]]
    to = modulemembership[[to_module]]
    nedges = length(to) * 0.5
    
    effect = modulenet[i,]$effect
    strength = modulenet[i, ]$strength
    
    net = expand.grid(from=from, to=to) %>% as_tibble() %>% group_by(to) %>% sample_n(nedges) %>% mutate(effect=effect, strength=strength, cooperativity=2)
    #tibble(from=sample(from, nedges, T), to=sample(to, nedges, T), effect=modulenet[i,]$effect, strength=modulenet[i, ]$strength, cooperativity=2)
  }) %>% bind_rows
  
  named_list(modulemembership, net)
}

# generate gene network
#' @import dplyr
modulenet_to_genenet = function(modulenet, modulenodes) {
  convertedmodulenet = modulenet_to_modules(modulenet, modulenodes)
  net = convertedmodulenet$net
  modulemembership = convertedmodulenet$modulemembership
  
  allgenes = ldtfs = sort(unique(union(net$from, net$to)))
  net = add_targets_shared(net)
  net$effect[net$effect == 0] = sample(c(1, -1), sum(net$effect == 0), replace=T)
  net = bind_rows(net, tibble(from=last(ldtfs), to=max(net$to)+1, effect=1, strength=1, cooperativity=2)) # add death marker
  
  net$from = paste0("G", net$from)
  net$to = paste0("G", net$to)
  modulemembership = map(modulemembership, ~paste0("G", .))
  
  allgenes = sort(unique(union(net$from, net$to)))
  geneinfo = tibble(gene=allgenes) %>% mutate(
    tf = gene %in% sort(unique(net$from)),
    ldtf = gene %in% paste0("G", ldtfs)
  )
  
  modulemembership = map(modulemembership, ~intersect(., allgenes)) # remove genes which are not regulated by any other gene
  
  named_list(net, geneinfo, modulemembership, allgenes)
}

# list genes
#' @import dplyr
list_genes = function(geneinfo, modulemembership, net, modulenodes) {
  allgenes = sort(unique(union(net$from, net$to)))
  gene2module = unlist(lapply(names(modulemembership), function(moduleid) setNames(rep(moduleid, length(modulemembership[[moduleid]])), modulemembership[[moduleid]])))
  gene2module = c(gene2module, net %>% group_by(to) %>% summarize(firstf = first(from)) %>% mutate(module=gene2module[as.character(firstf)]) %>% select(-firstf) %>% dplyr::rename(gene=to) %>% filter(!(gene %in% names(gene2module))) %>% {setNames(.$module, .$gene)}) # assign genes to its first parent module
  geneinfo$module = gene2module[geneinfo$gene]
  geneinfo = geneinfo %>% bind_rows(tibble(gene=allgenes[!(allgenes %in% geneinfo$gene)])) %>% # add genes not in one of the modules
    mutate(tf=gene %in% net$from) %>% 
    left_join(modulenodes, by="module")
  geneinfo$a0 = ifelse(geneinfo$ldtf, geneinfo$a0, NA) # target genes will adapt the a0, not use the a0 of the module
  
  geneinfo
}



#' # Plotting
plot_modulenet = function(model) {
  graph = igraph::graph_from_data_frame(model$modulenet, vertices = model$modulenodes)
  
  modulenames = unique(model$geneinfo$module)
  colors = rainbow(length(modulenames))
  V(graph)$color = colors[match(model$geneinfo$module[match(names(V(graph)), model$geneinfo$gene)], modulenames)]
  
  layout <- igraph::layout_with_fr(graph)
  #png(file.path(imagefolder, "net_consecutive_bifurcating.png"), pointsize = 30, width=1000, height=1000)
  igraph::plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = 20, edge.arrow.size=0.5, edge.loop.angle=0.1, vertex.color=c("#222222", "#662222")[model$modulenodes$a0+1], vertex.label.color="white", vertex.label.size=20)
  #dev.off()
}

plot_net = function(model) {
  graph = igraph::graph_from_data_frame(model$net)
  layout <- igraph::layout_with_fr(graph)
  ldtf_filter = as.numeric(factor(model$geneinfo$ldtf, levels = c(F, T)))
  igraph::plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "#7cd637")[as.numeric(factor(model$net$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = c(1,5)[ldtf_filter], edge.arrow.size=0.5, edge.loop.angle=0.1, vertex.color=c("black", "white")[ldtf_filter])
}

plot_net_overlaps = function(model) {
  jaccard = function(x, y) {length(intersect(x, y))/length(union(x,y))}
  pheatmap::pheatmap(sapply(model$geneinfo$gene, function(i) sapply(model$geneinfo$gene, function(j) jaccard(model$net$from[model$net$to==i], model$net$from[model$net$to==j]))))
}

plot_net_tfs = function(model) {
  goi = unique(model$net$from)
  subnet = model$net %>% filter(from %in% goi) %>% filter(to %in% goi)
  graph = graph_from_data_frame(subnet)
  
  modulenames = model$modulenodes$module
  
  colors = rainbow(length(modulenames))
  V(graph)$color = colors[match(model$geneinfo$module[match(names(V(graph)), model$geneinfo$gene)], modulenames)]
  E(graph)$color = subnet$effect %>% factor(levels=c(1, -1)) %>% as.numeric() %>% c("blue", "red")[.]
  plot(graph)
}

