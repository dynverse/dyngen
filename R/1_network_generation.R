# load module net
#' @importFrom readr read_tsv
#' @importFrom jsonlite read_json
load_modulenet <- function(modulenetname) {
  modulenodes <- read_tsv(paste0("data/modulenetworks/", modulenetname, "/modulenodes.tsv"), col_types=cols())
  if(!("celltype" %in% colnames(modulenodes))) {
    modulenodes$celltype <- 1
  }
  modulenet <- read_tsv(paste0("data/modulenetworks/", modulenetname, "/modulenet.tsv"), col_types=cols())
  celltypes <- read_tsv(paste0("data/modulenetworks/", modulenetname, "/celltypes.tsv"), col_types=cols())
  
  #statenet <- read_tsv(paste0("data/networks/", modulenetname, "/statenet.tsv"), col_types=cols())
  
  states <- jsonlite::read_json(paste0("data/modulenetworks/", modulenetname, "/states.json"), simplifyVector=TRUE)
  
  lst(modulenodes, modulenet, celltypes, modulenetname, states)
}

## add extra target genes to every tf
add_targets_individually <- function(net, geneinfo, tfs = geneinfo$gene[geneinfo$tf], genestart_id = max(geneinfo$gene), add_to_existing_net = TRUE) {
  allgenes <- tfs
  addnet <- tibble()
  for(tf in tfs) {
    nnewtargets <- sample(1:8, 1)
    if (nnewtargets > 10) {
      subnet <- dyngen::generate.ba.with.modules(nnewtargets, nnewtargets*2, 0.05, 0.05)$data.frame %>% rename(from=i, to=j) %>% mutate(from=paste0("G", from+genestart_id), to=to+genestart_id) %>% mutate(from=replace(from, from==paste0("G", genestart_id+1), tf))
    } else if (nnewtargets > 0) {
      subnet <- tibble(from=tf, to=(genestart_id+1):(genestart_id + nnewtargets + 1))
    } else {
      subnet <- tibble()
    }
    
    subnet <- subnet
    
    addnet <- bind_rows(addnet, subnet)
    genestart_id <- max(c(as.numeric(gsub("G([0-9]*)", "\\1", c(addnet$from, addnet$to))))+1, genestart_id)
  }
  
  addgeneinfo <- tibble(gene=sort(unique(union(addnet$from, addnet$to)))) %>% mutate(tf = gene %in% addnet$from, isgene=ifelse(gene %in% addnet$to, TRUE, FALSE), a0=NA, celltype=1)
  addnet <- addnet %>% mutate(effect=sample(c(-1, 1), length(to), TRUE), strength=1, cooperativity=1, randomize=TRUE)
  
  if(add_to_existing_net) {
    list(net=bind_rows(net, addnet), geneinfo=bind_rows(geneinfo, addgeneinfo) %>% group_by(gene) %>% filter(duplicated(gene) | n()==1))
  } else {
    list(net=addnet, geneinfo=addgeneinfo)
  }
  
}

## add extra target genes to every tf in a modular nature, where some target modules are regulated by multiple ldtf modules
add_targets_shared <- function(net, geneinfo, tfs = geneinfo$gene[geneinfo$tf], add_to_existing_net = TRUE) {
  allgenes <- geneinfo$gene
  addnet <- tibble()
  for(i in seq_along(tfs)) {
    nnewtargets <- sample(1:6, 1)
    newtargets <- (max(allgenes)+1):(max(allgenes) + nnewtargets + 1)
    
    subnet <- expand.grid(from=sample(tfs, size = sample(1:3, size=1, prob=c(0.4, 0.4, 0.2))), to=newtargets) %>% as.data.frame()
    
    subnet <- subnet %>% mutate(effect=0, strength=1, cooperativity=1)
    
    addnet <- bind_rows(addnet, subnet)
    allgenes <- sort(unique(union(subnet$from, subnet$to)))
  }
  
  addgeneinfo <- tibble(gene=sort(unique(union(addnet$from, addnet$to)))) %>% mutate(tf = gene %in% addnet$from, isgene=ifelse(gene %in% addnet$to, TRUE, FALSE), a0=NA, celltype=1, burn=TRUE)
  addnet <- addnet %>% mutate(effect=sample(c(-1, 1), length(to), TRUE), strength=1, cooperativity=1, randomize=TRUE)
  
  if(add_to_existing_net) {
    list(net=bind_rows(net, addnet), geneinfo=bind_rows(geneinfo, addgeneinfo) %>% group_by(gene) %>% filter(duplicated(gene) | n()==1))
  } else {
    list(net=addnet, geneinfo=addgeneinfo)
  }
}

# convert modulenet to modules
#' @param edge_retainment Function returning probabilities for the retainment of an edge (i ranging from 1 to `n`). Given `n` possible edges, what is the (relative) probability the i edges will be chosen? Choosing the default will retain on average half of the edges, taking into acount that always one edge should be chosen.
modulenet_to_modules <- function(
  modulenet, 
  modulenodes, 
  ngenes_per_module= function(n) sample(1:10, n, replace=TRUE), 
  gene_name_generator = function(i) paste0("G", i),
  edge_retainment = function(n) rep(1, n)
) {
  # generate tfs for each module, add to geneinfo
  geneinfo <- modulenodes %>% 
    mutate(
      ngenes = ngenes_per_module(n()),
      startgeneid = c(0, cumsum(ngenes)[-n()]),
      genes = map2(startgeneid, ngenes, ~gene_name_generator(seq(.x, .x+.y-1)))
    ) %>% 
    unnest(genes) %>% 
    rename(gene_id = genes) %>% 
    select(gene_id, module_id, a0, burn)
  
  # generate network between tfs
  # first generate the complete network for every from to every to
  completenet <- modulenet %>% 
    mutate(modulenet_edge_id=row_number()) %>% 
    left_join(geneinfo %>% select(gene_id, module_id) %>% rename(from_gene=gene_id), by=c("from"="module_id")) %>%
    left_join(geneinfo %>% select(gene_id, module_id) %>% rename(to_gene=gene_id), by=c("to"="module_id")) %>% 
    rename(from_module_id = from, to_module_id = to, from = from_gene, to = to_gene)
  
  # now subsample, making sure every "to" gene is regulated by at least one "from" gene for every module edge
  net <- completenet %>% 
    group_by(to, modulenet_edge_id) %>% 
    filter(row_number() %in% sample(seq_len(n()), sample(1:n(), prob=edge_retainment(n()))))
  
  lst(geneinfo, net)
}

# generate gene network
modulenet_to_genenet <- function(modulenet, modulenodes, genestart_id=0) {
  convertedmodulenet <- modulenet_to_modules(modulenet, modulenodes, genestart_id=genestart_id)
  net <- convertedmodulenet$net
  modulemembership <- convertedmodulenet$modulemembership
  
  allgenes <- tfs <- sort(unique(union(net$from, net$to)))
  net <- add_targets_shared(net, geneinfo=tibble(gene=allgenes, tf=TRUE))$net
  net$effect[net$effect == 0] <- sample(c(1, -1), sum(net$effect == 0), replace=T)
  net <- bind_rows(net, tibble(from=last(tfs), to=max(net$to)+1, effect=1, strength=1, cooperativity=2)) # add death marker
  
  allgenes <- sort(unique(union(net$from, net$to)))
  geneinfo <- tibble(gene=allgenes) %>% mutate(
    tf = gene %in% sort(unique(net$from)),
    isgene = TRUE
  )
  
  modulemembership <- map(modulemembership, ~intersect(., allgenes)) # remove genes which are not regulated by any other gene
  
  tibble::lst(net, geneinfo, modulemembership, allgenes)
}

# list genes
#' @importFrom stats setNames
list_genes <- function(geneinfo, modulemembership, net, modulenodes) {
  allgenes <- geneinfo$gene
  gene2module <- unlist(lapply(names(modulemembership), function(module_id) stats::setNames(rep(module_id, length(modulemembership[[module_id]])), modulemembership[[module_id]])))
  gene2module <- c( # assign genes to its first parent module
    gene2module, 
    net %>%
      group_by(to) %>% 
      summarize(firstf = first(from)) %>% 
      mutate(module=gene2module[as.character(firstf)]) %>% 
      select(-firstf) %>% dplyr::rename(gene=to) %>%
      filter(!(gene %in% names(gene2module))) %>% 
      {stats::setNames(.$module, .$gene)}
  ) 
  geneinfo$module <- gene2module[as.character(geneinfo$gene)]
  geneinfo <- geneinfo %>% bind_rows(tibble(gene=allgenes[!(allgenes %in% geneinfo$gene)])) %>% # add genes not in one of the modules
    mutate(tf = gene %in% net$from) %>% 
    left_join(modulenodes, by="module")
  
  geneinfo %<>% mutate(celltype=ifelse(is.na(celltype), 1, celltype))
  
  geneinfo
}



#' Plotting
#' 
#' @param model The model to plot
#' 
#' @importFrom igraph layout.graphopt graph_from_data_frame plot.igraph
plot_modulenet <- function(model) {
  graph <- igraph::graph_from_data_frame(model$modulenet, vertices = model$modulenodes)
  
  modulenames <- unique(model$geneinfo$module)
  colors <- rainbow(length(modulenames))
  V(graph)$color <- colors[match(names(V(graph)), modulenames)]
  
  layout <- igraph::layout.graphopt(graph, charge=0.01, niter=10000)
  #layout <- ForceAtlas2::layout.forceatlas2(graph, iterations=3000, plotstep=100)
  #png(file.path(imagefolder, "net_consecutive_bifurcating.png"), pointsize = 30, width=1000, height=1000)
  #igraph::plot.igraph(graph, edge.color = c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(-1,1, 0)))], layout=layout, vertex.size = 20, edge.arrow.size=0.5, edge.loop.angle=0.1, vertex.color=c("#222222", "#662222")[model$modulenodes$a0+1], vertex.label.color="white")
  
  igraph::plot.igraph(
    graph,
    edge.color = c("#d63737", "#3793d6", "green")[as.numeric(factor(model$modulenet$effect, levels = c(-1,1, 0)))], 
    layout = layout,
    vertex.size = 20, 
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1, 
    vertex.color = V(graph)$color, 
    vertex.label.color = "black"
  )
  #dev.off()
}

plot_net <- function(model) {
  graph <- igraph::graph_from_data_frame(model$net)
  layout <- igraph::layout_with_fr(graph)
  ldtf_filter <- as.numeric(factor(model$geneinfo$ldtf, levels = c(F, T)))
  igraph::plot.igraph(
    graph,
    edge.color = c("#d63737", "#3793d6", "#7cd637")[as.numeric(factor(model$net$effect, levels = c(-1,1, 0)))],
    layout = layout,
    vertex.size = c(1,5)[ldtf_filter], 
    edge.arrow.size = 0.5,
    edge.loop.angle = 0.1,
    vertex.color = c("black", "white")[ldtf_filter]
  )
}

#' @importFrom pheatmap pheatmap
plot_net_overlaps <- function(model) {
  jaccard <- function(x, y) {length(intersect(x, y))/length(union(x,y))}
  pheatmap::pheatmap(sapply(model$geneinfo$gene, function(i) sapply(model$geneinfo$gene, function(j) jaccard(model$net$from[model$net$to==i], model$net$from[model$net$to==j]))))
}

#' @importFrom igraph graph_from_data_frame V E plot.igraph
plot_net <- function(model) {
  goi <- unique(c(model$net$from, model$net$to))
  subnet <- model$net %>% filter(from %in% goi) %>% filter(to %in% goi)
  graph <- igraph::graph_from_data_frame(subnet, vertices=model$geneinfo$gene)
  
  modulenames <- model$modulenodes$module
  
  colors <- rainbow(length(modulenames))
  V(graph)$color <- colors[match(model$geneinfo$module[match(names(V(graph)), model$geneinfo$gene)], modulenames)]
  #V(graph)$label <- model$geneinfo$gene
  V(graph)$label <- ""
  E(graph)$color <- subnet$effect %>% factor(levels=c(1, -1)) %>% as.numeric() %>% c("blue", "red")[.]
  igraph::plot.igraph(graph, vertex.size=6, edge.arrow.size=0.5)
}

#' @importFrom igraph graph_from_data_frame V E plot.igraph
plot_net_tfs <- function(model) {
  goi <- unique(model$net$from)
  subnet <- model$net %>% filter(from %in% goi) %>% filter(to %in% goi)
  graph <- igraph::graph_from_data_frame(subnet)
  
  modulenames <- model$modulenodes$module
  
  colors <- rainbow(length(modulenames))
  V(graph)$color <- colors[match(model$geneinfo$module[match(names(V(graph)), model$geneinfo$gene)], modulenames)]
  E(graph)$color <- subnet$effect %>% factor(levels=c(1, -1)) %>% as.numeric() %>% c("blue", "red")[.]
  igraph::plot.igraph(graph)
}

