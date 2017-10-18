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
  states <- jsonlite::read_json(paste0("data/modulenetworks/", modulenetname, "/states.json"), simplifyVector=TRUE)
  
  lst(modulenodes, modulenet, celltypes, states)
}

# convert modulenet to modules regulating each other
#' @param edge_retainment Function returning probabilities for the retainment of an edge (i ranging from 1 to `n`). Given `n` possible edges, what is the (relative) probability the i edges will be chosen? Choosing the default will retain on average half of the edges, taking into acount that always one edge should be chosen.
modulenet_to_genenet <- function(
  modulenet, 
  modulenodes, 
  ngenes_per_module= function(n) sample(1:10, n, replace=TRUE), 
  gene_name_generator = function(i) paste0("GM", i),
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
    select(gene_id, module_id, a0, burn, celltype) # select only relevant columns
  
  # generate network between tfs
  # first generate the complete network for every from to every to
  completenet <- modulenet %>% 
    mutate(modulenet_edge_id=row_number()) %>% 
    left_join(geneinfo %>% select(gene_id, module_id) %>% rename(from_gene=gene_id), by=c("from"="module_id")) %>%
    left_join(geneinfo %>% select(gene_id, module_id) %>% rename(to_gene=gene_id), by=c("to"="module_id")) %>% 
    rename(from_module_id = from, to_module_id = to, from = from_gene, to = to_gene) %>% 
    drop_na()
  
  # now subsample, making sure every "to" gene is regulated by at least one "from" gene for every module edge
  net <- completenet %>% 
    group_by(to, modulenet_edge_id) %>% 
    filter(row_number() %in% sample(seq_len(n()), sample(1:n(), prob=edge_retainment(n())))) %>% 
    ungroup() %>% 
    select(-from_module_id, -to_module_id, -modulenet_edge_id)
  
  lst(geneinfo, net)
}





## add extra target genes to every tf
add_targets_individually <- function(
  net, 
  geneinfo, 
  tfs = geneinfo$gene_id, 
  add_to_existing_net = TRUE,
  gene_name_generator = function(i) paste0("GE", i),
) {
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