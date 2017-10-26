#' load module net
#' @param modulenetname A modulenetname as defined in data/modulenetworks or `tree` to generate a random tree
#' @importFrom readr read_tsv
#' @importFrom jsonlite read_json
load_modulenet <- function(modulenetname) {
  modulenodes <- read_tsv(paste0("data/modulenetworks/", modulenetname, "/modulenodes.tsv"), col_types=cols(module_id=col_character()))
  modulenet <- read_tsv(paste0("data/modulenetworks/", modulenetname, "/modulenet.tsv"), col_types=cols(from=col_character(), to=col_character()))
  if (!("randomize" %in% colnames(modulenet))) {
    modulenet$randomize <- TRUE
  }
  
  # checks
  if (!all(modulenet$from %in% modulenodes$module_id) | !all(modulenet$to %in% modulenodes$module_id)) stop("Not all nodes in modulenodes")
  
  if(file.exists(paste0("data/modulenetworks/", modulenetname, "/cells.tsv"))) {
    cells <- read_tsv(paste0("data/modulenetworks/", modulenetname, "/cells.tsv"), col_types=cols())
  } else {
    cells <- tibble(cell_id = 1, dies=FALSE)
    modulenodes$cell_id <- 1
  }
  
  states <- jsonlite::read_json(paste0("data/modulenetworks/", modulenetname, "/states.json"), simplifyVector=TRUE)
  
  lst(modulenodes, modulenet, cells, states)
}

#' Convert modulenet to modules regulating each other
#' @param modulenet Module network
#' @param modulenodes Module nodes
#' @param ngenes_per_modules Functions for sampling the number of genes per module
#' @param gene_name_generate FUnction for generating the name of a gene
#' @param edge_retainment Function for sampling the number of edges retained between tfs of module nodes
modulenet_to_genenet <- function(
  modulenet, 
  modulenodes, 
  ngenes_per_module= function(n) sample(1:10, n, replace=TRUE), 
  gene_name_generator = function(i) paste0("GM", i),
  edge_retainment = function(n) max(c(round(n/2), 1))
) {
  # generate tfs for each module, add to geneinfo
  geneinfo <- modulenodes %>% 
    mutate(
      ngenes = ngenes_per_module(n()),
      startgeneid = c(0, cumsum(ngenes)[-n()]),
      genes = map2(startgeneid, ngenes, ~gene_name_generator(seq(.x, .x+.y-1))),
      isgene = TRUE,
      main = TRUE
    ) %>% 
    unnest(genes) %>% 
    rename(gene_id = genes) %>% 
    select(gene_id, module_id, a0, burn, cell_id, isgene, main) # select only relevant columns
  
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
    filter(row_number() %in% sample(seq_len(n()), edge_retainment(n()))) %>% 
    ungroup() %>% 
    select(-from_module_id, -to_module_id, -modulenet_edge_id)
  
  lst(geneinfo, net)
}

#' Add targets based on a real network
#' 
#' @param net Network dataframe
#' @param geneinfo Geneinfo dataframe
#' @param realnet_name Name of the real network
#' 
#' @inheritParams extract_induced_subgraph_from_tf
#' 
#' @importFrom glue glue
add_targets_realnet <- function(
  net, 
  geneinfo, 
  realnet_name = "regulatorycircuits", 
  damping=0.05, 
  ntargets_sampler=function() {sample(20:100, 1)},
  gene_name_generator = function(i) paste0("G", i)
) {
  # get the real network
  realnet <- read_csv(glue::glue("data/realnetworks/{realnet_name}.csv"), col_types=cols(from=col_character(), to=col_character()))
  allgenes <- unique(c(realnet$from, realnet$to))
  realgene2gene <- set_names(gene_name_generator(seq_along(allgenes)), allgenes)
  realnet$from <- realgene2gene[realnet$from]
  realnet$to <- realgene2gene[realnet$to]
  
  tfs <- realnet$from %>% unique
  
  # map existing genes to real tfs
  oldtf_to_newtf_mapper <- geneinfo$gene_id %>% set_names(sample(tfs, nrow(geneinfo)), .)
  
  geneinfo$gene_id <- oldtf_to_newtf_mapper[geneinfo$gene_id]
  net$from <- oldtf_to_newtf_mapper[net$from]
  net$to <- oldtf_to_newtf_mapper[net$to]
  
  # extract the small induced subgraphs
  geneinfo <- geneinfo %>% 
    rowwise() %>% 
    mutate(target_net=list(extract_induced_subgraph_from_tf(realnet, gene_id, damping=damping, ngenesampler=ntargets_sampler))) %>% 
    ungroup()
  geneinfo$target <- map(geneinfo$target_net, ~unique(c(.$from, .$to)))
  
  # create the geneinfo of the targets, based on the geneinfo of its (indirect) tfs
  # this can contain duplicates, so these are removed
  added_geneinfo <- geneinfo %>%
    unnest(target) %>% 
    select(-gene_id) %>% 
    rename(gene_id=target) %>% 
    group_by(gene_id) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    mutate(
      isgene=TRUE,
      main=FALSE
    )
  
  # also combine the subnetworks
  added_net <- geneinfo$target_net %>% bind_rows() %>% group_by(from, to) %>% filter(row_number() == 1) %>% ungroup()
  
  # combine
  geneinfo = bind_rows(geneinfo, added_geneinfo) %>% group_by(gene_id) %>% filter(row_number() == 1) %>% ungroup()
  net = bind_rows(net, added_net) %>% group_by(from, to) %>% filter(row_number() == 1) %>% ungroup()
  
  list(
  )
}


#' Extract subnetwork
#' Given a network, extract an induced subgraph from from a central TF to (indirect) targets using Personalized PageRank
#' 
#' @param net Network as a data.frame or igraph object
#' @param damping Damping factor for personalized pagerank
#' @param ngenesampler Function for sampling the number of genes
#' 
#' @return Network as a data.frame
extract_induced_subgraph_from_tf <- function(net, tf_of_interest=NULL, damping=0.05, ngenesampler=function() {sample(20:100, 1)}) {
  if(!("igraph" %in% class(net))) net <- igraph::graph_from_data_frame(net)
  tfs <- net %>% igraph::degree(mode="out") %>% {which(.>0)} %>% names
  
  personalized <- rep(0, length(igraph::V(net))) %>% set_names(names(igraph::V(net)))
  
  if (is.null(tf_of_interest)) tf_of_interest <- sample(tfs, 1)
  personalized[tf_of_interest] <- 1
  p <- igraph::page_rank(net, personalized=personalized, directed=TRUE, damping=damping)
  
  n_genes <- ngenesampler()
  genes_of_interest <- tibble(score=p$vector, gene_id=names(p$vector)) %>% 
    mutate(score = score + runif(n(), 0, 1e-15)) %>%  # avoid ties
    sample_n(n_genes, weight = score) %>% 
    pull(gene_id) %>% c(tf_of_interest) %>% unique()
  
  # remove genes not connected
  real_genes_of_interest <- igraph::ego(igraph::induced_subgraph(net, genes_of_interest), 10, tf_of_interest, mode="out") %>% first %>% names
  igraph::induced_subgraph(net, real_genes_of_interest) %>% igraph::as_data_frame()
}

#' Randomize network parameters
#' Randomize network parameters: cooperativity, strengths and effect
#' 
#' @param net Network dataframe
randomize_network_parameters <- function(net) {
  net %>% mutate(
    effect = ifelse(!is.na(effect), effect, sample(c(-1, 1), n(), TRUE)),
    cooperativity = ifelse(!is.na(cooperativity), cooperativity, runif(n(), 0.5, 2)),
    strength = ifelse(!is.na(strength), strength, 1)
  )
}