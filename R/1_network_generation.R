#' Model generation from modulenet
#' 
#' Default parameters are defined in individual functions and are inherited in `zzz_param_inheritance.R`
#' 
#' @param target_adder_name The method for adding targets. Only "realnet" is currently supported.
#' @param platform The platform
#' @inheritParams load_modulenet
#' @inheritParams generate_system
#' @inheritParams add_targets_realnet
#' @inheritParams modulenet_to_genenet 
#' @inheritParams generate_random_tree
#' @param ngenes_per_module_sampler Function for number of genes per module, given n_genes and n_modules
#' @param ntargets_sampler Function for number of targets per regulator, given n_genes and n_regulators
#' @param main_targets_ratio Ratio of targets versus main tfs
#' @param verbose Whether or not to produce textual output
#' 
#' @export
generate_model_from_modulenet <- function(
  # module net --
  modulenet_name, 
  
  # reference
  platform,
  
  # params for modulenet_name == "tree"
  treeseed,
  decay,
  
  # module net to net --
  ngenes_per_module_sampler,
  edge_retainment,
  ntargets_sampler,
  main_targets_ratio = 0.05,
  
  # add targets --
  target_adder_name, 
  
  # params for target_adder_name == "realnet_name"
  realnet_name,
  damping,
  
  # params for system
  samplers,
  
  verbose=TRUE
) {
  # load modulenet
  if (modulenet_name == "tree") {
    stagenet <- generate_random_tree(treeseed)
    model <- from_stages_to_modulenet(stagenet)
  } else {
    model <- load_modulenet(modulenet_name)
  }
  
  # convert module network to gene network between modules
  if (verbose) print("Generating main network")
  # number of genes in main based on number of genes in platform
  n_genes_total <- round(platform$n_genes * platform$pct_changing)
  ngenes_in_main <- round(n_genes_total * main_targets_ratio)
  ngenes_in_targets <- n_genes_total - ngenes_in_main
  
  model <- modulenet_to_genenet(
    model$modulenet, 
    model$modulenodes, 
    ngenes_in_main,
    ngenes_per_module_sampler,
    edge_retainment = edge_retainment
  ) %>% c(model)
  
  # add some targets
  if (verbose) print("Sampling targets")
  if (target_adder_name == "realnet") {
    model <- add_targets_realnet(
      model$net, model$geneinfo,
      realnet_name = realnet_name,
      damping = damping,
      ntargets = ngenes_in_targets,
      ntargets_sampler = ntargets_sampler
    ) %>% dynutils::merge_lists(model, .)
  } else {
    warning("Target sampler not supported")
  }
  
  # randomize parameters 
  if (verbose) print("Randomizing network")
  model$net <- randomize_network_parameters(model$net)
  
  # generate thermodynamics formulae & kinetics
  if (verbose) print("Generating system")
  model$system <- generate_system(model$net, model$geneinfo, model$cells, samplers)
  
  model
}




#' Load module network from extra data
#' @param modulenet_name A modulenet_name as defined in data/modulenetworks or `tree` to generate a random tree
#' @importFrom readr read_tsv
#' @importFrom jsonlite read_json
load_modulenet <- function(modulenet_name) {
  folder <- paste0(find.package("dyngen"), "/ext_data/modulenetworks/", modulenet_name)
  
  modulenodes <- read_tsv(paste0(folder, "/modulenodes.tsv"), col_types=cols(module_id=col_character()))
  modulenet <- read_tsv(paste0(folder, "/modulenet.tsv"), col_types=cols(from=col_character(), to=col_character()))
  if (!("randomize" %in% colnames(modulenet))) {
    modulenet$randomize <- TRUE
  }
  
  # checks
  if (!all(modulenet$from %in% modulenodes$module_id) | !all(modulenet$to %in% modulenodes$module_id)) stop("Not all nodes in modulenodes")
  
  if(file.exists(paste0(folder, "/cells.tsv"))) {
    cells <- read_tsv(paste0(folder, "/cells.tsv"), col_types=cols())
  } else {
    cells <- tibble(cell_id = 1, dies=FALSE)
    modulenodes$cell_id <- 1
  }
  
  edge_operations <-  read_tsv(paste0(folder, "/edge_operations.tsv"), col_types=cols())
  
  lst(modulenodes, modulenet, cells, edge_operations)
}

#' Convert modulenet to modules regulating each other
#' @param modulenet Module network
#' @param modulenodes Module nodes
#' @param n_genes Number of genes to use
#' @param ngenes_per_module_sampler Functions for sampling the number of genes per module
#' @param gene_name_generator Function for generating the name of a gene
#' @param edge_retainment Function for sampling the number of edges retained between tfs of module nodes
modulenet_to_genenet <- function(
  modulenet, 
  modulenodes, 
  n_genes,
  ngenes_per_module_sampler = function(n_genes, n_modules) sample(1:10, n_modules, replace=TRUE), 
  gene_name_generator = function(i) paste0("GM", i),
  edge_retainment = function(n) max(c(round(n/2), 1))
) {
  # generate tfs for each module, add to geneinfo
  geneinfo <- modulenodes %>% 
    mutate(
      ngenes = ngenes_per_module_sampler(max(n_genes, n()), n()),
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
#' @param realnet_name Name of the real network.
#' @param damping A daping factor for personalized pagerank
#' @param ntargets Number of total targets
#' @param ntargets_sampler A sample function for the number of targets
#' @param gene_name_generator A function for determining the name of a gene
#' 
#' @inheritParams extract_induced_subgraph_from_tf
#' 
#' @importFrom glue glue
add_targets_realnet <- function(
  net, 
  geneinfo, 
  realnet_name = "regulatorycircuits", 
  damping=0.05, 
  ntargets,
  ntargets_sampler=function(n_genes, n_regulators) {sample(20:100, 1)},
  gene_name_generator = function(i) paste0("G", i)
) {
  # get the real network
  realnet <- read_csv(glue::glue(find.package("dyngen"), "/ext_data/realnetworks/{realnet_name}.csv"), col_types=cols(from=col_character(), to=col_character()))
  realnet <- bind_rows(realnet, realnet %>% mutate(from=paste0("2_", from), to=paste0("2_", to)) %>% sample_frac(0.5))
  allgenes <- unique(c(realnet$from, realnet$to))
  realgene2gene <- set_names(gene_name_generator(seq_along(allgenes)), allgenes)
  realnet$from <- realgene2gene[realnet$from]
  realnet$to <- realgene2gene[realnet$to]
  
  tfs <- realnet$from %>% unique
  
  # map existing tfs in network to real tfs, by randomly sampling real tfs for each gene
  if (length(tfs) < nrow(geneinfo)) {
    print(length(tfs))
    print(nrow(geneinfo))
    stop("Not enough tfs in real network")
  }
  
  oldtf_to_newtf_mapper <- geneinfo$gene_id %>% set_names(sample(tfs, nrow(geneinfo)), .)
  
  geneinfo$gene_id <- oldtf_to_newtf_mapper[geneinfo$gene_id]
  net$from <- oldtf_to_newtf_mapper[net$from]
  net$to <- oldtf_to_newtf_mapper[net$to]
  
  geneinfo$ntargets <- ntargets_sampler(max(nrow(geneinfo), ntargets), nrow(geneinfo))
  
  # extract the small induced subgraphs
  geneinfo <- geneinfo %>% 
    rowwise() %>% 
    mutate(target_net=list(extract_induced_subgraph_from_tf(realnet, gene_id, damping=damping, ngenes=ntargets))) %>% 
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
      main=FALSE,
      a0 = NA, # a0 is decided later based on regulation
      burn=FALSE # extra genes should be available during burn in
    )
  
  # also combine the subnetworks
  added_net <- geneinfo$target_net %>% bind_rows() %>% group_by(from, to) %>% filter(row_number() == 1) %>% ungroup()
  
  # remove connections between main tfs (to avoid ruining the given module network)
  added_net <- added_net %>% filter(!(from %in% geneinfo$gene_id & to %in% geneinfo$gene_id))
  
  # combine
  geneinfo = bind_rows(geneinfo, added_geneinfo) %>% group_by(gene_id) %>% filter(row_number() == 1) %>% ungroup()
  net = bind_rows(net, added_net) %>% group_by(from, to) %>% filter(row_number() == 1) %>% ungroup()
  
  lst(
    geneinfo, net
  )
}


#' Extract subnetwork
#' Given a network, extract an induced subgraph from from a central TF to (indirect) targets using Personalized PageRank
#' 
#' @param net Network as a data.frame or igraph object
#' @param tf_of_interest Specify a TF of interest
#' @param damping Damping factor for personalized pagerank
#' @param ngenes The number of genes
#' 
#' @return Network as a data.frame
#' 
#' @importFrom stats runif
extract_induced_subgraph_from_tf <- function(net, tf_of_interest=NULL, damping=0.05, ngenes=1) {
  if(!("igraph" %in% class(net))) net <- igraph::graph_from_data_frame(net)
  tfs <- net %>% igraph::degree(mode="out") %>% {which(.>0)} %>% names
  
  personalized <- rep(0, length(igraph::V(net))) %>% set_names(names(igraph::V(net)))
  
  if (is.null(tf_of_interest)) tf_of_interest <- sample(tfs, 1)
  personalized[tf_of_interest] <- 1
  p <- igraph::page_rank(net, personalized=personalized, directed=TRUE, damping=damping)
  
  genes_of_interest <- tibble(score=p$vector, gene_id=names(p$vector)) %>% 
    mutate(score = score + stats::runif(n(), 0, 1e-15)) %>%  # avoid ties
    sample_n(ngenes, weight = score) %>% 
    pull(gene_id) %>% c(tf_of_interest) %>% unique()
  
  # remove genes not connected
  real_genes_of_interest <- igraph::ego(igraph::induced_subgraph(net, genes_of_interest), 10, tf_of_interest, mode="out") %>% first %>% names
  igraph::induced_subgraph(net, real_genes_of_interest) %>% igraph::as_data_frame()
}

#' Randomize network parameters
#' 
#' @param net Network dataframe
#' 
#' @importFrom stats runif
randomize_network_parameters <- function(net) {
  net %>% mutate(
    effect = ifelse(!is.na(effect), effect, sample(c(-1, 1), n(), TRUE)),
    cooperativity = ifelse(!is.na(cooperativity), cooperativity, stats::runif(n(), 0.5, 2)),
    strength = ifelse(!is.na(strength), strength, 1)
  )
}