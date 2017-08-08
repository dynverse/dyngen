# genenetname = "tcells_shukla"
# genenetname = "myeloid_lymphoid_nakadai"
# genenetname = "bifurcating"
# genenetname = "mono_gran_olsson"
# genenetname = "B_plasma_singh"
# genenetname = "cell_cycle_arabidopsis"

load_circuit = function(genenetname) {
  geneinfo = read_tsv(file.path("data/genecircuits/", genenetname, "genes.tsv"))
  net = read_tsv(file.path("data/genecircuits/", genenetname, "net.tsv"))
  modulenodes = read_tsv(file.path("data/genecircuits/", genenetname, "modulenodes.tsv"))
  statenet = read_tsv(file.path("data/genecircuits/", genenetname, "statenet.tsv"))
  
  genemap = geneinfo$gene %>% {set_names(seq_along(.), .)}
  geneinfo$name = geneinfo$gene
  geneinfo$gene = genemap[geneinfo$gene]
  net$from = genemap[net$from]
  net$to = genemap[net$to]
  
  net = net %>% mutate(
    #strength=10^runif(length(strength)),
    strength=ifelse(is.na(strength), 10^runif(length(strength), 0, 2), strength),
    cooperativity = ifelse(is.na(cooperativity), 1.0123, cooperativity)
  )
  geneinfo = geneinfo %>% mutate(
    ldtf=T,
    #a0=ifelse(is.na(a0), runif(length(a0), 0, 1), a0)
    burn=as.logical(burn)
  )
  
  modulemembership = split(geneinfo$gene, geneinfo$module)
  
  celltypes = tibble(celltype=1, dies=F)
  
  model = tibble::lst(
    net, modulenodes, modulemembership, celltypes, statenet
  )
  
  model$geneinfo = list_genes(geneinfo, model$modulemembership, model$net, model$modulenodes)
  
  plot_net(model)
  
  model
}

