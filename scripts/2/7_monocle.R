counts = Efiltered %>% exprs %>% t
group.index = phenoData(Efiltered)$time %>% cut(5, labels=F)


verbose=F

#function(counts, group.index, num.genes, num.paths, verbose=F) {
  library(monocle)
  
  requireNamespace("dplyr")
  
  if (verbose) cat("Constructing data structures\n")
  counts.mat <- t(as.matrix(counts))
  pd <- new("AnnotatedDataFrame", data = data.frame(row.names=rownames(counts), num.genes=rowSums(counts!=min(counts))))
  fd <- new("AnnotatedDataFrame", data = data.frame(row.names=colnames(counts), times.expressed=colSums(counts!=min(counts))))
  cds <- newCellDataSet(counts.mat, phenoData = pd, featureData = fd)
  
  if (verbose) cat("Determining genes used for ordering\n")
  cds <- detectGenes(cds, min_expr = min(counts.mat)+1e-20)
  
  pct.cells.expressed <- fData(cds)$num_cells_expressed / nrow(counts)
  var <- apply(counts, 2, var)
  range <- apply(counts, 2, function(x) max(x) - min(x))
  score <- pct.cells.expressed * var * range
  
  ordering_genes <- names(score)[order(score, decreasing = T)[seq_len(min(num.genes, length(score)))]]
  
  # ordering_genes <- row.names(subset(fData(cds), num_cells_expressed >= num.cells.expressed))
  cds <- setOrderingFilter(cds, ordering_genes)
  
  if (verbose) cat("Reduce dimensionality\n")
  cds <- reduceDimension(cds)
  
  if (verbose) cat("Ordering the cells\n")
  cds <- orderCells(cds, num_paths=num.paths, reverse=F)
  
  space <- t(cds@reducedDimS)
  colnames(space) <- paste0("Comp", seq_len(ncol(space)))
  
  mst <- minSpanningTree(cds)
  diam <- as.data.frame(as.vector(V(mst)[get.diameter(mst, weights = NA)]$name))
  colnames(diam) <- c("sample_name")
  path <- space[diam$sample_name,]
  
  value <- dplyr::percent_rank(pData(cds)[,"Pseudotime"])
  names(value) <- rownames(counts)
  
  path <- path[order(value[rownames(path)]),]
  rownames(path) <- NULL
  
  list(space=space, path=path, value=value)
#}
  
  
monocle::plot_cell_trajectory(cds)
statemeantimes = tapply(phenoData(Efiltered)$time, phenoData(cds)$State, mean)
root_state = which.min(statemeantimes)
cds = orderCells(cds, root_state = root_state)

phenoData(cds)$time = phenoData(Efiltered)$time
plot_cell_trajectory(cds, color_by="time") + viridis::scale_color_viridis(option = "A", direction = -1, end=0.8) + theme_classic() + guides(color=FALSE)

plots = lapply(1:nrow(Emodules), function(moduleid) {
  phenoData(cds)$expression = apply(cds[modulemembership[[moduleid]],], 2, mean)
  plot_cell_trajectory(cds, color_by="expression") + viridis::scale_color_viridis(option="A", direction = -1) + theme(legend.position="none")
})
cowplot::plot_grid(plotlist=plots)
plots[[5]]

pheatmap(SCORPIUS::quant.scale(cds[,phenoData(cds)$time %>% order]), cluster_cols = F, scale="none", cluster_rows=F, labels_row = 1:nrow(cds))



plot(phenoData(Efiltered)$time, phenoData(cds)$Pseudotime)
