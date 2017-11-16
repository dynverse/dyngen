library(splatter)


data("sc_example_counts")
params <- splatEstimate(sc_example_counts)
sim <- splatSimulate(params, dropout.present = FALSE)



## preproc 10X
reference_expression = Matrix::readMM("data/realexpression/10x/pbmc8k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/GRCh38/matrix.mtx") %>% as.matrix()
write_rds(reference_expression, "data/realexpression/10x/pbmc8k_filtered_gene_bc_matrices/expression.rds", compress="gz")
reference_expression = read_rds("data/realexpression/10x/pbmc8k_filtered_gene_bc_matrices/expression.rds")
reference_expression = reference_expression[reference_expression %>% rowMeans() %>% order(decreasing = T) %>% .[1:2000],]
platform = list(
  params=splatEstimate(reference_expression), 
  originalGenes = nrow(reference_expression), 
  reference_dataset="10x/pbmc8k_filtered_gene_bc_matrices",
  name="10x_chromium_v2",
  id = "10x_chromium_v2_1"
)

## fluidigm
# warning from robrecht to wouter: do not use SCORPIUS::ginhoux anymore! See ?SCORPIUS::ginhoux for more details.
# I placed the ginhoux dataset in your data folder.
reference_expression = SCORPIUS::ginhoux$expression %>% t
reference_expression = reference_expression[reference_expression %>% rowMeans() %>% order(decreasing = T) %>% .[1:2000],]
platform = list(
  params=splatEstimate(reference_expression), 
  originalGenes = nrow(reference_expression), 
  reference_dataset="ginhoux",
  name="fluidigm_c1",
  id = "fluidigm_c1_ginhoux"
)
platform %>% save_platform()

##fluidigm (?) GSE66053
reference_expression = read_tsv("data/realexpression/preproc/GSE66053_singlecellseq_star_HTSeqCounts.tsv") %>% 
  dplyr::filter(substr(gene_id, 0, 2) != "__") %>% 
  dplyr::filter(substr(gene_id, 0, 4) != "ERCC")
reference_expression = reference_expression %>% reshape2::acast(sampleID~gene_id, value.var="counts") %>% t
reference_expression = reference_expression[reference_expression %>% rowMeans() %>% order(decreasing = T) %>% .[1:2000],]
reference_expression = reference_expression[, colSums(reference_expression) > 1000000]
platform = list(
  params=splatEstimate(reference_expression), 
  originalGenes = nrow(reference_expression), 
  reference_dataset="GSE66053",
  name="fluidigm_c1",
  id = "fluidigm_c1_GSE66053"
)
platform %>% save_platform()

## 


save_platform = function(platform) {
  write_rds(platform, paste0("data/platforms/", platform$id, ".rds"))
  read_rds("data/platforms/platforms.rds") %>% dplyr::filter(id != platform$id) %>% dplyr::bind_rows(platform %>% purrr::keep(is.character)) %>% write_rds("data/platforms/platforms.rds")
}

tibble() %>% write_rds("data/platforms/platforms.rds")



SCORPIUS::ginhoux$expression






rlnorm(1000, params@lib.loc, params@lib.scale) %>% hist



expression <- t(experiment$expression)
pheatmap::pheatmap(expression, cluster_cols=F, cluster_rows=T)


nGenes = estimate@nGenes = nrow(expression)
nCells = estimate@nCells = ncol(expression)
gene.names = rownames(expression)
cell.names = colnames(expression)
dummy.expression <- matrix(1, ncol = nCells, nrow = nGenes)
rownames(dummy.expression) <- gene.names
colnames(dummy.expression) <- cell.names
phenos <- new("AnnotatedDataFrame", data = data.frame(Cell = cell.names))
rownames(phenos) <- cell.names
features <- new("AnnotatedDataFrame", data = data.frame(Gene = gene.names))
rownames(features) <- gene.names
sim <- newSCESet(countData = dummy.expression, phenoData = phenos, featureData = features)
sim <- splatter:::splatSimLibSizes(sim, estimate)
sim <- splatter:::splatSimGeneMeans(sim, estimate)
sim <- splatter:::splatSimPathDE(sim, estimate)
sim <- splatter:::splatSimSingleCellMeans(sim, estimate)
sim <- splatter:::splatSimBCVMeans(sim, estimate)

phenoData(sim) %>% as("data.frame")
fData(sim) %>% as("data.frame") %>% str

libSizeScaler = nGenes
counts2 = lapply(seq_len(ncol(counts)), function(celli) {
  rmultinom(1, pData(sim)$ExpLibSize[[celli]]*libSizeScaler, abs(counts[, celli]))
}) %>% do.call(cbind, .)

# inspired by splatter:::splatSimDropout
drop.prob <- sapply(seq_len(ncol(counts2)), function(idx) {
  eta <- 
  return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
})
keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob), 
               nrow = nGenes, ncol = nCells)



pheatmap::pheatmap(counts2)
pheatmap::pheatmap(counts)


space = SCORPIUS::reduce_dimensionality(SCORPIUS::correlation_distance(t(counts2)))
SCORPIUS::draw_trajectory_plot(space)
space = SCORPIUS::reduce_dimensionality(SCORPIUS::correlation_distance(t(counts)))
SCORPIUS::draw_trajectory_plot(space)







expression



















function (params = newSplatParams(), method = c("single", "groups", 
                                                "paths"), verbose = TRUE, ...) 
{
  checkmate::assertClass(params, "SplatParams")
  method <- match.arg(method)
  if (verbose) {
    message("Getting parameters...")
  }
  params <- setParams(params, ...)
  params <- expandParams(params)
  validObject(params)
  seed <- getParam(params, "seed")
  set.seed(seed)
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nGroups <- getParam(params, "nGroups")
  group.cells <- getParam(params, "groupCells")
  if (nGroups == 1 && method == "groups") {
    warning("nGroups is 1, switching to single mode")
    method <- "single"
  }
  if (verbose) {
    message("Creating simulation object...")
  }
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  if (method == "groups") {
    group.names <- paste0("Group", seq_len(nGroups))
  }
  else if (method == "paths") {
    group.names <- paste0("Path", seq_len(nGroups))
  }
  dummy.counts <- matrix(1, ncol = nCells, nrow = nGenes)
  rownames(dummy.counts) <- gene.names
  colnames(dummy.counts) <- cell.names
  phenos <- new("AnnotatedDataFrame", data = data.frame(Cell = cell.names))
  rownames(phenos) <- cell.names
  features <- new("AnnotatedDataFrame", data = data.frame(Gene = gene.names))
  rownames(features) <- gene.names
  sim <- newSCESet(countData = dummy.counts, phenoData = phenos, 
                   featureData = features)
  if (method != "single") {
    groups <- lapply(seq_len(nGroups), function(i, g) {
      rep(i, g[i])
    }, g = group.cells)
    groups <- unlist(groups)
    pData(sim)$Group <- group.names[groups]
  }
  if (verbose) {
    message("Simulating library sizes...")
  }
  sim <- splatter:::splatSimLibSizes(sim, params)
  if (verbose) {
    message("Simulating gene means...")
  }
  sim <- splatter:::splatSimGeneMeans(sim, params)
  if (method == "single") {
    sim <- splatSimSingleCellMeans(sim, params)
  }
  else if (method == "groups") {
    if (verbose) {
      message("Simulating group DE...")
    }
    sim <- splatSimGroupDE(sim, params)
    if (verbose) {
      message("Simulating cell means...")
    }
    sim <- splatSimGroupCellMeans(sim, params)
  }
  else {
    if (verbose) {
      message("Simulating path endpoints...")
    }
    sim <- splatSimPathDE(sim, params)
    if (verbose) {
      message("Simulating path steps...")
    }
    sim <- splatSimPathCellMeans(sim, params)
  }
  if (verbose) {
    message("Simulating BCV...")
  }
  sim <- splatSimBCVMeans(sim, params)
  if (verbose) {
    message("Simulating counts..")
  }
  sim <- splatSimTrueCounts(sim, params)
  if (verbose) {
    message("Simulating dropout (if needed)...")
  }
  sim <- splatSimDropout(sim, params)
  if (verbose) {
    message("Creating final SCESet...")
  }
  sce <- newSCESet(countData = counts(sim), phenoData = new("AnnotatedDataFrame", 
                                                            data = pData(sim)), featureData = new("AnnotatedDataFrame", 
                                                                                                  data = fData(sim)))
  for (assay.name in names(assayData(sim))) {
    if (!(assay.name %in% names(assayData(sce)))) {
      set_exprs(sce, assay.name) <- get_exprs(sim, assay.name)
    }
  }
  if (verbose) {
    message("Done!")
  }
  return(sce)
}
