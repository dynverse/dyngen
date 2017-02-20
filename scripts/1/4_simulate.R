totaltime = 25
burntime = 2

##1: one cell at random time points
cell = simulate_cell(deterministic = T)
sampleids = sort(sample(length(cell$times), min(length(cell$times), 500)))
celltimes = cell$times[sampleids]
expression = cell$expression[sampleids,]
rownames(expression) = NULL

##2: multiple cells at random end points
celltimes = runif(50, 0, totaltime)
#cells = mclapply(celltimes, simulate_cell, mc.cores=8, deterministic=T)
cells = qsub.lapply(celltimes, function(celltime) {simulate_cell(celltime, deterministic=T)}, qsub.config = qsub.conf)
expression = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))

##3: combination
expressions = qsub.lapply(1:16, function(i) {
  cell = simulate_cell(deterministic = T)
  sampleids = sort(sample(length(cell$times), min(length(cell$times), 30)))
  celltimes = cell$times[sampleids]
  expression = cell$expression[sampleids,]
  rownames(expression) = NULL
  list(expression=expression, celltimes=celltimes, simulationid=i)
}, qsub.config=qsub.conf)

expression = do.call(rbind, map(expressions, ~.$expression))
celltimes = do.call(c, map(expressions, ~.$celltimes))

##
E = ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
Emrna = E[str_detect(featureNames(E), "x_")]
Eprot = E[str_detect(featureNames(E), "y_")]
Ebound = E[str_detect(featureNames(E), "b_")]

pheatmap(SCORPIUS::quant.scale(exprs(Emrna[,phenoData(Emrna)$time %>% order])), cluster_cols = F, scale="none")

Emrna2 = emdbook::rzinbinom(length(exprs(Emrna)), exprs(Emrna), 5, 0.2)
#Emrna2 = MASS::rnegbin(length(exprs(Emrna)), exprs(Emrna), 2)
Emrna2[is.na(Emrna2)] = 0
Emrna2 = matrix(Emrna2, nrow = nrow(Emrna), ncol=ncol(Emrna), dimnames = dimnames(Emrna))
Emrna2 = ExpressionSet(Emrna2, phenoData(Emrna), featureData(Emrna))

pheatmap(SCORPIUS::quant.scale(exprs(Emrna2[,phenoData(Emrna2)$time %>% order])), cluster_cols = F, scale="none")

saveRDS(Emrna2, file="results/linear.RDS")

Emrna2 = readRDS(file="results/linear.RDS")
