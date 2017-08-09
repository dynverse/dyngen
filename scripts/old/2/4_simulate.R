totaltime = 15
burntime = 2

estimate_convergence = function(nruns=8, cutoff=0.15, verbose=F) {
  convergences = mclapply(1:nruns, function(i) {
    E = expression_one_cell(burntime, 15)
    Emrna = E[str_detect(featureNames(E), "x_")]
    
    exprs(Emrna) %>% t %>% zoo::rollmean(10) %>% zoo::rollapply(100, sd, fill=0) %>% t %>% apply(2, quantile, probs=0.9)
  }, mc.cores=8) %>% do.call(rbind, .)
  
  timepoints = convergences %>% apply(1, function(convergences) {convergences %>% {.>=cutoff} %>% rev %>% cumsum %>% {. == 0} %>% rev %>% which() %>% first() %>% phenoData(Emrna)$time[.]}) # the first timepoint at which all following timepoints are below cutoff
  
  if (verbose) {
    plot = convergences %>% reshape2::melt(varnames=c("run", "time"), value.name="convergence") %>% mutate(run=factor(run)) %>% mutate(time=as(phenoData(Emrna), "data.frame")[,"time"][time]) %>% qplot(time, convergence, data=., group=run, color=run, geom="line") + geom_hline(aes(yintercept=cutoff)) + geom_vline(aes(xintercept = timepoint, color=run), data=tibble(timepoint=timepoints, run=factor(1:nruns)))
    print(plot)
  }
  timepoints
}
newtime = estimate_convergence(8, verbose=T) %>% max() %>% {.+1}

##1: one cell at random time points
expression_one_cell = function(burntime, totaltime) {
  cell = simulate_cell(deterministic = T, burngenes=burngenes, burntime=burntime, totaltime=totaltime)
  sampleids = sort(sample(length(cell$times), min(length(cell$times), 500)))
  celltimes = cell$times[sampleids]
  expression = cell$expression[sampleids,]
  rownames(expression) = NULL
  
  ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
}
E = expression_one_cell(2, newtime)

##2: multiple cells at random end points
expression_multiple_cells = function(burntime, totaltime) {
  celltimes = runif(500, 0, totaltime)
  #cells = mclapply(celltimes, simulate_cell, mc.cores=8, deterministic=T)
  cells = qsub_lapply(celltimes, function(celltime) {simulate_cell(celltime, deterministic=T, burngenes=burngenes)}, qsub_config = qsub_conf)
  expression = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
  
  ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
}
E = expression_multiple_cells(2, newtime)

##3: combination
expressions = qsub_lapply(1:16, function(i) {
  cell = simulate_cell(deterministic = T)
  sampleids = sort(sample(length(cell$times), min(length(cell$times), 30)))
  celltimes = cell$times[sampleids]
  expression = cell$expression[sampleids,]
  rownames(expression) = NULL
  list(expression=expression, celltimes=celltimes, simulationid=i)
}, qsub_config=qsub_conf)
expression = do.call(rbind, map(expressions, ~.$expression))
celltimes = do.call(c, map(expressions, ~.$celltimes))

##
Emrna = E[str_detect(featureNames(E), "x_")]
Eprot = E[str_detect(featureNames(E), "y_")]
Ebound = E[str_detect(featureNames(E), "b_")]
Erobrecht = E[str_detect(featureNames(E), "[xy]_")]

pheatmap(SCORPIUS::quant.scale(exprs(Emrna[,phenoData(Emrna)$time %>% order])), cluster_cols = F, cluster_rows=T,scale="none")
pheatmap(SCORPIUS::quant.scale(exprs(Eprot[,phenoData(Eprot)$time %>% order])), cluster_cols = F, scale="none", cluster_rows=F)
pheatmap(SCORPIUS::quant.scale(exprs(Erobrecht[rbind(seq(1, nrow(Erobrecht)/2), seq(nrow(Erobrecht)/2+1, nrow(Erobrecht))) %>% as.vector(),phenoData(Erobrecht)$time %>% order])), cluster_cols = F, scale="none", cluster_rows=F)

Emrna2 = emdbook::rzinbinom(length(exprs(Emrna)), exprs(Emrna), 5, 0.00)
#Emrna2 = MASS::rnegbin(length(exprs(Emrna)), exprs(Emrna), 2)
Emrna2[is.na(Emrna2)] = 0
Emrna2 = matrix(Emrna2, nrow = nrow(Emrna), ncol=ncol(Emrna), dimnames = dimnames(Emrna))
Emrna2 = ExpressionSet(Emrna2, phenoData(Emrna), featureData(Emrna))


## simulate sc RNA-seq
counts = Emrna %>% exprs %>% {.*100} %>% round() %>% abs()
amplify = function(a=0, steps=100, rate=0.05) {
  for (i in 1:steps) {
    a = a + rbinom(length(a), a, rate)
  }
  a
}
# result = lapply(c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75), function(n) {
#   map_dbl(rep(n, 1000), amplify) %>% tibble(a = ., n=n)
# }) %>% bind_rows()
# ggplot(result) + geom_density(aes(a, color=n, group=n))
libprep = function(curcounts, lysisrate = 0.6, capturerate = 0.1, amplify = T, amplifysteps = 100, amplifyrate = c(0.0001, 0.03), sequencerate=0.1) {
  curcounts = counts_lysed = apply(curcounts, 2, function(row) rmultinom(1, sum(row)*lysisrate, row)[,1])
  curcounts = counts_captured = apply(curcounts, 2, function(row) rmultinom(1, sum(row)*capturerate, row)[,1])
  
  # different rates per cell and per gene
  amplifyrates_cells = runif(ncol(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates_genes = runif(nrow(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates = sqrt(amplifyrates_cells %o% amplifyrates_genes) %>% t
  
  if(amplify) {
    curcounts = counts_amplified = lapply(1:ncol(curcounts), function(colid) {
      amplify(curcounts[, colid], amplifysteps, amplifyrates[,colid])
    }) %>% do.call(cbind, .)
  }
  
  curcounts = counts_sequences = apply(curcounts, 2, function(row) rmultinom(1, sum(row)*sequencerate, row)[,1])
  curcounts
}
Emrna2 = libprep(counts) %>% ExpressionSet(phenoData(Emrna))
pheatmap(exprs(Emrna2)[,phenoData(Emrna2)$time %>% order], cluster_cols = F, scale="none")













Emodules = lapply(modulemembership, function(module) apply(exprs(E)[intersect(featureNames(E), paste0("x_", module)),,drop=F], 2, mean)) %>% do.call(rbind, .)
pheatmap(SCORPIUS::quant.scale(t(Emodules[,phenoData(Emrna2)$time %>% order]), 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F, labels_row = 1:nrow(Emodules))

nodes = tibble(node=1:length(modules))
graph = graph_from_data_frame(modulenet, vertices = nodes)
layout <- layout_with_fr(graph)


expressionoi = Emodules[,138]
expressionoi = expressionoi/max(expressionoi)
V(graph)$color = viridis::viridis(100)[cut(expressionoi, seq(0, 1, length.out=100), labels=F)]
plot.igraph(graph,  edge.color = c("red", "blue")[as.numeric(factor(modulenet$effect, levels = c(-1,1)))], layout=layout, vertex.size = 20, edge.arrow.size=0.5, edge.loop.angle=0.1)

saveRDS(Emrna2, file="results/trifurcating.rds")
data = list(expression=Emrna2, net=net, modules=modules)
saveRDS(data, file="results/trifurcating_.rds")

data = readRDS(file="results/linear_.rds")
list2env(data, .GlobalEnv)


Emrna2 = readRDS(file="results/linear.rds")
