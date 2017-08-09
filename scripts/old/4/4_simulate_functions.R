simulate_cell = function(timeofsampling=NULL, deterministic=F, burngenes=c(), totaltime=10, burntime=2) {
  requireNamespace("fastgssa")
  requireNamespace("tidyr")
  requireNamespace("purrr")
  library(purrr)
  
  variables_burngenes = map(variables, "gene") %>% keep(~!is.null(.)) %>% unlist() %>% keep(~. %in% burngenes) %>% names
  formulae.nus.burn = formulae.nus
  formulae.nus.burn[setdiff(rownames(formulae.nus.burn), variables_burngenes),] = 0
  
  # burn in
  if (!deterministic) {
    out <- fastgssa::ssa(initial.state, formulae.strings, formulae.nus.burn, burntime, params, method=fastgssa::ssa.direct(), recalculate.all =F, stop.on.negstate = TRUE)
  } else {
    out <- fastgssa::ssa(initial.state, formulae.strings, formulae.nus.burn, burntime, params, method=fastgssa::ssa.em(), recalculate.all =F, stop.on.negstate = FALSE, stop.on.propensity=FALSE)
  }
  output = process_ssa(out)
  initial.state.burn = output$expression[nrow(output$expression), ] %>% abs
  
  # determine total time to simulate
  if (!is.null(timeofsampling)) {
    time = timeofsampling
  } else {
    time = totaltime
  }
  
  # actual simulation
  if (!deterministic) {
    out <- fastgssa::ssa(initial.state.burn, formulae.strings, formulae.nus, time, params, method=fastgssa::ssa.direct(), recalculate.all =F, stop.on.negstate = TRUE)
  } else {
    out <- fastgssa::ssa(initial.state.burn, formulae.strings, formulae.nus, time, params, method=fastgssa::ssa.em(), recalculate.all = FALSE, stop.on.negstate = FALSE, stop.on.propensity=FALSE)
  }
  output = process_ssa(out)
  expression = output$expression
  times = output$times
  
  # return either the full expression matrix, or the expression at timeofsampling
  if(is.null(timeofsampling)) {
    rownames(expression) = c(1:nrow(expression))
    return(list(expression=expression, times=times))
  } else {
    return(tail(expression, n=1)[1,])
  }
}

process_ssa <- function(out, starttime=0) {
  final.time <- out$args$final.time
  data <- out$timeseries
  lastrow <- tail(data, n=1)
  lastrow$t <- final.time
  data[nrow(data)+1,] <- lastrow
  data <- data[data$t<=final.time,]
  times <- data$t + starttime
  data <- as.matrix(data[,-1])
  
  return(list(times=times, expression=data))
}






# Emrna doesnt get updated, this function does not work
estimate_convergence = function(nruns=8, cutoff=0.15, verbose=F, totaltime=15) {
  convergences = mclapply(1:nruns, function(i) {
    E = expression_one_cell(burntime, totaltime)
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


expression_one_cell = function(burntime, totaltime, burngenes) {
  cell = simulate_cell(deterministic = T, burngenes=burngenes, burntime=burntime, totaltime=totaltime)
  sampleids = sort(sample(length(cell$times), min(length(cell$times), 500)))
  celltimes = cell$times[sampleids]
  expression = cell$expression[sampleids,]
  rownames(expression) = NULL
  
  ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
}


expression_multiple_cells = function(burntime, totaltime, burngenes, ncells=500) {
  celltimes = runif(ncells, 0, totaltime)
  #cells = mclapply(celltimes, simulate_cell, mc.cores=8, deterministic=T)
  cells = qsub_lapply(celltimes, function(celltime) {simulate_cell(celltime, deterministic=T, burngenes=burngenes, totaltime=totaltime, burntime=burntime)}, qsub_config = qsub_conf)
  expression = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
  
  ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
}


expression_multiple_cells_split = function(burntime, totaltime, burngenes, ncores=16, ncellspercore=30) {
  expressions = qsub_lapply(seq_len(ncores), function(i) {
    cell = simulate_cell(deterministic = T, burngenes=burngenes, totaltime=totaltime, burntime=burntime)
    sampleids = sort(sample(length(cell$times), min(length(cell$times), ncellspercore)))
    celltimes = cell$times[sampleids]
    expression = cell$expression[sampleids,]
    rownames(expression) = NULL
    list(expression=expression, celltimes=celltimes, simulationid=i)
  }, qsub_config=qsub_conf)
  expression = do.call(rbind, map(expressions, ~.$expression))
  celltimes = do.call(c, map(expressions, ~.$celltimes))
  
  ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
}

expression_multiple_cells_split_local = function(burntime, totaltime) {
  expressions = mclapply(1:16, function(i) {
    cell = simulate_cell(deterministic = T, burngenes=burngenes, totaltime=totaltime, burntime=burntime)
    sampleids = sort(sample(length(cell$times), min(length(cell$times), 30)))
    celltimes = cell$times[sampleids]
    expression = cell$expression[sampleids,]
    rownames(expression) = NULL
    list(expression=expression, celltimes=celltimes, simulationid=i)
  }, mc.cores = 8)
  expression = do.call(rbind, map(expressions, ~.$expression))
  celltimes = do.call(c, map(expressions, ~.$celltimes))
  
  ExpressionSet(t(expression), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
}

amplify = function(a=0, steps=100, rate=0.05) {
  for (i in 1:steps) {
    a = a + rbinom(length(a), a, rate)
  }
  a
}

libprep = function(counts, lysisrate = 0.6, capturerate = 0.1, amplify = T, amplifysteps = 100, amplifyrate = c(0.0001, 0.03), sequencerate=0.1, verbose=F, verbose_plot_cell=50, verbose_follow_gene="x_1") {
  curcounts = counts
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
  
  curcounts = counts_sequences = apply(curcounts, 2, function(row) rmultinom(1, sum(row)*sequencerate, row+0.000001)[,1]) # pseudocounts added because sum of probs cannot be zero
  
  
  if(verbose) {
    overview = tibble()
    overview = bind_rows(
      tibble(gene=rownames(counts), count = counts[, verbose_plot_cell], step="cell"),
      tibble(gene=rownames(counts), count = counts_lysed[, verbose_plot_cell], step="lysed"),
      tibble(gene=rownames(counts), count = counts_captured[, verbose_plot_cell], step="captured"),
      tibble(gene=rownames(counts), count = counts_amplified[, verbose_plot_cell], step="amplified"),
      tibble(gene=rownames(counts), count = counts_sequences[, verbose_plot_cell], step="sequenced")
    ) %>% mutate(step=factor(step, levels=unique(step)))
    plot = overview %>% ggplot() + geom_histogram(aes(count), bins=50) + facet_wrap(~step, ncol=1) + geom_vline(aes(xintercept=count, color=gene), data=overview %>% filter(gene %in% verbose_follow_gene))
    plot(plot)
  }
  curcounts
}