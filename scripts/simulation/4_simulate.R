simulate_cell = function(model, timeofsampling=NULL, deterministic=F, totaltime=10, burntime=2) {
  requireNamespace("fastgssa")
  requireNamespace("tidyr")
  requireNamespace("purrr")
  library(purrr)
  
  variables_burngenes = map(model$variables, "gene") %>% keep(~!is.null(.)) %>% unlist() %>% keep(~. %in% model$burngenes) %>% names
  formulae.nus.burn = model$formulae.nus
  formulae.nus.burn[setdiff(rownames(formulae.nus.burn), variables_burngenes),] = 0
  
  # burn in
  if (!deterministic) {
    out <- fastgssa::ssa(model$initial.state, model$formulae.strings, formulae.nus.burn, burntime, model$params, method=fastgssa::ssa.direct(), recalculate.all =F, stop.on.negstate = TRUE)
  } else {
    out <- fastgssa::ssa(model$initial.state, model$formulae.strings, formulae.nus.burn, burntime, model$params, method=fastgssa::ssa.em(), recalculate.all =F, stop.on.negstate = FALSE, stop.on.propensity=FALSE)
  }
  output = process_ssa(out)
  initial.state.burn = output$molecules[nrow(output$molecules), ] %>% abs
  
  # determine total time to simulate
  if (!is.null(timeofsampling)) {
    time = timeofsampling
  } else {
    time = totaltime
  }
  
  # actual simulation
  if (!deterministic) {
    out <- fastgssa::ssa(initial.state.burn, model$formulae.strings, model$formulae.nus, time, model$params, method=fastgssa::ssa.direct(), recalculate.all =F, stop.on.negstate = TRUE)
  } else {
    out <- fastgssa::ssa(initial.state.burn, model$formulae.strings, model$formulae.nus, time, model$params, method=fastgssa::ssa.em(), recalculate.all = FALSE, stop.on.negstate = FALSE, stop.on.propensity=FALSE)
  }
  output = process_ssa(out)
  molecules = output$molecules
  times = output$times
  
  # return either the full molecules matrix, or the molecules at timeofsampling
  if(is.null(timeofsampling)) {
    rownames(molecules) = c(1:nrow(molecules))
    return(list(molecules=molecules, times=times))
  } else {
    return(tail(molecules, n=1)[1,])
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
  
  return(list(times=times, molecules=data))
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


## Get the molecules matrix of multiple cells

simulate_one_cell = function(model, burntime, totaltime) {
  cell = simulate_cell(model, deterministic = T, burntime=burntime, totaltime=totaltime)
  sampleids = sort(sample(length(cell$times), min(length(cell$times), 500)))
  celltimes = cell$times[sampleids]
  molecules = cell$molecules[sampleids,]
  
  process_simulation(molecules, celltimes, 1)
}


simulate_multiple_cells = function(model, burntime, totaltime, ncells=500) {
  celltimes = runif(ncells, 0, totaltime)
  #cells = mclapply(celltimes, simulate_cell, mc.cores=8, deterministic=T)
  cells = qsub.lapply(celltimes, function(celltime) {simulate_cell(model, celltime, deterministic=T, totaltime=totaltime, burntime=burntime)}, qsub.config = qsub.conf)
  molecules = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
  
  process_simulation(molecules, celltimes, seq_len(nrow(molecules)))
}

simulate_multiple_cells_split = function(model, burntime, totaltime, nsimulations=16, ncellspersimulation=30, local=F) {
  if(!local) {
    multilapply = function(x, fun) {qsub.lapply(x, fun, qsub.config=qsub.conf)}
  } else {
    multilapply = function(x, fun) {mclapply(x, fun, mc.cores = 8)}
  }
  
  moleculess = multilapply(seq_len(nsimulations), function(i) {
    cell = simulate_cell(model, deterministic = T, totaltime=totaltime, burntime=burntime)
    sampleids = sort(sample(length(cell$times), min(length(cell$times), ncellspersimulation)))
    celltimes = cell$times[sampleids]
    molecules = cell$molecules[sampleids,]
    rownames(molecules) = NULL
    list(molecules=molecules, celltimes=celltimes, simulationids=rep(i, length(celltimes)))
  })
  
  process_simulation(
    do.call(rbind, map(moleculess, "molecules")), 
    do.call(c, map(moleculess, "celltimes")),
    do.call(c, map(moleculess, "simulationids"))
  )
}

process_simulation = function(molecules, celltimes, simulationids=1) {
  
  rownames(molecules) = paste0("C", seq_len(nrow(molecules)))
  cellinfo = tibble(cell=rownames(molecules), simulationtime=celltimes, simulationid=simulationids)
  
  molecules = molecules[order(cellinfo$simulationtime),]
  cellinfo = cellinfo[order(cellinfo$simulationtime),]
  expression = molecules[,str_detect(colnames(molecules), "x_")]
  colnames(expression) = gsub("x_(.*)", "\\1", colnames(expression))
  
  named.list(molecules, cellinfo, expression)
}

## scRNAseq functions
amplify = function(a=0, steps=100, rate=0.05) {
  for (i in 1:steps) {
    a = a + rbinom(length(a), a, rate)
  }
  a
}

libprep = function(counts, lysisrate = 0.6, capturerate = 0.1, amplify = T, amplifysteps = 100, amplifyrate = c(0.0001, 0.03), sequencerate=0.1, verbose=F, verbose_plot_cell=50, verbose_follow_gene="x_1") {
  curcounts = counts
  curcounts = counts_lysed = apply(curcounts, 1, function(col) rmultinom(1, sum(col)*lysisrate, col)[,1]) %>% t
  dimnames(curcounts) = dimnames(counts_lysed) = dimnames(counts)
  curcounts = counts_captured = apply(curcounts, 1, function(col) rmultinom(1, sum(col)*capturerate, col)[,1]) %>% t
  dimnames(curcounts) = dimnames(counts_captured) = dimnames(counts)
  
  # different rates per cell and per gene
  amplifyrates_cells = runif(nrow(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates_genes = runif(ncol(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates = sqrt(amplifyrates_cells %o% amplifyrates_genes)
  
  if(amplify) {
    curcounts = counts_amplified = lapply(1:nrow(curcounts), function(rowid) {
      amplify(curcounts[rowid, ], amplifysteps, amplifyrates[rowid, ])
    }) %>% do.call(rbind, .)
  }
  dimnames(curcounts) = dimnames(counts_amplified) = dimnames(counts)
  
  curcounts = counts_sequences = apply(curcounts, 1, function(col) rmultinom(1, sum(col)*sequencerate, col+0.000001)[,1]) %>% t # pseudocounts added because sum of probs cannot be zero
  dimnames(curcounts) = dimnames(counts_sequences) = dimnames(counts)
  
  if(verbose) {
    overview = tibble()
    overview = bind_rows(
      tibble(gene=colnames(counts), count = counts[verbose_plot_cell, ], step="cell"),
      tibble(gene=colnames(counts), count = counts_lysed[verbose_plot_cell, ], step="lysed"),
      tibble(gene=colnames(counts), count = counts_captured[verbose_plot_cell, ], step="captured"),
      tibble(gene=colnames(counts), count = counts_amplified[verbose_plot_cell, ], step="amplified"),
      tibble(gene=colnames(counts), count = counts_sequences[verbose_plot_cell, ], step="sequenced")
    ) %>% mutate(step=factor(step, levels=unique(step)))
    plot = overview %>% ggplot() + geom_histogram(aes(count), bins=50) + facet_wrap(~step, ncol=1) + geom_vline(aes(xintercept=count, color=gene), data=overview %>% filter(gene %in% verbose_follow_gene))
    plot(plot)
  }
  curcounts
}

simulate_scrnaseq = function(expression) {
  cellcounts = expression %>% {.*100} %>% round() %>% abs()
  libprep(cellcounts, amplifyrate = c(0.01, 0.03), verbose=T, verbose_plot_cell= "C1", verbose_follow_gene=c("G3","G10", "G50"))
}