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


simulate_multiple_cells = function(model, burntime, totaltime, ncells=500, qsub.conf=NULL) {
  celltimes = runif(ncells, 0, totaltime)
  #cells = mclapply(celltimes, simulate_cell, mc.cores=8, deterministic=T)
  cells = qsub.lapply(celltimes, function(celltime) {simulate_cell(model, celltime, deterministic=T, totaltime=totaltime, burntime=burntime)}, qsub.config = qsub.conf)
  molecules = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
  
  process_simulation(molecules, celltimes, seq_len(nrow(molecules)))
}

simulate_multiple_cells_split = function(model, burntime, totaltime, nsimulations=16, ncellspersimulation=30, local=F, qsub.conf=NULL) {
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