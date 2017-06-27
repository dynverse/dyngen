#' @import dplyr
#' @import purrr
#' @import fastgssa
simulate_cell = function(model, timeofsampling=NULL, deterministic=F, totaltime=10, burntime=2, ssa.algorithm = fastgssa::ssa.em(noise_strength=4)) {
  variables_burngenes = map(model$variables, "gene") %>% keep(~!is.null(.)) %>% unlist() %>% keep(~. %in% model$burngenes) %>% names
  variables_burngenes = c(variables_burngenes, model$vargroups$rg)
  formulae.nus.burn = model$formulae.nus
  formulae.nus.burn[setdiff(rownames(formulae.nus.burn), variables_burngenes),] = 0
  
  # burn in
  if (!deterministic) {
    out <- fastgssa::ssa(model$initial.state, model$formulae.strings, formulae.nus.burn, burntime, model$params, method=fastgssa::ssa.direct(), recalculate.all = TRUE, stop.on.negstate = TRUE)
  } else {
    out <- fastgssa::ssa(model$initial.state, model$formulae.strings, formulae.nus.burn, burntime, model$params, method=ssa.algorithm, recalculate.all = TRUE, stop.on.negstate = FALSE, stop.on.propensity=FALSE)
  }
  output = process_ssa(out)
  initial.state.burn = output$molecules[nrow(output$molecules), names(model$initial.state)]
  
  # determine total time to simulate
  if (!is.null(timeofsampling)) {
    time = timeofsampling
  } else {
    time = totaltime
  }
  
  # actual simulation
  if (!deterministic) {
    out <- fastgssa::ssa(initial.state.burn, model$formulae.strings, model$formulae.nus, time, model$params, method=fastgssa::ssa.direct(), recalculate.all = TRUE, stop.on.negstate = TRUE)
  } else {
    out <- fastgssa::ssa(initial.state.burn, model$formulae.strings, model$formulae.nus, time, model$params, method=ssa.algorithm, recalculate.all = TRUE, stop.on.negstate = FALSE, stop.on.propensity=FALSE)
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

#' @import dplyr
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
#' @import dplyr
estimate_convergence = function(nruns=8, cutoff=0.15, verbose=F, totaltime=15) {
  convergences = parallel::mclapply(1:nruns, function(i) {
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
#' @import dplyr
simulate_one_cell = function(model, burntime, totaltime, ncells=500, ssa.algorithm = fastgssa::ssa.em(noise_strength=4)) {
  cell = simulate_cell(model, deterministic = T, burntime=burntime, totaltime=totaltime, ssa.algorithm=ssa.algorithm)
  if(!is.null(ncells)) {
    sampleids = sort(sample(length(cell$times), min(length(cell$times), ncells)))
  } else {
    sampleids = seq_len(length(cell$times))
  }
  celltimes = cell$times[sampleids]
  molecules = cell$molecules[sampleids,]
  
  process_simulation(molecules, celltimes, 1, list(molecules))
}

#' @import dplyr
simulate_multiple_cells = function(model, burntime, totaltime, ncells=500, qsub_conf=NULL, ssa.algorithm = fastgssa::ssa.em(noise_strength=4)) {
  celltimes = runif(ncells, 0, totaltime)
  #cells = mclapply(celltimes, simulate_cell, mc.cores=8, deterministic=T)
  cells = PRISM::qsub_lapply(celltimes, function(celltime) {simulate_cell(model, celltime, deterministic=T, totaltime=totaltime, burntime=burntime, ssa.algorithm=ssa.algorithm)}, qsub_config = qsub_conf)
  molecules = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
  
  process_simulation(molecules, celltimes, seq_len(nrow(molecules)))
}

#' @import dplyr
#' @import purrr
simulate_multiple_cells_split = function(model, burntime, totaltime, nsimulations=16, ncellspersimulation=30, endonly=F, local=F, qsub_conf=NULL, ssa.algorithm = fastgssa::ssa.em(noise_strength=4)) {
  force(model) # force the evaluation of the model argument, as the qsub environment will be empty except for existing function arguments
  if(!local) {
    multilapply = function(x, fun) {PRISM::qsub_lapply(x, fun, qsub_config=qsub_conf, qsub_environment = formalArgs(simulate_multiple_cells_split))}
  } else {
    multilapply = function(x, fun) {parallel::mclapply(x, fun)}
  }
  
  moleculess = multilapply(seq_len(nsimulations), function(i) {
    cell = simulate_cell(model, deterministic = T, totaltime=totaltime, burntime=burntime, ssa.algorithm=ssa.algorithm)
    
    if(!endonly) {
      sampleids = sort(sample(length(cell$times), min(length(cell$times), ncellspersimulation)))
    } else {
      sampleids = length(cell$times) # only the last cell
    }
    celltimes = cell$times[sampleids]
    molecules = cell$molecules[sampleids,,drop=F]
    rownames(molecules) = NULL
    list(molecules=molecules, celltimes=celltimes, simulationids=rep(i, length(celltimes)), simulation=cell$molecules)
  })
  
  process_simulation(
    do.call(rbind, map(moleculess, "molecules")), 
    do.call(c, map(moleculess, "celltimes")),
    do.call(c, map(moleculess, "simulationids")),
    map(moleculess, "simulation")
  )
}

#' @import dplyr
#' @import ggplot2
#' @import purrr
process_simulation = function(molecules, celltimes, simulationids=1, simulations=NULL) {
  rownames(molecules) = paste0("C", seq_len(nrow(molecules)))
  cellinfo = tibble(cell=rownames(molecules), simulationtime=celltimes, simulationid=simulationids)
  
  molecules = molecules[order(cellinfo$simulationtime),]
  cellinfo = cellinfo[order(cellinfo$simulationtime),]
  expression = molecules[,str_detect(colnames(molecules), "x_")]
  expression[expression < 0] = 0
  expression = round(expression * 100)
  colnames(expression) = gsub("x_(.*)", "\\1", colnames(expression))
  
  simulations = map(simulations, function(molecules) {
    expression = molecules[,str_detect(colnames(molecules), "x_")]
    colnames(expression) = gsub("x_(.*)", "\\1", colnames(expression))
    list(molecules = molecules,expression = expression)
  })
  
  named_list(molecules, cellinfo, expression, simulations)
}

#' @import dplyr
#' @import ggplot2
#' @import purrr
plot_simulations = function(simulations, samplingrate=0.1) {
  for (i in seq_len(length(simulations))) {
    sample = sample(c(T, F), size=nrow(simulations[[i]]$expression), T, c(samplingrate, 1-samplingrate))
    
    sample[simulations[[i]]$expression[sample,] %>% apply(1, mean) %>% {. == 0}] = FALSE
    
    print(sum(sample))
    
    simulations[[i]]$subcellinfo = simulations[[i]]$cellinfo[sample,]
    simulations[[i]]$subexpression = simulations[[i]]$expression[sample,]
    
    simulations[[i]]$expression[sample,] %>% apply(1, sd) %>% {. == 0} %>% sum %>% print
    simulations[[i]]$subexpression %>% {apply(., 1, sd) == 0} %>% sum %>% print
    
    print("---")
  }
  
  overallexpression = map(simulations, "subexpression") %>% do.call(rbind, .)
  overallexpression = overallexpression %>% set_rownames(1:nrow(overallexpression))
  #overallcellinfo = map(simulations, ~.$subcellinfo) %>% bind_rows()
  overallcellinfo = assign_progression(overallexpression, reference)
  overallcellinfo$simulationid = map(seq_len(length(simulations)), ~rep(., nrow(simulations[[.]]$expression))) %>% unlist %>% factor()
  overallcellinfo$observed = T
  
  #overallexpression = rbind(overallexpression, reference_expression)
  #overallcellinfo = bind_rows(overallcellinfo, tibble(simulationid=rep(NA, nrow(reference_cellinfo)), observed=FALSE))
  
  space =  lmds(overallexpression, 3) %>% as.data.frame() %>% bind_cols(overallcellinfo)
  space %>% filter(observed) %>% ggplot() + geom_path(aes(Comp1, Comp2, group=simulationid, color=progression)) + geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3)
  space %>% filter(observed) %>% ggplot() + geom_path(aes(Comp1, Comp2, group=simulationid, color=progression)) + geom_path(aes(Comp1, Comp2), data=space %>% filter(!observed), size=3) + viridis::scale_color_viridis(option="A")
  space %>% filter(observed) %>% ggplot() + geom_path(aes(V1, V2, group=simulationid, color=piecestateid)) + geom_path(aes(V1, V2), data=space %>% filter(!observed), size=3)
  
  for (i in unique(as.numeric(space$simulationid))) rgl::lines3d(space %>% filter(simulationid == i), col = rainbow(length(unique(space$simulationid)))[[i]])
  
}




#' Add housekeeping genes
#' @export
add_housekeeping_poisson <- function(expression, geneinfo, ngenes=200, overallaverage = mean(expression)) {
  reference_expression <- (2^SCORPIUS::ginhoux$expression)-1
  meanpoissons <- colMeans(reference_expression) %>% {./mean(reference_expression)*overallaverage}
  
  additional_expression <- purrr::map(sample(meanpoissons, ngenes), ~rpois(nrow(expression), .)) %>% 
    invoke(cbind, .) %>% 
    magrittr::set_colnames(seq_len(ngenes)+ncol(expression))
  
  geneinfo <- dplyr::bind_rows(geneinfo %>% dplyr::mutate(housekeeping=F), tibble(gene=colnames(additional_expression) %>% as.numeric(), housekeeping=T))
  
  list(expression=cbind(expression, additional_expression), geneinfo=geneinfo)
}