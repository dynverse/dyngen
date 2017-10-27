#' Simulation one individual cells
#' 
#' @param model The model to simulation
#' @param timeofsampling `NULL` to return expression at all time steps, a double to return expression at a particular time step
#' @param totaltime The total simulation time
#' @param burntime A burnin time
#' @param ssa_algorithm Which SSA algorithm to use
#' 
#' @import fastgssa
simulate_cell <- function(model, timeofsampling = NULL, totaltime = 10, burntime = 2, ssa_algorithm = fastgssa::ssa.em(noise_strength = 4)) {
  if (burntime > 0) {
    nus_burn <- model$nus
    nus_burn[setdiff(rownames(nus_burn), model$burn_variables),] <- 0
    
    # burn in
    out <- fastgssa::ssa(
      model$initial_state, 
      model$formulae_strings,
      nus_burn,
      burntime, 
      model$params, 
      method = ssa_algorithm,
      recalculate.all = TRUE, 
      stop.on.negstate = FALSE,
      stop.on.propensity = FALSE
    )
    output <- process_ssa(out)
    initial_state_after_burn <- output$molecules[nrow(output$molecules), names(model$initial_state)]
  } else {
    initial_state_after_burn <- model$initial_state
  }

  # determine total time to simulate
  if (!is.null(timeofsampling)) {
    time <- timeofsampling
  } else {
    time <- totaltime
  }
  
  # actual simulation
  out <- fastgssa::ssa(
    initial_state_after_burn, 
    model$formulae_strings, 
    model$nus, time,
    model$params, 
    method = ssa_algorithm,
    recalculate.all = TRUE, 
    stop.on.negstate = FALSE,
    stop.on.propensity = FALSE
  )
  output <- process_ssa(out)
  molecules <- output$molecules
  times <- output$times
  
  # return either the full molecules matrix, or the molecules at timeofsampling
  if (is.null(timeofsampling)) {
    rownames(molecules) <- 1:nrow(molecules)
    list(molecules = molecules, times = times)
  } else {
    molecules %>% slice(n())
  }
}

#' @importFrom utils tail
process_ssa <- function(out, starttime=0) {
  final.time <- out$args$final.time
  data <- out$timeseries
  lastrow <- utils::tail(data, n=1)
  lastrow$t <- final.time
  data[nrow(data)+1,] <- lastrow
  data <- data[data$t<=final.time,]
  times <- data$t + starttime
  data <- as.matrix(data[,-1])
  
  return(list(times=times, molecules=data))
}

#' Simulate multiple cells
#' 
#' @param model The model to simulate
#' @param burntime The burn-in time before sampling
#' @param totaltime The total simulation time
#' @param nsimulations The number of simulations
#' @param local Whether or not to use \code{\link[PRISM]{qsub_lapply}}
#' @param ssa_algorithm Which GSSA algorithm to use
#' 
#' @importFrom PRISM qsub_lapply
#' @importFrom pbapply pblapply
#' @export
simulate_multiple <- function(model, burntime, totaltime, nsimulations = 16, local=FALSE, ssa_algorithm = fastgssa::ssa.em(noise_strength=4)) {
  force(model) # force the evaluation of the model argument, as the qsub environment will be empty except for existing function arguments
  if(!local) {
    multilapply = function(x, fun) {PRISM::qsub_lapply(x, fun, qsub_environment = list2env(list()))}
  } else {
    multilapply = function(x, fun) {pbapply::pblapply(x, fun, cl = 8)}
  }
  
  simulations = multilapply(seq_len(nsimulations), function(i) {
    cell = simulate_cell(model, totaltime=totaltime, burntime=burntime, ssa_algorithm=ssa_algorithm)
    
    rownames(cell$molecules) = paste0(i, "_", seq_len(nrow(cell$molecules)))
    
    expression = cell$molecules[,str_detect(colnames(cell$molecules), "x_")]
    colnames(expression) = gsub("x_(.*)", "\\1", colnames(expression))
    stepinfo = tibble(step_id = rownames(cell$molecules), step=seq_along(cell$times), simulationtime=cell$times, simulation_id=i)
    
    lst(molecules=cell$molecules, stepinfo=stepinfo, expression=expression)
  })
  
  molecules <- map(simulations, "molecules") %>% do.call(rbind, .)
  expression <- map(simulations, "expression") %>% do.call(rbind, .)
  stepinfo <- map(simulations, "stepinfo") %>% do.call(bind_rows, .)
  
  lst(molecules, expression, stepinfo)
}

# Emrna doesnt get updated, this function does not work
#' @importFrom zoo rollapply
#' @importFrom stats quantile sd
estimate_convergence = function(nruns=8, cutoff=0.15, verbose=F, totaltime=15) {
  convergences = parallel::mclapply(1:nruns, function(i) {
    E = expression_one_cell(burntime, totaltime)
    Emrna = E[str_detect(featureNames(E), "x_")]
    
    exprs(Emrna) %>% t %>% zoo::rollmean(10) %>% zoo::rollapply(100, stats::sd, fill=0) %>% t %>% apply(2, stats::quantile, probs=0.9)
  }, mc.cores=8) %>% do.call(rbind, .)
  
  timepoints = convergences %>% apply(1, function(convergences) {convergences %>% {.>=cutoff} %>% rev %>% cumsum %>% {. == 0} %>% rev %>% which() %>% first() %>% phenoData(Emrna)$time[.]}) # the first timepoint at which all following timepoints are below cutoff
  
  if (verbose) {
    plot = convergences %>%
      reshape2::melt(varnames=c("run", "time"), value.name="convergence") %>% 
      mutate(run=factor(run)) %>% 
      mutate(time=as(phenoData(Emrna), "data.frame")[,"time"][time]) %>%
      qplot(time, convergence, data=., group=run, color=run, geom="line") + 
      geom_hline(aes(yintercept=cutoff)) + 
      geom_vline(aes(xintercept = timepoint, color=run), data=tibble(timepoint=timepoints, run=factor(1:nruns)))
    print(plot)
  }
  timepoints
}