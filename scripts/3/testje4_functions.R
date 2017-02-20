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
    out <- ssa_deterministic(initial.state,formulae.strings, formulae.nus.burn, burntime, params)
  }
  output = process_ssa(out)
  initial.state.burn = output$expression[nrow(output$expression), ] %>% abs
  
  print(initial.state)
  print(initial.state.burn)
  
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
    out <- ssa_deterministic(initial.state.burn, formulae.strings, formulae.nus, time, params)
  }
  output = process_ssa(out)
  expression = output$expression
  times = output$times
  #burnin = times > burntime
  #expression = methods::rbind2(output$expression[max(1, findInterval(burntime, output$times)),], output$expression[burnin,])
  #times = c(0,times[burnin] - burntime)
  
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