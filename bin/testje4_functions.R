
simulate_cell = function(timeofsampling=NULL) {
  requireNamespace("fastgssa")
  
  if (!is.null(timeofsampling)) {
    time = timeofsampling
  } else {
    time = totaltime
  }
  
  out <- fastgssa::ssa(initial.state,formulae.strings,formulae.nus, burntime+time, params, method=fastgssa::ssa.direct(), recalculate.all =F)
  output = process_ssa(out)
  expression = output$expression
  times = output$times
  burnin = times > burntime
  expression = methods::rbind2(output$expression[max(1, findInterval(burntime, output$times)),], output$expression[burnin,])
  times = c(0,times[burnin] - burntime)
  
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