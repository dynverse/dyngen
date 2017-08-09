parse.propensity.functions <- function(propensity.funs, var.names) {
  names(propensity.funs) <- NULL
  lapply(propensity.funs, function(string) {
    aag <- strsplit(gsub("([A-Za-z][A-Za-z0-9_]*)", " \\1 ", string), " ")[[1]]
    aag.match <- match(aag, var.names)
    ix <- !is.na(aag.match)
    found.var.names <- unique(aag[ix])
    aag[ix] <- paste0("vec[[", aag.match[ix], "]]")
    new.aa <- paste0("function(vec) ", paste(aag, collapse=""))
    aa.fun <- eval(parse(text = new.aa))
    list(fun = aa.fun, var.names = found.var.names)
  })
}

process.parms <- function(parms) {
  if (is.null(parms)) {
    t(c(t = 0))
  } else if (is.vector(parms)) {
    t(c(t = 0, parms))
  } else if (is.data.frame(parms)) {
    as.matrix(parms)
  } else if (is.matrix(parms)) {
    parms
  } else {
    stop(sQuote("parms"), " needs to be NULL, a vector, a data frame or a matrix")
  }
}


#' Invoking the stochastic simulation algorithm
#'
#' Main interface function to the implemented \acronym{SSA} methods. Runs a
#' single realization of a predefined system.
#'
#' Although \code{ssa} can be invoked by only specifying the system arguments
#' (initial state vector \code{initial.state}, propensity vector \code{propensity.funs}, state-change
#' matrix \code{nu}), the final time (\code{final.time}), and the \acronym{SSA} method
#' to use, substantial improvements in speed and accuracy can be obtained by
#' adjusting the additional (and optional) \code{ssa} arguments. By default
#' \code{ssa} (tries to) use conservative default values for the these
#' arguments, prioritizing computational accuracy over computational speed.
#' These default values are, however, \strong{not} fool proof for the
#' approximate methods, and occasionally one will have to hand tweak them in
#' order for a stochastic model to run appropriately.
#'
#' @param initial.state numerical vector of initial states where the component elements
#'   must be named using the same notation as the corresponding state variable in
#'   the propensity vector, \code{a}.
#' @param propensity.funs character vector of propensity functions where state variables
#'   correspond to the names of the elements in \code{initial.state}.
#' @param nu numerical matrix of change if the number of individuals in each
#'   state (rows) caused by a single reaction of any given type (columns).
#' @param final.time final time.
#' @param parms named vector of model parameters.
#' @param method which method to use. Must be one of: \code{\link{ssa.direct}},
#'   \code{\link{ssa.btl}}, \code{\link{ssa.etl}}, or \code{\link{ssa.otl}}
#' @param stop.on.negstate boolean object indicating if negative state
#'   values should be ignored (this can occur in the \code{etl} method). If
#'   \code{stop.on.negstate=TRUE} the simulation finishes gracefully when
#'   encountering a negative population size (i.e. does not throw an error). If
#'   \code{stop.on.negstate=FALSE} the simulation stops with an error message
#'   when encountering a negative population size.
#' @param census.interval (approximate) interval between recording the state of
#'   the system. If \code{census.interval=0} \eqn{(t,x)} is recorded at each time
#'   step (or tau-leap). If \code{census.interval=Inf} only
#'   \eqn{(t_0,x_0)}{(t0,initial.state)} and \eqn{(t_f,x_t)}{(final.time,xf)} is recorded. Note, the
#'   size of the time step (or tau-leaps) ultimately limits the interval between
#'   subsequent recordings of the system state since the state of the system
#'   cannot be recorded at a finer time interval the size of the time steps (or
#'   tau-leaps).
#' @param verbose If \code{TRUE}, intermediary information pertaining to the simulation will be displayed at each step.
#'   If \code{verbose} is a numeric, it will only be displayed about every \code{verbose} seconds.
#'   \strong{Verbose runs drastically slows down simulations.}
#' @param max.duration maximum wall time duration (in seconds) that the
#'   simulation is allowed to run for before terminated. This option is usefull,
#'   in particular, for systems that can end up growing uncontrolably.
#'
#' @return Returns a list object with the following elements,
#'   \item{timeseries}{a data frame of the simulation time series where the first column is the time vector and subsequent columns are the state frequencies.}
#'   \item{stats}{a data frame containing several simulation statistics statistics.}
#'   \item{args}{the original parameters passed to \code{\link{ssa}}}
#'
#' @note Selecting the appropriate \acronym{SSA} method is a trade-off between
#'   computational speed, accuracy of the results, and which \acronym{SSA}
#'   actually works for a given scenario. This depends on the characteristics of
#'   the defined system (e.g. number of reaction channels, number of species, and
#'   the absolute and relative magnitude of the propensity functions).
#'   \strong{Not all methods are appropriate for all models.} When selecting a
#'   \acronym{SSA} method all of these factors have to be taken into
#'   consideration. The various tau-leap methods accept a number of additional
#'   arguments. While the default values of these arguments may work for some
#'   scenarios they may have to be adjusted for others. The default values for
#'   the tau-leap methods are conservative in terms of computational speed and
#'   substantial increase in efficiency may be gained by optimizing their values
#'   for a specific system.
#' @section Preparing a run: In order to invoke \acronym{SSA} the stochastic
#'   model needs at least four components, the initial state vector (\code{initial.state}),
#'   state-change matrix (\code{nu}), propensity vector (\code{a}), and the final
#'   time of the simulation (\code{final.time}). The initial state vector defines the
#'   population sizes in all the states at \eqn{t=0}, e.g. for a system with two
#'   species \code{X1} and \code{X2} where both have an initial population size
#'   of 1000 the initial state vector is defined as \code{initial.state <-
#'   c(X1=1000,X2=1000)}. The elements of the vector have to be labelled using
#'   the same notation as the state variables used in the propensity functions.
#'   The state-change matrix defines the change in the number of individuals in
#'   each state (rows) as caused by one reaction of a given type (columns). For
#'   example, the state-change matrix for system with the species \eqn{S_1}{S1}
#'   and \eqn{S_2}{S2} with two reactions \deqn{S_1
#'   \stackrel{c_1}{\longrightarrow} S_2}{S1 --c1--> S2} \deqn{S_2
#'   \stackrel{c_2}{\longrightarrow} 0}{S2 --c2--> 0}
#'
#'   is defined as \code{nu <- matrix(c(-1,0,+1,-1),nrow=2,byrow=TRUE)} where
#'   \eqn{c_1}{c1} and \eqn{c_2}{c2} are the per capita reaction probabilities.
#'   The propensity vector, \code{a}, defines the probabilities that a particular
#'   reaction will occur over the next infinitesimal time interval \eqn{\left[
#'   t,t+dt \right]}{[t,t+dt]}. For example, in the previous example the
#'   propensity vector is defined as \code{a <- c("c1*X1","c2*X2")}. The
#'   propensity vector consists of character elements of each reaction's
#'   propensity function where each state variable requires the corresponding
#'   named element label in the initial state vector (\code{initial.state}).
#'
#' @author Robrecht Cannoodt
#' @seealso \link{fastgssa-package}, \code{\link{ssa.direct}}, \code{\link{ssa.etl}}, \code{\link{ssa.btl}}, \code{\link{ssa.otl}}
#'
#' @keywords misc datagen ts
#' @examples
#'
#' ## Irreversible isomerization
#' ## Large initial population size (X=1000)
#' \dontrun{
#' parms <- c(c=0.5)
#' initial.state  <- c(X=10000)
#' a   <- c("c*X")
#' nu  <- matrix(-1)
#' out <- ssa(initial.state,a,nu,parms,final.time=10,simName="Irreversible isomerization") # Direct method
#' plot(out$data[,1],out$data[,2]/10000,col="red",cex=0.5,pch=19)
#' }
#'
#' ## Smaller initial population size (X=100)
#' \dontrun{
#' initial.state  <- c(X=100)
#' out <- ssa(initial.state,a,nu,parms,final.time=10) # Direct method
#' points(out$data[,1],out$data[,2]/100,col="green",cex=0.5,pch=19)
#' }
#'
#' ## Small initial population size (X=10)
#' \dontrun{
#' initial.state  <- c(X=10)
#' out <- ssa(initial.state,a,nu,parms,final.time=10) # Direct method
#' points(out$data[,1],out$data[,2]/10,col="blue",cex=0.5,pch=19)
#' }
#'
#' ## Logistic growth
#' \dontrun{
#' parms <- c(b=2, d=1, K=1000)
#' initial.state  <- c(N=500)
#' a   <- c("b*N", "(d+(b-d)*N/K)*N")
#' nu  <- matrix(c(+1,-1),ncol=2)
#' out <- ssa(initial.state,a,nu,parms,final.time=10,method="D",max.duration=5,simName="Logistic growth")
#' ssa.plot(out)
#' }
#'
#' ## Kermack-McKendrick SIR model
#' \dontrun{
#' parms <- c(beta=0.001, gamma=0.1)
#' initial.state  <- c(S=499,I=1,R=0)
#' a   <- c("beta*S*I","gamma*I")
#' nu  <- matrix(c(-1,0,+1,-1,0,+1),nrow=3,byrow=TRUE)
#' out <- ssa(initial.state,a,nu,parms,final.time=100,simName="SIR model")
#' ssa.plot(out)
#' }
#'
#' ## Lotka predator-prey model
#' \dontrun{
#' parms <- c(c1=10, c2=.01, c3=10)
#' initial.state  <- c(Y1=1000,Y2=1000)
#' a   <- c("c1*Y1","c2*Y1*Y2","c3*Y2")
#' nu  <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#' out <- ssa(initial.state,a,nu,parms,final.time=100,method="ETL",simName="Lotka predator-prey model")
#' ssa.plot(out)
#' }
#'
#' @export ssa
ssa_deterministic <- function(
  initial.state,
  propensity.funs,
  nu,
  final.time,
  parms = NULL,
  method = fastgssa::ssa.direct(),
  stop.on.negstate = TRUE,
  census.interval = 0,
  verbose = FALSE,
  max.duration = Inf,
  recalculate.all = F
) {
  # Take a snapshot of all the options so they can be saved later
  args <- list(initial.state = initial.state, propensity.funs = propensity.funs, nu = nu, final.time = final.time, parms = parms, method.name = method$name, method.args = method$params)
  
  # construct state space
  parms <- process.parms(parms)
  parms.time <- parms[, 1]
  parms <- parms[, -1, drop=F]
  parms.index <- 1
  
  x <- initial.state
  p <- parms[parms.index,]
  state.env <- c(x, p)
  varname.ix <- match(names(x), names(state.env))
  parname.ix <- match(names(p), names(state.env))
  
  # Initialize time-related counters
  if (is.numeric(verbose)) {
    console.interval <- verbose
    verbose <- T
  } else if (verbose) {
    console.interval <- 0
  } else {
    console.interval <- Inf
  }
  t <- 0
  t.next.census <- census.interval
  t.next.console <- 0
  
  # Initialise output
  timeseries.output <- vector('list', 1000)
  timeseries.output[[1]] <- c(t = t, x)
  timeseries.index <- 2
  step.sizes <- c()
  
  # Parse the propensity functions
  parsed.pfs <- parse.propensity.functions(propensity.funs, names(state.env))
  parsed.pf.funs <- lapply(parsed.pfs, function(pf) pf$fun)
  varname.map <- lapply(names(x), function(varname) {
    setNames(which(sapply(parsed.pfs, function(pf) varname %in% pf$var.names)), NULL)
  })
  
  # Evaluate initial transition rates
  a <- sapply(parsed.pf.funs, function(f) f(state.env))
  print(a)
  if (any(is.na(a))) (print(propensity.funs[is.na(a)]))
  if (any(a < 0)) stop("negative propensity function")
  
  #############################################################################
  # We are ready to roll, start the simulation loop...
  # Continue the simulation loop as long as we have not reached the end time,
  # all of the populations have not gone extincs, as long as no population is
  # negative (occasinal by-product of the ETL method), and as long as at
  # least one propensity function is larger than zero
  #############################################################################
  
  # Start the timer
  time.start <- Sys.time()
  elapsed.time <- 0
  
  if (verbose) {
    cat("Running ", method$name, " method with console output every ", console.interval, " time step\n", sep="")
    cat("Start wall time: ", format(time.start), "...\n" , sep = "")
    flush.console()
  }
  
  # Prepare method functions
  out.fun <- method$fun
  method.state <- method$initial.method.state
  method.params <- method$params
  
  while ( t < final.time  &&  elapsed.time <= max.duration )  {
    
    if (verbose && t.next.console <= t) {
      cat("(", elapsed.time, "s) t=", t, " : ", paste(x, collapse=","), "\n", sep="")
      flush.console()
      t.next.console <- t.next.console + console.interval
    }
    
    h = 0.001
    noise_strength = 5
    
    prevx = x
    preva = a
    
    #print(x)
    
    # euler method
    # https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method
    
    t = t + h
    x = x + rowSums(t(t(nu) * a)) * h + sqrt(abs(x)) * noise_strength * rnorm(length(x), 0, h) # this is the fastest way http://stackoverflow.com/questions/3643555/multiply-rows-of-matrix-by-vector
    
    if (any(is.na(x))) {print(propensity.funs[is.infinite(preva)]);print(preva);print(prevx);print(paste0("time", t));print("na ----");break}
    
    #out <- do.call(out.fun, c(list(x = x, a = a, nu = nu, method.state = method.state), method.params))
    
    # Record the simulation state
    if (t.next.census <= t) {
      timeseries.output[[timeseries.index]] <- c(t = t, x)
      timeseries.index <- timeseries.index + 1
      t.next.census <- t + census.interval
      
      # Expand the time series list if necessary
      if (timeseries.index > length(timeseries.output)) {
        new.ts <- vector(mode = "list", length = length(timeseries.output) * 2)
        new.ts[seq_along(timeseries.output)] <- timeseries.output
        timeseries.output <- new.ts
      }
    }
    
    # Update values for x
    state.env[varname.ix] <- x
    
    # Evaluate the transition rates for the next time step
    a <- sapply(parsed.pf.funs, function(f) f(state.env))
    
    time.end <- Sys.time()
    elapsed.time <- difftime(time.end, time.start, units = "secs")
  }
  
  # Display the last time step on the console
  if (verbose) {
    cat("t=", t, " : ", paste(x, collapse = ","), "\n", sep = "")
    flush.console()
  }
  
  # Record the final state of the system
  timeseries.output[[timeseries.index]] <- c(t = t, x)
  timeseries <- data.frame(do.call("rbind", timeseries.output[seq_len(timeseries.index)]))
  
  # Stop timer
  time.end <- Sys.time()
  
  # Calculate some stats for the used method
  stats <- data.frame(
    method             = method$name,
    final.time.reached = t >= final.time,
    extinction         = all(x == 0),
    negative.state     = any(x < 0),
    zero.prop          = all(a == 0),
    max.duration       = elapsed.time >= max.duration,
    start.time         = time.start,
    end.time           = time.end,
    elapsed.wall.time  = elapsed.time
  )
  if (verbose) {
    cat("final time = ", t, "\n", sep="")
    print(stats)
  }
  
  # Output results
  list(
    timeseries = timeseries,
    stats = stats,
    args  = args
  )
}

function(x, a, nu, method.state) {
  j    <- sample(seq_along(a), size=1, prob=a)
  nu_j <- nu[,j]
  tau  <- -log(runif(1))/sum(a)
  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method.state = method.state)
}



#ssa_deterministic(initial.state,formulae.strings,formulae.nus, 1, params)
