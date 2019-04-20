
fastgssa <- function(
  initial.state,
  propensity.funs,
  nu,
  final.time,
  params = NULL,
  method = ssa_direct(),
  stop.on.neg.state = TRUE,
  stop.on.neg.propensity = TRUE,
  verbose = FALSE,
  max.duration = Inf,
  time.next.console = 1
) {
  # check parameters
  assert_that(
    is.numeric(initial.state),
    !is.null(names(initial.state)),
    is.character(propensity.funs),
    dynutils::is_sparse(nu),
    is.numeric(final.time),
    is.numeric(time.next.console),
    is.logical(verbose),
    is(method, "fastgssa::ssa_method"),
    length(initial.state) == nrow(nu),
    length(propensity.funs) == ncol(nu),
    is.numeric(params),
    length(params) == 0 || !is.null(names(params)),
    !any(c("state", "time", "params") %in% names(params)),
    !any(c("state", "time", "params") %in% names(initial.state)),
    !any(duplicated(c(names(initial.state), names(params))))
  )
  
  state <- initial.state
  time <- 0
  
  # compile propensity funs
  parse_propensity_funs(propensity.funs, state, params, env = environment())
  
  transition_rates <- activation(state, params, time) %>% set_names(names(propensity.funs))
  
  # Initialise output
  output <- list()
  output[[length(output) + 1]] <- tibble(time, state = list(state), transition_rates = list(transition_rates))
  
  # # Start the timer
  time.next.console <- time.curr <- time.start <- Sys.time()
  
  if (verbose) {
    cat("Running ", method$name, " method with console output every ", console.interval, " time step\n", sep = "")
    cat("Start wall time: ", format(time.start), "\n" , sep = "")
    utils::flush.console()
  }
  
  while (time < final.time && (time.curr - time.start) <= max.duration)  {
    time.curr <- Sys.time()
    
    if (verbose && time.next.console <= time.curr) {
      cat(format(time.curr), " | time = ", time, " : ", paste(round(state, 2), collapse = ","), "\n", sep = "")
      utils::flush.console()
      time.next.console <- time.next.console + console.interval
    }
    
    method_out <- method$fun(state = state, transition_rates = transition_rates, nu = nu)
    
    state_prev <- state
    
    time <- time + method_out$dtime
    state <- state + method_out$dstate
    
    # Check that no states are negative (can occur in some tau-leaping methods)
    invalid_ix <- is.na(state) | state < 0
    if (any(invalid_ix)) {
      if (stop.on.neg.state) {
        stop("state vector contains negative values\n", paste(names(state)[invalid_ix], collapse = ", "))
      } else {
        state[invalid_ix] <- 0
      }
    }
    
    transition_rates <- activation(state, params, time) %>% set_names(names(propensity.funs))
    
    invalid_ix <- is.na(transition_rates) | transition_rates < 0
    if (any(invalid_ix)) {
      if (stop.on.neg.propensity) {
        stop("transition rate contains negative values\n", paste(names(transition_rates)[invalid_ix], collapse = ", "))
      } else {
        transition_rates[invalid_ix] <- 0
      }
    }
    
    output[[length(output) + 1]] <- tibble(time, state = list(state), transition_rates = list(transition_rates))
  }
  
  # Display the last time step on the console
  if (verbose) {
    cat("time = ", time, " : ", paste(round(state, 2), collapse = ","), "\n", sep = "")
    utils::flush.console()
  }
  
  # Stop timer
  time.end <- Sys.time()
  elapsed.time <- difftime(time.start, time.end, units = "secs")
  
  # Calculate some stats for the used method
  stats <- tibble(
    method             = method$name,
    final.time.reached = time >= final.time,
    extinction         = all(state == 0),
    negative.state     = any(state < 0),
    zero.prop          = all(transition_rates == 0),
    max.duration       = elapsed.time >= max.duration,
    start.time         = time.start,
    end.time           = time.end,
    elapsed.time       = elapsed.time
  )
  if (verbose) {
    print(stats)
  }
  
  
  # Output results
  list(
    timeseries = bind_rows(output),
    stats = stats
  )
  
}

#' @importFrom stringr str_count str_replace_all str_extract_all str_replace
parse_propensity_funs <- function(propensity.funs, state, params, env = parent.frame()) {
  buffer_size <- max(str_count(propensity.funs, "="))
  rcpp_prop_funs <- map_chr(
    propensity.funs,
    function(prop_fun) {
      buffer_names <- prop_fun %>% str_extract_all("[A-Za-z_0-9]* *=") %>% first() %>% str_replace_all("[ =]", "")
      prop_split <- gsub("([A-Za-z][A-Za-z0-9_]*)", " \\1 ", prop_fun) %>% strsplit(" ") %>% first() %>% discard(~ .=="")
      
      state_match <- match(prop_split, names(state))
      state_ix <- which(!is.na(state_match))
      prop_split[state_ix] <- paste0("state[", state_match[state_ix] - 1, "]")
      
      params_match <- match(prop_split, names(params))
      params_ix <- which(!is.na(params_match))
      prop_split[params_ix] <- paste0("params[", params_match[params_ix] - 1, "]")
      
      buffer_match <- match(prop_split, buffer_names)
      buffer_ix <- which(!is.na(buffer_match))
      prop_split[buffer_ix] <- paste0("buffer[", buffer_match[buffer_ix] - 1, "]")
      
      paste(prop_split, collapse = "")
    }
  )
  
  buffers <- rcpp_prop_funs %>% str_replace(";[^=]*$", ";") %>% {ifelse(grepl("=", .), ., "")} %>% str_replace_all("([^;]*;)", "  \\1\n")
  calculations <- rcpp_prop_funs %>% str_replace("^.*;", "") %>% {paste0("  out[", seq_along(.)-1, "] = ", ., ";\n")}
  
  rcpp_code <- paste0(
    "NumericVector activation(NumericVector state, NumericVector params, double time) {\n",
    ifelse(buffer_size == 0, "", paste0("  NumericVector buffer = no_init(", buffer_size, ");\n")),
    "  NumericVector out = no_init(", length(propensity.funs), ");\n",
    paste(paste0(buffers, calculations), collapse = ""),
    "  return out;\n",
    "}\n"
  )
  
  Rcpp::cppFunction(rcpp_code, env = env)
}


ssa_method <- function(name, fun) {
  l <- lst(
    name,
    fun
  )
  class(l) <- "fastgssa::ssa_method"
  l
}

#' Euler-Maruyama method (EM)
#'
#' Euler-Maruyama method implementation
#'
#' @param h h parameter
#' @param noise_strength noise_strength parameter
#'
#' @importFrom stats rnorm
#'
# @export
ssa_em <- function(h = 0.01, noise_strength = 2) {
  ssa_method(
    name = "em",
    fun = function(state, transition_rates, nu) {
      dtime <- h
      
      dstate <- (nu %*% transition_rates)[,1] * h + sqrt(abs(state)) * noise_strength * stats::rnorm(length(state), 0, h)
      
      list(
        dtime = dtime,
        dstate = dstate,
        j = which(dstate != 0)
      )
    }
  )
}


#' Direct method (D)
#'
#' Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#'
#' @return an object of to be used by \code{\link{ssa}}.
# @export
ssa_direct <- function() {
  ssa_method(
    name = "direct",
    fun = function(state, transition_rates, nu) {
      
      j <- sample.int(length(transition_rates), size = 1, prob = transition_rates)
      dstate <- nu[,j]
      
      dtime <- -log(stats::runif(1)) / sum(transition_rates)
      
      list(
        dtime = dtime,
        dstate = dstate,
        j = which(dstate != 0)
      )
    }
  )
}
