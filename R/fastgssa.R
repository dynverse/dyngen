# @useDynLib dyngen, .registration = TRUE
#' @useDynLib dyngen
#' 
#' @importFrom utils flush.console
#' @importFrom RcppXPtrUtils cppXPtr
fastgssa <- function(
  initial_state,
  propensity_funs,
  nu,
  final_time,
  params = NULL,
  method = ssa_direct(),
  stop_on_neg_state = TRUE,
  stop_on_neg_propensity = TRUE,
  verbose = FALSE,
  max_duration = Inf,
  console_interval = 1
) {
  # check parameters
  assert_that(
    is.numeric(initial_state),
    !is.null(names(initial_state)),
    is.character(propensity_funs),
    dynutils::is_sparse(nu),
    is.numeric(final_time),
    is.numeric(console_interval),
    is.logical(verbose),
    is(method, "fastgssa::ssa_method"),
    length(initial_state) == nrow(nu),
    length(propensity_funs) == ncol(nu),
    is.numeric(params),
    length(params) == 0 || !is.null(names(params)),
    !any(c("state", "time", "params") %in% names(params)),
    !any(c("state", "time", "params") %in% names(initial_state)),
    !any(duplicated(c(names(initial_state), names(params))))
  )
  
  state <- initial_state
  time <- 0
  
  prop_fun_comp <- compile_propensity_functions(propensity_funs, state, params, env = environment())
  
  output <- dyngen:::simulate(
    transition_fun = prop_fun_comp,
    ssa_alg = dyngen:::make_ssa_em(h = .01, noise_strength = 2),
    initial_state = initial_state,
    params = params,
    nu = as.matrix(nu),
    final_time = final_time,
    max_duration = max_duration,
    stop_on_neg_state = stop_on_neg_state,
    stop_on_neg_propensity = stop_on_neg_propensity,
    verbose = verbose,
    console_interval = console_interval
  )
  
  
}

#' @importFrom stringr str_count str_replace_all str_extract_all str_replace
#' @importFrom Rcpp sourceCpp
compile_propensity_functions <- function(propensity_funs, state, params, env = parent.frame()) {
  buffer_size <- max(str_count(propensity_funs, "="))
  rcpp_prop_funs <- map_chr(
    propensity_funs,
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
  
  buffers <- rcpp_prop_funs %>% str_replace(";[^=]*$", ";") %>% {ifelse(grepl("=", .), ., "")} %>% str_replace_all("([^;]*;)", "    \\1\n")
  calculations <- rcpp_prop_funs %>% str_replace("^.*;", "") %>% {paste0("  transition_rates[", seq_along(.)-1, "] = ", ., ";\n")}

  rcpp_code <- paste0("void calculate_transition_rates(
  const NumericVector& state,
  const NumericVector& params,
  const NumericVector& time,
  NumericVector& transition_rates
) {
",
  ifelse(buffer_size == 0, "", paste0("    NumericVector buffer = no_init(", buffer_size, ");\n")),
  paste(paste0(buffers, calculations), collapse = ""),
"}
")
  tmpdir <- dynutils::safe_tempdir("fastgssa")
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))
  
  transition_functon <- RcppXPtrUtils::cppXPtr(rcpp_code, cacheDir = tmpdir)

  transition_functon
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
#' @export
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
#' @export
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
