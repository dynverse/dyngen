#' @export
tf_random <- function(
  main_targets_ratio = 0.05,
  ngenes_per_module_sampler = function(n_features, n_modules) sample(1:10, n_modules, replace = TRUE),
  edge_retainment = function(n) max(c(round(n/2), 1))
) {
  lst(
    main_targets_ratio,
    ngenes_per_module_sampler,
    edge_retainment
  )
}

#' ... todo
#' @export
target_realnet <- function(
  realnet_name = realnets$name,
  damping = 0.05,
  ntargets_sampler = function(n_features, n_regulators) sample(20:100, 1)
) {
  realnet_name <- match.arg(realnet_name)
  
  data(realnets, package = "dyngen")
  assert_that(realnet_name %all_in% realnets$name)
  
  realnet_url <- realnets$url[[match(realnet_name, realnets$name)]]
  
  lst(
    realnet_name,
    realnet_url,
    damping,
    ntargets_sampler
  )
}
