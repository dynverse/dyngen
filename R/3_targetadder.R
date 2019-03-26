#' ... todo
#' @export
targetadder_realnet <- function(
  realnet_name = realnets$name,
  damping = 0.05,
  ntargets_sampler = function(n_features, n_regulators) sample(20:100, 1),
  gene_name_generator = function(i) paste0("G", i)
) {
  realnet_name <- match.arg(realnet_name)
  
  data(realnets, package = "dyngen")
  assert_that(realnet_name %all_in% realnets$name)
  
  realnet_url <- realnets$url[[match(realnet_name, realnets$name)]]
  
  lst(
    realnet_name,
    realnet_url,
    damping,
    ntargets_sampler,
    gene_name_generator
  )
}