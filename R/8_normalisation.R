#' @inheritParams dynnormaliser::normalise_filter_counts
#' 
#' @export
normalisation_default <- function(
  filter_cells = TRUE,
  filter_features = TRUE,
  filter_hvg = TRUE,
  normalisation = "scran_size_features",
  has_spike = FALSE,
  nmads = 3,
  min_ave_expression = 0.02,
  hvg_fdr = 0.05,
  hvg_bio = .5,
  min_variable_fraction = 0.15
) {
  lst(
    filter_cells,
    filter_features,
    filter_hvg,
    normalisation,
    has_spike,
    nmads,
    min_ave_expression,
    hvg_fdr,
    hvg_bio,
    min_variable_fraction
  )
}
