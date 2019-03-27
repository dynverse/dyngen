#' @export
goldstandard_default <- function(
  max_path_length = 10,
  reference_length = 100,
  smooth_window = 50
) {
  lst(
    max_path_length,
    reference_length,
    smooth_window
  )
}