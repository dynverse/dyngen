#' @export
generate_toy_milestone_network <- function(ti_type=c("linear", "bifurcating", "cycle")) {
  if(!(ti_type %in% c("linear", "bifurcating", "cycle"))) stop("invalid ti_type")
  
  if(identical(ti_type, "linear")) {
    tibble(from="M1", to="M2", length=1)
  } else if (identical(ti_type, "bifurcating")) {
    tibble(from=c("M1", "M2", "M2"), to=c("M2", "M3", "M4"), length=1)
  } else if (identical(ti_type, "cycle")) {
    tibble(from=c("M1", "M2", "M3"), to=c("M2", "M3", "M1"), length=1)
  }
}