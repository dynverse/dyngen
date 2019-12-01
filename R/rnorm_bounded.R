#' A bounded version of rnorm
#' 
#' @inheritParams stats::rnorm
#' @param min lower limits of the distribution.
#' @param max upper limits of the distribution.
#' 
#' @export 
rnorm_bounded <- function(n, mean = 0, sd = 1, min = -Inf, max = Inf) {
  unif_min <- pnorm(min, mean = mean, sd = sd)
  unif_max <- pnorm(max, mean = mean, sd = sd)
  quan <- runif(n, unif_min, unif_max)
  qnorm(quan, mean = mean, sd = sd)
}