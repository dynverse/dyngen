#' A bounded version of rnorm
#' 
#' @inheritParams stats::rnorm
#' @param min lower limits of the distribution.
#' @param max upper limits of the distribution.
#'
#' @importFrom stats pnorm qnorm
#' 
#' @return Generates values with rnorm, bounded by \[min, max\]
#' 
#' @export 
#' 
#' @examples
#' rnorm_bounded(10)
rnorm_bounded <- function(n, mean = 0, sd = 1, min = -Inf, max = Inf) {
  unif_min <- pnorm(min, mean = mean, sd = sd)
  unif_max <- pnorm(max, mean = mean, sd = sd)
  quan <- runif(n, unif_min, unif_max)
  qnorm(quan, mean = mean, sd = sd)
}

#' A subrange version of runif
#' 
#' Will generate numbers from a random subrange within the given range. 
#' For example, if \[min, max\]` is set to \[0, 10\], this function
#' could decide to generate `n` numbers between 2 and 6.
#' 
#' @param n Number of observations
#' @param min Lower limits of the distribution.
#' @param max Upper limits of the distribution.
#' 
#' @importFrom stats runif
#' 
#' @return Generates values with runif, bounded by a range drawn from `sort(runif(2, min, max))`.
#' 
#' @export
#' 
#' @examples 
#' runif_subrange(20, 0, 10)
runif_subrange <- function(n, min, max) {
  range <- sort(runif(2, min, max))
  runif(n, range[[1]], range[[2]])
}