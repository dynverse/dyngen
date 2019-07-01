# why does rlang::`%|%` not work like this?
`%|||%` <- function(x, y) {
  ifelse(is.na(x), y, x)
}