colour_brighten <- function(cols, pcts) {
  rgb_mat <- 1 - grDevices::col2rgb(cols) / 255
  mat <- (1 - sweep(rgb_mat, 2, pcts, "*"))
  grDevices::rgb(mat[1,], mat[2,], mat[3,])
}